#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从已保存100³热端场重启冷却，并在Ag有液相时限制最大步长。"""
from __future__ import annotations

import csv
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
PHASE = PROJECT / "phase_change_validation"
SOLVER = PROJECT / "solver_optimization"
for _path in (str(SOLVER), str(PHASE), str(PROJECT)):
    if _path in sys.path:
        sys.path.remove(_path)
    sys.path.insert(0, _path)
sys.path.insert(0, "/tmp/solver-opt-deps")

from checkpointing import CheckpointStore
from diagnostics import phase_metrics
from implicit_enthalpy_fvm import BoundarySet
from material_database_v2.fast_materials import FastPankratzMaterials
from material_model import CONFIG
from optimized_solver import NumbaEnthalpyFVM
from profiled_advance import advance_profiled


def main() -> int:
    import argparse
    from numba import set_num_threads

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--liquid-step-cap-ps", type=float, default=10.0)
    parser.add_argument("--wall-limit-s", type=float, default=1800.0)
    parser.add_argument("--stop-at-freeze", action="store_true",
                        help="只保存首次完全凝固场，用于三维图，不积分到2 μs")
    args = parser.parse_args()
    if args.liquid_step_cap_ps <= 0.0:
        raise ValueError("液相存在期间的步长上限必须为正数")
    set_num_threads(4)

    n = 100
    base_case = HERE / "results" / "pankratz_v2_n100_dt1"
    hot_path = base_case / "heated_state.npz"
    if not hot_path.is_file():
        raise FileNotFoundError(f"缺少100³原步长加热终点场：{hot_path}")
    with np.load(hot_path, allow_pickle=False) as data:
        temperature = np.asarray(data["temperature_K"], dtype=np.float64)
        ag = np.ascontiguousarray(data["ag_fraction"], dtype=np.float64)
        dx = float(data["spacing_m"][0])
        heat_end_time = float(data["time_s"])

    numerical = CONFIG["numerical_heating_test"]
    heat_duration = float(numerical["heat_duration_s_at_10um"])
    cool_duration = float(numerical["cool_duration_s_at_10um"])
    source = np.ascontiguousarray(
        float(numerical["source_rate_W_m3_at_10um"]) * ag
        * _sphere_mask(n, dx, float(numerical["source_region_radius_fraction_of_edge"])))
    material = FastPankratzMaterials(1.0, profile=True)
    sim = NumbaEnthalpyFVM(
        ag, (dx, dx, dx), materials=material, interface_resistance_m2_K_W=0.0,
        initial_temperature_K=300.0, nonlinear_relative_tolerance=1e-8,
        linear_relative_tolerance=1e-8, maximum_nonlinear_iterations=48,
        maximum_linear_iterations=400)

    # 将保存的热端温度映射回焓场，并保留300 K初始焓作为总能量基准。
    sim.temperature_K[...] = temperature
    sim.enthalpy_density_J_m3[...] = material.enthalpy_density_J_m3(temperature, ag)
    sim.time_s = heat_end_time
    sim.source_energy_J = float(np.sum(source) * dx**3 * heat_duration)
    sim.boundary_energy_out_J = 0.0

    cap_s = args.liquid_step_cap_ps * 1e-12
    cooling_nominal_dt = cool_duration / 40.0
    config = {
        "case": "100³_Pankratz_v2_freezing_event_refinement",
        "restarted_from": str(hot_path.relative_to(HERE)),
        "resolution": n, "domain_edge_m": n * dx, "Ag_fraction_mean": float(ag.mean()),
        "Ag_mushy_width_K": 1.0, "initial_cooling_time_s": heat_end_time,
        "source_energy_from_heating_J": sim.source_energy_J,
        "cooling_duration_s": cool_duration, "cooling_nominal_dt_s": cooling_nominal_dt,
        "liquid_present_dt_cap_s": cap_s,
        "boundary_temperature_K": 300.0, "interface_resistance_m2_K_W": 0.0,
        "solver": "Numba + EW-safe + FFT", "numba_threads": 4,
        "nonlinear_rtol": 1e-8, "linear_rtol": 1e-8,
        "source_is_numerical_validation_pulse": True,
    }
    cap_tag = f"{args.liquid_step_cap_ps:g}ps"
    suffix = "_freeze_event_only" if args.stop_at_freeze else ""
    case_dir = HERE / "results" / f"pankratz_v2_n100_freezing_event_cap_{cap_tag}{suffix}"
    case_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = case_dir / "checkpoint"
    store = CheckpointStore(checkpoint_dir, sim, source, config)
    store.save("cooling", cap_s, {"initialized_from_heat_end_field": True}, [],
               {"phase": "cooling", "phase_elapsed_s": 0.0, "accepted_total": 0})
    history_path = case_dir / "accepted_steps.csv"
    retries_path = case_dir / "retries.jsonl"
    retries_path.write_text("", encoding="utf-8")
    sink = BoundarySet(dirichlet_K={f"{axis}_{side}": 300.0
                                    for axis in "xyz" for side in ("low", "high")})
    fields = (
        "time_s", "dt_s", "max_temperature_K", "Ag_bulk_liquid_fraction",
        "Ag_max_liquid_fraction", "Ag_liquid_volume_m3", "stored_energy_change_J",
        "boundary_energy_out_J", "normalized_cumulative_energy_residual",
        "nonlinear_iterations", "pcg_iterations", "step_wall_seconds", "maxrss_mib",
        "retries_before_accept", "liquid_step_cap_active",
    )
    history: list[dict] = []
    retries = 0
    next_dt = cap_s
    started = time.perf_counter()
    stop_reason = None
    cooling_end = heat_duration + cool_duration
    while sim.time_s < cooling_end - max(1e-18, cool_duration * 1e-14):
        if time.perf_counter() - started >= args.wall_limit_s:
            stop_reason = "wall_limit"
            break
        pre_liquid = material.liquid_fraction(sim.temperature_K)
        liquid_present = bool(np.any(pre_liquid > 0.0))
        dt_limit = cap_s if liquid_present else cooling_nominal_dt
        trial = min(next_dt, dt_limit, cooling_end - sim.time_s)
        failures = 0
        while True:
            try:
                record, profiles = advance_profiled(
                    sim, trial, 0.0, sink, "ew_safe",
                    fft_workers=4, preconditioner_method="fft")
                break
            except RuntimeError as exc:
                failures += 1
                retries += 1
                trial *= 0.5
                with retries_path.open("a", encoding="utf-8") as stream:
                    stream.write(json.dumps({
                        "time_s": sim.time_s, "next_dt_s": trial,
                        "failure_count": failures, "reason": str(exc),
                    }, allow_nan=False) + "\n")
                if failures >= 24 or trial < dt_limit / 4096.0:
                    raise RuntimeError(f"冷却事件细化步长重试耗尽：{exc}") from exc
        liquid = material.liquid_fraction(sim.temperature_K)
        stored = float(np.sum(sim.enthalpy_density_J_m3) * sim.volume_m3
                       - sim.initial_stored_enthalpy_J)
        scale = max(abs(sim.source_energy_J), abs(sim.boundary_energy_out_J), abs(stored), 1e-30)
        row = {
            "time_s": float(sim.time_s), "dt_s": float(trial),
            "max_temperature_K": float(sim.temperature_K.max()),
            "Ag_bulk_liquid_fraction": float(np.sum(ag * liquid) / np.sum(ag)),
            "Ag_max_liquid_fraction": float(np.max(liquid[ag > 1e-14], initial=0.0)),
            "Ag_liquid_volume_m3": float(np.sum(ag * liquid) * dx**3),
            "stored_energy_change_J": stored,
            "boundary_energy_out_J": float(sim.boundary_energy_out_J),
            "normalized_cumulative_energy_residual": abs(float(record["energy_residual_J"])) / scale,
            "nonlinear_iterations": int(record["nonlinear_iterations"]),
            "pcg_iterations": int(record["pcg_iterations_total"]),
            "step_wall_seconds": float(record["step_wall_seconds"]),
            "maxrss_mib": float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0,
            "retries_before_accept": failures,
            "liquid_step_cap_active": bool(liquid_present),
        }
        history.append(row)
        with history_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(history)
        store.append_log(record, profiles)
        after_liquid = bool(np.any(liquid > 0.0))
        if liquid_present and not after_liquid:
            np.savez_compressed(case_dir / "first_frozen_state.npz",
                                temperature_K=sim.temperature_K.astype(np.float64),
                                ag_fraction=ag, liquid_fraction=liquid,
                                spacing_m=np.array([dx] * 3), time_s=np.array(sim.time_s))
            if args.stop_at_freeze:
                stop_reason = "complete_freeze_event"
                break
        next_dt = min(cooling_nominal_dt, trial * 1.5) if after_liquid else cooling_nominal_dt
        store.save("cooling", next_dt, record, profiles,
                   {"phase": "cooling", "phase_elapsed_s": sim.time_s - heat_end_time,
                    "accepted_total": len(history)})
        if len(history) % 10 == 0 or not after_liquid:
            np.savez_compressed(case_dir / "last_cooling_state.npz",
                                temperature_K=sim.temperature_K.astype(np.float64),
                                ag_fraction=ag, liquid_fraction=liquid,
                                spacing_m=np.array([dx] * 3), time_s=np.array(sim.time_s))

    metrics = phase_metrics(sim)
    cooling_rows = history
    frozen = next((row for row in cooling_rows if row["Ag_bulk_liquid_fraction"] <= 1e-8), None)
    summary = {
        "status": "COMPLETE" if sim.time_s >= cooling_end - max(1e-18, cool_duration * 1e-14)
                  else "PARTIAL",
        "resolution": n, "domain_edge_m": n * dx,
        "first_complete_freeze_time_s": float(frozen["time_s"]) if frozen else None,
        "cooling_time_to_complete_freeze_s": (float(frozen["time_s"]) - heat_end_time)
                                               if frozen else None,
        "final_ag_bulk_liquid_fraction": metrics["Ag_bulk_liquid_fraction"],
        "final_max_temperature_K": float(sim.temperature_K.max()),
        "final_normalized_energy_residual": abs(metrics["energy_residual_J"]) / max(
            abs(metrics["source_energy_J"]), abs(metrics["boundary_energy_out_J"]),
            abs(metrics["stored_energy_change_J"]), 1e-30),
        "maximum_normalized_energy_residual": max(
            (row["normalized_cumulative_energy_residual"] for row in history), default=0.0),
        "accepted_cooling_steps": len(history), "retries": retries,
        "minimum_dt_s": min((row["dt_s"] for row in history), default=None),
        "maximum_dt_while_liquid_s": max(
            (row["dt_s"] for row in history if row["liquid_step_cap_active"]), default=None),
        "run_wall_seconds": time.perf_counter() - started,
        "solver_step_wall_seconds": sum(row["step_wall_seconds"] for row in history),
        "peak_rss_mib": float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0,
        "stop_reason": stop_reason, "run_configuration": config,
        "reinitialized_from_saved_heat_end_field": True,
        "stopped_at_first_freeze_event": args.stop_at_freeze,
    }
    (case_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


def _sphere_mask(n: int, dx: float, radius_fraction: float) -> np.ndarray:
    # 热源区域按基准RVE中心球定义，与原始50³/100³输入保持一致。
    coordinates = (np.arange(n, dtype=np.float64) + 0.5 - 0.5 * n) * dx
    radius = radius_fraction * n * dx
    return (coordinates[:, None, None]**2 + coordinates[None, :, None]**2
            + coordinates[None, None, :]**2) <= radius**2


if __name__ == "__main__":
    raise SystemExit(main())
