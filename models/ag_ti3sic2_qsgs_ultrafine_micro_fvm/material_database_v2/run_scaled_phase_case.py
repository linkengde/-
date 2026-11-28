#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""运行50³ v2完整相变算例：微米尺度相似性或单一探索性热源。"""
from __future__ import annotations

import argparse
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
from run_3d_phase_change import load_seed_ag, make_ag_fraction, source_field


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True, help="结果目录标签，只允许字母数字、横线和下划线")
    parser.add_argument("--edge-scale", type=float, default=1.0,
                        help="相对10 μm基准边长；微米模型使用20")
    parser.add_argument("--source-scale", type=float, default=1.0,
                        help="独立于相似缩放的体热源功率倍率")
    parser.add_argument("--wall-limit-s", type=float, default=1800.0)
    args = parser.parse_args()
    if not args.label.replace("-", "").replace("_", "").isalnum():
        raise ValueError("label只能由字母、数字、横线和下划线组成")
    if args.edge_scale <= 0.0 or args.source_scale <= 0.0:
        raise ValueError("尺度和热源倍率必须为正数")

    from numba import set_num_threads
    set_num_threads(4)
    n = 50
    numerical = CONFIG["numerical_heating_test"]
    base_edge = float(CONFIG["domain_edge_m"]["ultrafine_500nm"])
    edge = base_edge * args.edge_scale
    scale2 = args.edge_scale**2
    seed = load_seed_ag()
    # 同一seed与同一网格保持无量纲QGS几何一致，只改变真实空间尺度。
    ag = np.ascontiguousarray(make_ag_fraction("bicontinuous_qsgs_seed1105", n, seed))
    dx = edge / n
    rate = (float(numerical["source_rate_W_m3_at_10um"]) / scale2
            * args.source_scale)
    source = np.ascontiguousarray(source_field(
        ag, dx, rate, float(numerical["source_region_radius_fraction_of_edge"])))
    heat_duration = float(numerical["heat_duration_s_at_10um"]) * scale2
    cool_duration = float(numerical["cool_duration_s_at_10um"]) * scale2
    # 保持与10 μm基准相同的无量纲时间步；dt随L²缩放。
    heat_dt = heat_duration / 10.0 * 0.25
    cool_dt = cool_duration / 20.0 * 0.25
    material = FastPankratzMaterials(1.0, profile=True)
    sim = NumbaEnthalpyFVM(
        ag, (dx, dx, dx), materials=material, interface_resistance_m2_K_W=0.0,
        initial_temperature_K=300.0, nonlinear_relative_tolerance=1e-8,
        linear_relative_tolerance=1e-8, maximum_nonlinear_iterations=48,
        maximum_linear_iterations=400)
    config = {
        "label": args.label, "material_set": "Pankratz_v2", "resolution": n,
        "base_edge_m": base_edge, "domain_edge_m": edge, "edge_scale": args.edge_scale,
        "grid_spacing_m": dx, "Ag_fraction_mean": float(ag.mean()),
        "source_rate_W_m3": rate, "source_scale_independent_of_similarity": args.source_scale,
        "source_radius_fraction": float(numerical["source_region_radius_fraction_of_edge"]),
        "heating_duration_s": heat_duration, "cooling_duration_s": cool_duration,
        "heat_dt_s": heat_dt, "cool_dt_s": cool_dt, "dt_scale_relation": "L^2",
        "initial_temperature_K": 300.0, "boundary_temperature_K": 300.0,
        "interface_resistance_m2_K_W": 0.0, "Ag_mushy_width_K": 1.0,
        "solver": "Numba + EW-safe + FFT", "numba_threads": 4,
        "nonlinear_rtol": 1e-8, "linear_rtol": 1e-8,
        "source_is_numerical_validation_pulse": True,
    }
    case_dir = HERE / "results" / f"scaled_n50_{args.label}"
    case_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = case_dir / "checkpoint"
    history_path = case_dir / "accepted_steps.csv"
    retry_path = case_dir / "retries.jsonl"
    store = CheckpointStore(checkpoint_dir, sim, source, config)
    store.save("heating", heat_dt, {"initial": True}, [],
               {"phase": "heating", "phase_elapsed_s": 0.0, "accepted_total": 0})
    retry_path.write_text("", encoding="utf-8")
    sink = BoundarySet(dirichlet_K={f"{axis}_{side}": 300.0
                                    for axis in "xyz" for side in ("low", "high")})
    phases = (("heating", heat_duration, heat_dt, source, None, 0.0),
              ("cooling", cool_duration, cool_dt, 0.0, sink, heat_duration))
    fields = (
        "phase", "time_s", "dt_s", "max_temperature_K", "Ag_mean_temperature_K",
        "Ti3SiC2_mean_temperature_K", "Ag_max_liquid_fraction", "Ag_bulk_liquid_fraction",
        "Ag_liquid_volume_m3", "source_energy_J", "boundary_energy_out_J",
        "stored_energy_change_J", "energy_residual_J", "normalized_energy_residual",
        "nonlinear_iterations", "pcg_iterations", "step_wall_seconds", "maxrss_mib",
    )
    history: list[dict] = []
    retries = 0
    started = time.perf_counter()
    stop_reason = None
    for phase, duration, nominal_dt, phase_source, boundary, phase_origin in phases:
        elapsed = 0.0
        phase_rows: list[dict] = []
        while elapsed < duration - max(1e-18, duration * 1e-14):
            if time.perf_counter() - started >= args.wall_limit_s:
                stop_reason = "wall_limit"
                break
            trial = min(nominal_dt, duration - elapsed)
            failures = 0
            while True:
                try:
                    record, profiles = advance_profiled(
                        sim, trial, phase_source, boundary, "ew_safe",
                        fft_workers=4, preconditioner_method="fft")
                    break
                except RuntimeError as exc:
                    failures += 1
                    retries += 1
                    trial *= 0.5
                    with retry_path.open("a", encoding="utf-8") as stream:
                        stream.write(json.dumps({
                            "phase": phase, "time_s": sim.time_s, "next_dt_s": trial,
                            "failure_count": failures, "reason": str(exc),
                        }, allow_nan=False) + "\n")
                    if failures >= 24 or trial < nominal_dt / 4096.0:
                        raise RuntimeError(f"EW-safe时间步重试耗尽：{exc}") from exc
            elapsed += trial
            liquid = material.liquid_fraction(sim.temperature_K)
            stored = float(np.sum(sim.enthalpy_density_J_m3) * sim.volume_m3
                           - sim.initial_stored_enthalpy_J)
            scale = max(abs(sim.source_energy_J), abs(sim.boundary_energy_out_J), abs(stored), 1e-30)
            row = {
                "phase": phase, "time_s": float(sim.time_s), "dt_s": float(trial),
                "max_temperature_K": float(sim.temperature_K.max()),
                "Ag_mean_temperature_K": float(np.sum(sim.temperature_K * ag) / np.sum(ag)),
                "Ti3SiC2_mean_temperature_K": float(np.sum(sim.temperature_K * (1.0-ag)) / np.sum(1.0-ag)),
                "Ag_max_liquid_fraction": float(np.max(liquid[ag > 1e-14], initial=0.0)),
                "Ag_bulk_liquid_fraction": float(np.sum(ag * liquid) / np.sum(ag)),
                "Ag_liquid_volume_m3": float(np.sum(ag * liquid) * dx**3),
                "source_energy_J": float(sim.source_energy_J),
                "boundary_energy_out_J": float(sim.boundary_energy_out_J),
                "stored_energy_change_J": stored,
                "energy_residual_J": float(record["energy_residual_J"]),
                "normalized_energy_residual": abs(float(record["energy_residual_J"])) / scale,
                "nonlinear_iterations": int(record["nonlinear_iterations"]),
                "pcg_iterations": int(record["pcg_iterations_total"]),
                "step_wall_seconds": float(record["step_wall_seconds"]),
                "maxrss_mib": float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0,
            }
            history.append(row)
            phase_rows.append(row)
            with history_path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
                writer.writeheader()
                writer.writerows(history)
            store.append_log(record, profiles)
            done = elapsed >= duration - max(1e-18, duration * 1e-14)
            next_phase = "cooling" if phase == "heating" and done else (
                "complete" if phase == "cooling" and done else phase)
            store.save(next_phase, min(nominal_dt, trial * 1.5), record, profiles,
                       {"phase": next_phase,
                        "phase_elapsed_s": 0.0 if next_phase == "cooling" else elapsed,
                        "accepted_total": len(history)})
            if phase == "heating" and done:
                np.savez_compressed(case_dir / "heated_state.npz",
                                    temperature_K=sim.temperature_K.astype(np.float64),
                                    ag_fraction=ag, liquid_fraction=liquid,
                                    spacing_m=np.array([dx] * 3), time_s=np.array(sim.time_s))
        if stop_reason:
            break

    liquid = material.liquid_fraction(sim.temperature_K)
    metrics = phase_metrics(sim)
    heating = [row for row in history if row["phase"] == "heating"]
    cooling = [row for row in history if row["phase"] == "cooling"]
    first_melt = next((row for row in heating if row["Ag_max_liquid_fraction"] > 0.0), None)
    frozen = next((row for row in cooling if row["Ag_bulk_liquid_fraction"] <= 1e-8), None)
    peak = max(history, key=lambda row: row["max_temperature_K"]) if history else None
    complete = sim.time_s >= heat_duration + cool_duration - max(1e-18, cool_duration * 1e-14)
    summary = {
        "status": "COMPLETE" if complete else "PARTIAL",
        "label": args.label, "material_set": "Pankratz_v2", "resolution": n,
        "domain_edge_m": edge, "edge_scale": args.edge_scale,
        "first_melt_time_s": float(first_melt["time_s"]) if first_melt else None,
        "first_complete_freeze_time_s": float(frozen["time_s"]) if frozen else None,
        "cooling_time_to_complete_freeze_s": float(frozen["time_s"]) - heat_duration if frozen else None,
        "peak_temperature_K": float(peak["max_temperature_K"]) if peak else float(sim.temperature_K.max()),
        "heat_end_max_temperature_K": float(heating[-1]["max_temperature_K"]) if heating else None,
        "heat_end_ag_max_liquid_fraction": float(heating[-1]["Ag_max_liquid_fraction"]) if heating else None,
        "heat_end_ag_bulk_liquid_fraction": float(heating[-1]["Ag_bulk_liquid_fraction"]) if heating else None,
        "heat_end_ag_liquid_volume_m3": float(heating[-1]["Ag_liquid_volume_m3"]) if heating else None,
        "final_ag_bulk_liquid_fraction": metrics["Ag_bulk_liquid_fraction"],
        "final_normalized_energy_residual": abs(metrics["energy_residual_J"]) / max(
            abs(metrics["source_energy_J"]), abs(metrics["boundary_energy_out_J"]),
            abs(metrics["stored_energy_change_J"]), 1e-30),
        "max_normalized_energy_residual": max((row["normalized_energy_residual"] for row in history), default=0.0),
        "accepted_steps": len(history), "retries": retries,
        "newton_iterations": sum(row["nonlinear_iterations"] for row in history),
        "pcg_iterations": sum(row["pcg_iterations"] for row in history),
        "solver_step_wall_seconds": sum(row["step_wall_seconds"] for row in history),
        "run_wall_seconds": time.perf_counter() - started,
        "peak_rss_mib": float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0,
        "stop_reason": stop_reason, "run_configuration": config,
        "phase_metrics": metrics,
    }
    np.savez_compressed(case_dir / "final_state.npz",
                        temperature_K=sim.temperature_K.astype(np.float64),
                        ag_fraction=ag, liquid_fraction=liquid,
                        spacing_m=np.array([dx] * 3), time_s=np.array(sim.time_s))
    (case_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
