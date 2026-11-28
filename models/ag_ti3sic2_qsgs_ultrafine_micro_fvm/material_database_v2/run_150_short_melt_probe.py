#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用v2物性和EW-safe求解器做150³首次熔化短探针，不运行冷却全循环。"""
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
from material_database_v2.fast_materials import FastPankratzMaterials
from material_model import CONFIG
from optimized_solver import NumbaEnthalpyFVM
from profiled_advance import advance_profiled
from run_3d_phase_change import load_seed_ag, make_ag_fraction, source_field


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-steps", type=int, default=28)
    parser.add_argument("--wall-limit-s", type=float, default=2400.0)
    args = parser.parse_args()
    from numba import set_num_threads

    set_num_threads(4)
    n = 150
    numerical = CONFIG["numerical_heating_test"]
    edge = float(CONFIG["domain_edge_m"]["ultrafine_500nm"])
    seed = load_seed_ag()
    ag = np.ascontiguousarray(make_ag_fraction("bicontinuous_qsgs_seed1105", n, seed))
    dx = edge / n
    source = np.ascontiguousarray(source_field(
        ag, dx, float(numerical["source_rate_W_m3_at_10um"]),
        float(numerical["source_region_radius_fraction_of_edge"])))
    heat_duration = float(numerical["heat_duration_s_at_10um"])
    nominal_dt = heat_duration / 20.0
    material = FastPankratzMaterials(1.0, profile=True)
    warm = np.full((4, 4, 4), 300.0, dtype=np.float64)
    material.enthalpy_density_J_m3(warm, np.full_like(warm, 0.7))
    sim = NumbaEnthalpyFVM(
        ag, (dx, dx, dx), materials=material, interface_resistance_m2_K_W=0.0,
        nonlinear_relative_tolerance=1e-8, linear_relative_tolerance=1e-8,
        maximum_nonlinear_iterations=48, maximum_linear_iterations=400,
    )

    # 输出与检查点单独保存；探针仅从室温加热到首次越过熔点。
    case_dir = HERE / "results" / "pankratz_v2_n150_first_melt_probe"
    case_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = case_dir / "checkpoint"
    history_path = case_dir / "accepted_steps.csv"
    retries_path = case_dir / "retries.jsonl"
    config = {
        "case": "bicontinuous_qsgs_seed1105", "material_set": "Pankratz_v2",
        "resolution": n, "domain_edge_m": edge, "grid_spacing_m": dx,
        "Ag_fraction_mean": float(ag.mean()),
        "source_rate_W_m3": float(numerical["source_rate_W_m3_at_10um"]),
        "source_radius_fraction": float(numerical["source_region_radius_fraction_of_edge"]),
        "heating_duration_cap_s": heat_duration, "nominal_dt_s": nominal_dt,
        "initial_temperature_K": 300.0, "interface_resistance_m2_K_W": 0.0,
        "solver": "Numba + EW-safe + FFT", "numba_threads": 4,
        "nonlinear_rtol": 1e-8, "linear_rtol": 1e-8,
        "Ag_mushy_width_K": 1.0, "probe_only_until_first_melt": True,
    }
    store = CheckpointStore(checkpoint_dir, sim, source, config)
    history: list[dict] = []
    retries = 0
    next_dt = nominal_dt
    store.save("heating", next_dt, {"initial": True}, [],
               {"phase": "heating", "phase_elapsed_s": 0.0, "accepted_total": 0})
    retries_path.write_text("", encoding="utf-8")
    fields = (
        "time_s", "dt_s", "max_temperature_K", "Ag_max_liquid_fraction",
        "Ag_bulk_liquid_fraction", "Ag_liquid_volume_m3", "source_energy_J",
        "normalized_cumulative_energy_residual", "nonlinear_iterations",
        "pcg_iterations", "step_wall_seconds", "maxrss_mib", "retries_before_accept",
    )
    started = time.perf_counter()
    cpu_started = time.process_time()
    status = "STEP_CAP"
    while len(history) < args.max_steps and sim.time_s < heat_duration:
        if time.perf_counter() - started >= args.wall_limit_s:
            status = "WALL_LIMIT"
            break
        trial_dt = min(next_dt, heat_duration - sim.time_s)
        failures = 0
        while True:
            try:
                record, profiles = advance_profiled(
                    sim, trial_dt, source, None, "ew_safe",
                    fft_workers=4, preconditioner_method="fft")
                break
            except RuntimeError as exc:
                failures += 1
                retries += 1
                trial_dt *= 0.5
                with retries_path.open("a", encoding="utf-8") as stream:
                    stream.write(json.dumps({
                        "time_s": sim.time_s, "next_dt_s": trial_dt,
                        "failure_count": failures, "reason": str(exc),
                    }, allow_nan=False) + "\n")
                if failures >= 24 or trial_dt < nominal_dt / 4096.0:
                    raise RuntimeError(f"短探针步长重试耗尽：{exc}") from exc

        liquid = material.liquid_fraction(sim.temperature_K)
        ag_liquid = float(np.sum(ag * liquid) / np.sum(ag))
        liquid_volume = float(np.sum(ag * liquid) * dx**3)
        stored = float(np.sum(sim.enthalpy_density_J_m3) * sim.volume_m3
                       - sim.initial_stored_enthalpy_J)
        scale = max(abs(sim.source_energy_J), abs(sim.boundary_energy_out_J),
                    abs(stored), 1e-30)
        row = {
            "time_s": float(sim.time_s), "dt_s": float(trial_dt),
            "max_temperature_K": float(sim.temperature_K.max()),
            "Ag_max_liquid_fraction": float(np.max(liquid[ag > 1e-14])),
            "Ag_bulk_liquid_fraction": ag_liquid, "Ag_liquid_volume_m3": liquid_volume,
            "source_energy_J": float(sim.source_energy_J),
            "normalized_cumulative_energy_residual": abs(float(record["energy_residual_J"])) / scale,
            "nonlinear_iterations": int(record["nonlinear_iterations"]),
            "pcg_iterations": int(record["pcg_iterations_total"]),
            "step_wall_seconds": float(record["step_wall_seconds"]),
            "maxrss_mib": float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0,
            "retries_before_accept": failures,
        }
        history.append(row)
        with history_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(history)
        store.append_log(record, profiles)
        store.save("heating", min(nominal_dt, trial_dt * 1.5), record, profiles,
                   {"phase": "heating", "phase_elapsed_s": sim.time_s,
                    "accepted_total": len(history)})
        next_dt = min(nominal_dt, trial_dt * 1.5)
        if ag_liquid > 0.0:
            status = "FIRST_MELT_REACHED"
            break

    liquid = material.liquid_fraction(sim.temperature_K)
    metrics = phase_metrics(sim)
    np.savez_compressed(case_dir / "first_melt_state.npz",
                        temperature_K=sim.temperature_K.astype(np.float64),
                        ag_fraction=ag, liquid_fraction=liquid,
                        spacing_m=np.array([dx] * 3), time_s=np.array(sim.time_s))
    summary = {
        "status": status, "resolution": n, "domain_edge_m": edge,
        "grid_spacing_m": dx, "material_set": "Pankratz_v2",
        "time_s": float(sim.time_s), "accepted_steps": len(history),
        "retries": retries, "peak_temperature_K": float(sim.temperature_K.max()),
        "Ag_max_liquid_fraction": float(np.max(liquid[ag > 1e-14])),
        "Ag_bulk_liquid_fraction": float(np.sum(ag * liquid) / np.sum(ag)),
        "Ag_liquid_volume_m3": float(np.sum(ag * liquid) * dx**3),
        "normalized_cumulative_energy_residual": abs(metrics["energy_residual_J"]) / max(
            abs(metrics["source_energy_J"]), abs(metrics["boundary_energy_out_J"]),
            abs(metrics["stored_energy_change_J"]), 1e-30),
        "maxrss_mib": float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0,
        "wall_seconds": time.perf_counter() - started,
        "cpu_seconds": time.process_time() - cpu_started,
        "run_configuration": config,
        "source_is_numerical_validation_pulse": True,
        "full_heat_cycle_started": False,
    }
    (case_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
