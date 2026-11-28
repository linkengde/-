#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""100³ Pankratz v2 Numba/EW-safe完整热循环与时间步收敛。"""
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

OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)


def maxrss_mib() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0


def read_history(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_history(path: Path, rows: list[dict]):
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_case(dt_factor: float, wall_limit_s: float, max_steps: int, resume: bool) -> dict:
    from numba import set_num_threads
    set_num_threads(4)
    n = 100
    numerical = CONFIG["numerical_heating_test"]
    edge = float(CONFIG["domain_edge_m"]["ultrafine_500nm"])
    seed = load_seed_ag()
    ag = np.ascontiguousarray(make_ag_fraction("bicontinuous_qsgs_seed1105", n, seed))
    dx = edge / n
    source = np.ascontiguousarray(source_field(
        ag, dx, float(numerical["source_rate_W_m3_at_10um"]),
        float(numerical["source_region_radius_fraction_of_edge"])))
    heat_duration = float(numerical["heat_duration_s_at_10um"])
    cool_duration = float(numerical["cool_duration_s_at_10um"])
    heat_dt = heat_duration / 20.0 * dt_factor
    cool_dt = cool_duration / 40.0 * dt_factor
    material = FastPankratzMaterials(1.0, profile=True)
    # 先预编译物性内核；冷启动到首次有效步进的耗时单独记录。
    warm_temperature = np.full((4, 4, 4), 300.0, dtype=np.float64)
    material.enthalpy_density_J_m3(warm_temperature, np.full_like(warm_temperature, 0.7))
    options = dict(nonlinear_relative_tolerance=1e-8, linear_relative_tolerance=1e-8,
                   maximum_nonlinear_iterations=48, maximum_linear_iterations=400)
    sim = NumbaEnthalpyFVM(ag, (dx, dx, dx), materials=material,
                           interface_resistance_m2_K_W=0.0, **options)
    config = {
        "case": "bicontinuous_qsgs_seed1105", "material_set": "Pankratz_v2",
        "resolution": n, "domain_edge_m": edge, "Ag_fraction_mean": float(ag.mean()),
        "source_rate_W_m3": float(numerical["source_rate_W_m3_at_10um"]),
        "source_radius_fraction": float(numerical["source_region_radius_fraction_of_edge"]),
        "heating_duration_s": heat_duration, "cooling_duration_s": cool_duration,
        "heat_nominal_dt_s": heat_dt, "cool_nominal_dt_s": cool_dt,
        "dt_factor": dt_factor, "initial_temperature_K": 300.0,
        "boundary_temperature_K": 300.0, "interface_resistance_m2_K_W": 0.0,
        "solver": "Numba + EW-safe + FFT", "numba_threads": 4,
        "nonlinear_rtol": 1e-8, "linear_rtol": 1e-8,
        "max_newton": 48, "max_pcg": 400, "Ag_mushy_width_K": 1.0,
    }
    tag = f"pankratz_v2_n100_dt{dt_factor:g}"
    case_dir = OUT / tag
    case_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = case_dir / "checkpoint"
    history_path = case_dir / "accepted_steps.csv"
    retry_path = case_dir / "retries.jsonl"
    store = CheckpointStore(checkpoint_dir, sim, source, config)
    progress_phase = "heating"
    phase_elapsed = 0.0
    next_dt = heat_dt
    if resume and (checkpoint_dir / "latest_checkpoint.npz").exists():
        meta = CheckpointStore.restore(checkpoint_dir / "latest_checkpoint.npz", sim, source)
        progress = meta.get("progress", {})
        progress_phase = progress.get("phase", "heating")
        phase_elapsed = float(progress.get("phase_elapsed_s", 0.0))
        next_dt = float(meta["next_dt_s"])
        history = read_history(history_path)
        history = [row for row in history if float(row["time_s"]) <= sim.time_s + 1e-18]
        write_history(history_path, history)
    else:
        history = []
        (checkpoint_dir / "step_log.jsonl").unlink(missing_ok=True)
        retry_path.write_text("", encoding="utf-8")
        store.save("heating", heat_dt, {"initial": True}, [],
                   {"phase": "heating", "phase_elapsed_s": 0.0, "accepted_total": 0})

    sink = BoundarySet(dirichlet_K={f"{axis}_{side}": 300.0
                                    for axis in "xyz" for side in ("low", "high")})
    phases = (("heating", heat_duration, heat_dt, source, None, 0.0),
              ("cooling", cool_duration, cool_dt, 0.0, sink, heat_duration))
    started = time.perf_counter()
    started_cpu = time.process_time()
    warm_wall = None
    warm_cpu = None
    retries_total = sum(1 for _ in retry_path.open(encoding="utf-8")) if retry_path.exists() else 0
    stop_reason = None
    for phase_name, duration, nominal_dt, phase_source, boundary, phase_origin in phases:
        phase_end = phase_origin + duration
        if sim.time_s >= phase_end - max(1e-18, duration * 1e-14):
            continue
        if progress_phase == "cooling" and phase_name == "heating":
            continue
        local_elapsed = max(0.0, sim.time_s - phase_origin)
        if progress_phase == phase_name:
            local_elapsed = phase_elapsed
        elif phase_name == "cooling":
            local_elapsed = 0.0
            next_dt = nominal_dt
        phase_records = []
        while local_elapsed < duration - max(1e-18, duration * 1e-14):
            if time.perf_counter() - started >= wall_limit_s:
                stop_reason = "wall_limit"
                break
            if max_steps and len(history) >= max_steps:
                stop_reason = "max_steps"
                break
            trial_dt = min(next_dt, duration - local_elapsed)
            failure_count = 0
            while True:
                try:
                    record, profiles = advance_profiled(
                        sim, trial_dt, phase_source, boundary, "ew_safe",
                        fft_workers=4, preconditioner_method="fft")
                    break
                except RuntimeError as exc:
                    failure_count += 1
                    retries_total += 1
                    trial_dt *= 0.5
                    retry_row = {"phase": phase_name, "time_s": sim.time_s,
                                 "rejected_dt_s": 2.0 * trial_dt, "next_dt_s": trial_dt,
                                 "failure_count": failure_count, "reason": str(exc)}
                    with retry_path.open("a", encoding="utf-8") as stream:
                        stream.write(json.dumps(retry_row, allow_nan=False) + "\n")
                    if failure_count >= 24 or trial_dt < nominal_dt / 4096.0:
                        raise RuntimeError(f"EW-safe dt retries exhausted: {exc}") from exc
                if failure_count >= 24:
                    raise RuntimeError("EW-safe dt retries exhausted")
            if warm_wall is None:
                warm_wall = time.perf_counter() - started
                warm_cpu = time.process_time() - started_cpu
            local_elapsed += trial_dt
            liquid = material.liquid_fraction(sim.temperature_K)
            total_storage = float(np.sum(sim.enthalpy_density_J_m3) * sim.volume_m3
                                  - sim.initial_stored_enthalpy_J)
            energy_scale = max(abs(sim.source_energy_J), abs(sim.boundary_energy_out_J),
                               abs(total_storage), 1e-30)
            row = {
                "phase": phase_name, "time_s": sim.time_s, "phase_elapsed_s": local_elapsed,
                "dt_s": trial_dt, "retries_before_accept": failure_count,
                "max_temperature_K": float(sim.temperature_K.max()),
                "ag_mean_temperature_K": float(np.sum(sim.temperature_K * ag) / np.sum(ag)),
                "ti_mean_temperature_K": float(np.sum(sim.temperature_K * (1.0-ag)) / np.sum(1.0-ag)),
                "ag_max_liquid_fraction": float(np.max(liquid[ag > 1e-14])),
                "ag_bulk_liquid_fraction": float(np.sum(ag * liquid) / np.sum(ag)),
                "ag_liquid_volume_m3": float(np.sum(ag * liquid) * dx**3),
                "source_energy_J": float(sim.source_energy_J),
                "boundary_energy_out_J": float(sim.boundary_energy_out_J),
                "stored_energy_change_J": total_storage,
                "energy_residual_J": float(record["energy_residual_J"]),
                "normalized_cumulative_energy_residual": abs(float(record["energy_residual_J"])) / energy_scale,
                "nonlinear_iterations": int(record["nonlinear_iterations"]),
                "pcg_iterations": int(record["pcg_iterations_total"]),
                "step_wall_seconds": float(record["step_wall_seconds"]),
                "maxrss_mib": maxrss_mib(),
            }
            history.append(row)
            phase_records.append(row)
            write_history(history_path, history)
            store.append_log(record, profiles)
            done = local_elapsed >= duration - max(1e-18, duration * 1e-14)
            progress_phase = ("cooling" if phase_name == "heating" and done else
                              "complete" if phase_name == "cooling" and done else phase_name)
            store.save(progress_phase, min(nominal_dt, trial_dt * 1.5), record, profiles,
                       {"phase": progress_phase,
                        "phase_elapsed_s": 0.0 if progress_phase == "cooling" else local_elapsed,
                        "accepted_total": len(history),
                        "cumulative_wall_s": time.perf_counter() - started})
            next_dt = min(nominal_dt, trial_dt * 1.5)
            if phase_name == "heating" and done:
                np.savez_compressed(case_dir / "heated_state.npz",
                                    temperature_K=sim.temperature_K.astype(np.float64),
                                    ag_fraction=ag, liquid_fraction=liquid,
                                    spacing_m=np.array([dx] * 3), time_s=np.array(sim.time_s))
        if stop_reason:
            break

    liquid = material.liquid_fraction(sim.temperature_K)
    metric = phase_metrics(sim)
    heat_rows = [row for row in history if row["phase"] == "heating"]
    cool_rows = [row for row in history if row["phase"] == "cooling"]
    heat_end = heat_rows[-1] if heat_rows else None
    peak = max(history, key=lambda row: float(row["max_temperature_K"])) if history else None
    melt = next((row for row in heat_rows if float(row["ag_max_liquid_fraction"]) > 0.0), None)
    frozen = next((row for row in cool_rows if float(row["ag_bulk_liquid_fraction"]) <= 1e-8), None)
    latent_heat_end = (float(heat_end["ag_liquid_volume_m3"]) * material.rho_ag_solid_kg_m3
                       * material.Lf_J_kg) if heat_end else None
    stored_heat_end = float(heat_end["stored_energy_change_J"]) if heat_end else None
    complete = sim.time_s >= heat_duration + cool_duration - max(1e-18, cool_duration * 1e-14)
    summary = {
        "status": "COMPLETE" if complete else "PARTIAL", "dt_factor": dt_factor,
        "material_set": "Pankratz_v2", "resolution": n, "domain_edge_m": edge,
        "grid_spacing_m": dx, "Ag_fraction_mean": float(ag.mean()),
        "first_melt_time_s": float(melt["time_s"]) if melt else None,
        "first_complete_freeze_time_s": float(frozen["time_s"]) if frozen else None,
        "cooling_time_to_complete_freeze_s": float(frozen["time_s"]) - heat_duration if frozen else None,
        "peak_temperature_K": float(peak["max_temperature_K"]) if peak else float(sim.temperature_K.max()),
        "peak_temperature_time_s": float(peak["time_s"]) if peak else sim.time_s,
        "heat_end_max_temperature_K": float(heat_end["max_temperature_K"]) if heat_end else None,
        "heat_end_ag_mean_temperature_K": float(heat_end["ag_mean_temperature_K"]) if heat_end else None,
        "heat_end_ti_mean_temperature_K": float(heat_end["ti_mean_temperature_K"]) if heat_end else None,
        "heat_end_ag_max_liquid_fraction": float(heat_end["ag_max_liquid_fraction"]) if heat_end else None,
        "heat_end_ag_bulk_liquid_fraction": float(heat_end["ag_bulk_liquid_fraction"]) if heat_end else None,
        "heat_end_ag_liquid_volume_m3": float(heat_end["ag_liquid_volume_m3"]) if heat_end else None,
        "heat_end_latent_energy_J": latent_heat_end,
        "heat_end_sensible_energy_J": stored_heat_end - latent_heat_end if stored_heat_end is not None else None,
        "final_ag_bulk_liquid_fraction": metric["Ag_bulk_liquid_fraction"],
        "final_energy_residual_J": metric["energy_residual_J"],
        "final_normalized_energy_residual": abs(metric["energy_residual_J"]) / max(
            abs(metric["source_energy_J"]), abs(metric["boundary_energy_out_J"]),
            abs(metric["stored_energy_change_J"]), 1e-30),
        "max_normalized_cumulative_energy_residual": max(
            (float(row["normalized_cumulative_energy_residual"]) for row in history), default=0.0),
        "accepted_steps": len(history), "retries": retries_total,
        "newton_iterations": sum(int(row["nonlinear_iterations"]) for row in history),
        "pcg_iterations": sum(int(row["pcg_iterations"]) for row in history),
        "minimum_accepted_dt_s": min((float(row["dt_s"]) for row in history), default=None),
        "maximum_accepted_dt_s": max((float(row["dt_s"]) for row in history), default=None),
        "solver_step_wall_seconds": sum(float(row["step_wall_seconds"]) for row in history),
        "run_wall_seconds": time.perf_counter() - started,
        "cold_start_through_first_accepted_step_wall_s": warm_wall,
        "cold_start_through_first_accepted_step_cpu_s": warm_cpu,
        "peak_rss_mib": maxrss_mib(), "stop_reason": stop_reason,
        "run_configuration": config,
        "source_is_numerical_validation_pulse": True,
    }
    np.savez_compressed(case_dir / "final_state.npz",
                        temperature_K=sim.temperature_K.astype(np.float64),
                        ag_fraction=ag, liquid_fraction=liquid,
                        spacing_m=np.array([dx] * 3), time_s=np.array(sim.time_s))
    (case_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--factors", nargs="+", type=float, default=[1.0, 0.5])
    parser.add_argument("--wall-limit-s", type=float, default=3600.0)
    parser.add_argument("--max-steps", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--summarize-only", action="store_true",
                        help="只比较已保存完整算例，不启动新的FVM计算")
    args = parser.parse_args()
    # 支持分批完成步长组，并汇总已经完整保存在结果目录中的算例。
    by_factor = {}
    for case_dir in OUT.glob("pankratz_v2_n100_dt*"):
        summary_path = case_dir / "summary.json"
        if summary_path.exists():
            row = json.loads(summary_path.read_text(encoding="utf-8"))
            by_factor[float(row["dt_factor"])] = row
    if not args.summarize_only:
        for factor in args.factors:
            by_factor[float(factor)] = run_case(factor, args.wall_limit_s, args.max_steps, args.resume)
    summaries = [by_factor[factor] for factor in sorted(by_factor, reverse=True)]
    comparison = {f"dt_factor_{row['dt_factor']:g}": row for row in summaries}
    metrics = ("first_melt_time_s", "peak_temperature_K",
               "heat_end_ag_max_liquid_fraction", "heat_end_ag_bulk_liquid_fraction",
               "heat_end_ag_liquid_volume_m3", "cooling_time_to_complete_freeze_s")
    adjacent_changes = {}
    for coarse, fine in zip(summaries, summaries[1:]):
        adjacent_changes[f"dt_factor_{fine['dt_factor']:g}_vs_{coarse['dt_factor']:g}"] = {
            metric: abs(float(fine[metric]) - float(coarse[metric]))
                    / max(abs(float(coarse[metric])), 1e-30)
            if fine.get(metric) is not None and coarse.get(metric) is not None else None
            for metric in metrics
        }
    reference = summaries[0] if summaries else None
    relative_to_original = {}
    if reference is not None:
        for case in summaries[1:]:
            relative_to_original[f"dt_factor_{case['dt_factor']:g}_vs_{reference['dt_factor']:g}"] = {
                metric: abs(float(case[metric]) - float(reference[metric]))
                        / max(abs(float(reference[metric])), 1e-30)
                if case.get(metric) is not None and reference.get(metric) is not None else None
                for metric in metrics
            }
    payload = {
        "cases": comparison, "adjacent_relative_changes": adjacent_changes,
        "relative_to_original_dt": relative_to_original,
        "time_step_convergence_3percent": len(summaries) >= 2 and all(
            value is not None and value <= 0.03
            for group in adjacent_changes.values() for value in group.values()),
        "interpretation": "Pankratz v2; identical geometry, pulse, initial/boundary conditions and EW-safe solver; not an arc model.",
    }
    (OUT / "time_refinement_n100_comparison.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    fields = ("dt_factor", "status", "first_melt_time_s", "peak_temperature_K",
              "heat_end_max_temperature_K", "heat_end_ag_mean_temperature_K",
              "heat_end_ti_mean_temperature_K", "heat_end_ag_max_liquid_fraction",
              "heat_end_ag_bulk_liquid_fraction", "heat_end_ag_liquid_volume_m3",
              "heat_end_latent_energy_J", "heat_end_sensible_energy_J",
              "cooling_time_to_complete_freeze_s", "final_ag_bulk_liquid_fraction",
              "final_normalized_energy_residual", "max_normalized_cumulative_energy_residual",
              "accepted_steps", "retries",
              "newton_iterations", "pcg_iterations", "solver_step_wall_seconds",
              "run_wall_seconds", "peak_rss_mib")
    with (OUT / "time_refinement_n100_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: row.get(key) for key in fields} for row in summaries)
    print(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
