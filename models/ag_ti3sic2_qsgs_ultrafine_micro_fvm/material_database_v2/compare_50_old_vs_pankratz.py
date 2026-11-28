#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""50³ seed1105旧/新物性EW-safe同热源完整加热—冷却对照。"""
from __future__ import annotations

import argparse
import csv
import json
import os
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
from optimized_solver import FastCandidateMaterials, NumbaEnthalpyFVM
from profiled_advance import advance_profiled
from run_3d_phase_change import load_seed_ag, make_ag_fraction, source_field

OUT = HERE / "results" / "same_load_50"
OUT.mkdir(parents=True, exist_ok=True)


def maxrss_mib() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0


def load_rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def append_row(path: Path, row: dict):
    rows = load_rows(path)
    rows.append(row)
    keys = list(dict.fromkeys(k for item in rows for k in item))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_material_set(name: str, new_database: bool, resume: bool, wall_limit_s: float,
                     mushy_width_K: float = 1.0) -> dict:
    from numba import set_num_threads
    set_num_threads(4)
    numerical = CONFIG["numerical_heating_test"]
    edge = float(CONFIG["domain_edge_m"]["ultrafine_500nm"])
    n = 50
    seed = load_seed_ag()
    ag = np.ascontiguousarray(make_ag_fraction("bicontinuous_qsgs_seed1105", n, seed))
    dx = edge / n
    source = np.ascontiguousarray(source_field(
        ag, dx, float(numerical["source_rate_W_m3_at_10um"]),
        float(numerical["source_region_radius_fraction_of_edge"])))
    heat_duration = float(numerical["heat_duration_s_at_10um"])
    cool_duration = float(numerical["cool_duration_s_at_10um"])
    dt_factor = 0.25
    heat_dt = heat_duration / 10.0 * dt_factor
    cool_dt = cool_duration / 20.0 * dt_factor
    material = (FastPankratzMaterials(mushy_width_K, profile=True) if new_database
                else FastCandidateMaterials(profile=True))
    sim = NumbaEnthalpyFVM(
        ag, (dx, dx, dx), materials=material, interface_resistance_m2_K_W=0.0,
        initial_temperature_K=300.0, nonlinear_relative_tolerance=1e-8,
        linear_relative_tolerance=1e-8, maximum_nonlinear_iterations=48,
        maximum_linear_iterations=400)
    config = {
        "material_set": name, "n": n, "geometry": "QGS seed 1105, conservative remap",
        "edge_m": edge, "Ag_volume_fraction": float(ag.mean()),
        "source_rate_W_m3": float(numerical["source_rate_W_m3_at_10um"]),
        "source_radius_fraction": float(numerical["source_region_radius_fraction_of_edge"]),
        "heating_duration_s": heat_duration, "cooling_duration_s": cool_duration,
        "heat_dt_s": heat_dt, "cool_dt_s": cool_dt, "dt_factor": dt_factor,
        "Ag_mushy_interval_width_K": mushy_width_K,
        "initial_temperature_K": 300.0, "boundary_temperature_K": 300.0,
        "interface_resistance_m2_K_W": 0.0, "solver": "Numba + EW-safe + FFT",
        "numba_threads": 4, "nonlinear_rtol": 1e-8, "linear_rtol": 1e-8,
        "max_newton": 48, "max_pcg": 400,
    }
    tag = "pankratz_v2" if new_database else "old_candidate"
    if new_database and mushy_width_K != 1.0:
        tag += f"_w{str(mushy_width_K).replace('.', 'p')}K"
    checkpoint = OUT / f"checkpoint_{tag}"
    store = CheckpointStore(checkpoint, sim, source, config)
    history_path = OUT / f"{tag}_step_history.csv"
    previous_summary_path = OUT / f"{tag}_summary.json"
    previous_summary = (json.loads(previous_summary_path.read_text(encoding="utf-8"))
                        if previous_summary_path.exists() else {})
    progress_phase = "heating"
    if resume and (checkpoint / "latest_checkpoint.npz").is_file():
        meta = CheckpointStore.restore(checkpoint / "latest_checkpoint.npz", sim, source)
        progress = meta.get("progress", {})
        progress_phase = progress.get("phase", "heating")
        history = load_rows(history_path)
        history = [row for row in history if float(row["time_s"]) <= sim.time_s + 1e-18]
        if history:
            with history_path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(history[0]), lineterminator="\n")
                writer.writeheader(); writer.writerows(history)
    else:
        history = []
        progress_phase = "heating"
        (checkpoint / "step_log.jsonl").unlink(missing_ok=True)
        store.save(progress_phase, heat_dt, {"initial": True}, [], {"phase": progress_phase})

    sink = BoundarySet(dirichlet_K={f"{axis}_{side}": 300.0
                                    for axis in "xyz" for side in ("low", "high")})
    phases = (("heating", heat_duration, heat_dt, source, None),
              ("cooling", cool_duration, cool_dt, 0.0, sink))
    started = time.perf_counter()
    melt_rows = [row for row in history if row["phase"] == "heating"
                 and float(row["ag_max_liquid_fraction"]) > 0.0]
    first_melt_s = float(melt_rows[0]["time_s"]) if melt_rows else None
    max_energy_norm = max((float(row["normalized_energy_residual"]) for row in history), default=0.0)
    for phase, duration, nominal_dt, phase_source, boundary in phases:
        phase_start_s = 0.0 if phase == "heating" else heat_duration
        phase_end_s = phase_start_s + duration
        if progress_phase == "cooling" and phase == "heating":
            continue
        if sim.time_s >= phase_end_s - max(1e-18, duration * 1e-14):
            continue
        start_elapsed = max(0.0, sim.time_s - phase_start_s)
        local_elapsed = start_elapsed
        while local_elapsed < duration - max(1e-18, duration * 1e-14):
            if time.perf_counter() - started >= wall_limit_s:
                break
            trial_dt = min(nominal_dt, duration - local_elapsed)
            record, newton_profile = advance_profiled(
                sim, trial_dt, phase_source, boundary, "ew_safe",
                fft_workers=4, preconditioner_method="fft")
            local_elapsed += trial_dt
            step_energy_scale = max(abs(sim.source_energy_J), abs(sim.boundary_energy_out_J),
                                    abs(float(np.sum(sim.enthalpy_density_J_m3)) * sim.volume_m3
                                        - sim.initial_stored_enthalpy_J), 1e-30)
            energy_norm = abs(float(record["energy_residual_J"])) / step_energy_scale
            max_energy_norm = max(max_energy_norm, energy_norm)
            liquid = material.liquid_fraction(sim.temperature_K)
            if phase == "heating" and first_melt_s is None and float(liquid[ag > 1e-14].max()) > 0.0:
                first_melt_s = sim.time_s
            row = {
                "phase": phase, "time_s": sim.time_s, "phase_elapsed_s": local_elapsed,
                "dt_s": trial_dt, "max_temperature_K": float(sim.temperature_K.max()),
                "ag_mean_temperature_K": float(np.sum(sim.temperature_K * ag) / np.sum(ag)),
                "ti_mean_temperature_K": float(np.sum(sim.temperature_K * (1.0 - ag)) / np.sum(1.0 - ag)),
                "ag_max_liquid_fraction": float(np.max(liquid[ag > 1e-14])),
                "ag_bulk_liquid_fraction": float(np.sum(ag * liquid) / np.sum(ag)),
                "ag_liquid_volume_m3": float(np.sum(ag * liquid) * dx**3),
                "source_energy_J": float(sim.source_energy_J),
                "boundary_energy_out_J": float(sim.boundary_energy_out_J),
                "stored_energy_change_J": float(np.sum(sim.enthalpy_density_J_m3) * sim.volume_m3
                                                 - sim.initial_stored_enthalpy_J),
                "energy_residual_J": float(record["energy_residual_J"]),
                "normalized_energy_residual": energy_norm,
                "nonlinear_iterations": int(record["nonlinear_iterations"]),
                "pcg_iterations": int(record["pcg_iterations_total"]),
                "step_wall_seconds": float(record["step_wall_seconds"]),
                "retry_count": int(record.get("adaptive_retries", 0)),
                "maxrss_mib": maxrss_mib(),
            }
            history.append(row)
            with history_path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(history[0]), lineterminator="\n")
                writer.writeheader(); writer.writerows(history)
            store.append_log(record, newton_profile)
            if local_elapsed >= duration - max(1e-18, duration * 1e-14):
                progress_phase = "cooling" if phase == "heating" else "complete"
            else:
                progress_phase = phase
            store.save(progress_phase, nominal_dt, record, newton_profile,
                       {"phase": progress_phase, "phase_elapsed_s": local_elapsed,
                        "accepted_steps": len(history), "first_melt_time_s": first_melt_s})
            if phase == "heating" and local_elapsed >= duration - max(1e-18, duration * 1e-14):
                np.savez_compressed(OUT / f"{tag}_n50_heated_state.npz",
                                    temperature_K=sim.temperature_K.astype(np.float64),
                                    ag_fraction=ag.astype(np.float64),
                                    liquid_fraction=liquid.astype(np.float64),
                                    spacing_m=np.array([dx, dx, dx]),
                                    material_set=np.array(name), time_s=np.array(sim.time_s))
        if time.perf_counter() - started >= wall_limit_s:
            break

    liquid = material.liquid_fraction(sim.temperature_K)
    metric = phase_metrics(sim)
    sensible_and_latent = {
        "latent_energy_stored_J": float(np.sum(ag * material.rho_ag_solid_kg_m3
                                                 * material.Lf_J_kg * liquid) * dx**3),
        "total_stored_energy_change_J": float(np.sum(sim.enthalpy_density_J_m3) * sim.volume_m3
                                                - sim.initial_stored_enthalpy_J),
    }
    sensible_and_latent["sensible_energy_stored_J"] = (
        sensible_and_latent["total_stored_energy_change_J"] - sensible_and_latent["latent_energy_stored_J"])
    complete = sim.time_s >= heat_duration + cool_duration - max(1e-18, cool_duration * 1e-14)
    heating_rows = [row for row in history if row["phase"] == "heating"]
    cooling_rows = [row for row in history if row["phase"] == "cooling"]
    heated_endpoint = heating_rows[-1] if heating_rows else None
    peak_row = max(history, key=lambda row: float(row["max_temperature_K"])) if history else None
    first_frozen_row = next((row for row in cooling_rows
                             if float(row["ag_bulk_liquid_fraction"]) <= 1.0e-8), None)
    latent_at_heat_end = None
    sensible_at_heat_end = None
    if heated_endpoint is not None:
        latent_at_heat_end = (float(heated_endpoint["ag_liquid_volume_m3"])
                              * material.rho_ag_solid_kg_m3 * material.Lf_J_kg)
        sensible_at_heat_end = float(heated_endpoint["stored_energy_change_J"]) - latent_at_heat_end
    observed_wall = (float(previous_summary.get("wall_seconds", 0.0))
                     if resume and previous_summary.get("status") == "COMPLETE"
                     else time.perf_counter() - started)
    summary = {
        "status": "COMPLETE" if complete else "PARTIAL", "material_set": name,
        "Ag_mushy_interval_width_K": mushy_width_K,
        "first_melt_time_s": first_melt_s, "physical_time_s": sim.time_s,
        "heating_duration_s": heat_duration, "cooling_duration_s": cool_duration,
        "maximum_temperature_K": float(peak_row["max_temperature_K"]) if peak_row else float(sim.temperature_K.max()),
        "maximum_temperature_time_s": float(peak_row["time_s"]) if peak_row else sim.time_s,
        "Ag_mean_temperature_K": float(np.sum(sim.temperature_K * ag) / np.sum(ag)),
        "Ti3SiC2_mean_temperature_K": float(np.sum(sim.temperature_K * (1.0-ag)) / np.sum(1.0-ag)),
        "Ag_mean_temperature_at_heating_end_K": float(heated_endpoint["ag_mean_temperature_K"]) if heated_endpoint else None,
        "Ti3SiC2_mean_temperature_at_heating_end_K": float(heated_endpoint["ti_mean_temperature_K"]) if heated_endpoint else None,
        "Ag_max_liquid_fraction_at_heating_end": float(heated_endpoint["ag_max_liquid_fraction"]) if heated_endpoint else None,
        "Ag_bulk_liquid_fraction_at_heating_end": float(heated_endpoint["ag_bulk_liquid_fraction"]) if heated_endpoint else None,
        "Ag_liquid_volume_at_heating_end_m3": float(heated_endpoint["ag_liquid_volume_m3"]) if heated_endpoint else None,
        "first_complete_freeze_time_s": float(first_frozen_row["time_s"]) if first_frozen_row else None,
        "cooling_time_to_complete_freeze_s": (float(first_frozen_row["time_s"]) - heat_duration)
                                               if first_frozen_row else None,
        "latent_energy_at_heating_end_J": latent_at_heat_end,
        "sensible_energy_at_heating_end_J": sensible_at_heat_end,
        "stored_energy_at_heating_end_J": float(heated_endpoint["stored_energy_change_J"]) if heated_endpoint else None,
        "Ag_max_liquid_fraction": metric["Ag_max_liquid_fraction"],
        "Ag_bulk_liquid_fraction": metric["Ag_bulk_liquid_fraction"],
        "Ag_liquid_volume_m3": metric["Ag_liquid_weighted_volume_m3"],
        "phase_metrics": metric, **sensible_and_latent,
        "normalized_energy_residual_final": abs(metric["energy_residual_J"]) / max(
            abs(metric["source_energy_J"]), abs(metric["boundary_energy_out_J"]),
            abs(metric["stored_energy_change_J"]), 1e-30),
        "maximum_step_energy_residual_norm": max_energy_norm,
        "step_count": int(sim.step_count),
        "total_newton_iterations": int(sum(int(row["nonlinear_iterations"]) for row in history)),
        "total_pcg_iterations": int(sum(int(row["pcg_iterations"]) for row in history)),
        "wall_seconds": observed_wall, "resume_reporting_wall_seconds": time.perf_counter() - started,
        "sum_solver_step_wall_seconds": float(sum(float(row["step_wall_seconds"]) for row in history)),
        "peak_rss_mib": maxrss_mib(),
        "Ag_solid_k_at_Tm_W_m_K": float(material.ag_conductivity_W_m_K(material.Ts_K)),
        "Ag_liquid_k_at_Tm_W_m_K": float(material.ag_conductivity_W_m_K(material.Tl_K)),
        "run_configuration": config, "source_is_numerical_only": True,
        "phase_progress": progress_phase,
    }
    np.savez_compressed(OUT / f"{tag}_n50_final_state.npz",
                        temperature_K=sim.temperature_K.astype(np.float64),
                        ag_fraction=ag.astype(np.float64), liquid_fraction=liquid.astype(np.float64),
                        spacing_m=np.array([dx, dx, dx]),
                        material_set=np.array(name), time_s=np.array(sim.time_s))
    summary_path = OUT / f"{tag}_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("old", "pankratz_v2", "both"), default="both")
    parser.add_argument("--wall-limit-s", type=float, default=600.0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--mushy-width-K", type=float, default=1.0)
    args = parser.parse_args()
    which = ("old", "pankratz_v2") if args.case == "both" else (args.case,)
    summaries = []
    for case in which:
        set_name = (f"Pankratz_v2_{args.mushy_width_K:g}K" if case == "pankratz_v2" and args.mushy_width_K != 1.0
                    else case)
        summaries.append(run_material_set(set_name, case == "pankratz_v2", args.resume,
                                          args.wall_limit_s, args.mushy_width_K))
    comparison = {row["material_set"]: row for row in summaries}
    (OUT / "old_vs_pankratz_v2_comparison.json").write_text(
        json.dumps(comparison, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    csv_path = OUT / "old_vs_pankratz_v2_comparison.csv"
    fields = ("material_set", "status", "maximum_temperature_K", "Ag_mean_temperature_K",
              "Ti3SiC2_mean_temperature_K", "first_melt_time_s", "Ag_max_liquid_fraction",
              "Ag_bulk_liquid_fraction", "Ag_liquid_volume_m3", "latent_energy_stored_J",
              "sensible_energy_stored_J", "normalized_energy_residual_final", "step_count",
              "total_newton_iterations", "total_pcg_iterations", "wall_seconds", "peak_rss_mib")
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: row.get(k) for k in fields} for row in summaries)
    print(json.dumps(comparison, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
