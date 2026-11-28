#!/usr/bin/env python3
"""Repeatable phase-change performance run with per-step logs and atomic restart."""

from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
PHASE = PROJECT / "phase_change_validation"
sys.path.insert(0, str(PHASE))
sys.path.insert(0, str(HERE))

import numpy as np

from checkpointing import CheckpointStore
from diagnostics import mapping_diagnostics, phase_metrics, threshold_sensitivity, write_csv
from implicit_enthalpy_fvm import BoundarySet, ImplicitEnthalpyFVM
from material_model import CandidateMaterials, CONFIG
from optimized_solver import FastCandidateMaterials, NumbaEnthalpyFVM
from profiled_advance import advance_profiled
from run_3d_phase_change import load_seed_ag, make_ag_fraction, source_field


OUT = HERE / "results"
OUT.mkdir(exist_ok=True)


def maxrss_mib():
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024.0 if sys.platform != "darwin" else value / (1024.0 * 1024.0)


def atomic_json(path: Path, payload: dict):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    os.replace(temp, path)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("pure_ag", "bilayer", "bicontinuous_qsgs_seed1105"),
                        default="bicontinuous_qsgs_seed1105")
    parser.add_argument("--n", type=int, choices=(50, 100, 150), default=50)
    parser.add_argument("--solver", choices=("baseline", "numba"), default="numba")
    parser.add_argument("--linear-policy", choices=("fixed", "ew", "ew_safe"), default="fixed")
    parser.add_argument("--preconditioner", choices=("fft", "amg"), default="fft")
    parser.add_argument("--threads", type=int, choices=(1, 2, 4), default=1)
    parser.add_argument("--dt-factor", type=float, default=1.0)
    parser.add_argument("--wall-limit-s", type=float, default=900.0)
    parser.add_argument("--max-steps", type=int, default=0,
                        help="stop after this many accepted steps; zero runs the full heat/freeze cycle")
    parser.add_argument("--checkpoint-dir", type=Path, default=None)
    parser.add_argument("--resume", action="store_true")
    return parser


def warm_numba():
    from numba_kernels import apply_stencil_3d, conductance_faces_3d, material_fields_3d
    from optimized_solver import FastCandidateMaterials

    materials = FastCandidateMaterials()
    shape = (4, 4, 4)
    temperature = np.full(shape, 300.0, dtype=np.float64)
    ag = np.full(shape, 0.7, dtype=np.float64)
    h, c, k = material_fields_3d(temperature, ag, materials._kernel_parameters)
    periodic = np.zeros(3, dtype=np.bool_)
    faces = conductance_faces_3d(k, ag, 1.0e-7, 0.0, periodic)
    apply_stencil_3d(temperature, faces[0], faces[1], faces[2], periodic)


def main():
    args = build_parser().parse_args()
    if args.threads not in (1, 2, 4):
        raise ValueError("threads must be one of 1, 2, or 4")
    if args.solver == "numba":
        from numba import set_num_threads
        set_num_threads(args.threads)
        jit_started = time.perf_counter()
        warm_numba()
        jit_warmup_seconds = time.perf_counter() - jit_started
    else:
        jit_warmup_seconds = 0.0

    config = CONFIG["numerical_heating_test"]
    edge_m = float(CONFIG["domain_edge_m"]["ultrafine_500nm"])
    heat_duration = float(config["heat_duration_s_at_10um"])
    cool_duration = float(config["cool_duration_s_at_10um"])
    seed = load_seed_ag()
    ag = np.ascontiguousarray(make_ag_fraction(args.case, args.n, seed))
    dx = edge_m / args.n
    radius_fraction = float(config["source_region_radius_fraction_of_edge"])
    source = np.ascontiguousarray(source_field(ag, dx,
                                               float(config["source_rate_W_m3_at_10um"]),
                                               radius_fraction))
    if args.solver == "numba":
        materials = FastCandidateMaterials(profile=True)
        solver_type = NumbaEnthalpyFVM
    else:
        materials = CandidateMaterials()
        solver_type = ImplicitEnthalpyFVM
    solver_options = {
        "nonlinear_rtol": 1.0e-8, "linear_rtol": 1.0e-8,
        "max_newton": 48, "max_pcg": 400,
        "interface_resistance_m2_K_W": 0.0,
    }
    sim = solver_type(
        ag, (dx, dx, dx), materials=materials,
        nonlinear_relative_tolerance=solver_options["nonlinear_rtol"],
        linear_relative_tolerance=solver_options["linear_rtol"],
        maximum_nonlinear_iterations=solver_options["max_newton"],
        maximum_linear_iterations=solver_options["max_pcg"],
        interface_resistance_m2_K_W=solver_options["interface_resistance_m2_K_W"],
    )
    checkpoint_dir = args.checkpoint_dir or (OUT / "checkpoints" / (
        f"{args.case}_n{args.n}_{args.solver}_{args.linear_policy}_pc{args.preconditioner}_"
        f"dt{args.dt_factor:g}_t{args.threads}"
    ))
    run_configuration = {
        "case": args.case, "resolution": args.n, "domain_edge_m": edge_m,
        "source_rate_W_m3_on_Ag_volume": float(config["source_rate_W_m3_at_10um"]),
        "source_region_radius_fraction_of_edge": radius_fraction,
        "heat_duration_s": heat_duration, "cool_duration_s": cool_duration,
        "dt_factor": args.dt_factor, "solver": args.solver,
        "linear_policy": args.linear_policy,
        "preconditioner": args.preconditioner,
        "numba_threads": args.threads if args.solver == "numba" else 1,
        "fft_workers": args.threads,
        "solver_parameters": solver_options,
    }
    store = CheckpointStore(checkpoint_dir, sim, source, run_configuration)
    if args.resume:
        state_path = checkpoint_dir / "latest_checkpoint.npz"
        if not state_path.exists():
            raise FileNotFoundError(f"no checkpoint to resume: {state_path}")
        checkpoint_meta = CheckpointStore.restore(state_path, sim, source)
        if checkpoint_meta.get("run_configuration") != run_configuration:
            raise ValueError("checkpoint run configuration differs from current command line")
        progress = checkpoint_meta.get("progress", {})
        phase = progress.get("phase", "heating")
        phase_elapsed = float(progress.get("phase_elapsed_s", 0.0))
        next_dt = float(checkpoint_meta["next_dt_s"])
        accepted_total = int(progress.get("accepted_total", 0))
        heated_endpoint_metrics = progress.get("heated_endpoint_metrics")
        heated_endpoint_thresholds = progress.get("heated_endpoint_thresholds")
        prior_wall_seconds = float(progress.get("cumulative_wall_seconds", 0.0))
        prior_cpu_seconds = float(progress.get("cumulative_cpu_seconds", 0.0))
        prior_peak_rss_mib = float(progress.get("peak_rss_mib", 0.0))
        print(json.dumps({"resume_from_s": sim.time_s, "phase": phase,
                          "phase_elapsed_s": phase_elapsed, "next_dt_s": next_dt}), flush=True)
    else:
        phase, phase_elapsed, accepted_total = "heating", 0.0, 0
        heated_endpoint_metrics = None
        heated_endpoint_thresholds = None
        prior_wall_seconds, prior_cpu_seconds, prior_peak_rss_mib = 0.0, 0.0, 0.0
        next_dt = heat_duration / (20.0 if args.n >= 100 else 10.0) * args.dt_factor
        initial_record = {"time_s": sim.time_s, "dt_s": next_dt, "initial": True,
                          "source_energy_J": sim.source_energy_J,
                          "boundary_energy_out_J": sim.boundary_energy_out_J}
        initial_progress = {"phase": phase, "phase_elapsed_s": phase_elapsed,
                            "accepted_total": accepted_total, "max_steps": args.max_steps}
        store.save(phase, next_dt, initial_record, [], initial_progress)

    sink = BoundarySet(dirichlet_K={
        "x_low": 300.0, "x_high": 300.0, "y_low": 300.0,
        "y_high": 300.0, "z_low": 300.0, "z_high": 300.0,
    })
    phase_duration = {"heating": heat_duration, "cooling": cool_duration}
    nominal_dt = {
        "heating": heat_duration / (20.0 if args.n >= 100 else 10.0) * args.dt_factor,
        "cooling": cool_duration / (40.0 if args.n >= 100 else 20.0) * args.dt_factor,
    }
    phase_order = ("heating", "cooling")
    if phase not in phase_order:
        raise ValueError(f"invalid checkpoint phase: {phase}")
    phase_index = phase_order.index(phase)
    started_wall, started_cpu = time.perf_counter(), time.process_time()
    accepted_this_invocation = 0
    retries_total = 0
    while phase_index < len(phase_order):
        phase = phase_order[phase_index]
        duration = phase_duration[phase]
        remaining = max(0.0, duration - phase_elapsed)
        phase_source = source if phase == "heating" else 0.0
        phase_boundary = None if phase == "heating" else sink
        if phase == "cooling" and phase_elapsed == 0.0 and phase_index > 0:
            next_dt = nominal_dt["cooling"]
        failures = 0
        while remaining > max(1.0e-18, duration * 1.0e-14):
            elapsed_wall = time.perf_counter() - started_wall
            if elapsed_wall >= args.wall_limit_s:
                break
            if args.max_steps and accepted_this_invocation >= args.max_steps:
                break
            trial_dt = min(next_dt, remaining)
            step_cpu_start = time.process_time()
            try:
                record, newton_profile = advance_profiled(
                    sim, trial_dt, phase_source, phase_boundary, args.linear_policy,
                    fft_workers=args.threads, preconditioner_method=args.preconditioner,
                )
            except RuntimeError as exc:
                failures += 1
                retries_total += 1
                next_dt = trial_dt * 0.5
                retry_payload = {
                    "status": "RETRY", "phase": phase, "time_s": sim.time_s,
                    "failed_trial_dt_s": trial_dt, "next_dt_s": next_dt,
                    "reason": str(exc), "wall_seconds": time.perf_counter() - started_wall,
                }
                with (checkpoint_dir / "retry_log.jsonl").open("a", encoding="utf-8") as stream:
                    stream.write(json.dumps(retry_payload, allow_nan=False) + "\n")
                    stream.flush()
                    os.fsync(stream.fileno())
                if next_dt < nominal_dt[phase] / 256.0:
                    raise RuntimeError(f"adaptive step fell below minimum after retry: {exc}") from exc
                continue

            failures_in_step = failures
            failures = 0
            remaining -= trial_dt
            phase_elapsed += trial_dt
            accepted_total += 1
            accepted_this_invocation += 1
            next_dt = min(nominal_dt[phase], trial_dt * 1.5)
            record.update({"phase": phase, "adaptive_retries": failures_in_step,
                           "maxrss_mib": maxrss_mib(),
                           "step_cpu_seconds": time.process_time() - step_cpu_start,
                           "cpu_seconds": prior_cpu_seconds + time.process_time() - started_cpu,
                           "elapsed_wall_seconds": prior_wall_seconds + time.perf_counter() - started_wall})
            current_peak_rss = max(prior_peak_rss_mib, maxrss_mib())
            progress = {"phase": phase, "phase_elapsed_s": phase_elapsed,
                        "accepted_total": accepted_total, "max_steps": args.max_steps,
                        "heated_endpoint_metrics": heated_endpoint_metrics,
                        "heated_endpoint_thresholds": heated_endpoint_thresholds,
                        "cumulative_wall_seconds": prior_wall_seconds + time.perf_counter() - started_wall,
                        "cumulative_cpu_seconds": prior_cpu_seconds + time.process_time() - started_cpu,
                        "peak_rss_mib": current_peak_rss}
            store.append_log(record, newton_profile)
            checkpoint_path = store.save(phase, next_dt, record, newton_profile, progress)
            print(json.dumps({"case": args.case, "n": args.n, "solver": args.solver,
                              "policy": args.linear_policy, "threads": args.threads,
                              "step": {key: record[key] for key in (
                                  "phase", "time_s", "dt_s", "nonlinear_iterations",
                                  "pcg_iterations_total", "max_temperature_K",
                                  "ag_liquid_volume_fraction", "energy_residual_J",
                                  "step_wall_seconds", "adaptive_retries")},
                              "peak_rss_mib": maxrss_mib(),
                              "checkpoint": str(checkpoint_path)},
                             ensure_ascii=False, allow_nan=False), flush=True)
        if remaining > max(1.0e-18, duration * 1.0e-14):
            break
        phase_index += 1
        if phase_index < len(phase_order):
            if phase == "heating":
                heated_endpoint_metrics = phase_metrics(sim)
                heated_endpoint_thresholds = threshold_sensitivity(sim)
            phase = phase_order[phase_index]
            phase_elapsed = 0.0
            next_dt = nominal_dt[phase]
            progress = {"phase": phase, "phase_elapsed_s": phase_elapsed,
                        "accepted_total": accepted_total, "max_steps": args.max_steps,
                        "heated_endpoint_metrics": heated_endpoint_metrics,
                        "heated_endpoint_thresholds": heated_endpoint_thresholds,
                        "cumulative_wall_seconds": prior_wall_seconds + time.perf_counter() - started_wall,
                        "cumulative_cpu_seconds": prior_cpu_seconds + time.process_time() - started_cpu,
                        "peak_rss_mib": max(prior_peak_rss_mib, maxrss_mib())}
            store.save(phase, next_dt, {"phase_transition": True, "time_s": sim.time_s}, [], progress)

    complete = phase_index == len(phase_order)
    final_metrics = phase_metrics(sim)
    summary = {
        "status": "COMPLETE" if complete else "LIMIT_REACHED",
        "whole_run_timing_reliable": True,
        "case": args.case, "resolution": args.n, "solver": args.solver,
        "linear_policy": args.linear_policy, "preconditioner": args.preconditioner,
        "threads": args.threads,
        "run_configuration": run_configuration,
        "phase": phase, "phase_elapsed_s": phase_elapsed,
        "physical_time_s": sim.time_s, "accepted_steps_total": accepted_total,
        "accepted_steps_this_invocation": accepted_this_invocation,
        "retry_count_this_invocation": retries_total,
        "wall_seconds_this_invocation": time.perf_counter() - started_wall,
        "cpu_seconds_this_invocation": time.process_time() - started_cpu,
        "wall_seconds_total_estimated": prior_wall_seconds + time.perf_counter() - started_wall,
        "cpu_seconds_total_estimated": prior_cpu_seconds + time.process_time() - started_cpu,
        "average_cpu_cores_used": (prior_cpu_seconds + time.process_time() - started_cpu)
        / max(prior_wall_seconds + time.perf_counter() - started_wall, 1.0e-30),
        "peak_rss_mib": max(prior_peak_rss_mib, maxrss_mib()), "energy_residual_J": final_metrics["energy_residual_J"],
        "normalized_energy_residual": abs(final_metrics["energy_residual_J"])
        / max(abs(sim.source_energy_J), abs(sim.boundary_energy_out_J),
              abs(final_metrics["stored_energy_change_J"]), 1.0e-30),
        **final_metrics,
        "heated_endpoint_metrics": heated_endpoint_metrics,
        "jit_warmup_seconds": jit_warmup_seconds,
        "source_interpretation": "numerical localized volumetric heat pulse; not an arc boundary",
    }
    atomic_json(checkpoint_dir / "run_summary.json", summary)
    sens = threshold_sensitivity(sim)
    write_csv(checkpoint_dir / "liquid_threshold_sensitivity.csv", sens)
    if heated_endpoint_thresholds is not None:
        write_csv(checkpoint_dir / "liquid_threshold_sensitivity_heated.csv", heated_endpoint_thresholds)
    print(json.dumps(summary, ensure_ascii=False, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
