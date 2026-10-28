#!/usr/bin/env python3
"""Short numerical acceptance test and restartable Case A run for Ag 60.7 wt%."""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import signal
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT_DIR = HERE.parents[2]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = HERE / "nano_geometry_ag607wt_ti3sic2_393wt.npz"
COMPOSITION = HERE / "composition_geometry.json"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid, write_csv, write_cluster_csv  # noqa: E402
from qsgs_liquid_tracking import LiquidClusterTracker  # noqa: E402


STOP_REQUESTED = False


def _handle_stop(signum: int, frame: object) -> None:
    global STOP_REQUESTED
    STOP_REQUESTED = True


signal.signal(signal.SIGTERM, _handle_stop)
signal.signal(signal.SIGINT, _handle_stop)


def atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(value, indent=2, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def append_jsonl(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(value, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def tag_ps(value: float) -> str:
    return f"{value:.9f}".rstrip("0").rstrip(".").replace(".", "p") or "0"


def model_for(refinement: int) -> QSGSEnthalpyGrid:
    return QSGSEnthalpyGrid(
        None,
        "A",
        native_geometry_npz=GEOMETRY,
        source_geometry_npz=GEOMETRY,
        mesh_refinement=refinement,
    )


def tracker_for(model: QSGSEnthalpyGrid) -> LiquidClusterTracker:
    return LiquidClusterTracker(model.spacing_A, model.size_A, liquid_threshold=0.5)


def _relative_energy_error(model: QSGSEnthalpyGrid) -> tuple[float, float]:
    total = float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3)
    error = total - model.injected_energy_J
    relative = error / model.injected_energy_J if model.injected_energy_J else 0.0
    return error, relative


def assess_stable_pool(
    tracker: LiquidClusterTracker,
    end_ps: float,
    minimum_track_ps: float = 3.0,
    window_ps: float = 2.0,
    maximum_range_fraction: float = 0.05,
) -> dict[str, Any]:
    """Require a persistent main liquid track and a 5% volume plateau."""
    main_records = [record for record in tracker.history if record.get("is_main_cluster")]
    current = [record for record in main_records if abs(float(record["time_ps"]) - end_ps) <= 1.0e-8]
    if not current:
        return {"stable": False, "reason": "no main liquid cluster at current time"}
    record = max(current, key=lambda item: int(item["voxel_count"]))
    track_id = int(record["track_id"])
    duration = float(record["continuous_duration_ps"])
    if duration < minimum_track_ps - 1.0e-8:
        return {
            "stable": False,
            "reason": "main liquid track has not persisted for the required duration",
            "track_id": track_id,
            "continuous_duration_ps": duration,
        }
    start_ps = end_ps - window_ps
    records = sorted(
        [
            item for item in main_records
            if int(item["track_id"]) == track_id
            and start_ps - 1.0e-9 <= float(item["time_ps"]) <= end_ps + 1.0e-9
        ],
        key=lambda item: float(item["time_ps"]),
    )
    if len(records) < 2 or float(records[0]["time_ps"]) > start_ps + 0.02:
        return {
            "stable": False,
            "reason": "main track does not cover the full stability window",
            "track_id": track_id,
            "samples": len(records),
        }
    times = [float(item["time_ps"]) for item in records]
    gaps = [right - left for left, right in zip(times, times[1:])]
    max_gap = max(gaps, default=0.0)
    if max_gap > 0.025:
        return {
            "stable": False,
            "reason": "sampling gap in the stability window",
            "track_id": track_id,
            "max_gap_ps": max_gap,
        }
    volumes = [float(item["volume_A3"]) for item in records]
    mean_volume = float(np.mean(volumes))
    relative_range = (max(volumes) - min(volumes)) / mean_volume if mean_volume else float("inf")
    return {
        "stable": relative_range <= maximum_range_fraction,
        "reason": "stability thresholds met" if relative_range <= maximum_range_fraction else "pool volume range exceeds threshold",
        "track_id": track_id,
        "continuous_duration_ps": duration,
        "window_start_ps": start_ps,
        "window_end_ps": end_ps,
        "samples": len(records),
        "volume_mean_A3": mean_volume,
        "volume_min_A3": min(volumes),
        "volume_max_A3": max(volumes),
        "relative_range": relative_range,
        "maximum_range_fraction": maximum_range_fraction,
    }


def validate_checkpoint(
    model: QSGSEnthalpyGrid,
    tracker: LiquidClusterTracker,
    checkpoint: Path,
    *,
    phase: str,
    stability: dict[str, Any] | None = None,
) -> dict[str, Any]:
    restored = model.load_checkpoint(checkpoint, tracker)
    energy_error, relative_error = _relative_energy_error(model)
    peak = float(model.state.temperature_K.max())
    if not np.isfinite(peak) or not np.all(np.isfinite(model.enthalpy_density_J_m3)):
        raise FloatingPointError(f"Checkpoint restore produced a non-finite state: {checkpoint}")
    if np.any(model.state.ag_liquid_fraction < 0.0) or np.any(model.state.ag_liquid_fraction > 1.0):
        raise ArithmeticError(f"Checkpoint restore has Ag liquid fraction outside [0, 1]: {checkpoint}")
    if abs(relative_error) > 1.0e-8:
        raise ArithmeticError(f"Checkpoint energy error exceeds 1e-8: {relative_error:.3e}")
    return {
        "status": "CHECKPOINT_RESTORE_VALIDATED",
        "phase": phase,
        "checkpoint": str(checkpoint),
        "time_ps": float(restored["time_ps"]),
        "peak_temperature_K": peak,
        "max_ag_liquid_fraction": float(model.state.ag_liquid_fraction.max()),
        "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
        "source_power_W": model.actual_source_power_W,
        "injected_energy_J": model.injected_energy_J,
        "total_enthalpy_J": float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3),
        "energy_balance_error_J": energy_error,
        "relative_energy_balance_error": relative_error,
        "first_liquid_time_ps": tracker.first_liquid_time_ps,
        "first_connected_liquid_time_ps": tracker.first_connected_liquid_time_ps,
        "first_multivoxel_pool_time_ps": tracker.first_multivoxel_pool_time_ps,
        "stable_pool_assessment": stability,
    }


def run_smoke_test(refinement: int, output_dir: Path, smoke_ps: float) -> dict[str, Any]:
    """Exercise exact source, stable dt, energy closure, determinism, and restart."""
    output_dir.mkdir(parents=True, exist_ok=True)
    with COMPOSITION.open(encoding="utf-8") as handle:
        composition = json.load(handle)
    target_ag = float(composition["target_phase_volume_fraction"]["Ag"])
    repeats: list[QSGSEnthalpyGrid] = []
    repeat_trackers: list[LiquidClusterTracker] = []
    durations: list[float] = []
    for _ in range(2):
        model = model_for(refinement)
        tracker = tracker_for(model)
        initial_temperature = model.state.temperature_K.copy()
        if not np.all(initial_temperature == 300.0):
            raise AssertionError("Fresh smoke-test model is not initialized at 300 K")
        if abs(float(model.ag_fraction.mean()) - target_ag) > 1.0 / model.ag_fraction.size:
            raise AssertionError("Mesh Ag fraction does not match the pore-free target")
        if float(model.source_cell_power_W[model.ag_fraction == 0.0].sum()) != 0.0:
            raise AssertionError("External source deposits power outside Ag voxels")
        tracker.update(0.0, model.state.ag_liquid_fraction, model.ag_fraction > 0.0)
        started = time.perf_counter()
        model.run(
            smoke_ps,
            source_on=True,
            sample_interval_ps=smoke_ps,
            track_liquid_connectivity=True,
            liquid_threshold=0.5,
            tracker=tracker,
            progress=False,
        )
        durations.append(time.perf_counter() - started)
        error, relative = _relative_energy_error(model)
        expected = model.actual_source_power_W * model.time_s
        if not np.isclose(model.injected_energy_J, expected, rtol=1.0e-13, atol=1.0e-32):
            raise AssertionError("Injected energy does not equal exact source power times duration")
        if abs(relative) > 1.0e-8:
            raise AssertionError(f"Smoke-test energy balance failed: {relative:.3e}")
        if float(model.state.temperature_K.max()) <= 300.0:
            raise AssertionError("Heat source did not produce a measurable temperature increase")
        if np.any(model.state.ag_liquid_fraction != 0.0):
            raise AssertionError("1 fs smoke test unexpectedly created liquid")
        repeats.append(model)
        repeat_trackers.append(tracker)

    first, second = repeats
    if not np.array_equal(first.enthalpy_density_J_m3, second.enthalpy_density_J_m3):
        raise AssertionError("Independent smoke-test enthalpy fields are not bitwise reproducible")
    if not np.array_equal(first.state.temperature_K, second.state.temperature_K):
        raise AssertionError("Independent smoke-test temperature fields are not bitwise reproducible")
    if first.actual_source_power_W != second.actual_source_power_W:
        raise AssertionError("Independent smoke-test source powers differ")

    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = checkpoint_dir / f"caseA_ag607_{refinement}x_smoke_{tag_ps(smoke_ps)}ps.npz"
    first.save_checkpoint(checkpoint, repeat_trackers[0])
    restored_model = model_for(refinement)
    restored_tracker = tracker_for(restored_model)
    receipt = validate_checkpoint(restored_model, restored_tracker, checkpoint, phase="smoke")
    if not np.array_equal(restored_model.enthalpy_density_J_m3, first.enthalpy_density_J_m3):
        raise AssertionError("Smoke-test checkpoint restore changed enthalpy")
    if not np.array_equal(restored_model.state.temperature_K, first.state.temperature_K):
        raise AssertionError("Smoke-test checkpoint restore changed temperature")
    receipt["checkpoint_sha256"] = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    receipt["grid_shape_xyz"] = list(first.temperature_shape)
    receipt["ag_volume_fraction"] = float(first.ag_fraction.mean())
    receipt["scaffold_volume_fraction"] = float(first.sc_fraction.mean())
    receipt["initial_stable_timestep_fs"] = first.stable_timestep_s() * 1.0e15
    receipt["smoke_duration_ps"] = smoke_ps
    receipt["reproducibility_bitwise_equal"] = True
    receipt["repeat_wall_time_s"] = durations
    atomic_json(output_dir / "preflight.json", receipt)
    return receipt


def _write_segment_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if rows:
        write_csv(path, rows)


def run_full(args: argparse.Namespace) -> dict[str, Any]:
    global STOP_REQUESTED
    STOP_REQUESTED = False
    refinement = args.refinement
    output_dir: Path = args.output_dir.resolve()
    checkpoint_dir = output_dir / "checkpoints"
    segment_dir = output_dir / "segments"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    segment_dir.mkdir(parents=True, exist_ok=True)
    lock_handle = (output_dir / ".runner.lock").open("a+", encoding="utf-8")
    try:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        raise RuntimeError(f"Another solver already owns this output directory: {output_dir}") from exc
    lock_handle.seek(0)
    lock_handle.truncate()
    lock_handle.write(f"pid={os.getpid()}\n")
    lock_handle.flush()
    state_path = output_dir / "run_state.json"

    model = model_for(refinement)
    tracker = tracker_for(model)
    run_state: dict[str, Any] = {}
    if args.resume:
        if not state_path.exists():
            raise FileNotFoundError(f"Cannot resume: missing {state_path}")
        run_state = json.loads(state_path.read_text(encoding="utf-8"))
        if int(run_state["refinement"]) != refinement:
            raise ValueError("Resume refinement differs from run_state.json")
        checkpoint = (output_dir / run_state["latest_checkpoint"]).resolve()
        validate_checkpoint(model, tracker, checkpoint, phase=run_state["phase"])
        phase = str(run_state["phase"])
        heat_off_time_ps = run_state.get("heat_off_time_ps")
        heating_stop_reason = run_state.get("heating_stop_reason")
        heat_off_stability_assessment = run_state.get("heat_off_stability_assessment")
        start_wall = time.perf_counter()
    elif args.restart_from is not None:
        if state_path.exists():
            raise FileExistsError("Use --resume for an output directory that already has run_state.json")
        restored = validate_checkpoint(model, tracker, args.restart_from.resolve(), phase=args.restart_phase)
        phase = args.restart_phase
        heat_off_time_ps = None
        heating_stop_reason = None
        heat_off_stability_assessment = None
        start_wall = time.perf_counter()
        run_state = {"restart_checkpoint_validation": restored, "latest_checkpoint": None}
    else:
        if state_path.exists() or any(checkpoint_dir.glob("*.npz")):
            raise FileExistsError("Output already contains a run; use --resume or a fresh output directory")
        tracker.update(0.0, model.state.ag_liquid_fraction, model.ag_fraction > 0.0)
        initial_checkpoint = checkpoint_dir / f"caseA_ag607_{refinement}x_initial_0ps.npz"
        model.save_checkpoint(initial_checkpoint, tracker)
        initial_receipt = validate_checkpoint(model, tracker, initial_checkpoint, phase="heating")
        atomic_json(checkpoint_dir / (initial_checkpoint.stem + ".validation.json"), initial_receipt)
        phase = "heating"
        heat_off_time_ps = None
        heating_stop_reason = None
        heat_off_stability_assessment = None
        start_wall = time.perf_counter()
        run_state = {"latest_checkpoint": str(initial_checkpoint.relative_to(output_dir))}

    if phase == "complete":
        summary_path = output_dir / "summary.json"
        if summary_path.exists():
            return json.loads(summary_path.read_text(encoding="utf-8"))

    geometry_sha = QSGSEnthalpyGrid._sha256_file(GEOMETRY)
    milestones: dict[str, Any] = dict(run_state.get("milestones", {}))
    segment_number = int(run_state.get("segment_number", 0))
    last_checkpoint = run_state.get("latest_checkpoint")
    status = "running"
    stop_reason: str | None = None
    stability_assessment = assess_stable_pool(
        tracker,
        model.time_s * 1.0e12,
        args.stability_duration_ps,
        args.stability_window_ps,
        args.stability_range_fraction,
    )

    def save_state(current_phase: str, current_status: str = "running") -> None:
        nonlocal last_checkpoint
        payload = {
            "status": current_status,
            "refinement": refinement,
            "phase": current_phase,
            "latest_checkpoint": last_checkpoint,
            "time_ps": model.time_s * 1.0e12,
            "heat_off_time_ps": heat_off_time_ps,
            "heating_stop_reason": heating_stop_reason,
            "heat_off_stability_assessment": heat_off_stability_assessment,
            "segment_number": segment_number,
            "peak_temperature_K": float(model.state.temperature_K.max()),
            "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
            "first_liquid_time_ps": tracker.first_liquid_time_ps,
            "first_connected_liquid_time_ps": tracker.first_connected_liquid_time_ps,
            "first_multivoxel_pool_time_ps": tracker.first_multivoxel_pool_time_ps,
            "stable_pool_assessment": stability_assessment,
            "milestones": milestones,
            "geometry_sha256": geometry_sha,
            "process_parameters": {
                "boundary_case": "A",
                "source_peak_flux_eV_ps_A2": 5.0,
                "source_gaussian_radius_A": 30.0,
                "source_absorption_depth_A": 20.0,
                "maximum_temperature_K": args.max_temperature_K,
                "maximum_heating_time_ps": args.max_heating_ps,
                "maximum_cooling_duration_ps": args.max_cooling_ps,
                "stability_duration_ps": args.stability_duration_ps,
                "stability_window_ps": args.stability_window_ps,
                "stability_volume_range_fraction": args.stability_range_fraction,
                "sample_interval_ps": args.sample_ps,
                "segment_interval_ps": args.segment_ps,
            },
        }
        atomic_json(state_path, payload)

    def save_segment_checkpoint(current_phase: str, stability: dict[str, Any]) -> Path:
        nonlocal last_checkpoint
        current_ps = model.time_s * 1.0e12
        prefix = f"caseA_ag607_{refinement}x_{current_phase}"
        cp = checkpoint_dir / f"{prefix}_{tag_ps(current_ps)}ps.npz"
        model.save_checkpoint(cp, tracker)
        receipt = validate_checkpoint(model, tracker, cp, phase=current_phase, stability=stability)
        receipt["checkpoint_sha256"] = hashlib.sha256(cp.read_bytes()).hexdigest()
        atomic_json(checkpoint_dir / (cp.stem + ".validation.json"), receipt)
        last_checkpoint = str(cp.relative_to(output_dir))
        return cp

    def run_segment(current_phase: str, end_ps: float) -> dict[str, Any]:
        nonlocal segment_number, stability_assessment
        current_ps = model.time_s * 1.0e12
        duration = end_ps - current_ps
        if duration <= 1.0e-12:
            raise ValueError(f"Non-positive segment duration: {current_ps} to {end_ps} ps")
        segment_number += 1
        started = time.perf_counter()
        rows = model.run(
            duration,
            source_on=(current_phase == "heating"),
            sample_interval_ps=min(args.sample_ps, duration),
            track_liquid_connectivity=True,
            liquid_threshold=0.5,
            progress=True,
            tracker=tracker,
        )
        _write_segment_csv(
            segment_dir / f"{current_phase}_{segment_number:05d}_{tag_ps(model.time_s * 1.0e12)}ps.csv",
            rows,
        )
        stability_assessment = assess_stable_pool(
            tracker,
            model.time_s * 1.0e12,
            args.stability_duration_ps,
            args.stability_window_ps,
            args.stability_range_fraction,
        )
        checkpoint = save_segment_checkpoint(current_phase, stability_assessment)
        peak = float(model.state.temperature_K.max())
        energy_error, relative_error = _relative_energy_error(model)
        if not np.isfinite(peak) or not np.isfinite(relative_error):
            raise FloatingPointError("Non-finite state or energy ledger")
        if np.any(model.state.ag_liquid_fraction < 0.0) or np.any(model.state.ag_liquid_fraction > 1.0):
            raise ArithmeticError("Ag liquid fraction left [0, 1]")
        if abs(relative_error) > 1.0e-8:
            raise ArithmeticError(f"Energy balance relative error {relative_error:.3e} exceeds 1e-8")
        last_row = rows[-1] if rows else {}
        event_names: list[str] = []
        for key, milestone in (
            ("first_liquid_time_ps", tracker.first_liquid_time_ps),
            ("first_connected_liquid_time_ps", tracker.first_connected_liquid_time_ps),
            ("first_multivoxel_pool_time_ps", tracker.first_multivoxel_pool_time_ps),
        ):
            if milestone is not None and key not in milestones:
                milestones[key] = float(milestone)
                milestones[key + "_checkpoint"] = str(checkpoint.relative_to(output_dir))
                event_names.append(key)
        record = {
            "segment": segment_number,
            "phase": current_phase,
            "start_time_ps": current_ps,
            "end_time_ps": model.time_s * 1.0e12,
            "wall_time_s": time.perf_counter() - started,
            "peak_temperature_K": peak,
            "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
            "max_ag_liquid_fraction": float(model.state.ag_liquid_fraction.max()),
            "connected_cluster_count": int(last_row.get("connected_liquid_cluster_count", 0)),
            "largest_cluster_voxels": int(last_row.get("largest_liquid_cluster_voxels", 0)),
            "largest_cluster_volume_A3": float(last_row.get("largest_liquid_cluster_volume_A3", 0.0)),
            "first_liquid_time_ps": tracker.first_liquid_time_ps,
            "first_connected_liquid_time_ps": tracker.first_connected_liquid_time_ps,
            "first_multivoxel_pool_time_ps": tracker.first_multivoxel_pool_time_ps,
            "stable_pool_assessment": stability_assessment,
            "injected_energy_J": model.injected_energy_J,
            "total_enthalpy_J": float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3),
            "energy_balance_error_J": energy_error,
            "relative_energy_balance_error": relative_error,
            "checkpoint": str(checkpoint.relative_to(output_dir)),
            "events": event_names,
        }
        append_jsonl(output_dir / "segments.jsonl", record)
        return record

    if not args.resume and args.restart_from is None:
        save_state("heating", "running")

    # Continue the heating phase until the predefined stable-pool rule, the
    # heating time cap, or the temperature protection limit is reached.
    if phase == "heating":
        while model.time_s * 1.0e12 < args.max_heating_ps - 1.0e-10:
            current_ps = model.time_s * 1.0e12
            peak = float(model.state.temperature_K.max())
            if peak >= args.max_temperature_K:
                status, stop_reason = "safety_stop", "temperature_protection_reached"
                save_state("complete", status)
                phase = "complete"
                break
            if STOP_REQUESTED:
                status, stop_reason = "paused", "safe_stop_requested"
                save_state("heating", status)
                break
            interval = min(args.segment_ps, args.max_heating_ps - current_ps)
            if peak >= args.max_temperature_K - 50.0:
                interval = min(interval, 0.002)
            run_segment("heating", current_ps + interval)
            save_state("heating")
            if float(model.state.temperature_K.max()) >= args.max_temperature_K:
                status, stop_reason = "safety_stop", "temperature_protection_reached"
                save_state("complete", status)
                phase = "complete"
                break
            if stability_assessment.get("stable"):
                heat_off_time_ps = model.time_s * 1.0e12
                heating_stop_reason = "stable_pool_thresholds_met"
                heat_off_stability_assessment = dict(stability_assessment)
                milestones["stable_pool_time_ps"] = heat_off_time_ps
                milestones["stable_pool_checkpoint"] = last_checkpoint
                phase = "cooling"
                status = "running"
                save_state(phase, status)
                break
            if STOP_REQUESTED:
                status, stop_reason = "paused", "safe_stop_requested"
                save_state("heating", status)
                break
        if phase == "heating" and status == "running":
            heat_off_time_ps = model.time_s * 1.0e12
            heating_stop_reason = "maximum_heating_time_reached_without_stable_pool"
            heat_off_stability_assessment = dict(stability_assessment)
            phase = "cooling"
            save_state(phase)

    # Once heating has stopped, retain the insulated Case A domain and let
    # internal Ag conduction redistribute energy until all Ag is solid.
    if phase == "cooling":
        if heat_off_time_ps is None:
            heat_off_time_ps = model.time_s * 1.0e12
        cooling_end_ps = float(heat_off_time_ps) + args.max_cooling_ps
        while model.time_s * 1.0e12 < cooling_end_ps - 1.0e-10:
            current_ps = model.time_s * 1.0e12
            if not np.any(model.state.ag_liquid_fraction > 1.0e-12):
                stop_reason = "fully_solid"
                status = "complete" if (heat_off_stability_assessment or {}).get("stable") else "complete_without_stable_pool"
                phase = "complete"
                save_state(phase, status)
                break
            if STOP_REQUESTED:
                status, stop_reason = "paused", "safe_stop_requested"
                save_state("cooling", status)
                break
            interval = min(args.segment_ps, cooling_end_ps - current_ps)
            run_segment("cooling", current_ps + interval)
            save_state("cooling")
            if float(model.state.temperature_K.max()) >= args.max_temperature_K:
                status, stop_reason = "safety_stop", "temperature_protection_reached_during_cooling"
                phase = "complete"
                save_state(phase, status)
                break
        if phase == "cooling" and status == "running":
            if not np.any(model.state.ag_liquid_fraction > 1.0e-12):
                stop_reason = "fully_solid"
                status = "complete" if (heat_off_stability_assessment or {}).get("stable") else "complete_without_stable_pool"
            else:
                stop_reason = "cooling_time_cap_reached_before_full_solidification"
                status = "complete_partial_cooling"
            phase = "complete"
            save_state(phase, status)

    # A safety stop is deliberately terminal because vaporization is absent.
    if phase == "complete" and stop_reason is None:
        stop_reason = "completed"
    if phase == "heating" and status == "paused":
        stop_reason = "safe_stop_requested"
    final_ps = model.time_s * 1.0e12
    if tracker.history:
        write_cluster_csv(output_dir / "heating_and_cooling.clusters.csv", tracker.history)
    energy_error, relative_error = _relative_energy_error(model)
    stable_final = assess_stable_pool(
        tracker,
        final_ps if phase == "heating" else (float(heat_off_time_ps) if heat_off_time_ps is not None else final_ps),
        args.stability_duration_ps,
        args.stability_window_ps,
        args.stability_range_fraction,
    )
    solidified = not bool(np.any(model.state.ag_liquid_fraction > 1.0e-12))
    summary = {
        "status": status,
        "refinement": refinement,
        "grid_shape_xyz": list(model.temperature_shape),
        "geometry": str(GEOMETRY),
        "geometry_sha256": geometry_sha,
        "composition": json.loads(COMPOSITION.read_text(encoding="utf-8")),
        "boundary_case": "A: zero Ti3SiC2 heat storage/conductivity; insulated external faces",
        "source": {
            "active_phase": "Ag only",
            "spatial_form": "lateral Gaussian times depth exponential",
            "peak_flux_eV_ps_A2": 5.0,
            "gaussian_radius_A": 30.0,
            "absorption_depth_A": 20.0,
            "integration": model.source_integration_method,
            "power_W": model.actual_source_power_W,
        },
        "initial_temperature_K": 300.0,
        "melting_temperature_reference_K": model.melting_temperature_K,
        "heating_stop_reason": heating_stop_reason,
        "stop_reason": stop_reason,
        "heat_off_time_ps": heat_off_time_ps,
        "final_time_ps": final_ps,
        "final_peak_temperature_K": float(model.state.temperature_K.max()),
        "maximum_temperature_seen_K": model.max_temperature_seen_K,
        "safety_temperature_limit_K": args.max_temperature_K,
        "first_liquid_time_ps": tracker.first_liquid_time_ps,
        "first_connected_liquid_time_ps": tracker.first_connected_liquid_time_ps,
        "first_multivoxel_pool_time_ps": tracker.first_multivoxel_pool_time_ps,
        "stable_pool_assessment_at_heat_off": heat_off_stability_assessment,
        "stable_pool_assessment_final": stable_final,
        "ag_liquid_volume_fraction_final": model._liquid_ag_volume_fraction(),
        "max_ag_liquid_fraction_final": float(model.state.ag_liquid_fraction.max()),
        "fully_solid_final": solidified,
        "injected_energy_J": model.injected_energy_J,
        "total_enthalpy_J": float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3),
        "energy_balance_error_J": energy_error,
        "relative_energy_balance_error": relative_error,
        "elapsed_wall_time_s": time.perf_counter() - start_wall,
        "checkpoint_count": len(list(checkpoint_dir.glob("*.npz"))),
        "last_checkpoint": last_checkpoint,
        "milestones": milestones,
        "stability_rule": {
            "same_main_cluster_persists_ps": args.stability_duration_ps,
            "latest_volume_window_ps": args.stability_window_ps,
            "maximum_volume_range_fraction": args.stability_range_fraction,
            "cluster_sampling_ps": args.sample_ps,
        },
        "physical_limitations": [
            "Case A assigns zero thermal capacity and conductivity to Ti3SiC2 and insulates all external faces.",
            "Cooling is internal energy redistribution and does not model environmental heat loss.",
            "The phase model has no melt flow, surface deformation, vaporization, or atomistic dynamics.",
            "Thermophysical values remain candidate values pending source-table verification.",
        ],
    }
    atomic_json(output_dir / "summary.json", summary)
    if status in {"complete", "complete_without_stable_pool", "complete_partial_cooling", "safety_stop"}:
        save_state("complete", status)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refinement", type=int, choices=(1, 2), required=True)
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--smoke-ps", type=float, default=0.001)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--restart-from", type=Path, default=None)
    parser.add_argument("--restart-phase", choices=("heating", "cooling"), default="heating")
    parser.add_argument("--max-heating-ps", type=float, default=15.0)
    parser.add_argument("--max-cooling-ps", type=float, default=5.0)
    parser.add_argument("--max-temperature-K", type=float, default=2300.0)
    parser.add_argument("--segment-ps", type=float, default=0.05)
    parser.add_argument("--sample-ps", type=float, default=0.01)
    parser.add_argument("--stability-duration-ps", type=float, default=3.0)
    parser.add_argument("--stability-window-ps", type=float, default=2.0)
    parser.add_argument("--stability-range-fraction", type=float, default=0.05)
    args = parser.parse_args()
    if not GEOMETRY.exists():
        parser.error(f"Missing new-composition geometry: {GEOMETRY}")
    if args.smoke_test:
        if args.smoke_ps <= 0.0:
            parser.error("--smoke-ps must be positive")
        result = run_smoke_test(args.refinement, args.output_dir.resolve(), args.smoke_ps)
    else:
        if args.max_heating_ps <= 0.0 or args.max_cooling_ps < 0.0:
            parser.error("Heating time must be positive and cooling time nonnegative")
        if args.max_temperature_K <= 300.0 or args.segment_ps <= 0.0 or args.sample_ps <= 0.0:
            parser.error("Temperature and sampling/segment limits must be positive")
        if args.stability_duration_ps <= 0.0 or args.stability_window_ps <= 0.0:
            parser.error("Stable-pool duration and window must be positive")
        if not 0.0 < args.stability_range_fraction <= 1.0:
            parser.error("Stability volume range must be in (0,1]")
        result = run_full(args)
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
