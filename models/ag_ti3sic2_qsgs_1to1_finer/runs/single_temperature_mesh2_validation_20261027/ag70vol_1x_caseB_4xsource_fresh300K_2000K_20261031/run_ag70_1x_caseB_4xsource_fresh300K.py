#!/usr/bin/env python3
"""Preflight and run the Ag70 Case B model on the 1x grid."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


OUT_DIR = Path(__file__).resolve().parent
RUNS_BASE = OUT_DIR.parent
PROJECT_DIR = OUT_DIR.parents[2]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = RUNS_BASE / "qsgs_morphology_ag70vol.npz"
NUMBA_RUNTIME = Path("/workspace/.numba-runtime")
if NUMBA_RUNTIME.is_dir():
    sys.path.insert(0, str(NUMBA_RUNTIME))
sys.path.insert(0, str(PHASE_DIR))

try:
    import numba
except ImportError as exc:  # pragma: no cover - this run requires the accelerator
    raise SystemExit("Numba is required for this 1x Case B long run; install numba and llvmlite first") from exc

from qsgs_enthalpy_solver import (  # noqa: E402
    EV_PER_PS_ANGSTROM2_TO_W_M2,
    QSGSEnthalpyGrid,
    write_cluster_csv,
)
from qsgs_liquid_tracking import LiquidClusterTracker  # noqa: E402


PEAK_FLUX_EV_PS_A2 = 5.3953185
BASE_PEAK_FLUX_EV_PS_A2 = PEAK_FLUX_EV_PS_A2 / 3.0
SOURCE_SWITCH_MULTIPLIER = 3.5
SOURCE_SWITCH_GAP_MW = 0.01
GAUSSIAN_RADIUS_A = 50.0
ABSORPTION_DEPTH_A = 20.0
SOURCE_CENTER_XY_A = (149.8175, 150.9287)
SINK_TEMPERATURE_K = 300.0
SCAFFOLD_LIMIT_K = 2000.0
# Segment lengths keep each checkpointing loop bounded without imposing a
# total heating or cooling time limit.
HEAT_SEGMENT_PS = 100.0
COOL_SEGMENT_PS = 100.0
SAMPLE_INTERVAL_PS = 0.1
PROGRESS_INTERVAL_PS = 0.01
CHECKPOINT_INTERVAL_PS = 0.5
LIQUID_THRESHOLD = 0.5


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def build_model(source_multiplier: float = 3.0) -> QSGSEnthalpyGrid:
    if not GEOMETRY.is_file():
        raise FileNotFoundError(f"Missing Ag70 QSGS geometry: {GEOMETRY}")
    model = QSGSEnthalpyGrid(
        None,
        "B",
        peak_flux_eV_ps_A2=PEAK_FLUX_EV_PS_A2,
        gaussian_radius_A=GAUSSIAN_RADIUS_A,
        absorption_depth_A=ABSORPTION_DEPTH_A,
        source_center_xy_A=SOURCE_CENTER_XY_A,
        native_geometry_npz=GEOMETRY,
        source_geometry_npz=GEOMETRY,
        mesh_refinement=1,
        bottom_sink_temperature_K=SINK_TEMPERATURE_K,
    )
    model.source_multiplier = 3.0
    model.source_schedule_events = []
    requested_multiplier = float(source_multiplier)
    if requested_multiplier > 3.0 and requested_multiplier != SOURCE_SWITCH_MULTIPLIER:
        # Preserve the same floating-point source-scaling path used by staged
        # heating transitions (3x -> 3.5x -> requested multiplier).
        set_source_multiplier(model, SOURCE_SWITCH_MULTIPLIER)
        set_source_multiplier(model, requested_multiplier)
    elif requested_multiplier != 3.0:
        set_source_multiplier(model, requested_multiplier)
    return model


def set_source_multiplier(model: QSGSEnthalpyGrid, source_multiplier: float) -> None:
    """Scale the fixed normalized Ag source without changing its spatial shape."""
    old_multiplier = float(getattr(model, "source_multiplier", 3.0))
    if old_multiplier <= 0.0 or source_multiplier <= 0.0:
        raise ValueError("Source multipliers must be positive")
    scale = float(source_multiplier) / old_multiplier
    model.peak_flux_W_m2 = (
        BASE_PEAK_FLUX_EV_PS_A2
        * float(source_multiplier)
        * EV_PER_PS_ANGSTROM2_TO_W_M2
    )
    model.source_W_m3 = model.source_W_m3 * scale
    model.source_cell_power_W = model.source_cell_power_W * scale
    model.source_power_in_roundoff_only_cells_W *= scale
    model.actual_source_power_W = float(model.source_cell_power_W.sum(dtype=np.float64))
    model.source_multiplier = float(source_multiplier)


def volume_integral(field: np.ndarray, model: QSGSEnthalpyGrid) -> float:
    return float(np.sum(field, dtype=np.float64) * model.cell_volume_m3)


def phase_energy_breakdown(model: QSGSEnthalpyGrid) -> dict[str, float]:
    temperature = model.state.temperature_K
    liquid = model.state.ag_liquid_fraction
    ag_mask = model.ag_fraction > 0.0
    sc_mask = model.sc_fraction > 0.0
    ag_temperature = np.where(ag_mask, temperature, 300.0)
    ag_fraction = np.where(ag_mask, liquid, 0.0)
    ag_h = np.asarray(model.ag.specific_enthalpy_J_kg(ag_temperature, ag_fraction))
    ag_density = model.ag_fraction * model.ag.reference_density_kg_m3 * ag_h
    ag_latent_density = (
        model.ag_fraction * model.ag.reference_density_kg_m3 * model.ag.latent_heat_J_kg * liquid
    )
    sc_temperature = np.where(sc_mask, temperature, 300.0)
    sc_h = np.asarray(model.sc.specific_enthalpy_J_kg(sc_temperature))
    sc_density = model.sc_fraction * model.sc.density_kg_m3 * sc_h
    ag_total = volume_integral(ag_density, model)
    ag_latent = volume_integral(ag_latent_density, model)
    sc_total = volume_integral(sc_density, model)
    system_total = volume_integral(model.enthalpy_density_J_m3, model)
    return {
        "ag_sensible_enthalpy_J": ag_total - ag_latent,
        "ag_latent_enthalpy_J": ag_latent,
        "ag_total_enthalpy_J": ag_total,
        "ti3sic2_stored_enthalpy_J": sc_total,
        "system_total_enthalpy_J": system_total,
        "cumulative_source_energy_J": model.injected_energy_J,
        "cumulative_bottom_sink_energy_J": model.removed_energy_J,
        "energy_balance_residual_J": system_total - model.injected_energy_J + model.removed_energy_J,
    }


def metrics(
    model: QSGSEnthalpyGrid,
    tracker_summary: dict[str, Any] | None = None,
    source_on: bool = True,
) -> dict[str, Any]:
    temperature = model.state.temperature_K
    liquid = model.state.ag_liquid_fraction
    ag_mask = model.ag_fraction > 0.0
    sc_mask = model.sc_fraction > 0.0
    ag_count = float(model.ag_fraction.sum())
    sc_count = float(model.sc_fraction.sum())
    row: dict[str, Any] = {
        "time_ps": model.time_s * 1.0e12,
        "dt_fs": model.last_run_step_limit_s * 1.0e15 if model.last_run_step_limit_s else model.stable_timestep_s() * 1.0e15,
        "ag_peak_temperature_K": float(temperature[ag_mask].max()),
        "ag_mean_temperature_K": float(np.sum(model.ag_fraction * temperature) / ag_count),
        "ti3sic2_peak_temperature_K": float(temperature[sc_mask].max()),
        "ti3sic2_mean_temperature_K": float(np.sum(model.sc_fraction * temperature) / sc_count),
        "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
        "ag_peak_liquid_fraction": float(liquid[ag_mask].max()),
        "bottom_sink_power_W": float(
            model.bottom_sink_power_density_W_m3().sum(dtype=np.float64) * model.cell_volume_m3
        ),
        "source_power_W": model.actual_source_power_W if source_on else 0.0,
        "cumulative_source_energy_J": model.injected_energy_J,
        "cumulative_bottom_sink_energy_J": model.removed_energy_J,
        "total_enthalpy_J": volume_integral(model.enthalpy_density_J_m3, model),
    }
    if tracker_summary:
        row.update(tracker_summary)
        row["largest_liquid_cluster_volume_nm3"] = (
            float(row.get("largest_liquid_cluster_volume_A3", 0.0)) / 1000.0
        )
    row.update(phase_energy_breakdown(model))
    return row


def append_csv(path: Path, row: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and path.stat().st_size > 0
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), extrasaction="ignore", lineterminator="\n")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def runtime_metadata(model: QSGSEnthalpyGrid | None = None) -> dict[str, Any]:
    try:
        branch = subprocess.check_output(
            ["git", "-C", str(PROJECT_DIR), "branch", "--show-current"], text=True
        ).strip()
        commit = subprocess.check_output(
            ["git", "-C", str(PROJECT_DIR), "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception:
        branch, commit = "unknown", "unknown"
    multiplier = float(getattr(model, "source_multiplier", 3.0))
    active_flux = BASE_PEAK_FLUX_EV_PS_A2 * multiplier
    active_power_mW = (
        float(model.actual_source_power_W) * 1.0e3
        if model is not None
        else 4.654075634 * multiplier / 3.0
    )
    schedule_events = list(getattr(model, "source_schedule_events", []))
    schedule_multipliers = {3.0, multiplier}
    for event in schedule_events:
        for key in ("from_multiplier", "to_multiplier"):
            if event.get(key) is not None:
                schedule_multipliers.add(float(event[key]))
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_branch": branch,
        "git_commit": commit,
        "numba_version": numba.__version__,
        "numba_threads": numba.get_num_threads(),
        "source_multiplier_vs_baseline": multiplier,
        "source_schedule_multipliers": sorted(schedule_multipliers),
        "source_switch_gap_mW": SOURCE_SWITCH_GAP_MW,
        "source_schedule_events": schedule_events,
        "source_peak_flux_eV_ps_A2": active_flux,
        "source_peak_flux_W_m2": active_flux * EV_PER_PS_ANGSTROM2_TO_W_M2,
        "source_gaussian_radius_A": GAUSSIAN_RADIUS_A,
        "source_absorption_depth_A": ABSORPTION_DEPTH_A,
        "source_center_xy_A": list(SOURCE_CENTER_XY_A),
        "requested_source_power_mW": active_power_mW,
        "sink_temperature_K": SINK_TEMPERATURE_K,
        "scaffold_property_limit_K": SCAFFOLD_LIMIT_K,
        "thermal_model": "single-temperature enthalpy solver; no TTM",
        "outer_boundary": "x/y periodic; z insulated except bottom Ag-contact ideal isothermal sink",
        "temperature_range_scope": "Barsoum Cp/k correlations; values beyond source ranges are extrapolations",
    }


def preflight() -> dict[str, Any]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    model = build_model()
    assert model.temperature_shape == (100, 100, 50)
    assert model.mesh_refinement == 1
    assert abs(float(model.ag_fraction.mean()) - 0.699996) < 2.0e-8
    assert abs(float(model.sc_fraction.mean()) - 0.300004) < 2.0e-8
    assert model.bottom_sink_contact_voxels == 8249
    assert np.isclose(model.bottom_sink_contact_area_m2 * 1.0e18, 746.0978590068572, atol=1.0e-8)
    assert model.actual_source_power_W > 0.0
    assert np.isclose(model.actual_source_power_W * 1.0e3, 4.654075634, rtol=2.0e-4)
    assert np.all(model.source_cell_power_W[model.sc_fraction > 0.0] == 0.0)
    assert model.sc.maximum_temperature_K == 2000.0
    assert np.isfinite(model.sc.cp_J_kg_K(2000.0))
    assert np.isfinite(model.sc.conductivity_W_m_K(2000.0))
    assert np.isfinite(model.sc.specific_enthalpy_J_kg(2000.0))
    with np.testing.assert_raises_regex(ValueError, "2000.00 K"):
        model.sc.cp_J_kg_K(2000.01)
    stable_dt_s = model.stable_timestep_s()
    assert np.isfinite(stable_dt_s) and stable_dt_s > 0.0

    deterministic = build_model()
    for _ in range(20):
        model.advance(stable_dt_s, source_on=True, enforce_stability=False, collect_metrics=False)
        deterministic.advance(stable_dt_s, source_on=True, enforce_stability=False, collect_metrics=False)
    assert np.array_equal(model.enthalpy_density_J_m3, deterministic.enthalpy_density_J_m3)
    assert np.array_equal(model.state.temperature_K, deterministic.state.temperature_K)
    assert np.array_equal(model.state.ag_liquid_fraction, deterministic.state.ag_liquid_fraction)

    short_runs: list[dict[str, Any]] = []
    for duration_ps in (0.1, 1.0):
        current = build_model()
        start = time.monotonic()
        rows: list[dict[str, Any]] = []
        next_progress = PROGRESS_INTERVAL_PS
        next_sample = SAMPLE_INTERVAL_PS
        target_time_s = duration_ps * 1.0e-12
        step_limit_s = current.stable_timestep_s()
        current.last_run_step_limit_s = step_limit_s
        while current.time_s < target_time_s - 1.0e-24:
            next_progress_s = min(target_time_s, next_progress * 1.0e-12)
            next_sample_s = min(target_time_s, next_sample * 1.0e-12)
            next_target_s = min(next_progress_s, next_sample_s)
            dt = min(step_limit_s, next_target_s - current.time_s)
            is_sample = next_sample_s <= next_progress_s + 1.0e-24
            temperature_event = (
                SCAFFOLD_LIMIT_K
                if current.max_temperature_seen_K >= SCAFFOLD_LIMIT_K - 100.0
                else None
            )
            step_row = current.advance(
                dt,
                source_on=True,
                enforce_stability=False,
                stop_at_sc_temperature_K=temperature_event,
                collect_metrics=(
                    is_sample
                    or next_target_s >= target_time_s - 1.0e-24
                    or temperature_event is not None
                ),
            )
            reached_event = bool(step_row.get("sc_temperature_limit_reached", 0.0))
            if current.time_s >= next_sample_s - 1.0e-24 or reached_event:
                sample = metrics(current)
                relative_error = abs(sample["energy_balance_residual_J"]) / max(
                    abs(sample["cumulative_source_energy_J"]), 1.0e-300
                )
                assert relative_error < 1.0e-8, f"Energy residual too large: {relative_error:g}"
                assert sample["ag_peak_temperature_K"] > 300.0
                assert sample["ti3sic2_peak_temperature_K"] >= 300.0
                assert sample["ti3sic2_peak_temperature_K"] <= SCAFFOLD_LIMIT_K + 1.0e-7
                rows.append(sample)
                if current.time_s >= next_sample_s - 1.0e-24:
                    next_sample += SAMPLE_INTERVAL_PS
            if current.time_s >= next_progress_s - 1.0e-24:
                print(
                    f"preflight {duration_ps:g} ps: t={current.time_s*1e12:.2f} ps, "
                    f"peak={step_row['max_temperature_K']:.2f} K, "
                    f"scaffold={step_row['max_scaffold_temperature_K']:.2f} K",
                    file=sys.stderr,
                    flush=True,
                )
                next_progress += PROGRESS_INTERVAL_PS
            if reached_event:
                print(
                    f"preflight {duration_ps:g} ps: stopped at Ti3SiC2 {SCAFFOLD_LIMIT_K:.0f} K event",
                    file=sys.stderr,
                    flush=True,
                )
                break
        requested_tag = str(duration_ps).replace('.', 'p')
        actual_time_ps = current.time_s * 1.0e12
        actual_tag = f"{actual_time_ps:.3f}".replace('.', 'p')
        checkpoint_dir = OUT_DIR / "preflight" / f"{requested_tag}ps" / "checkpoints"
        checkpoint_path = checkpoint_dir / (
            f"caseB_ag70_1x_3xsource_{actual_tag}ps_"
            f"{'2000K_event' if reached_event else 'duration_complete'}.npz"
        )
        current.save_checkpoint(checkpoint_path)
        checkpoint_validation = current.load_checkpoint(checkpoint_path)
        short_runs.append(
            {
                "requested_duration_ps": duration_ps,
                "actual_duration_ps": actual_time_ps,
                "stop_reason": "ti3sic2_reached_2000K" if reached_event else "requested_duration_completed",
                "wall_time_s": time.monotonic() - start,
                "rows": rows,
                "max_temperature_seen_K": current.max_temperature_seen_K,
                "max_ti3sic2_temperature_K": float(current.state.temperature_K[current.sc_fraction > 0.0].max()),
                "final_energy_balance_residual_J": float(rows[-1]["energy_balance_residual_J"]),
                "checkpoint": str(checkpoint_path.relative_to(OUT_DIR)),
                "checkpoint_reload": checkpoint_validation,
            }
        )

    report = {
        "status": "CASE_B_1X_3X_HEAT_SOURCE_PREFLIGHT_PASSED",
        "configuration": runtime_metadata(),
        "geometry": str(GEOMETRY.relative_to(PROJECT_DIR)),
        "geometry_sha256": QSGSEnthalpyGrid._sha256_file(GEOMETRY),
        "grid_shape_xyz": list(model.temperature_shape),
        "mesh_refinement": model.mesh_refinement,
        "spacing_A": model.spacing_A.tolist(),
        "ag_volume_fraction": float(model.ag_fraction.mean()),
        "ti3sic2_volume_fraction": float(model.sc_fraction.mean()),
        "source_power_W": model.actual_source_power_W,
        "source_power_mW": model.actual_source_power_W * 1.0e3,
        "source_only_on_ag": True,
        "bottom_sink_contact_voxels": model.bottom_sink_contact_voxels,
        "bottom_sink_contact_area_nm2": model.bottom_sink_contact_area_m2 * 1.0e18,
        "stable_timestep_at_300K_fs": stable_dt_s * 1.0e15,
        "high_temperature_properties": {
            "Ti3SiC2_limit_K": model.sc.maximum_temperature_K,
            "Cp_2000K_J_kg_K": model.sc.cp_J_kg_K(2000.0),
            "k_2000K_W_m_K": model.sc.conductivity_W_m_K(2000.0),
            "h_2000K_J_kg": model.sc.specific_enthalpy_J_kg(2000.0),
            "2000K_values_are_extrapolated": True,
        },
        "determinism_20_steps": "passed; enthalpy, temperature, and liquid fields are bitwise equal",
        "short_runs": short_runs,
        "limitations": [
            "Ti3SiC2 Barsoum Cp and k values above reported evidence ranges are correlation extrapolations, not measurements.",
            "Case B assumes ideal Ag/Ti3SiC2 contact; no finite interface resistance is included.",
            "The 300 K sink is a zero-resistance Dirichlet reservoir on bottom Ag contacts; it is not an electrode model.",
        ],
    }
    atomic_json(OUT_DIR / "caseB_preflight.json", report)
    return report


def write_cluster_history(path: Path, tracker: LiquidClusterTracker) -> None:
    write_cluster_csv(path, tracker.history)


def append_series_row(
    path: Path,
    model: QSGSEnthalpyGrid,
    tracker: LiquidClusterTracker | None,
    source_on: bool = True,
) -> None:
    tracker_summary = None
    if tracker is not None:
        tracker_summary, _ = tracker.update(
            model.time_s * 1.0e12,
            model.state.ag_liquid_fraction,
            model.ag_fraction > 0.0,
        )
    row = metrics(model, tracker_summary, source_on=source_on)
    append_csv(path, row)


def save_run_state(path: Path, phase: str, reason: str, model: QSGSEnthalpyGrid, checkpoint: Path) -> None:
    peak_ag = float(model.state.temperature_K[model.ag_fraction > 0.0].max())
    peak_sc = float(model.state.temperature_K[model.sc_fraction > 0.0].max())
    atomic_json(
        path,
        {
            "phase": phase,
            "reason": reason,
            "time_ps": model.time_s * 1.0e12,
            "ag_peak_temperature_K": peak_ag,
            "ti3sic2_peak_temperature_K": peak_sc,
            "max_temperature_seen_K": model.max_temperature_seen_K,
            "source_energy_J": model.injected_energy_J,
            "bottom_sink_energy_J": model.removed_energy_J,
            "source_multiplier": float(getattr(model, "source_multiplier", 3.0)),
            "source_on": phase == "heating",
            "active_source_power_W": model.actual_source_power_W if phase == "heating" else 0.0,
            "configured_source_power_W": model.actual_source_power_W,
            "source_schedule_events": list(getattr(model, "source_schedule_events", [])),
            "pending_source_transition": getattr(model, "pending_source_transition", None),
            "checkpoint": str(checkpoint),
            "updated_utc": datetime.now(timezone.utc).isoformat(),
        },
    )


def run_phase(
    model: QSGSEnthalpyGrid,
    tracker: LiquidClusterTracker,
    source_on: bool,
    duration_limit_ps: float,
    phase: str,
    resume_phase: bool = False,
) -> tuple[str, str]:
    output_csv = OUT_DIR / f"{phase}_timeseries.csv"
    checkpoint_dir = OUT_DIR / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    latest_checkpoint = checkpoint_dir / f"latest_{phase}.npz"
    state_path = OUT_DIR / "run_state.json"
    step_limit_s = model.stable_timestep_s()
    model.last_run_step_limit_s = step_limit_s
    next_progress_ps = (np.floor(model.time_s * 1.0e12 / PROGRESS_INTERVAL_PS + 1.0e-10) + 1.0) * PROGRESS_INTERVAL_PS
    next_sample_ps = (np.floor(model.time_s * 1.0e12 / SAMPLE_INTERVAL_PS + 1.0e-10) + 1.0) * SAMPLE_INTERVAL_PS
    next_checkpoint_ps = (
        np.floor(model.time_s * 1.0e12 / CHECKPOINT_INTERVAL_PS + 1.0e-10) + 1.0
    ) * CHECKPOINT_INTERVAL_PS
    end_time_s = model.time_s + duration_limit_ps * 1.0e-12
    stop_reason = "max_time_reached"

    while model.time_s < end_time_s - 1.0e-24:
        next_target_s = min(
            end_time_s,
            next_progress_ps * 1.0e-12,
            next_sample_ps * 1.0e-12,
            next_checkpoint_ps * 1.0e-12,
        )
        dt = min(step_limit_s, next_target_s - model.time_s)
        if dt <= 1.0e-25:
            if next_target_s * 1.0e12 >= next_progress_ps - 1.0e-12:
                next_progress_ps += PROGRESS_INTERVAL_PS
            if next_target_s * 1.0e12 >= next_sample_ps - 1.0e-12:
                next_sample_ps += SAMPLE_INTERVAL_PS
            if next_target_s * 1.0e12 >= next_checkpoint_ps - 1.0e-12:
                next_checkpoint_ps += CHECKPOINT_INTERVAL_PS
            continue
        current_ps = model.time_s * 1.0e12
        temperature_event = SCAFFOLD_LIMIT_K if source_on else None
        if temperature_event is not None and model.max_temperature_seen_K < temperature_event - 100.0:
            temperature_event = None
        collect_metrics = (
            next_target_s * 1.0e12 >= next_sample_ps - 1.0e-12
            or next_target_s * 1.0e12 >= next_checkpoint_ps - 1.0e-12
            or next_target_s >= end_time_s - 1.0e-24
            or temperature_event is not None
        )
        row = model.advance(
            dt,
            source_on=source_on,
            enforce_stability=False,
            stop_at_sc_temperature_K=temperature_event,
            collect_metrics=collect_metrics,
        )
        time_ps = model.time_s * 1.0e12
        reached_event = bool(row.get("sc_temperature_limit_reached", 0.0))
        if time_ps >= next_progress_ps - 1.0e-10 or reached_event:
            print(
                f"{phase}: t={time_ps:.4f} ps, Ag peak={row['max_temperature_K']:.2f} K, "
                f"scaffold peak={row['max_scaffold_temperature_K']:.2f} K, "
                f"bottom sink={row['bottom_sink_power_W']*1.0e3:.4g} mW",
                file=sys.stderr,
                flush=True,
            )
            while next_progress_ps <= time_ps + 1.0e-10:
                next_progress_ps += PROGRESS_INTERVAL_PS
        sample_written = False
        if time_ps >= next_sample_ps - 1.0e-10 or reached_event:
            append_series_row(output_csv, model, tracker, source_on=source_on)
            sample_written = True
            while next_sample_ps <= time_ps + 1.0e-10:
                next_sample_ps += SAMPLE_INTERVAL_PS
        if time_ps >= next_checkpoint_ps - 1.0e-10 or reached_event:
            model.save_checkpoint(latest_checkpoint, tracker)
            save_run_state(state_path, phase, "temperature_limit" if reached_event else "running", model, latest_checkpoint)
            while next_checkpoint_ps <= time_ps + 1.0e-10:
                next_checkpoint_ps += CHECKPOINT_INTERVAL_PS
        if reached_event:
            stop_reason = "ti3sic2_reached_2000K"
            model.save_checkpoint(checkpoint_dir / "heating_complete_2000K.npz", tracker)
            save_run_state(state_path, "cooling", stop_reason, model, checkpoint_dir / "heating_complete_2000K.npz")
            return "cooling", stop_reason
        if source_on and float(getattr(model, "source_multiplier", 3.0)) < SOURCE_SWITCH_MULTIPLIER:
            bottom_sink_power_W = float(row["bottom_sink_power_W"])
            source_sink_gap_mW = (model.actual_source_power_W - bottom_sink_power_W) * 1.0e3
            if source_sink_gap_mW <= SOURCE_SWITCH_GAP_MW:
                if not sample_written:
                    append_series_row(output_csv, model, tracker, source_on=source_on)
                pending_transition = {
                    "time_ps": float(model.time_s * 1.0e12),
                    "from_multiplier": float(model.source_multiplier),
                    "to_multiplier": SOURCE_SWITCH_MULTIPLIER,
                    "source_sink_gap_mW_at_switch": float(source_sink_gap_mW),
                    "source_power_before_W": float(model.actual_source_power_W),
                    "reason": "3x source-sink gap reached 0.01 mW; pause before raising source to 3.5x",
                }
                model.pending_source_transition = pending_transition
                pause_checkpoint = checkpoint_dir / f"heating_3x_balance_pause_{time_ps:.6f}ps.npz"
                model.save_checkpoint(pause_checkpoint, tracker)
                save_run_state(
                    state_path,
                    "source_switch_pending",
                    "source_balance_gap_reached_0p01mW_pause_before_3p5x",
                    model,
                    pause_checkpoint,
                )
                print(
                    f"source pause: t={time_ps:.4f} ps, 3x source paused at "
                    f"gap={source_sink_gap_mW:.5f} mW; checkpoint={pause_checkpoint.name}",
                    file=sys.stderr,
                    flush=True,
                )
                return "source_switch_pending", "source_balance_gap_reached_0p01mW_pause_before_3p5x"
        if not source_on:
            ag_mask = model.ag_fraction > 0.0
            max_remaining_liquid = float(np.max(model.state.ag_liquid_fraction[ag_mask]))
            if max_remaining_liquid <= 1.0e-8:
                stop_reason = "all_ag_solidified"
                if not sample_written:
                    append_series_row(output_csv, model, tracker, source_on=source_on)
                solid_checkpoint = checkpoint_dir / "cooling_solidification_complete.npz"
                model.save_checkpoint(solid_checkpoint, tracker)
                save_run_state(state_path, "complete", stop_reason, model, solid_checkpoint)
                return "complete", stop_reason
        if time_ps >= end_time_s * 1.0e12 - 1.0e-10:
            break

    model.save_checkpoint(latest_checkpoint, tracker)
    if source_on:
        stop_reason = "heating_segment_complete_before_2000K"
        segment_checkpoint = checkpoint_dir / f"heating_segment_{time_ps:.6f}ps.npz"
        model.save_checkpoint(segment_checkpoint, tracker)
        save_run_state(state_path, "heating", stop_reason, model, segment_checkpoint)
        return "heating", stop_reason
    stop_reason = "cooling_segment_complete_before_full_solidification"
    segment_checkpoint = checkpoint_dir / f"cooling_segment_{time_ps:.6f}ps.npz"
    model.save_checkpoint(segment_checkpoint, tracker)
    save_run_state(state_path, "cooling", stop_reason, model, segment_checkpoint)
    return "cooling", stop_reason


def run_full(
    resume: bool = False,
    source_multiplier_override: float | None = None,
) -> dict[str, Any]:
    preflight_path = OUT_DIR / "caseB_preflight.json"
    if not preflight_path.is_file():
        raise RuntimeError("Run and pass the 0.1 ps and 1 ps preflight before starting the full heating calculation")
    preflight_report = json.loads(preflight_path.read_text())
    if preflight_report.get("status") != "CASE_B_1X_3X_HEAT_SOURCE_PREFLIGHT_PASSED":
        raise RuntimeError("Preflight report is not in the passed state")
    checkpoint_dir = OUT_DIR / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    if source_multiplier_override is not None and float(source_multiplier_override) <= 0.0:
        raise ValueError("Source multiplier must be positive")
    state = json.loads((OUT_DIR / "run_state.json").read_text()) if resume else {}
    source_multiplier = (
        float(state.get("source_multiplier", 3.0))
        if resume
        else float(source_multiplier_override if source_multiplier_override is not None else 3.0)
    )
    model = build_model(source_multiplier=source_multiplier)
    tracker = LiquidClusterTracker(model.spacing_A, model.size_A, liquid_threshold=LIQUID_THRESHOLD)
    if resume:
        checkpoint_path = Path(state["checkpoint"])
        model.load_checkpoint(checkpoint_path, tracker)
        model.source_schedule_events = list(state.get("source_schedule_events", []))
        model.pending_source_transition = state.get("pending_source_transition")
        phase = str(state["phase"])
        stop_reason = str(state.get("reason", "resumed"))
    else:
        phase = "heating"
        stop_reason = "fresh_run"
        if source_multiplier_override is not None:
            stop_reason = f"fresh_run_from_300K_at_{source_multiplier:g}x"
            model.source_schedule_events.append(
                {
                    "time_ps": 0.0,
                    "to_multiplier": source_multiplier,
                    "source_power_after_W": float(model.actual_source_power_W),
                    "reason": "fresh run initialized at 300 K with the user-requested source multiplier",
                }
            )
        model.save_checkpoint(checkpoint_dir / "initial_300K.npz", tracker)
        save_run_state(OUT_DIR / "run_state.json", phase, stop_reason, model, checkpoint_dir / "initial_300K.npz")
        if model.source_schedule_events:
            atomic_json(OUT_DIR / "source_schedule_events.json", {"events": model.source_schedule_events})

    if source_multiplier_override is not None and resume:
        if phase != "heating":
            raise RuntimeError("A source multiplier override is only allowed when resuming a heating run")
        target_multiplier = float(source_multiplier_override)
        previous_multiplier = float(model.source_multiplier)
        if not np.isclose(target_multiplier, previous_multiplier, rtol=0.0, atol=1.0e-12):
            previous_power_W = float(model.actual_source_power_W)
            set_source_multiplier(model, target_multiplier)
            event = {
                "time_ps": float(model.time_s * 1.0e12),
                "from_multiplier": previous_multiplier,
                "to_multiplier": target_multiplier,
                "source_power_before_W": previous_power_W,
                "source_power_after_W": float(model.actual_source_power_W),
                "reason": "user-requested source increase after a validated checkpoint",
            }
            model.source_schedule_events.append(event)
            model.pending_source_transition = None
            tag = f"{target_multiplier:g}".replace(".", "p")
            start_checkpoint = checkpoint_dir / (
                f"heating_{tag}x_start_{model.time_s * 1.0e12:.6f}ps.npz"
            )
            model.save_checkpoint(start_checkpoint, tracker)
            model.save_checkpoint(checkpoint_dir / "latest_heating.npz", tracker)
            save_run_state(
                OUT_DIR / "run_state.json",
                "heating",
                f"reheated_at_{target_multiplier:g}x_after_validated_checkpoint",
                model,
                start_checkpoint,
            )
            atomic_json(OUT_DIR / "source_schedule_events.json", {"events": model.source_schedule_events})
            print(
                f"source restart: t={model.time_s * 1.0e12:.4f} ps, "
                f"{target_multiplier:g}x source={model.actual_source_power_W * 1.0e3:.6f} mW, "
                f"from {previous_multiplier:g}x after validated checkpoint",
                file=sys.stderr,
                flush=True,
            )
            stop_reason = f"reheated_at_{target_multiplier:g}x_after_validated_checkpoint"

    while phase in {"heating", "cooling", "source_switch_pending"}:
        if phase == "source_switch_pending":
            pending = getattr(model, "pending_source_transition", None)
            if pending is None:
                raise RuntimeError("Source-switch pause checkpoint has no pending transition metadata")
            previous_power_W = float(model.actual_source_power_W)
            set_source_multiplier(model, SOURCE_SWITCH_MULTIPLIER)
            event = {
                **pending,
                "source_power_before_W": previous_power_W,
                "source_power_after_W": float(model.actual_source_power_W),
                "reason": "3x source-sink gap reached 0.01 mW; source raised to 3.5x after a validated pause checkpoint",
            }
            model.source_schedule_events.append(event)
            model.pending_source_transition = None
            start_checkpoint = checkpoint_dir / f"heating_3p5x_start_{model.time_s*1.0e12:.6f}ps.npz"
            model.save_checkpoint(start_checkpoint, tracker)
            save_run_state(
                OUT_DIR / "run_state.json",
                "heating",
                "reheated_at_3p5x_after_0p01mW_balance_pause",
                model,
                start_checkpoint,
            )
            atomic_json(OUT_DIR / "source_schedule_events.json", {"events": model.source_schedule_events})
            print(
                f"source restart: t={model.time_s*1.0e12:.4f} ps, "
                f"3.5x source={model.actual_source_power_W*1.0e3:.6f} mW, "
                f"3x pause gap={pending['source_sink_gap_mW_at_switch']:.5f} mW",
                file=sys.stderr,
                flush=True,
            )
            phase = "heating"
            stop_reason = "reheated_at_3p5x_after_0p01mW_balance_pause"
            resume = False
            continue
        is_heating = phase == "heating"
        limit_ps = HEAT_SEGMENT_PS if is_heating else COOL_SEGMENT_PS
        phase, stop_reason = run_phase(
            model,
            tracker,
            source_on=is_heating,
            duration_limit_ps=limit_ps,
            phase="heating" if is_heating else "cooling",
            resume_phase=resume,
        )
        resume = False

    final_state = json.loads((OUT_DIR / "run_state.json").read_text())
    final_checkpoint = Path(final_state["checkpoint"])
    write_cluster_history(OUT_DIR / "liquid_clusters.csv", tracker)
    summary = {
        "status": "AG70_CASE_B_1X_STAGED_SOURCE_COMPLETE" if phase == "complete" else "AG70_CASE_B_1X_STAGED_SOURCE_INCOMPLETE",
        "stop_reason": stop_reason,
        "configuration": runtime_metadata(model),
        "geometry_sha256": QSGSEnthalpyGrid._sha256_file(GEOMETRY),
        "grid_shape_xyz": list(model.temperature_shape),
        "mesh_refinement": model.mesh_refinement,
        "time_ps": model.time_s * 1.0e12,
        "ag_peak_temperature_K": float(model.state.temperature_K[model.ag_fraction > 0.0].max()),
        "ti3sic2_peak_temperature_K": float(model.state.temperature_K[model.sc_fraction > 0.0].max()),
        "max_temperature_seen_K": model.max_temperature_seen_K,
        "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
        "all_ag_solidified": bool(
            float(np.max(model.state.ag_liquid_fraction[model.ag_fraction > 0.0])) <= 1.0e-8
        ),
        "bottom_sink_contact_voxels": model.bottom_sink_contact_voxels,
        "bottom_sink_contact_area_nm2": model.bottom_sink_contact_area_m2 * 1.0e18,
        "actual_source_power_W": model.actual_source_power_W,
        "cumulative_source_energy_J": model.injected_energy_J,
        "cumulative_bottom_sink_energy_J": model.removed_energy_J,
        "energy_breakdown": phase_energy_breakdown(model),
        "checkpoint": str(final_checkpoint.relative_to(OUT_DIR)),
        "limitations": [
            "This is an idealized single-temperature continuum Case B calculation, not TTM or MD.",
            "The ideal bottom sink is a 300 K zero-resistance reservoir on the verified Ag contact mask.",
            "Ti3SiC2 properties above the reported experimental ranges are Barsoum-correlation extrapolations.",
            "Cooling includes external energy removal through the ideal bottom sink; this differs from the prior insulated Case B run.",
        ],
    }
    atomic_json(OUT_DIR / "caseB_run_summary.json", summary)
    write_cluster_history(OUT_DIR / "liquid_clusters.csv", tracker)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("preflight", "run", "resume"))
    parser.add_argument(
        "--source-multiplier",
        type=float,
        default=None,
        help="set the multiplier for a fresh run or change it after loading a heating checkpoint",
    )
    args = parser.parse_args()
    if args.mode == "preflight":
        print(json.dumps(preflight(), indent=2, sort_keys=True))
    else:
        summary = run_full(
            resume=(args.mode == "resume"),
            source_multiplier_override=args.source_multiplier,
        )
        print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
