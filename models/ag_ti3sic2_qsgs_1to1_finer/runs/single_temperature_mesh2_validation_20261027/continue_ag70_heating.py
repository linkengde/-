#!/usr/bin/env python3
"""Safely continue the Ag 70 vol% Case A melt-window run from a checkpoint.

The driver keeps each segment short, writes a restart after every segment, and
stops when a multi-voxel liquid cluster has grown for 0.02 ps beyond its first
appearance or when the configured safety/time limit is reached.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np


RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
INITIAL_1X = RUN_DIR / "ag70vol_melt_window/1x_0p2ps/checkpoints/caseA_ag70_0p2ps_0p2ps.npz"
INITIAL_2X = RUN_DIR / "ag70vol_melt_window/2x_0p2ps/checkpoints/caseA_ag70_0p2ps_0p2ps.npz"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid, write_cluster_csv, write_csv  # noqa: E402
from qsgs_liquid_tracking import LiquidClusterTracker  # noqa: E402


def tag(time_ps: float) -> str:
    return f"{time_ps:.6f}".rstrip("0").rstrip(".").replace(".", "p") + "ps"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refinement", type=int, choices=(1, 2), default=1)
    parser.add_argument("--restart-from", type=Path, default=None)
    parser.add_argument("--max-time-ps", type=float, default=1.0)
    parser.add_argument("--max-temperature-K", type=float, default=2300.0)
    parser.add_argument("--pool-followup-ps", type=float, default=0.02)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--restore-only", action="store_true", help="Validate a checkpoint and exit")
    args = parser.parse_args()
    if args.max_time_ps <= 0 or args.max_temperature_K <= 0 or args.pool_followup_ps < 0:
        parser.error("time/temperature limits must be positive and pool follow-up nonnegative")
    if not GEOMETRY.exists():
        parser.error(f"Missing target 70 vol% geometry: {GEOMETRY}")
    initial = INITIAL_1X if args.refinement == 1 else INITIAL_2X
    restart = (args.restart_from or initial).resolve()
    if not restart.exists():
        parser.error(f"Missing restart checkpoint: {restart}")
    output_dir = (args.output_dir or RUN_DIR / f"ag70vol_resume_{args.refinement}x").resolve()
    checkpoint_dir = output_dir / "checkpoints"
    segment_dir = output_dir / "segments"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    segment_dir.mkdir(parents=True, exist_ok=True)

    model = QSGSEnthalpyGrid(
        None, "A", native_geometry_npz=GEOMETRY, source_geometry_npz=GEOMETRY,
        mesh_refinement=args.refinement,
    )
    tracker = LiquidClusterTracker(model.spacing_A, model.size_A, liquid_threshold=0.5)
    restored = model.load_checkpoint(restart, tracker)
    if model.time_s * 1e12 >= args.max_time_ps:
        raise ValueError("Restart checkpoint is already at or beyond --max-time-ps")
    if args.restore_only:
        print(json.dumps({
            "status": "CHECKPOINT_RESTORE_VALIDATED",
            "checkpoint": str(restart),
            "time_ps": model.time_s * 1e12,
            "peak_temperature_K": float(model.state.temperature_K.max()),
            "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
            "injected_energy_J": model.injected_energy_J,
            "energy_balance_error_J": restored["energy_balance_error_J"],
            "tracker_previous_time_ps": tracker.previous_time_ps,
            "tracker_history_rows": len(tracker.history),
        }, indent=2))
        return

    start_time = time.perf_counter()
    stage_number = 0
    first_pool_ps = tracker.first_multivoxel_pool_time_ps
    maximum_observed_K = float(model.state.temperature_K.max())
    stopped_reason = "maximum_time_reached"
    while model.time_s * 1e12 < args.max_time_ps - 1e-12:
        current_ps = model.time_s * 1e12
        peak_K = float(model.state.temperature_K.max())
        maximum_observed_K = max(maximum_observed_K, peak_K)
        if not np.isfinite(peak_K) or not np.all(np.isfinite(model.enthalpy_density_J_m3)):
            raise FloatingPointError("Non-finite state detected before heating segment")
        if np.any(model.state.ag_liquid_fraction < 0.0) or np.any(model.state.ag_liquid_fraction > 1.0):
            raise FloatingPointError("Ag liquid fraction left [0, 1]")
        if peak_K >= args.max_temperature_K:
            stopped_reason = "temperature_protection_reached"
            break

        # Coarse intervals are capped at 0.05 ps. Near Tm they contract to
        # 0.02 ps and 0.01 ps so onset and connected-pool time are bracketed.
        if peak_K < model.melting_temperature_K - 180.0:
            interval_ps = 0.05
        elif peak_K < model.melting_temperature_K - 50.0:
            interval_ps = 0.02
        else:
            interval_ps = 0.01
        if peak_K >= args.max_temperature_K - 50.0:
            interval_ps = min(interval_ps, 0.002)
        end_ps = min(current_ps + interval_ps, args.max_time_ps)
        if first_pool_ps is not None:
            end_ps = min(end_ps, first_pool_ps + args.pool_followup_ps)
            if current_ps >= first_pool_ps + args.pool_followup_ps - 1e-10:
                stopped_reason = "multi_voxel_pool_followup_completed"
                break
        if end_ps <= current_ps + 1e-12:
            stopped_reason = "target_time_reached"
            break

        stage_number += 1
        end_tag = tag(end_ps)
        stage_rows = model.run(
            end_ps - current_ps,
            source_on=True,
            sample_interval_ps=min(0.01, end_ps - current_ps),
            track_liquid_connectivity=True,
            liquid_threshold=0.5,
            progress=True,
            tracker=tracker,
            checkpoint_times_ps=(end_ps,),
            checkpoint_dir=checkpoint_dir,
            checkpoint_prefix=f"caseA_ag70_{args.refinement}x_heat",
        )
        peak_K = float(model.state.temperature_K.max())
        maximum_observed_K = max(maximum_observed_K, peak_K)
        total_energy = float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3)
        energy_error = total_energy - model.injected_energy_J
        relative_error = energy_error / model.injected_energy_J if model.injected_energy_J else 0.0
        if not np.isfinite(peak_K) or not np.isfinite(relative_error):
            raise FloatingPointError("Non-finite state or energy balance after heating segment")
        if np.any(model.state.ag_liquid_fraction < 0.0) or np.any(model.state.ag_liquid_fraction > 1.0):
            raise FloatingPointError("Ag liquid fraction left [0, 1]")
        if abs(relative_error) > 1.0e-8:
            raise ArithmeticError(f"Energy balance relative error {relative_error:.3e} exceeds 1e-8")
        checkpoint = checkpoint_dir / f"caseA_ag70_{args.refinement}x_heat_{end_tag}.npz"
        if not checkpoint.exists():
            raise FileNotFoundError(f"Expected segment checkpoint was not written: {checkpoint}")
        # Every segment is independently auditable, even if the driver is stopped.
        write_csv(segment_dir / f"heating_{end_tag}.csv", stage_rows)
        segment = {
            "segment": stage_number,
            "start_time_ps": current_ps,
            "end_time_ps": end_ps,
            "interval_ps": end_ps - current_ps,
            "peak_temperature_K": peak_K,
            "max_ag_liquid_fraction": float(model.state.ag_liquid_fraction.max()),
            "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
            "connected_cluster_count": int(stage_rows[-1].get("connected_liquid_cluster_count", 0)),
            "largest_cluster_voxels": int(stage_rows[-1].get("largest_liquid_cluster_voxels", 0)),
            "largest_cluster_volume_A3": float(stage_rows[-1].get("largest_liquid_cluster_volume_A3", 0.0)),
            "first_liquid_time_ps": tracker.first_liquid_time_ps,
            "first_connected_liquid_time_ps": tracker.first_connected_liquid_time_ps,
            "first_multivoxel_pool_time_ps": tracker.first_multivoxel_pool_time_ps,
            "injected_energy_J": model.injected_energy_J,
            "total_enthalpy_J": total_energy,
            "energy_balance_error_J": energy_error,
            "relative_energy_balance_error": relative_error,
            "checkpoint": str(checkpoint.relative_to(output_dir)),
        }
        with (output_dir / "segments.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(segment, sort_keys=True) + "\n")
        print(json.dumps(segment, sort_keys=True), flush=True)
        if tracker.first_multivoxel_pool_time_ps is not None and first_pool_ps is None:
            first_pool_ps = tracker.first_multivoxel_pool_time_ps
        if peak_K >= args.max_temperature_K:
            stopped_reason = "temperature_protection_reached"
            break
        if first_pool_ps is not None and model.time_s * 1e12 >= first_pool_ps + args.pool_followup_ps - 1e-10:
            stopped_reason = "multi_voxel_pool_followup_completed"
            break

    elapsed_s = time.perf_counter() - start_time
    # Persist full tracker lineage and a final restart point at the actual stop time.
    final_ps = model.time_s * 1e12
    final_checkpoint = checkpoint_dir / f"caseA_ag70_{args.refinement}x_final_{tag(final_ps)}.npz"
    model.save_checkpoint(final_checkpoint, tracker)
    if tracker.history:
        write_cluster_csv(output_dir / "heating.clusters.csv", tracker.history)
    summary = {
        "status": "AG70VOL_CASE_A_HEATING_CONTINUED",
        "refinement": args.refinement,
        "geometry_sha256": QSGSEnthalpyGrid._sha256_file(GEOMETRY),
        "restart_checkpoint": str(restart),
        "final_checkpoint": str(final_checkpoint.relative_to(output_dir)),
        "boundary_case": "A: zero Ti3SiC2 heat storage/conductivity; insulated external faces",
        "source_on": True,
        "time_start_ps": restored["time_ps"],
        "time_end_ps": final_ps,
        "max_temperature_K": float(model.state.temperature_K.max()),
        "maximum_temperature_seen_K": maximum_observed_K,
        "melting_temperature_reference_K": model.melting_temperature_K,
        "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
        "max_ag_liquid_fraction": float(model.state.ag_liquid_fraction.max()),
        "first_liquid_time_ps": tracker.first_liquid_time_ps,
        "first_connected_liquid_time_ps": tracker.first_connected_liquid_time_ps,
        "first_multivoxel_pool_time_ps": tracker.first_multivoxel_pool_time_ps,
        "tracker_event_rows": len(tracker.history),
        "largest_track_continuous_duration_ps": max(
            (float(row["continuous_duration_ps"]) for row in tracker.history if row["is_main_cluster"]),
            default=0.0,
        ),
        "injected_energy_J": model.injected_energy_J,
        "total_enthalpy_J": float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3),
        "energy_balance_error_J": float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3 - model.injected_energy_J),
        "relative_energy_balance_error": float((model.enthalpy_density_J_m3.sum() * model.cell_volume_m3 - model.injected_energy_J) / model.injected_energy_J) if model.injected_energy_J else 0.0,
        "elapsed_wall_time_s": elapsed_s,
        "stopped_reason": stopped_reason,
        "safety_temperature_limit_K": args.max_temperature_K,
        "sampled_cluster_tracking_is_not_a_stability_proof": True,
        "physical_limitations": [
            "candidate thermophysical properties remain to be source-verified",
            "Case A has zero scaffold heat capacity/conductivity and insulated external faces",
            "vaporization is not modeled; the run is stopped at the conservative 2300 K numerical protection line",
        ],
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
