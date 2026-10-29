#!/usr/bin/env python3
"""Safely continue the Ag 70 vol% Case A melt-window run from a checkpoint.

The driver keeps each segment short, writes a restart after every segment, and
continues until a persistent, volume-stable liquid pool forms or the configured
temperature/time safety limit is reached.
"""

from __future__ import annotations

import argparse
import csv
import json
import signal
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

STOP_REQUESTED = False


def _handle_stop(signum: int, frame: object) -> None:
    global STOP_REQUESTED
    STOP_REQUESTED = True


signal.signal(signal.SIGTERM, _handle_stop)
signal.signal(signal.SIGINT, _handle_stop)


def tag(time_ps: float) -> str:
    return f"{time_ps:.6f}".rstrip("0").rstrip(".").replace(".", "p") + "ps"


def assess_stable_pool(
    tracker: LiquidClusterTracker,
    end_ps: float,
    minimum_track_ps: float,
    window_ps: float,
    maximum_range_fraction: float,
) -> dict[str, object]:
    """Require a persistent main liquid track and a volume plateau."""
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
    max_gap = max((right - left for left, right in zip(times, times[1:])), default=0.0)
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
        "volume_range_fraction": relative_range,
        "maximum_range_fraction": maximum_range_fraction,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refinement", type=int, choices=(1, 2), default=1)
    parser.add_argument("--restart-from", type=Path, default=None)
    parser.add_argument(
        "--max-time-ps",
        type=float,
        default=None,
        help="optional tranche endpoint; omit to continue until temperature protection or stable-pool criteria",
    )
    parser.add_argument("--max-temperature-K", type=float, default=2300.0)
    parser.add_argument("--stable-track-ps", type=float, default=3.0)
    parser.add_argument("--stability-window-ps", type=float, default=2.0)
    parser.add_argument("--stable-volume-range", type=float, default=0.05)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--restore-only", action="store_true", help="Validate a checkpoint and exit")
    args = parser.parse_args()
    if (
        (args.max_time_ps is not None and args.max_time_ps <= 0)
        or args.max_temperature_K <= 0
        or args.stable_track_ps <= 0
        or args.stability_window_ps <= 0
        or not 0 < args.stable_volume_range <= 1
    ):
        parser.error("time, temperature, and stability limits must be positive; stable volume range must be in (0, 1]")
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
    if args.max_time_ps is not None and model.time_s * 1e12 >= args.max_time_ps:
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
    segment_log = output_dir / "segments.jsonl"
    stage_number = sum(1 for _ in segment_log.open(encoding="utf-8")) if segment_log.exists() else 0
    maximum_observed_K = float(model.state.temperature_K.max())
    stability_assessment: dict[str, object] = {"stable": False, "reason": "stable pool not yet detected"}
    stopped_reason = "maximum_time_reached" if args.max_time_ps is not None else "temperature_protection_or_stable_pool"
    while args.max_time_ps is None or model.time_s * 1e12 < args.max_time_ps - 1e-12:
        if STOP_REQUESTED:
            stopped_reason = "safe_stop_requested"
            break
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
        end_ps = current_ps + interval_ps
        if args.max_time_ps is not None:
            end_ps = min(end_ps, args.max_time_ps)
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
        stability_assessment = assess_stable_pool(
            tracker,
            model.time_s * 1e12,
            args.stable_track_ps,
            args.stability_window_ps,
            args.stable_volume_range,
        )
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
            "stable_pool_assessment": stability_assessment,
            "injected_energy_J": model.injected_energy_J,
            "total_enthalpy_J": total_energy,
            "energy_balance_error_J": energy_error,
            "relative_energy_balance_error": relative_error,
            "checkpoint": str(checkpoint.relative_to(output_dir)),
        }
        with (output_dir / "segments.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(segment, sort_keys=True) + "\n")
        print(json.dumps(segment, sort_keys=True), flush=True)
        if peak_K >= args.max_temperature_K:
            stopped_reason = "temperature_protection_reached"
            break
        if stability_assessment.get("stable"):
            stopped_reason = "stable_pool_thresholds_met"
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
        "time_limit_ps": args.max_time_ps,
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
        "stable_pool_thresholds": {
            "minimum_track_duration_ps": args.stable_track_ps,
            "volume_window_ps": args.stability_window_ps,
            "maximum_volume_range_fraction": args.stable_volume_range,
        },
        "stable_pool_assessment": stability_assessment,
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
