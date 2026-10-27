#!/usr/bin/env python3
"""Run the Case A heating window on the isolated Ag 70 vol% QSGS mask."""

from __future__ import annotations

import argparse
import gc
import json
import sys
import time
from pathlib import Path

import numpy as np


RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid, write_cluster_csv, write_csv  # noqa: E402
from qsgs_liquid_tracking import LiquidClusterTracker  # noqa: E402


def time_tag(time_ps: float) -> str:
    return f"{time_ps:.9f}".rstrip("0").rstrip(".").replace(".", "p") + "ps"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refinement", type=int, choices=(1, 2), required=True)
    parser.add_argument("--duration-ps", type=float, default=0.2)
    parser.add_argument("--sample-ps", type=float, default=0.002)
    args = parser.parse_args()
    if args.duration_ps <= 0.0 or args.sample_ps <= 0.0:
        parser.error("--duration-ps and --sample-ps must be positive")
    if not GEOMETRY.exists():
        parser.error(f"Missing target 70 vol% geometry: {GEOMETRY}")

    output_dir = RUN_DIR / "ag70vol_melt_window" / f"{args.refinement}x_{time_tag(args.duration_ps)}"
    checkpoint_dir = output_dir / "checkpoints"
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    with np.load(GEOMETRY, allow_pickle=False) as archive:
        base_shape = tuple(np.asarray(archive["skeleton"]).shape)
        size_A = np.asarray(archive["size_A"], dtype=float)
    model = QSGSEnthalpyGrid(
        None,
        "A",
        native_geometry_npz=GEOMETRY,
        source_geometry_npz=GEOMETRY,
        tile_xy=1,
        mesh_refinement=args.refinement,
    )
    tracker = LiquidClusterTracker(model.spacing_A, model.size_A, liquid_threshold=0.5)
    stable_dt_s = model.stable_timestep_s()
    start = time.perf_counter()
    rows = model.run(
        args.duration_ps,
        source_on=True,
        sample_interval_ps=args.sample_ps,
        track_liquid_connectivity=True,
        liquid_threshold=0.5,
        progress=True,
        tracker=tracker,
        checkpoint_times_ps=(args.duration_ps,),
        checkpoint_dir=checkpoint_dir,
        checkpoint_prefix=f"caseA_ag70_{time_tag(args.duration_ps)}",
    )
    elapsed = time.perf_counter() - start

    csv_path = output_dir / "heating.csv"
    clusters_path = output_dir / "heating.clusters.csv"
    write_csv(csv_path, rows)
    write_cluster_csv(clusters_path, tracker.history)
    if not clusters_path.exists():
        clusters_path.write_text(
            "time_ps,track_id,parent_track_ids,event,is_main_cluster,continuous_duration_ps,component_id,voxel_count,volume_A3,centroid_A,depth_from_top_A,z_min_A,z_max_A,lateral_size_x_A,lateral_size_y_A\n",
            encoding="utf-8",
        )
    checkpoint_path = checkpoint_dir / f"caseA_ag70_{time_tag(args.duration_ps)}_{time_tag(args.duration_ps)}.npz"
    # run() uses the checkpoint prefix and appends its target time; handle the
    # common case without relying on implementation-specific file naming.
    if not checkpoint_path.exists():
        candidates = sorted(checkpoint_dir.glob("*.npz"))
        if len(candidates) != 1:
            raise RuntimeError(f"Expected exactly one final checkpoint, found {len(candidates)}")
        checkpoint_path = candidates[0]

    final = rows[-1]
    initial_field_digest = model._sha256_array(model.enthalpy_density_J_m3)
    restored = model.load_checkpoint(checkpoint_path, tracker)
    restart_round_trip = bool(
        model._sha256_array(model.enthalpy_density_J_m3) == initial_field_digest
        and np.isclose(model.time_s * 1.0e12, args.duration_ps, rtol=0.0, atol=1.0e-10)
        and np.isclose(float(restored["injected_energy_J"]), float(final["injected_energy_J"]), rtol=1.0e-12, atol=1.0e-30)
    )
    if not restart_round_trip:
        raise AssertionError("Target-composition checkpoint did not round-trip to the final heating state")

    track_duration: dict[int, float] = {}
    for event in tracker.history:
        track_id = int(event["track_id"])
        track_duration[track_id] = max(track_duration.get(track_id, 0.0), float(event["continuous_duration_ps"]))
    main_track_duration = max(
        (float(event["continuous_duration_ps"]) for event in tracker.history if bool(event["is_main_cluster"])),
        default=0.0,
    )
    final_cluster_events = [
        event for event in tracker.history
        if np.isclose(float(event["time_ps"]), args.duration_ps, rtol=0.0, atol=1.0e-10)
        and bool(event["is_main_cluster"])
    ]
    final_main_cluster = max(final_cluster_events, key=lambda event: int(event["voxel_count"]), default=None)
    total_energy = float(final["total_energy_J"])
    energy_error = float(final["energy_balance_error_J"])
    report = {
        "status": "AG70VOL_CASE_A_HEATING_WINDOW_COMPLETED",
        "geometry_path": str(GEOMETRY.relative_to(PROJECT_DIR)),
        "geometry_sha256": QSGSEnthalpyGrid._sha256_file(GEOMETRY),
        "parent_geometry_sha256": json.loads((RUN_DIR / "composition_70vol_geometry.json").read_text())["construction"]["parent_geometry_sha256"],
        "composition_target": "Ag:Ti3SiC2 = 70:30 by volume",
        "boundary_case": "A: zero Ti3SiC2 heat storage/conductivity; insulated external faces",
        "initial_state": "uniform 300 K, zero enthalpy relative to 300 K",
        "refinement": args.refinement,
        "native_parent_grid_xyz": list(base_shape),
        "grid_shape_xyz": list(model.temperature_shape),
        "spacing_A": model.spacing_A.tolist(),
        "domain_size_A": size_A.tolist(),
        "source_integration_method": model.source_integration_method,
        "source_power_W": float(model.actual_source_power_W),
        "peak_flux_eV_ps_A2": float(model.peak_flux_W_m2 / (1.602176634e-19 / 1.0e-12 / 1.0e-20)),
        "gaussian_radius_A": model.gaussian_radius_A,
        "absorption_depth_A": model.absorption_depth_A,
        "heat_duration_ps": args.duration_ps,
        "sample_interval_ps": args.sample_ps,
        "stable_dt_fs_initial": stable_dt_s * 1.0e15,
        "approximate_steps": int(np.ceil(args.duration_ps * 1000.0 / (stable_dt_s * 1.0e15))),
        "wall_time_s": elapsed,
        "max_temperature_K": float(final["max_temperature_K"]),
        "max_temperature_seen_K": float(model.max_temperature_seen_K),
        "mean_ag_temperature_K": float(final["mean_ag_temperature_K"]),
        "ag_liquid_volume_fraction": float(final["ag_liquid_volume_fraction"]),
        "max_ag_liquid_fraction": float(final["max_ag_liquid_fraction"]),
        "connected_liquid_cluster_count_final": int(final.get("connected_liquid_cluster_count", 0)),
        "tracked_cluster_count": len(track_duration),
        "cluster_event_rows": len(tracker.history),
        "main_track_max_sampled_duration_ps": main_track_duration,
        "first_liquid_time_ps": tracker.first_liquid_time_ps,
        "first_connected_liquid_time_ps": tracker.first_connected_liquid_time_ps,
        "first_multivoxel_pool_time_ps": tracker.first_multivoxel_pool_time_ps,
        "largest_final_cluster": (
            {
                "track_id": int(final_main_cluster["track_id"]),
                "voxel_count": int(final_main_cluster["voxel_count"]),
                "volume_A3": float(final_main_cluster["volume_A3"]),
                "centroid_A": final_main_cluster["centroid_A"],
                "depth_from_top_A": float(final_main_cluster["depth_from_top_A"]),
                "lateral_size_x_A": float(final_main_cluster["lateral_size_x_A"]),
                "lateral_size_y_A": float(final_main_cluster["lateral_size_y_A"]),
            }
            if final_main_cluster is not None
            else None
        ),
        "injected_energy_J": float(final["injected_energy_J"]),
        "total_enthalpy_J": total_energy,
        "energy_balance_error_J": energy_error,
        "relative_energy_balance_error": energy_error / total_energy if total_energy else 0.0,
        "boiling_reference_K": 2435.0,
        "peak_exceeds_boiling_reference": bool(model.max_temperature_seen_K > 2435.0),
        "checkpoint": str(checkpoint_path.relative_to(RUN_DIR)),
        "checkpoint_restart_round_trip_passed": restart_round_trip,
        "timeseries_csv": str(csv_path.relative_to(RUN_DIR)),
        "cluster_history_csv": str(clusters_path.relative_to(RUN_DIR)),
        "limitations": [
            "This is a candidate-property exploratory run; thermal-property references are not yet archived.",
            "Boundary A treats the ceramic skeleton as a zero-heat-capacity, zero-conductivity adiabatic limit.",
            "External faces are insulated; this is internal redistribution, not environmental cooling.",
            "The model has no vaporization; temperatures above the Ag normal-boiling reference are outside its reliable regime.",
            "A single 1x/2x comparison is a mesh-sensitivity check, not proof of fully converged pool dimensions or lifetime.",
        ],
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    del model, tracker, rows
    gc.collect()
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
