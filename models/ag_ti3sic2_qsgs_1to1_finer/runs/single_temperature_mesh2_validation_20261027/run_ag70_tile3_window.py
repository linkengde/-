#!/usr/bin/env python3
"""Optional 70:30 3x3 Case A central-source short run from uniform 300 K.

Run only after the source/energy preflight and only when the remaining wall
budget covers this job plus the required final reporting reserve.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
PHASE = PROJECT / "phase_field_qsgs"
GEOMETRY = HERE / "qsgs_morphology_ag70vol.npz"
sys.path.insert(0, str(PHASE))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid, write_cluster_csv, write_csv  # noqa: E402
from qsgs_liquid_connectivity import liquid_cluster_labels  # noqa: E402
from qsgs_liquid_tracking import LiquidClusterTracker  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration-ps", type=float, default=0.2)
    parser.add_argument("--sample-ps", type=float, default=0.01)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.duration_ps <= 0 or args.sample_ps <= 0:
        parser.error("duration and sampling interval must be positive")
    out = args.output_dir.resolve()
    checkpoint_dir = out / "checkpoints"
    out.mkdir(parents=True, exist_ok=True)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        parent_ag = ~np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
    base_shape = np.asarray(parent_ag.shape)
    model = QSGSEnthalpyGrid(
        None,
        "A",
        native_geometry_npz=GEOMETRY,
        source_geometry_npz=GEOMETRY,
        tile_xy=3,
        mesh_refinement=1,
    )
    if model.source_center_xy_A != tuple((model.size_A[:2] / 2.0).tolist()):
        raise AssertionError("The source is not centered in the middle repeated unit cell")
    tracker = LiquidClusterTracker(model.spacing_A, model.size_A, liquid_threshold=0.5)
    started = time.perf_counter()
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
        checkpoint_prefix="caseA_ag70_tile3_heat",
    )
    elapsed = time.perf_counter() - started
    write_csv(out / "heating.csv", rows)
    if tracker.history:
        write_cluster_csv(out / "heating.clusters.csv", tracker.history)
    else:
        (out / "heating.clusters.csv").write_text(
            "time_ps,track_id,parent_track_ids,event,is_main_cluster,continuous_duration_ps,component_id,voxel_count,volume_A3,centroid_A,depth_from_top_A,z_min_A,z_max_A,lateral_size_x_A,lateral_size_y_A\n",
            encoding="utf-8",
        )

    n = base_shape
    center = (slice(int(n[0]), int(2 * n[0])), slice(int(n[1]), int(2 * n[1])), slice(None))
    center_ag = model.ag_fraction[center] > 0.0
    center_temp = model.state.temperature_K[center]
    center_liquid = model.state.ag_liquid_fraction[center]
    center_peak_flat = np.where(center_ag, center_temp, -np.inf).argmax()
    center_peak_index = np.unravel_index(center_peak_flat, center_temp.shape)
    global_mask = (model.ag_fraction > 0.0) & (model.state.ag_liquid_fraction >= 0.5)
    labels, global_count = liquid_cluster_labels(global_mask)
    center_labels = np.unique(labels[center])
    center_labels = center_labels[center_labels > 0]
    cluster_sizes = np.bincount(labels.ravel())
    center_overlapping_clusters = [int(label) for label in center_labels]
    largest_center_overlapping_cluster = max(
        center_overlapping_clusters, key=lambda label: int(cluster_sizes[label]), default=None
    )
    global_volume_A3 = float(np.prod(model.size_A))
    global_voxel_volume_A3 = float(np.prod(model.spacing_A))
    global_energy = float(model.enthalpy_density_J_m3.sum() * model.cell_volume_m3)
    energy_error = global_energy - model.injected_energy_J
    energy_relative_error = energy_error / model.injected_energy_J if model.injected_energy_J else 0.0
    if not np.all(np.isfinite(model.state.temperature_K)) or not np.all(np.isfinite(model.state.ag_liquid_fraction)):
        raise FloatingPointError("Non-finite 3x3 field detected")
    if np.any(model.state.ag_liquid_fraction < 0.0) or np.any(model.state.ag_liquid_fraction > 1.0):
        raise FloatingPointError("3x3 Ag liquid fraction left [0, 1]")
    if abs(energy_relative_error) > 1.0e-8:
        raise ArithmeticError(f"3x3 energy-balance relative error {energy_relative_error:.3e} exceeds 1e-8")
    checkpoint = checkpoint_dir / f"caseA_ag70_tile3_heat_{args.duration_ps:g}ps.npz"
    if not checkpoint.exists():
        candidates = list(checkpoint_dir.glob("*.npz"))
        if len(candidates) != 1:
            raise FileNotFoundError("Expected exactly one final 3x3 checkpoint")
        checkpoint = candidates[0]
    checkpoint_digest = model._sha256_array(model.enthalpy_density_J_m3)
    restored = QSGSEnthalpyGrid(
        None,
        "A",
        native_geometry_npz=GEOMETRY,
        source_geometry_npz=GEOMETRY,
        tile_xy=3,
        mesh_refinement=1,
    )
    restored_tracker = LiquidClusterTracker(restored.spacing_A, restored.size_A, liquid_threshold=0.5)
    checkpoint_info = restored.load_checkpoint(checkpoint, restored_tracker)
    restart_ok = restored._sha256_array(restored.enthalpy_density_J_m3) == checkpoint_digest
    if not restart_ok:
        raise AssertionError("3x3 final checkpoint did not round-trip")

    report = {
        "status": "AG70_TILE3_CASE_A_HEATING_COMPLETED",
        "geometry_sha256": QSGSEnthalpyGrid._sha256_file(GEOMETRY),
        "composition_target": "Ag:Ti3SiC2 = 70:30 by volume",
        "boundary_case": "A: zero Ti3SiC2 heat storage/conductivity; insulated external faces",
        "grid_shape_xyz": list(model.temperature_shape),
        "base_cell_shape_xyz": list(base_shape),
        "tile_xy": 3,
        "source_count": 1,
        "source_center_xy_A": list(model.source_center_xy_A),
        "base_cell_size_A": size_A.tolist(),
        "source_power_W": model.actual_source_power_W,
        "time_ps": model.time_s * 1.0e12,
        "peak_temperature_global_K": float(model.state.temperature_K.max()),
        "peak_temperature_center_cell_K": float(center_temp[center_peak_index]),
        "peak_position_center_cell_local_A": ((np.asarray(center_peak_index) + 0.5) * model.spacing_A).tolist(),
        "mean_ag_temperature_center_cell_K": float(np.sum(center_temp * center_ag) / np.sum(center_ag)),
        "ag_liquid_volume_fraction_center_cell": float(np.sum(center_liquid * center_ag) / np.sum(center_ag)),
        "max_ag_liquid_fraction_center_cell": float(center_liquid.max()),
        "connected_cluster_count_global_at_threshold_0p5": int(global_count),
        "connected_clusters_overlapping_center_cell": center_overlapping_clusters,
        "largest_center_overlapping_cluster_volume_A3": (
            None if largest_center_overlapping_cluster is None
            else float(cluster_sizes[largest_center_overlapping_cluster] * global_voxel_volume_A3)
        ),
        "global_ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
        "injected_energy_J": model.injected_energy_J,
        "total_enthalpy_J": global_energy,
        "energy_balance_error_J": energy_error,
        "relative_energy_balance_error": energy_relative_error,
        "elapsed_wall_time_s": elapsed,
        "checkpoint": str(checkpoint.relative_to(out)),
        "checkpoint_restart_round_trip_passed": restart_ok,
        "checkpoint_reload_energy_error_J": checkpoint_info["energy_balance_error_J"],
        "interpretation": "The central unit cell has one source; outer repeated cells are unheated. Only central-cell local metrics should be compared to the single-cell calculation.",
    }
    (out / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
