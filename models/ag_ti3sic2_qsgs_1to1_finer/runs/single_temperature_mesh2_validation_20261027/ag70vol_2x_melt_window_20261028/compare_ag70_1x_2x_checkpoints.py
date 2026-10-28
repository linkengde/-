#!/usr/bin/env python3
"""Compare matched 1x/2x checkpoint fields and physical liquid-cluster metrics."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RUN_ROOT = HERE.parent
ONE_X_DIR = RUN_ROOT / "ag70vol_resume_1x_from1ps_20261028" / "checkpoints"
GEOMETRY = RUN_ROOT / "qsgs_morphology_ag70vol.npz"
DEFAULT_TWO_X_DIR = HERE / "heating_0p36_to_pool" / "checkpoints"
FIELDS = ("temperature_K", "enthalpy_density_J_m3", "ag_liquid_fraction", "tracker_track_labels")


def tag(time_ps: float) -> str:
    return f"{time_ps:.6f}".rstrip("0").rstrip(".").replace(".", "p") + "ps"


def checkpoint(directory: Path, refinement: int, time_ps: float) -> Path:
    name = f"caseA_ag70_{refinement}x_heat_{tag(time_ps)}.npz"
    path = directory / name
    if not path.exists() and refinement == 1 and abs(time_ps - 4.52) < 1e-8:
        path = directory / f"caseA_ag70_1x_final_{tag(time_ps)}.npz"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def read_checkpoint(path: Path) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: np.asarray(archive[name]) for name in FIELDS}
        metadata = json.loads(str(archive["metadata_json"].item()))
    return arrays, metadata


def metrics(
    arrays: dict[str, np.ndarray],
    metadata: dict[str, object],
    ag_mask: np.ndarray,
    size_A: np.ndarray,
    time_ps: float,
) -> dict[str, object]:
    temperature = arrays["temperature_K"]
    liquid = arrays["ag_liquid_fraction"]
    factor = tuple(int(n // coarse) for n, coarse in zip(temperature.shape, ag_mask.shape, strict=True))
    if factor == (1, 1, 1):
        projected_temperature = temperature
    else:
        projected_temperature = temperature.reshape(
            ag_mask.shape[0], factor[0], ag_mask.shape[1], factor[1], ag_mask.shape[2], factor[2]
        ).mean(axis=(1, 3, 5))
    ag_cells = int(ag_mask.sum())
    spacing_A = size_A / np.asarray(ag_mask.shape, dtype=float)
    cell_volume_A3 = float(np.prod(spacing_A))
    fine_ag_mask = ag_mask if factor == (1, 1, 1) else np.repeat(np.repeat(np.repeat(ag_mask, factor[0], 0), factor[1], 1), factor[2], 2)
    ag_liquid_fraction = float(np.sum(liquid * fine_ag_mask) / (ag_cells * np.prod(factor)))
    from qsgs_liquid_tracking import LiquidClusterTracker  # imported after CLI adds project path

    fine_spacing_A = size_A / np.asarray(liquid.shape, dtype=float)
    snapshot_tracker = LiquidClusterTracker(fine_spacing_A, size_A, liquid_threshold=0.5)
    cluster_summary, cluster_records = snapshot_tracker.update(time_ps, liquid, fine_ag_mask)
    largest_record = max(cluster_records, key=lambda row: int(row["voxel_count"]), default={})
    tracker_metadata = metadata.get("tracker_metadata") or {}
    stored_main_rows = [
        row
        for row in tracker_metadata.get("history", [])
        if abs(float(row.get("time_ps", -1.0)) - time_ps) <= 1.0e-8 and row.get("is_main_cluster")
    ]
    stored_main = max(stored_main_rows, key=lambda row: int(row.get("voxel_count", 0)), default={})
    largest_voxels = int(cluster_summary["largest_liquid_cluster_voxels"])
    largest_volume_A3 = largest_voxels * cell_volume_A3 / np.prod(factor)
    return {
        "time_ps": time_ps,
        "checkpoint_time_ps": float(metadata["time_s"]) * 1.0e12,
        "peak_temperature_K": float(temperature.max()),
        "ag_peak_temperature_K": float(temperature[fine_ag_mask].max()),
        "mean_ag_temperature_K": float(temperature[fine_ag_mask].mean()),
        "block_averaged_temperature_rmse_vs_1x_K": None,
        "ag_liquid_volume_fraction": ag_liquid_fraction,
        "total_ag_liquid_volume_A3": float(np.sum(liquid * fine_ag_mask) * cell_volume_A3 / np.prod(factor)),
        "connected_cluster_count": int(cluster_summary["connected_liquid_cluster_count"]),
        "liquid_voxel_count": int(cluster_summary["liquid_voxel_count"]),
        "largest_cluster_voxels": largest_voxels,
        "largest_cluster_volume_A3": float(largest_volume_A3),
        "largest_cluster_centroid_A": largest_record.get("centroid_A"),
        "largest_cluster_depth_from_top_A": largest_record.get("depth_from_top_A"),
        "largest_cluster_lateral_size_x_A": largest_record.get("lateral_size_x_A"),
        "largest_cluster_lateral_size_y_A": largest_record.get("lateral_size_y_A"),
        "largest_cluster_track_id": stored_main.get("track_id"),
        "largest_cluster_track_duration_ps": stored_main.get("continuous_duration_ps"),
        "relative_energy_balance_error": float(metadata["energy_balance_error_J"]) / float(metadata["injected_energy_J"]),
        "projected_temperature_K": projected_temperature,
        "ag_mask": ag_mask,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--one-x-checkpoint-dir", type=Path, default=ONE_X_DIR)
    parser.add_argument("--two-x-checkpoint-dir", type=Path, default=DEFAULT_TWO_X_DIR)
    parser.add_argument("--times-ps", type=float, nargs="+", default=(4.12, 4.22, 4.52))
    parser.add_argument("--output", type=Path, default=HERE / "MESH_COMPARISON_RESULTS.csv")
    parser.add_argument("--json-output", type=Path, default=HERE / "MESH_COMPARISON_RESULTS.json")
    args = parser.parse_args()

    import sys

    sys.path.insert(0, str(RUN_ROOT.parent.parent / "phase_field_qsgs"))
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        ag_mask = ~np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)

    output_rows = []
    for time_ps in args.times_ps:
        one_arrays, one_meta = read_checkpoint(checkpoint(args.one_x_checkpoint_dir, 1, time_ps))
        two_arrays, two_meta = read_checkpoint(checkpoint(args.two_x_checkpoint_dir, 2, time_ps))
        one = metrics(one_arrays, one_meta, ag_mask, size_A, time_ps)
        two = metrics(two_arrays, two_meta, ag_mask, size_A, time_ps)
        temperature_difference = two["projected_temperature_K"] - one_arrays["temperature_K"]
        temperature_difference_ag = temperature_difference[ag_mask]
        two["block_averaged_temperature_rmse_vs_1x_K"] = float(
            np.sqrt(np.mean(temperature_difference_ag**2))
        )
        row = {
            "time_ps": time_ps,
            "one_x_peak_temperature_K": one["peak_temperature_K"],
            "two_x_peak_temperature_K": two["peak_temperature_K"],
            "one_x_ag_peak_temperature_K": one["ag_peak_temperature_K"],
            "two_x_ag_peak_temperature_K": two["ag_peak_temperature_K"],
            "ag_peak_temperature_difference_K": two["ag_peak_temperature_K"] - one["ag_peak_temperature_K"],
            "one_x_ag_liquid_volume_fraction": one["ag_liquid_volume_fraction"],
            "two_x_ag_liquid_volume_fraction": two["ag_liquid_volume_fraction"],
            "one_x_total_ag_liquid_volume_A3": one["total_ag_liquid_volume_A3"],
            "two_x_total_ag_liquid_volume_A3": two["total_ag_liquid_volume_A3"],
            "one_x_largest_cluster_voxels": one["largest_cluster_voxels"],
            "two_x_largest_cluster_voxels": two["largest_cluster_voxels"],
            "one_x_largest_cluster_volume_A3": one["largest_cluster_volume_A3"],
            "two_x_largest_cluster_volume_A3": two["largest_cluster_volume_A3"],
            "one_x_largest_cluster_centroid_A": one["largest_cluster_centroid_A"],
            "two_x_largest_cluster_centroid_A": two["largest_cluster_centroid_A"],
            "one_x_largest_cluster_depth_from_top_A": one["largest_cluster_depth_from_top_A"],
            "two_x_largest_cluster_depth_from_top_A": two["largest_cluster_depth_from_top_A"],
            "one_x_largest_cluster_lateral_size_x_A": one["largest_cluster_lateral_size_x_A"],
            "two_x_largest_cluster_lateral_size_x_A": two["largest_cluster_lateral_size_x_A"],
            "one_x_largest_cluster_lateral_size_y_A": one["largest_cluster_lateral_size_y_A"],
            "two_x_largest_cluster_lateral_size_y_A": two["largest_cluster_lateral_size_y_A"],
            "one_x_connected_cluster_count": one["connected_cluster_count"],
            "two_x_connected_cluster_count": two["connected_cluster_count"],
            "block_averaged_ag_temperature_rmse_K": two["block_averaged_temperature_rmse_vs_1x_K"],
            "one_x_relative_energy_residual": one["relative_energy_balance_error"],
            "two_x_relative_energy_residual": two["relative_energy_balance_error"],
        }
        output_rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    args.json_output.write_text(json.dumps(output_rows, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output_rows, indent=2))


if __name__ == "__main__":
    main()
