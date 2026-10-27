#!/usr/bin/env python3
"""Compare a same-time native 1x checkpoint and strict-subdivision 2x state."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT / "phase_field_qsgs"))
from qsgs_liquid_connectivity import liquid_cluster_labels  # noqa: E402


def read(path: Path) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: np.asarray(archive[name]).copy() for name in (
            "temperature_K", "ag_liquid_fraction", "enthalpy_density_J_m3"
        )}
        metadata = json.loads(str(archive["metadata_json"].item()))
    return arrays, metadata


def liquid_metrics(liquid: np.ndarray, ag: np.ndarray, voxel_volume_A3: float) -> dict[str, object]:
    liquid_any = ag & (liquid > 1.0e-12)
    liquid_mask = ag & (liquid >= 0.5)
    labels, count = liquid_cluster_labels(liquid_mask)
    sizes = np.bincount(labels.ravel())[1:]
    max_voxels = int(sizes.max()) if sizes.size else 0
    return {
        "ag_liquid_volume_fraction": float(np.sum(liquid * ag) / ag.sum()),
        "ag_volume_fraction_any_liquid": float(liquid_any.sum() / ag.sum()),
        "max_liquid_fraction": float(liquid.max()),
        "liquid_voxels_at_threshold_0p5": int(liquid_mask.sum()),
        "connected_cluster_count_at_threshold_0p5": int(count),
        "largest_cluster_voxels": max_voxels,
        "largest_cluster_volume_A3": float(max_voxels * voxel_volume_A3),
    }


def checkpoint_main_cluster(metadata: dict[str, object]) -> dict[str, object] | None:
    tracker = metadata.get("tracker_metadata")
    if not isinstance(tracker, dict):
        return None
    target = float(metadata["time_ps"])
    rows = [row for row in tracker.get("history", [])
            if abs(float(row["time_ps"]) - target) <= 1.0e-9 and bool(row["is_main_cluster"])]
    if not rows:
        return None
    return max(rows, key=lambda row: int(row["voxel_count"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint_1x", type=Path)
    parser.add_argument("checkpoint_2x", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    one, meta1 = read(args.checkpoint_1x)
    two, meta2 = read(args.checkpoint_2x)
    t1 = float(meta1["time_ps"]); t2 = float(meta2["time_ps"])
    if not np.isclose(t1, t2, rtol=0.0, atol=1.0e-9):
        raise ValueError(f"Checkpoint times differ: {t1} vs {t2} ps")
    t2c = two["temperature_K"].reshape(one["temperature_K"].shape[0], 2, one["temperature_K"].shape[1], 2, one["temperature_K"].shape[2], 2).mean(axis=(1, 3, 5))
    l2c = two["ag_liquid_fraction"].reshape(one["temperature_K"].shape[0], 2, one["temperature_K"].shape[1], 2, one["temperature_K"].shape[2], 2).mean(axis=(1, 3, 5))
    with np.load(HERE / "qsgs_morphology_ag70vol.npz", allow_pickle=False) as archive:
        ag = ~np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
    if two["temperature_K"].shape != tuple(np.asarray(ag.shape) * 2):
        raise ValueError("2x field is not a strict subdivision of the 1x Ag70 QSGS grid")
    ag2 = np.repeat(np.repeat(np.repeat(ag, 2, axis=0), 2, axis=1), 2, axis=2)
    mask1 = ag & (one["ag_liquid_fraction"] >= 0.5)
    mask2_native = ag2 & (two["ag_liquid_fraction"] >= 0.5)
    mask2_coarse = mask2_native.reshape(ag.shape[0], 2, ag.shape[1], 2, ag.shape[2], 2).mean(axis=(1, 3, 5)) >= 0.5
    union = int(np.count_nonzero(mask1 | mask2_coarse))
    intersection = int(np.count_nonzero(mask1 & mask2_coarse))
    diff = t2c - one["temperature_K"]
    ag_weights = ag.astype(float)
    mae = float(np.sum(np.abs(diff) * ag_weights) / np.sum(ag_weights))
    rmse = float(np.sqrt(np.sum(diff**2 * ag_weights) / np.sum(ag_weights)))
    rise = float(np.sum((one["temperature_K"] - 300.0) * ag_weights) / np.sum(ag_weights))
    idx1 = np.unravel_index(np.argmax(one["temperature_K"]), one["temperature_K"].shape)
    idx2 = np.unravel_index(np.argmax(two["temperature_K"]), two["temperature_K"].shape)
    spacing1 = size_A / np.asarray(ag.shape)
    spacing2 = spacing1 / 2.0
    result = {
        "status": "AG70_MATCHED_TIME_1X_2X_COMPARISON",
        "time_ps": t1,
        "checkpoint_1x": str(args.checkpoint_1x),
        "checkpoint_2x": str(args.checkpoint_2x),
        "peak_temperature_1x_K": float(one["temperature_K"].max()),
        "peak_temperature_2x_K": float(two["temperature_K"].max()),
        "peak_temperature_difference_2x_minus_1x_K": float(two["temperature_K"].max() - one["temperature_K"].max()),
        "peak_position_1x_A": ((np.asarray(idx1) + 0.5) * spacing1).tolist(),
        "peak_position_2x_A": ((np.asarray(idx2) + 0.5) * spacing2).tolist(),
        "peak_lateral_offset_from_source_1x_A": float(np.linalg.norm(((np.asarray(idx1[:2]) + 0.5) * spacing1[:2]) - size_A[:2] / 2.0)),
        "peak_lateral_offset_from_source_2x_A": float(np.linalg.norm(((np.asarray(idx2[:2]) + 0.5) * spacing2[:2]) - size_A[:2] / 2.0)),
        "block_averaged_2x_minus_1x_ag_temperature_K": {
            "mean_signed": float(np.sum(diff * ag_weights) / np.sum(ag_weights)),
            "mean_absolute": mae,
            "rmse": rmse,
            "maximum_absolute": float(np.max(np.abs(diff[ag]))),
            "l1_over_1x_ag_mean_temperature_rise": float(np.sum(np.abs(diff) * ag_weights) / max(np.sum((one["temperature_K"] - 300.0) * ag_weights), 1.0e-300)),
        },
        "liquid_1x": liquid_metrics(one["ag_liquid_fraction"], ag, float(np.prod(spacing1))),
        "liquid_2x_native": liquid_metrics(two["ag_liquid_fraction"], ag2, float(np.prod(spacing2))),
        "coarsened_liquid_mask_iou_threshold_0p5": None if union == 0 else intersection / union,
        "coarsened_liquid_mask_intersection_voxels": intersection,
        "coarsened_liquid_mask_union_voxels": union,
        "largest_tracked_cluster_1x": checkpoint_main_cluster(meta1),
        "largest_tracked_cluster_2x": checkpoint_main_cluster(meta2),
        "injected_energy_1x_J": float(meta1["injected_energy_J"]),
        "injected_energy_2x_J": float(meta2["injected_energy_J"]),
        "total_enthalpy_1x_J": float(meta1["total_enthalpy_J"]),
        "total_enthalpy_2x_J": float(meta2["total_enthalpy_J"]),
        "energy_balance_error_1x_J": float(meta1["energy_balance_error_J"]),
        "energy_balance_error_2x_J": float(meta2["energy_balance_error_J"]),
        "interpretation": "Mesh-sensitivity comparison at this time only. A single 1x/2x pair is not a mesh-convergence proof.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
