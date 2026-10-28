#!/usr/bin/env python3
"""Compare matched Ag 60.7 wt% 1x/2x checkpoint fields and melt metrics."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
GEOMETRY = HERE / "nano_geometry_ag607wt_ti3sic2_393wt.npz"


def read_run(output: Path) -> tuple[dict[str, Any], dict[float, dict[str, Any]]]:
    summary_path = output / "summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing completed run summary: {summary_path}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    checkpoints: dict[float, dict[str, Any]] = {}
    for receipt_path in sorted((output / "checkpoints").glob("*.validation.json")):
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt.get("status") != "CHECKPOINT_RESTORE_VALIDATED":
            continue
        checkpoint = receipt_path.with_name(receipt_path.name.removesuffix(".validation.json") + ".npz")
        if not checkpoint.exists():
            continue
        key = round(float(receipt["time_ps"]), 9)
        checkpoints[key] = {
            "time_ps": float(receipt["time_ps"]),
            "phase": receipt["phase"],
            "checkpoint": checkpoint,
            "receipt": receipt_path,
        }
    if not checkpoints:
        raise RuntimeError(f"No validated checkpoints found in {output}")
    return summary, checkpoints


def cluster_at_checkpoint(checkpoint: Path, time_ps: float) -> dict[str, Any]:
    with np.load(checkpoint, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata_json"].item()))
    tracker = metadata.get("tracker_metadata") or {}
    records = [
        item for item in tracker.get("history", [])
        if item.get("is_main_cluster") and abs(float(item["time_ps"]) - time_ps) <= 1.0e-8
    ]
    if not records:
        return {"main_cluster_voxels": 0, "main_cluster_volume_A3": 0.0, "main_track_id": None, "track_duration_ps": 0.0}
    item = max(records, key=lambda row: int(row["voxel_count"]))
    return {
        "main_cluster_voxels": int(item["voxel_count"]),
        "main_cluster_volume_A3": float(item["volume_A3"]),
        "main_track_id": int(item["track_id"]),
        "track_duration_ps": float(item["continuous_duration_ps"]),
    }


def load_fields(checkpoint: Path, refinement: int) -> tuple[np.ndarray, np.ndarray]:
    with np.load(checkpoint, allow_pickle=False) as archive:
        temperature = np.asarray(archive["temperature_K"], dtype=np.float64).copy()
        liquid = np.asarray(archive["ag_liquid_fraction"], dtype=np.float64).copy()
    if refinement == 1:
        return temperature, liquid
    if refinement != 2:
        raise ValueError("Only 1x and 2x checkpoint grids are supported")
    expected = tuple(2 * value for value in (100, 100, 50))
    if temperature.shape != expected or liquid.shape != expected:
        raise ValueError(f"Unexpected 2x field shape: {temperature.shape}")
    return (
        temperature.reshape(100, 2, 100, 2, 50, 2).mean(axis=(1, 3, 5)),
        liquid.reshape(100, 2, 100, 2, 50, 2).mean(axis=(1, 3, 5)),
    )


def select_times(
    summary_1x: dict[str, Any],
    summary_2x: dict[str, Any],
    checkpoints_1x: dict[float, dict[str, Any]],
    checkpoints_2x: dict[float, dict[str, Any]],
) -> dict[float, list[str]]:
    common = sorted(set(checkpoints_1x) & set(checkpoints_2x))
    if not common:
        raise RuntimeError("No exact common checkpoint times exist between the two runs")
    selected: dict[float, list[str]] = {}

    def add_nearest(target: float, label: str, tolerance: float = 0.026) -> None:
        actual = min(common, key=lambda value: abs(value - target))
        if abs(actual - target) <= tolerance:
            selected.setdefault(actual, []).append(label)

    selected.setdefault(common[0], []).append("initial")
    selected.setdefault(common[-1], []).append("latest_shared")
    for time_ps in common:
        if abs(time_ps * 2.0 - round(time_ps * 2.0)) <= 1.0e-8:
            selected.setdefault(time_ps, []).append("periodic_0p5ps")

    for label, summary in (("1x", summary_1x), ("2x", summary_2x)):
        milestones = summary.get("milestones") or {}
        for key, value in milestones.items():
            if key.endswith("_checkpoint") or not key.endswith("_time_ps") or value is None:
                continue
            add_nearest(float(value), f"{label}_{key}")
        heat_off = summary.get("heat_off_time_ps")
        if heat_off is not None:
            add_nearest(float(heat_off), f"{label}_heat_off")
        final_time = summary.get("final_time_ps")
        if final_time is not None:
            add_nearest(float(final_time), f"{label}_final")
    return selected


def compare(
    output_1x: Path,
    output_2x: Path,
    report_dir: Path,
) -> dict[str, Any]:
    summary_1x, cp_1x = read_run(output_1x)
    summary_2x, cp_2x = read_run(output_2x)
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        ag_mask = ~np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
    voxel_volume_A3 = float(np.prod(size_A / np.asarray(ag_mask.shape)))
    ag_count = int(ag_mask.sum())
    selected = select_times(summary_1x, summary_2x, cp_1x, cp_2x)
    results: list[dict[str, Any]] = []
    for time_key in sorted(selected):
        left = cp_1x[time_key]
        right = cp_2x[time_key]
        temp_1, liquid_1 = load_fields(left["checkpoint"], 1)
        temp_2, liquid_2 = load_fields(right["checkpoint"], 2)
        if temp_1.shape != ag_mask.shape or temp_2.shape != ag_mask.shape:
            raise ValueError("Coarsened temperature fields do not match the native geometry")
        diff = temp_2 - temp_1
        ag_diff = diff[ag_mask]
        mean_1 = float(temp_1[ag_mask].mean())
        mean_2 = float(temp_2[ag_mask].mean())
        liquid_1_total = float(np.sum(liquid_1[ag_mask]))
        liquid_2_total = float(np.sum(liquid_2[ag_mask]))
        cluster_1 = cluster_at_checkpoint(left["checkpoint"], left["time_ps"])
        cluster_2 = cluster_at_checkpoint(right["checkpoint"], right["time_ps"])
        row = {
            "time_1x_ps": left["time_ps"],
            "time_2x_ps": right["time_ps"],
            "selected_labels": ";".join(sorted(set(selected[time_key]))),
            "phase_1x": left["phase"],
            "phase_2x": right["phase"],
            "phase_matches": left["phase"] == right["phase"],
            "peak_1x_K": float(temp_1.max()),
            "peak_2x_K": float(temp_2.max()),
            "peak_difference_2x_minus_1x_K": float(temp_2.max() - temp_1.max()),
            "ag_mean_1x_K": mean_1,
            "ag_mean_2x_K": mean_2,
            "ag_field_mae_K": float(np.mean(np.abs(ag_diff))),
            "ag_field_rmse_K": float(np.sqrt(np.mean(ag_diff**2))),
            "ag_field_max_abs_difference_K": float(np.max(np.abs(ag_diff))),
            "ag_temperature_rise_rmse_fraction_1x": float(
                np.sqrt(np.mean(ag_diff**2)) / max(mean_1 - 300.0, 1.0e-12)
            ),
            "ag_liquid_volume_fraction_1x": liquid_1_total / ag_count,
            "ag_liquid_volume_fraction_2x": liquid_2_total / ag_count,
            "liquid_volume_1x_A3": liquid_1_total * voxel_volume_A3,
            "liquid_volume_2x_A3": liquid_2_total * voxel_volume_A3,
            **{f"{key}_1x": value for key, value in cluster_1.items()},
            **{f"{key}_2x": value for key, value in cluster_2.items()},
            "checkpoint_1x": str(left["checkpoint"].relative_to(output_1x)),
            "checkpoint_2x": str(right["checkpoint"].relative_to(output_2x)),
        }
        results.append(row)

    report_dir.mkdir(parents=True, exist_ok=True)
    csv_path = report_dir / "mesh_comparison.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)
    report = {
        "status": "MATCHED_CHECKPOINT_FIELD_COMPARISON_COMPLETED",
        "composition": "Ag 60.7 wt% / Ti3SiC2 39.3 wt%; pore-free target 39.9585798/60.0414202 vol%",
        "geometry_sha256": summary_1x["geometry_sha256"],
        "refinement_1x": summary_1x["grid_shape_xyz"],
        "refinement_2x": summary_2x["grid_shape_xyz"],
        "common_checkpoint_count": len(set(cp_1x) & set(cp_2x)),
        "comparison_row_count": len(results),
        "comparison_csv": csv_path.name,
        "1x_summary_status": summary_1x["status"],
        "2x_summary_status": summary_2x["status"],
        "first_liquid_time_ps": {
            "1x": summary_1x.get("first_liquid_time_ps"),
            "2x": summary_2x.get("first_liquid_time_ps"),
        },
        "first_multivoxel_pool_time_ps": {
            "1x": summary_1x.get("first_multivoxel_pool_time_ps"),
            "2x": summary_2x.get("first_multivoxel_pool_time_ps"),
        },
        "heat_off_time_ps": {
            "1x": summary_1x.get("heat_off_time_ps"),
            "2x": summary_2x.get("heat_off_time_ps"),
        },
        "fully_solid_final": {
            "1x": summary_1x.get("fully_solid_final"),
            "2x": summary_2x.get("fully_solid_final"),
        },
        "interpretation": [
            "2x fields are block-averaged over each 2x2x2 subvoxel block before comparison with 1x.",
            "Temperature errors are reported on the same native Ag mask and are not alone a grid-convergence claim.",
            "Rows with phase_matches=false compare different heating/cooling states because each grid closes its heat source when its own stable-pool criterion is met.",
            "Case A cooling is internal energy redistribution in an insulated domain, not external heat loss.",
        ],
        "rows": results,
    }
    (report_dir / "mesh_comparison.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-1x", type=Path, default=HERE / "runs" / "full_1x")
    parser.add_argument("--output-2x", type=Path, default=HERE / "runs" / "full_2x")
    parser.add_argument("--report-dir", type=Path, default=HERE / "runs")
    args = parser.parse_args()
    result = compare(args.output_1x.resolve(), args.output_2x.resolve(), args.report_dir.resolve())
    print(json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
