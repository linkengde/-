#!/usr/bin/env python3
"""在固定小网格筛选参数下筛选 150³ 随机种子。"""

from __future__ import annotations

import csv
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import time

import numpy as np
from scipy import ndimage

from reference_generate_geometry import qsgs_four_parameter


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation"
SHAPE = (150, 150, 150)
P0 = 0.05
GROWTH = (0.11, 0.11, 0.11)
TARGET_FRACTION = 0.30
SIGMA = 4.0
CONNECTIVITY = ndimage.generate_binary_structure(3, 1)


def evaluate_seed(seed: int) -> dict[str, object]:
    raw, generations = qsgs_four_parameter(SHAPE, P0, GROWTH, TARGET_FRACTION, seed)
    score = ndimage.gaussian_filter(raw.astype(np.float64), SIGMA, mode="reflect")
    count = int(round(raw.size * TARGET_FRACTION))
    chosen = np.argpartition(score.ravel(), -count)[-count:]
    skeleton = np.zeros(SHAPE, dtype=bool)
    skeleton.ravel()[chosen] = True
    del score, raw, chosen

    labels, component_count = ndimage.label(skeleton, structure=CONNECTIVITY)
    component_sizes = np.bincount(labels.ravel(), minlength=component_count + 1)
    component_sizes[0] = 0
    largest_id = int(np.argmax(component_sizes))
    skeleton_largest_fraction = float(component_sizes[largest_id] / skeleton.sum())
    skeleton_spans = [
        bool(
            np.any(np.take(labels, 0, axis=axis) == largest_id)
            and np.any(np.take(labels, -1, axis=axis) == largest_id)
        )
        for axis in range(3)
    ]
    del labels, component_sizes

    ag = ~skeleton
    ag_labels, ag_component_count = ndimage.label(ag, structure=CONNECTIVITY)
    ag_sizes = np.bincount(ag_labels.ravel(), minlength=ag_component_count + 1)
    ag_sizes[0] = 0
    ag_largest_id = int(np.argmax(ag_sizes))
    ag_largest_fraction = float(ag_sizes[ag_largest_id] / ag.sum())
    ag_spans = [
        bool(
            np.any(np.take(ag_labels, 0, axis=axis) == ag_largest_id)
            and np.any(np.take(ag_labels, -1, axis=axis) == ag_largest_id)
        )
        for axis in range(3)
    ]
    row: dict[str, object] = {
        "seed": seed,
        "p0": P0,
        "px_py_pz": GROWTH[0],
        "target_skeleton_fraction": TARGET_FRACTION,
        "correlation_sigma_voxels": SIGMA,
        "qsgs_growth_generations": generations,
        "skeleton_fraction": float(skeleton.mean()),
        "ag_fraction": float(ag.mean()),
        "skeleton_component_count_6face": int(component_count),
        "skeleton_largest_component_fraction": skeleton_largest_fraction,
        "skeleton_spans_xyz": skeleton_spans,
        "ag_component_count_6face": int(ag_component_count),
        "ag_largest_component_fraction": ag_largest_fraction,
        "ag_spans_xyz": ag_spans,
    }

    if (
        skeleton_largest_fraction >= 0.99
        and ag_largest_fraction >= 0.99
        and all(skeleton_spans)
        and all(ag_spans)
    ):
        distance = ndimage.distance_transform_edt(skeleton)
        maxima = skeleton & (distance == ndimage.maximum_filter(distance, size=3, mode="constant"))
        sphere_diameters = 2.0 * distance[maxima]
        row.update(
            {
                "maximal_inscribed_sphere_centers": int(sphere_diameters.size),
                "skeleton_diameter_p10_voxels": float(np.percentile(sphere_diameters, 10)),
                "skeleton_diameter_p50_voxels": float(np.percentile(sphere_diameters, 50)),
                "skeleton_diameter_p90_voxels": float(np.percentile(sphere_diameters, 90)),
            }
        )
    else:
        row.update(
            {
                "maximal_inscribed_sphere_centers": None,
                "skeleton_diameter_p10_voxels": None,
                "skeleton_diameter_p50_voxels": None,
                "skeleton_diameter_p90_voxels": None,
            }
        )
    return row


def main() -> None:
    started = time.perf_counter()
    seeds = range(1100, 1251)
    with ProcessPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(evaluate_seed, seeds))
    columns = sorted({key for row in rows for key in row})
    with (OUT / "150cube_seed_screen.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    passed = [
        row
        for row in rows
        if row["skeleton_largest_component_fraction"] >= 0.99
        and row["ag_largest_component_fraction"] >= 0.99
        and all(row["skeleton_spans_xyz"])
        and all(row["ag_spans_xyz"])
        and 6.0 <= row["skeleton_diameter_p50_voxels"] <= 9.0
    ]
    summary = {
        "shape_xyz": SHAPE,
        "p0": P0,
        "px_py_pz": GROWTH,
        "target_skeleton_fraction": TARGET_FRACTION,
        "correlation_sigma_voxels": SIGMA,
        "seed_range_inclusive": [1100, 1250],
        "tested_seed_count": len(rows),
        "passing_seed_count": len(passed),
        "passing_seeds": passed,
        "wall_seconds": time.perf_counter() - started,
    }
    (OUT / "150cube_seed_screen.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
