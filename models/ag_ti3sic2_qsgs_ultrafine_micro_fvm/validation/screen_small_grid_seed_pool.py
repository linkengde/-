#!/usr/bin/env python3
"""用独立种子池检查50³/100³网格下QSGS算法的统计敏感性。"""

from __future__ import annotations

import csv
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "geometry_generation"))
from reference_generate_geometry import qsgs_four_parameter  # noqa: E402

OUT = ROOT / "validation"
SEEDS = range(1100, 1251)
SIZES = (50, 100)
P0 = 0.05
GROWTH = (0.11, 0.11, 0.11)
FRACTION = 0.30
SIGMA = 4.0
STRUCTURE = ndimage.generate_binary_structure(3, 1)


def summarize(mask: np.ndarray) -> tuple[int, float, list[bool]]:
    labels, count = ndimage.label(mask, structure=STRUCTURE)
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    sizes[0] = 0
    largest_id = int(np.argmax(sizes)) if count else 0
    largest_fraction = float(sizes[largest_id] / mask.sum()) if mask.any() else 0.0
    spans = [bool(largest_id and np.any(np.take(labels, 0, axis=a) == largest_id)
                  and np.any(np.take(labels, -1, axis=a) == largest_id)) for a in range(3)]
    return int(count), largest_fraction, spans


def evaluate(task: tuple[int, int]) -> dict[str, object]:
    n, seed = task
    shape = (n, n, n)
    raw, generations = qsgs_four_parameter(shape, P0, GROWTH, FRACTION, seed)
    score = ndimage.gaussian_filter(raw.astype(np.float64), sigma=SIGMA, mode="reflect")
    target = int(round(raw.size * FRACTION))
    chosen = np.argpartition(score.ravel(), -target)[-target:]
    skeleton = np.zeros(shape, dtype=bool)
    skeleton.ravel()[chosen] = True
    ag = ~skeleton
    sc_count, sc_fraction, sc_spans = summarize(skeleton)
    ag_count, ag_fraction, ag_spans = summarize(ag)
    p50 = None
    p10 = None
    p90 = None
    if sc_fraction >= 0.99 and all(sc_spans):
        d = ndimage.distance_transform_edt(skeleton)
        maxima = skeleton & (d >= ndimage.maximum_filter(d, size=3, mode="constant") - 1e-12)
        diameters = 2.0 * d[maxima]
        p10, p50, p90 = (float(np.percentile(diameters, q)) for q in (10, 50, 90))
    return {
        "grid_n": n, "seed": seed, "p0": P0, "growth_probability_xyz": GROWTH,
        "target_skeleton_fraction": FRACTION, "correlation_sigma_voxels": SIGMA,
        "qsgs_growth_generations": generations, "skeleton_fraction": float(skeleton.mean()),
        "ag_fraction": float(ag.mean()), "skeleton_component_count_6face": sc_count,
        "skeleton_largest_component_fraction": sc_fraction, "skeleton_spans_xyz": sc_spans,
        "ag_component_count_6face": ag_count, "ag_largest_component_fraction": ag_fraction,
        "ag_spans_xyz": ag_spans, "skeleton_p10_diameter_voxels": p10,
        "skeleton_p50_diameter_voxels": p50, "skeleton_p90_diameter_voxels": p90,
        "both_phases_pass_99pct_spanning": bool(sc_fraction >= 0.99 and ag_fraction >= 0.99
                                                   and all(sc_spans) and all(ag_spans)),
    }


def main() -> None:
    tasks = [(n, seed) for n in SIZES for seed in SEEDS]
    with ProcessPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(evaluate, tasks, chunksize=1))
    columns = sorted({k for row in rows for k in row})
    with (OUT / "small_grid_seed_screen.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    passing = {}
    for n in SIZES:
        passing[str(n)] = [int(r["seed"]) for r in rows if r["grid_n"] == n
                           and r["both_phases_pass_99pct_spanning"]]
    summary = {
        "parameter_values": {"p0": P0, "growth_probability_xyz": GROWTH,
                             "target_skeleton_fraction": FRACTION,
                             "correlation_sigma_voxels": SIGMA},
        "seed_range_inclusive": [min(SEEDS), max(SEEDS)],
        "grid_sizes": list(SIZES), "tested_case_count": len(rows),
        "passing_seed_count_by_grid": {n: len(seeds) for n, seeds in passing.items()},
        "passing_seeds_by_grid": passing,
        "warning": "Uses the same fixed voxel correlation length; component statistics vary with grid size. This is an algorithm and size-sensitivity screen, not an SEM calibration.",
    }
    (OUT / "small_grid_seed_screen.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
