#!/usr/bin/env python3
"""对最终150³种子补做50³/100³几何尺寸敏感性筛查。"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "geometry_generation"))
from reference_generate_geometry import qsgs_four_parameter


OUT = ROOT / "validation"
SEEDS = (1105, 1136, 1182)
SIZES = (50, 100)
P0 = 0.05
GROWTH = (0.11, 0.11, 0.11)
FRACTION = 0.30
SIGMA = 4.0
STRUCTURE = ndimage.generate_binary_structure(3, 1)


def phase_stats(mask: np.ndarray) -> dict[str, object]:
    labels, count = ndimage.label(mask, structure=STRUCTURE)
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    sizes[0] = 0
    largest_id = int(np.argmax(sizes)) if count else 0
    largest = float(sizes[largest_id] / mask.sum()) if mask.any() else 0.0
    spans = [
        bool(largest_id and np.any(np.take(labels, 0, axis=a) == largest_id)
             and np.any(np.take(labels, -1, axis=a) == largest_id))
        for a in range(3)
    ]
    distance = ndimage.distance_transform_edt(mask)
    maxima = mask & (distance >= ndimage.maximum_filter(distance, size=3, mode="constant") - 1e-12)
    diameters = 2.0 * distance[maxima]
    return {
        "component_count_6face": int(count),
        "largest_component_fraction": largest,
        "largest_component_spans_xyz": spans,
        "sphere_center_diameter_p10_voxels": float(np.percentile(diameters, 10)),
        "sphere_center_diameter_p50_voxels": float(np.percentile(diameters, 50)),
        "sphere_center_diameter_p90_voxels": float(np.percentile(diameters, 90)),
        "all_components_ge_99pct": largest >= 0.99 and all(spans),
    }


def main() -> None:
    rows: list[dict[str, object]] = []
    for n in SIZES:
        shape = (n, n, n)
        for seed in SEEDS:
            raw, generations = qsgs_four_parameter(shape, P0, GROWTH, FRACTION, seed)
            score = ndimage.gaussian_filter(raw.astype(np.float64), sigma=SIGMA, mode="reflect")
            count = int(round(raw.size * FRACTION))
            chosen = np.argpartition(score.ravel(), -count)[-count:]
            skeleton = np.zeros(shape, dtype=bool)
            skeleton.ravel()[chosen] = True
            ag = ~skeleton
            sc_stats = phase_stats(skeleton)
            ag_stats = phase_stats(ag)
            rows.append({
                "grid_n": n,
                "seed": seed,
                "p0": P0,
                "growth_xyz": json.dumps(GROWTH),
                "correlation_sigma_voxels": SIGMA,
                "growth_generations": generations,
                "skeleton_fraction": float(skeleton.mean()),
                "ag_fraction": float(ag.mean()),
                "skeleton_stats": json.dumps(sc_stats, separators=(",", ":")),
                "ag_stats": json.dumps(ag_stats, separators=(",", ":")),
                "both_phases_pass_99pct_spanning_screen": bool(
                    sc_stats["all_components_ge_99pct"] and ag_stats["all_components_ge_99pct"]
                ),
            })
    with (OUT / "final_seed_small_grid_checks.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "grid_sizes": list(SIZES),
        "seeds": list(SEEDS),
        "parameter_values": {"p0": P0, "growth_xyz": GROWTH, "target_skeleton_fraction": FRACTION,
                             "correlation_sigma_voxels": SIGMA},
        "warning": "Seed qualification was performed at 150³. These lower-resolution runs are separate scale-sensitivity checks; pass/fail can differ with grid size.",
        "cases": rows,
    }
    (OUT / "final_seed_small_grid_checks.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"cases": len(rows), "passes": sum(r["both_phases_pass_99pct_spanning_screen"] for r in rows)}), flush=True)


if __name__ == "__main__":
    main()
