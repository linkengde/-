#!/usr/bin/env python3
"""建立200³/300³网格的固定物理相关尺度几何，不改写150³基准。"""

from __future__ import annotations

import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
from scipy import ndimage


HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
PROJECT = STUDY.parent
sys.path.insert(0, str(PROJECT / "geometry_generation"))
from reference_generate_geometry import qsgs_four_parameter  # noqa: E402


GRID_SIZES = (200, 300)
SEEDS = (1105, 1136, 1182)
BASE_N = 150
BASE_P0 = 0.05
BASE_GROWTH = 0.11
BASE_SIGMA_VOXELS = 4.0
SKELETON_FRACTION = 0.30
MODEL_SIZE_M = {"ultrafine_500nm": 10.0e-6, "micro_10um": 200.0e-6}


def scale_probabilities(n: int) -> tuple[float, float, float]:
    """按细分体素数换算成核密度和粗体素跨度的等效概率。"""
    ratio = BASE_N / n
    p0 = 1.0 - (1.0 - BASE_P0) ** (ratio**3)
    growth = 1.0 - (1.0 - BASE_GROWTH) ** ratio
    return p0, growth, growth


def correlated_exact_fraction(raw: np.ndarray, sigma_voxels: float) -> np.ndarray:
    score = ndimage.gaussian_filter(raw.astype(np.float64), sigma=sigma_voxels, mode="reflect")
    count = int(round(raw.size * SKELETON_FRACTION))
    selected = np.argpartition(score.ravel(), -count)[-count:]
    skeleton = np.zeros(raw.shape, dtype=bool)
    skeleton.ravel()[selected] = True
    return skeleton


def main() -> None:
    all_rows = []
    geometry_dir = STUDY / "geometries"
    for n in GRID_SIZES:
        shape = (n, n, n)
        sigma = BASE_SIGMA_VOXELS * n / BASE_N
        p0, growth_x, growth_y = scale_probabilities(n)
        growth = (growth_x, growth_y, growth_y)
        for seed in SEEDS:
            started = time.perf_counter()
            print(json.dumps({"status": "generating", "n": n, "seed": seed,
                              "p0": p0, "growth_xyz": growth, "sigma_voxels": sigma}), flush=True)
            raw, generations = qsgs_four_parameter(shape, p0, growth, SKELETON_FRACTION, seed)
            qsgs_seconds = time.perf_counter() - started
            skeleton = correlated_exact_fraction(raw, sigma)
            del raw
            output_dir = geometry_dir / f"n{n}"
            output_dir.mkdir(parents=True, exist_ok=True)
            output = output_dir / f"qsgs_seed_{seed}.npz"
            np.savez_compressed(
                output,
                skeleton=skeleton,
                shape_xyz=np.asarray(shape, dtype=np.int32),
                seed=np.asarray(seed, dtype=np.int32),
                p0=np.asarray(p0, dtype=np.float64),
                growth_xyz=np.asarray(growth, dtype=np.float64),
                sigma_voxels=np.asarray(sigma, dtype=np.float64),
                sigma_fraction_of_domain=np.asarray(BASE_SIGMA_VOXELS / BASE_N, dtype=np.float64),
                target_skeleton_fraction=np.asarray(SKELETON_FRACTION, dtype=np.float64),
                growth_generations=np.asarray(generations, dtype=np.int32),
                model_size_m_by_name=np.asarray(list(MODEL_SIZE_M.values()), dtype=np.float64),
                sigma_m_by_model=np.asarray(
                    [BASE_SIGMA_VOXELS / BASE_N * length for length in MODEL_SIZE_M.values()],
                    dtype=np.float64,
                ),
            )
            row = {
                "n": n,
                "shape_xyz": shape,
                "seed": seed,
                "p0": p0,
                "growth_xyz": growth,
                "qsgs_generations": generations,
                "sigma_voxels": sigma,
                "sigma_fraction_of_domain": BASE_SIGMA_VOXELS / BASE_N,
                "skeleton_fraction": float(skeleton.mean()),
                "qsgs_wall_seconds": qsgs_seconds,
                "total_wall_seconds": time.perf_counter() - started,
                "process_peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
                "geometry_file": str(output.relative_to(STUDY)),
                "geometry_bytes": output.stat().st_size,
            }
            all_rows.append(row)
            print(json.dumps(row), flush=True)
    summary = {
        "status": "RESOLUTION_GEOMETRIES_GENERATED",
        "baseline_grid_n": BASE_N,
        "grid_sizes": list(GRID_SIZES),
        "seeds": list(SEEDS),
        "base_parameters": {
            "p0": BASE_P0,
            "growth_xyz": [BASE_GROWTH] * 3,
            "sigma_voxels": BASE_SIGMA_VOXELS,
            "sigma_fraction_of_domain": BASE_SIGMA_VOXELS / BASE_N,
            "skeleton_fraction": SKELETON_FRACTION,
        },
        "probability_scaling": {
            "p0_n": "1 - (1 - p0_150)^((150/n)^3), maintaining expected nucleus occupancy per physical volume",
            "growth_n": "1 - (1 - growth_150)^(150/n), matching coarse-cell growth success probability",
            "sigma_n": "4 * n / 150 voxels, keeping Gaussian correlation length fixed as a fraction of physical RVE length",
        },
        "physical_sigma_m_by_model": {
            name: BASE_SIGMA_VOXELS / BASE_N * length for name, length in MODEL_SIZE_M.items()
        },
        "model_sizes_m": MODEL_SIZE_M,
        "caveat": "These are scaled starting probabilities; geometry metrics decide whether further parameter scans are needed. They do not constitute mesh-convergence evidence by themselves.",
        "cases": all_rows,
    }
    (STUDY / "geometry_generation" / "generation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
