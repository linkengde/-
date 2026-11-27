#!/usr/bin/env python3
"""生成 150³ Ag/Ti₃SiC₂ 双连续随机几何对照。"""

from __future__ import annotations

import json
from pathlib import Path
import resource
import time

import numpy as np
from scipy import ndimage

from reference_generate_geometry import qsgs_four_parameter


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "geometry_generation"
SHAPE = (150, 150, 150)
SEEDS = (1105, 1136, 1182)
P0 = 0.05
GROWTH = (0.11, 0.11, 0.11)
TARGET_SKELETON_FRACTION = 0.30
CORRELATION_SIGMA_VOXELS = 4.0


def qsgs_growth_rank(
    shape: tuple[int, int, int],
    p0: float,
    growth: tuple[float, float, float],
    target_fraction: float,
    seed: int,
) -> tuple[np.ndarray, int]:
    """复现仓库四参数生长随机流，同时记录每格首次加入的生长代数。"""
    rng = np.random.default_rng(seed)
    occupied = rng.random(shape) < p0
    rank = np.full(shape, -1, dtype=np.int16)
    rank[occupied] = 0
    target_count = int(round(target_fraction * occupied.size))
    if not occupied.any():
        occupied.flat[int(rng.integers(occupied.size))] = True
        rank.flat[int(np.flatnonzero(occupied)[0])] = 0

    generations = 0
    while int(occupied.sum()) < target_count:
        frontier_accept = np.zeros(shape, dtype=bool)
        for axis, probability in enumerate(growth):
            for direction in (-1, 1):
                neighbor = np.roll(occupied, direction, axis=axis)
                edge = [slice(None)] * 3
                edge[axis] = 0 if direction == 1 else -1
                neighbor[tuple(edge)] = False
                frontier_accept |= neighbor & (rng.random(shape) < probability)

        frontier = np.flatnonzero(frontier_accept & ~occupied)
        needed = target_count - int(occupied.sum())
        if frontier.size == 0:
            raise RuntimeError("QSGS 生长前沿在达到目标相分数前消失")
        if frontier.size > needed:
            frontier = rng.choice(frontier, size=needed, replace=False)
        occupied.flat[frontier] = True
        rank.flat[frontier] = generations + 1
        generations += 1

    reference_mask, reference_generations = qsgs_four_parameter(
        shape, p0, growth, target_fraction, seed
    )
    if generations != reference_generations or not np.array_equal(occupied, reference_mask):
        raise RuntimeError("记录生长代数时偏离仓库 QSGS 随机流")
    return rank, generations


def correlated_exact_fraction(
    qsgs_skeleton: np.ndarray,
    sigma_voxels: float,
) -> tuple[np.ndarray, np.ndarray]:
    """平滑 QSGS 指示场以设定相关尺度，再按秩选择精确总体积分数。"""
    score = ndimage.gaussian_filter(
        qsgs_skeleton.astype(np.float64), sigma=sigma_voxels, mode="reflect"
    )
    skeleton_count = int(round(qsgs_skeleton.size * TARGET_SKELETON_FRACTION))
    selected = np.argpartition(score.ravel(), -skeleton_count)[-skeleton_count:]
    skeleton = np.zeros(qsgs_skeleton.shape, dtype=bool)
    skeleton.ravel()[selected] = True
    return skeleton, score.astype(np.float32)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    run_started = time.perf_counter()
    summaries = []
    for seed in SEEDS:
        started = time.perf_counter()
        rank, generations = qsgs_growth_rank(
            SHAPE, P0, GROWTH, TARGET_SKELETON_FRACTION, seed
        )
        qsgs_mask = rank >= 0
        after_growth = time.perf_counter()
        skeleton, score = correlated_exact_fraction(qsgs_mask, CORRELATION_SIGMA_VOXELS)
        after_correlation = time.perf_counter()
        ag = ~skeleton
        destination = OUTPUT / f"qsgs_seed_{seed}.npz"
        np.savez_compressed(
            destination,
            qsgs_growth_rank=rank,
            qsgs_skeleton=qsgs_mask,
            correlation_score=score,
            skeleton=skeleton,
            ag=ag,
            shape_xyz=np.asarray(SHAPE, dtype=np.int32),
            seed=np.asarray(seed, dtype=np.int32),
            qsgs_p0=np.asarray(P0, dtype=np.float64),
            qsgs_growth=np.asarray(GROWTH, dtype=np.float64),
            qsgs_generations=np.asarray(generations, dtype=np.int32),
            correlation_sigma_voxels=np.asarray(CORRELATION_SIGMA_VOXELS, dtype=np.float64),
        )
        summaries.append(
            {
                "seed": seed,
                "shape_xyz": SHAPE,
                "p0": P0,
                "px_py_pz": GROWTH,
                "qsgs_generations": generations,
                "qsgs_initial_fraction": float(qsgs_mask.mean()),
                "final_skeleton_fraction": float(skeleton.mean()),
                "final_ag_fraction": float(ag.mean()),
                "correlation_sigma_voxels": CORRELATION_SIGMA_VOXELS,
                "archive": str(destination.relative_to(ROOT)),
                "archive_bytes": destination.stat().st_size,
                "growth_seconds": after_growth - started,
                "correlation_seconds": after_correlation - after_growth,
                "process_peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
            }
        )
        print(json.dumps(summaries[-1]), flush=True)

    payload = {
        "algorithm": "repository four-parameter QSGS growth rank, followed by a Gaussian-correlated rank threshold",
        "screened_parameters": {
            "p0": P0,
            "px_py_pz": GROWTH,
            "target_skeleton_fraction": TARGET_SKELETON_FRACTION,
            "correlation_sigma_voxels": CORRELATION_SIGMA_VOXELS,
            "seed_selection": "selected from a 151-seed 150³ screen after phase connectivity and P50 thickness filters; 50³/100³ outcomes are recorded separately and are not selection criteria",
        },
        "models": summaries,
        "total_wall_seconds": time.perf_counter() - run_started,
    }
    (OUTPUT / "generation_summary.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
