#!/usr/bin/env python3
"""Build an isolated Ag 70 vol% / Ti3SiC2 30 vol% thermal mask.

This reuses the committed QSGS growth-rank realization and the original
periodic face-connectivity cleanup. It never edits the parent mask or
regenerates the stochastic growth history.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage


RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
PARENT_GEOMETRY = PROJECT_DIR / "qsgs_morphology.npz"
GROWTH_RANK = PROJECT_DIR / "qsgs_growth_rank.npy"
OUTPUT_GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
OUTPUT_REPORT = RUN_DIR / "composition_70vol_geometry.json"
TARGET_SKELETON_FRACTION = 0.30


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def periodic_labels(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Six-face connected labels with x/y seam unions and open z faces."""
    labels, count = ndimage.label(mask)
    parent = np.arange(count + 1)

    def root(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    for axis in (0, 1):
        low = np.take(labels, 0, axis=axis).ravel()
        high = np.take(labels, -1, axis=axis).ravel()
        for first, second in np.unique(np.stack((low, high), axis=1), axis=0):
            if first and second:
                first_root, second_root = root(int(first)), root(int(second))
                if first_root != second_root:
                    parent[second_root] = first_root

    roots = np.fromiter((root(index) for index in range(count + 1)), dtype=np.int64)
    labels = roots[labels]
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    sizes[0] = 0
    return labels, sizes


def clean_mask(mask: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """Match build_model.py's alternating largest-component cleanup."""
    changes: list[int] = []
    mask = np.asarray(mask, dtype=bool).copy()
    for _ in range(8):
        previous = mask.copy()
        labels, sizes = periodic_labels(mask)
        mask = labels == int(np.argmax(sizes))
        labels, sizes = periodic_labels(~mask)
        mask = ~(labels == int(np.argmax(sizes)))
        changes.append(int(np.count_nonzero(mask != previous)))
        if changes[-1] == 0:
            break
    return mask, changes


def largest_component_metrics(mask: np.ndarray) -> dict[str, object]:
    labels, sizes = periodic_labels(mask)
    largest = int(np.argmax(sizes))
    largest_mask = labels == largest
    coordinates = np.where(largest_mask)
    seam_edges = [
        int(np.count_nonzero(np.take(largest_mask, 0, axis=axis) & np.take(largest_mask, -1, axis=axis)))
        for axis in (0, 1)
    ]
    return {
        "periodic_face_connected_components": int(np.count_nonzero(sizes)),
        "largest_component_voxel_fraction_of_phase": float(sizes[largest] / np.count_nonzero(mask)),
        "largest_component_seam_edges_xy": seam_edges,
        "largest_component_open_z_extent_voxels": int(coordinates[2].max() - coordinates[2].min() + 1),
    }


def select_mask(rank: np.ndarray) -> tuple[np.ndarray, int, int, list[int]]:
    """Binary-search the existing realization's rank threshold after cleanup."""
    lower = int(rank.size * 0.25)
    upper = int(rank.size * 0.40)
    best: tuple[float, int, np.ndarray, int, list[int]] | None = None
    candidates: dict[int, tuple[np.ndarray, int, list[int]]] = {}
    for _ in range(22):
        threshold = (lower + upper) // 2
        if threshold not in candidates:
            cleaned, changes = clean_mask(rank < threshold)
            candidates[threshold] = (cleaned, int(np.count_nonzero(rank < threshold)), changes)
        cleaned, raw_count, changes = candidates[threshold]
        error = abs(float(cleaned.mean()) - TARGET_SKELETON_FRACTION)
        if best is None or error < best[0]:
            best = (error, threshold, cleaned.copy(), raw_count, changes)
        if cleaned.mean() < TARGET_SKELETON_FRACTION:
            lower = threshold + 1
        else:
            upper = threshold
        if lower >= upper:
            for final_threshold in {lower - 1, lower, upper, upper + 1}:
                if final_threshold < 1 or final_threshold in candidates:
                    continue
                cleaned, changes = clean_mask(rank < final_threshold)
                candidates[final_threshold] = (
                    cleaned,
                    int(np.count_nonzero(rank < final_threshold)),
                    changes,
                )
                error = abs(float(cleaned.mean()) - TARGET_SKELETON_FRACTION)
                if best is None or error < best[0]:
                    best = (error, final_threshold, cleaned.copy(), candidates[final_threshold][1], changes)
            break
    if best is None:
        raise RuntimeError("No candidate QSGS rank threshold was evaluated")
    _, threshold, mask, raw_count, changes = best
    return mask, threshold, raw_count, changes


def main() -> None:
    with np.load(PARENT_GEOMETRY, allow_pickle=False) as archive:
        parent_mask = np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
    rank = np.load(GROWTH_RANK, mmap_mode="r")
    if rank.shape != parent_mask.shape:
        raise ValueError(f"Growth rank shape {rank.shape} does not match parent mask {parent_mask.shape}")

    parent_hash = sha256(PARENT_GEOMETRY)
    rank_hash = sha256(GROWTH_RANK)
    mask, threshold, raw_count, cleanup_changes = select_mask(rank)
    if not np.isclose(mask.mean() + (~mask).mean(), 1.0, rtol=0.0, atol=1.0e-15):
        raise AssertionError("The binary phase masks do not partition the domain")
    ceramic_metrics = largest_component_metrics(mask)
    silver_metrics = largest_component_metrics(~mask)
    if ceramic_metrics["periodic_face_connected_components"] != 1:
        raise AssertionError("Ti3SiC2 phase is not one periodic face-connected component")
    if silver_metrics["periodic_face_connected_components"] != 1:
        raise AssertionError("Ag phase is not one periodic face-connected component")
    if any(value == 0 for value in ceramic_metrics["largest_component_seam_edges_xy"]):
        raise AssertionError("Ti3SiC2 does not connect across both periodic lateral seams")
    if any(value == 0 for value in silver_metrics["largest_component_seam_edges_xy"]):
        raise AssertionError("Ag does not connect across both periodic lateral seams")

    np.savez_compressed(
        OUTPUT_GEOMETRY,
        skeleton=mask,
        size_A=size_A,
        seed=np.asarray(20261024, dtype=np.int64),
        composition_basis=np.asarray("Ag:Ti3SiC2 volume fraction = 0.70:0.30"),
        parent_geometry_sha256=np.asarray(parent_hash),
        growth_rank_sha256=np.asarray(rank_hash),
    )
    with np.load(OUTPUT_GEOMETRY, allow_pickle=False) as saved:
        if not np.array_equal(saved["skeleton"], mask) or not np.array_equal(saved["size_A"], size_A):
            raise AssertionError("Written composition geometry did not round-trip")

    report = {
        "status": "INDEPENDENT_70_30_VOLUME_THERMAL_GEOMETRY_CREATED",
        "target_phase_volume_fraction": {"Ag": 0.70, "Ti3SiC2": 0.30},
        "achieved_phase_volume_fraction": {"Ag": float((~mask).mean()), "Ti3SiC2": float(mask.mean())},
        "absolute_skeleton_fraction_error": float(abs(mask.mean() - TARGET_SKELETON_FRACTION)),
        "geometry_path": str(OUTPUT_GEOMETRY.relative_to(PROJECT_DIR)),
        "geometry_sha256": sha256(OUTPUT_GEOMETRY),
        "shape_xyz": list(mask.shape),
        "domain_size_A": size_A.tolist(),
        "voxel_size_A": (size_A / np.asarray(mask.shape)).tolist(),
        "construction": {
            "method": "threshold the already generated QSGS growth-rank realization, then apply the original alternating largest face-connected component cleanup",
            "stochastic_growth_regenerated": False,
            "rank_threshold_exclusive": threshold,
            "raw_skeleton_voxel_count": raw_count,
            "cleanup_voxel_changes_by_iteration": cleanup_changes,
            "seed": 20261024,
            "parent_geometry_sha256": parent_hash,
            "growth_rank_sha256": rank_hash,
            "parent_geometry_edited": False,
        },
        "phase_connectivity": {"Ti3SiC2": ceramic_metrics, "Ag": silver_metrics},
        "mass_ratio_from_candidate_bulk_densities_Ag_over_Ti3SiC2": (
            0.70 * 10490.0 / (0.30 * 4520.0)
        ),
        "mass_fraction_from_candidate_bulk_densities": {
            "Ag": 0.70 * 10490.0 / (0.70 * 10490.0 + 0.30 * 4520.0),
            "Ti3SiC2": 0.30 * 4520.0 / (0.70 * 10490.0 + 0.30 * 4520.0),
        },
        "limitations": [
            "This is a new continuum thermal mask derived from the existing QSGS growth history; it is not an independent random realization or a regenerated atomic model.",
            "The exact experimental density/porosity mapping remains unverified; mass conversion uses candidate densities pending literature audit.",
            "Connectivity checks establish periodic face connectivity, not minimum neck thickness or statistical representativeness.",
        ],
    }
    OUTPUT_REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if sha256(PARENT_GEOMETRY) != parent_hash:
        raise AssertionError("Parent QSGS geometry changed during candidate creation")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
