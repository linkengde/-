#!/usr/bin/env python3
"""Audit thresholded phase connectivity on projected QSGS fractions.

This is a resolution diagnostic only. The thermal model should use the
unthresholded phase fractions unless a separate connectivity-preserving
binary geometry has been justified.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.ndimage import label


HERE = Path(__file__).resolve().parent


def _component_metrics(mask: np.ndarray) -> dict[str, object]:
    mask = np.asarray(mask, dtype=bool)
    structure = np.ones((3, 3, 3), dtype=bool)  # 26-neighbor connectivity
    labels, local_count = label(mask, structure=structure)
    parent = np.arange(local_count + 1, dtype=np.int64)

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    def union(a: int, b: int) -> None:
        if not a or not b:
            return
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[root_b] = root_a

    nx, ny, nz = mask.shape
    # Merge 26-neighbor contacts across the periodic x seam; z remains open.
    for dy in (-1, 0, 1):
        for dz in (-1, 0, 1):
            for y in range(ny):
                z0 = max(0, -dz)
                z1 = min(nz, nz - dz)
                for z in range(z0, z1):
                    union(int(labels[0, y, z]), int(labels[-1, (y + dy) % ny, z + dz]))
    # Merge contacts across the periodic y seam, including x/y corners.
    for dx in (-1, 0, 1):
        for dz in (-1, 0, 1):
            for x in range(nx):
                z0 = max(0, -dz)
                z1 = min(nz, nz - dz)
                for z in range(z0, z1):
                    union(int(labels[x, 0, z]), int(labels[(x + dx) % nx, -1, z + dz]))

    local_sizes = np.bincount(labels.ravel(), minlength=local_count + 1)
    component_sizes: dict[int, int] = {}
    for local_id in range(1, local_count + 1):
        root = find(local_id)
        component_sizes[root] = component_sizes.get(root, 0) + int(local_sizes[local_id])
    occupied = int(mask.sum())
    if not component_sizes:
        return {
            "occupied_cells": 0,
            "component_count": 0,
            "largest_component_cells": 0,
            "largest_component_fraction": 0.0,
            "largest_component_touches_z_min": False,
            "largest_component_touches_z_max": False,
        }

    largest_root = max(component_sizes, key=component_sizes.get)
    roots = np.zeros(local_count + 1, dtype=np.int64)
    for local_id in range(1, local_count + 1):
        roots[local_id] = find(local_id)
    largest_mask = (labels > 0) & (roots[labels] == largest_root)
    z_indices = np.flatnonzero(largest_mask.any(axis=(0, 1)))
    largest_cells = component_sizes[largest_root]
    return {
        "occupied_cells": occupied,
        "component_count": len(component_sizes),
        "largest_component_cells": largest_cells,
        "largest_component_fraction": largest_cells / occupied if occupied else 0.0,
        "largest_component_touches_z_min": bool(z_indices.size and z_indices[0] == 0),
        "largest_component_touches_z_max": bool(z_indices.size and z_indices[-1] == nz - 1),
    }


def main() -> None:
    records = []
    for tag in ("5A", "10A"):
        with np.load(HERE / f"qsgs_phase_fractions_{tag}.npz", allow_pickle=False) as data:
            ag_fraction = data["ag_fraction"]
            skeleton_fraction = data["skeleton_fraction"]
            spacing_A = data["spacing_A"].tolist()
        records.append({
            "requested_spacing": tag,
            "spacing_A_xyz": spacing_A,
            "threshold_fraction": 0.5,
            "connectivity": "26-neighbor; periodic x/y; open z",
            "Ag": _component_metrics(ag_fraction > 0.5),
            "Ti3SiC2_skeleton": _component_metrics(skeleton_fraction > 0.5),
        })
    report = {
        "status": "THRESHOLDED_COARSE_GRID_TOPOLOGY_DIAGNOSTIC_ONLY",
        "records": records,
        "limitations": [
            "Thresholded connectivity is not used to construct or modify the conservative volume-fraction fields.",
            "This diagnostic reports connected components and z-surface contact, not periodic winding rank.",
            "It cannot establish material-scale percolation or any Ag/scaffold interface physics.",
        ],
    }
    output = HERE / "qsgs_projected_connectivity_report.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
