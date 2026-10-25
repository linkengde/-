#!/usr/bin/env python3
"""Audit Ag face links introduced by conservative coarse-grid averaging.

The thermal solver assigns a nonzero Ag conductivity to every neighboring
pair of cells whose Ag volume fractions are both positive. This diagnostic
compares those coarse-grid links with exact face contact in the saved binary
QSGS geometry. It reports local face links only; it does not infer global
thermal error or replace the unthresholded finite-volume model.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL_DIR = HERE.parent


def _breaks(lo: float, hi: float, pitch: float) -> np.ndarray:
    """Exact tangential integration breaks for a piecewise-constant mask."""
    first = int(np.floor(lo / pitch)) + 1
    last = int(np.ceil(hi / pitch))
    internal = np.arange(first, last, dtype=float) * pitch
    internal = internal[(internal > lo) & (internal < hi)]
    return np.concatenate(([lo], internal, [hi]))


def _source_index(position: float, pitch: float, count: int) -> int:
    return min(max(int(np.floor(position / pitch)), 0), count - 1)


def audit_face_connectivity(
    source_ag: np.ndarray,
    ag_fraction: np.ndarray,
    size_A: np.ndarray,
    spacing_A: np.ndarray,
    positive_tolerance: float = 1.0e-12,
) -> dict[str, object]:
    """Compare nonzero coarse links to exact Ag contact on their shared faces.

    Source voxels are piecewise constant. For each coarse face the tangential
    face is partitioned at source-voxel boundaries, so the reported Ag contact
    area is exact for the stored binary QSGS field and axis-aligned grids.
    x/y coarse faces include the periodic seam; z has open outer boundaries.
    """
    source_ag = np.asarray(source_ag, dtype=bool)
    coarse_ag = np.asarray(ag_fraction, dtype=float)
    size = np.asarray(size_A, dtype=float)
    spacing = np.asarray(spacing_A, dtype=float)
    if source_ag.ndim != 3 or coarse_ag.ndim != 3:
        raise ValueError("Source and projected Ag fields must be 3D")
    if size.shape != (3,) or spacing.shape != (3,):
        raise ValueError("size_A and spacing_A must each have three components")
    if np.any(size <= 0.0) or np.any(spacing <= 0.0):
        raise ValueError("Grid dimensions must be positive")
    if not np.all(np.isfinite(coarse_ag)) or np.any(coarse_ag < 0.0):
        raise ValueError("Projected Ag fractions must be finite and non-negative")
    source_pitch = size / np.asarray(source_ag.shape, dtype=float)
    if not np.allclose(coarse_ag.shape * spacing, size, rtol=0.0, atol=1.0e-9):
        raise ValueError("Projected grid dimensions do not match the material box")

    axes: list[dict[str, object]] = []
    n_cells = coarse_ag.size
    active_ids = np.flatnonzero(coarse_ag.ravel() > positive_tolerance)
    coarse_parent = np.arange(n_cells, dtype=np.int64)
    contact_parent = np.arange(n_cells, dtype=np.int64)

    def find(parent: np.ndarray, node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = int(parent[node])
        return node

    def union(parent: np.ndarray, first: int, second: int) -> None:
        root_first = find(parent, first)
        root_second = find(parent, second)
        if root_first != root_second:
            parent[root_second] = root_first

    for axis in range(3):
        tangential = [d for d in range(3) if d != axis]
        cell_count = coarse_ag.shape[axis]
        face_ids = list(range(1, cell_count))
        periodic = axis < 2
        if periodic:
            face_ids.append(0)
        positive_faces = 0
        actual_contact_faces = 0
        false_positive_faces = 0
        contact_fractions: list[float] = []

        for face_id in face_ids:
            left_cell_axis = (face_id - 1) % cell_count
            right_cell_axis = face_id % cell_count
            face_position = face_id * spacing[axis]
            if face_id == 0:
                left_position = size[axis] - source_pitch[axis] * 1.0e-9
                right_position = source_pitch[axis] * 1.0e-9
            else:
                left_position = face_position - source_pitch[axis] * 1.0e-9
                right_position = face_position + source_pitch[axis] * 1.0e-9
            source_left_axis = _source_index(left_position, source_pitch[axis], source_ag.shape[axis])
            source_right_axis = _source_index(right_position, source_pitch[axis], source_ag.shape[axis])

            for first_tangent in range(coarse_ag.shape[tangential[0]]):
                lo0 = first_tangent * spacing[tangential[0]]
                hi0 = (first_tangent + 1) * spacing[tangential[0]]
                breaks0 = _breaks(lo0, hi0, source_pitch[tangential[0]])
                mid0 = 0.5 * (breaks0[:-1] + breaks0[1:])
                index0 = np.minimum((mid0 / source_pitch[tangential[0]]).astype(int), source_ag.shape[tangential[0]] - 1)

                for second_tangent in range(coarse_ag.shape[tangential[1]]):
                    lo1 = second_tangent * spacing[tangential[1]]
                    hi1 = (second_tangent + 1) * spacing[tangential[1]]
                    breaks1 = _breaks(lo1, hi1, source_pitch[tangential[1]])
                    mid1 = 0.5 * (breaks1[:-1] + breaks1[1:])
                    index1 = np.minimum((mid1 / source_pitch[tangential[1]]).astype(int), source_ag.shape[tangential[1]] - 1)

                    left = [0, 0, 0]
                    right = [0, 0, 0]
                    left[axis] = source_left_axis
                    right[axis] = source_right_axis
                    coarse_left = [0, 0, 0]
                    coarse_right = [0, 0, 0]
                    coarse_left[axis] = left_cell_axis
                    coarse_right[axis] = right_cell_axis
                    coarse_left[tangential[0]] = coarse_right[tangential[0]] = first_tangent
                    coarse_left[tangential[1]] = coarse_right[tangential[1]] = second_tangent
                    if coarse_ag[tuple(coarse_left)] <= positive_tolerance or coarse_ag[tuple(coarse_right)] <= positive_tolerance:
                        continue
                    positive_faces += 1
                    left_id = int(np.ravel_multi_index(tuple(coarse_left), coarse_ag.shape))
                    right_id = int(np.ravel_multi_index(tuple(coarse_right), coarse_ag.shape))
                    union(coarse_parent, left_id, right_id)

                    contact_area = 0.0
                    for i0, s0 in enumerate(index0):
                        left[tangential[0]] = right[tangential[0]] = int(s0)
                        width0 = breaks0[i0 + 1] - breaks0[i0]
                        for i1, s1 in enumerate(index1):
                            left[tangential[1]] = right[tangential[1]] = int(s1)
                            if source_ag[tuple(left)] and source_ag[tuple(right)]:
                                contact_area += width0 * (breaks1[i1 + 1] - breaks1[i1])
                    face_area = (hi0 - lo0) * (hi1 - lo1)
                    contact_fraction = contact_area / face_area
                    contact_fractions.append(contact_fraction)
                    if contact_area > 0.0:
                        actual_contact_faces += 1
                        union(contact_parent, left_id, right_id)
                    else:
                        false_positive_faces += 1

        axes.append({
            "axis": "xyz"[axis],
            "coarse_positive_conductance_faces": positive_faces,
            "faces_with_actual_qsgs_ag_contact": actual_contact_faces,
            "coarse_positive_faces_without_actual_ag_contact": false_positive_faces,
            "false_positive_fraction": false_positive_faces / positive_faces if positive_faces else 0.0,
            "mean_actual_ag_contact_area_fraction": float(np.mean(contact_fractions)) if contact_fractions else 0.0,
        })

    def component_summary(parent: np.ndarray) -> dict[str, float | int]:
        roots = [find(parent, int(node)) for node in active_ids]
        if not roots:
            return {"component_count": 0, "largest_component_fraction": 0.0}
        _, counts = np.unique(roots, return_counts=True)
        return {
            "component_count": int(counts.size),
            "largest_component_fraction": float(counts.max() / active_ids.size),
        }

    return {
        "status": "LOCAL_FACE_CONNECTIVITY_DIAGNOSTIC",
        "source_ag_shape_xyz": list(source_ag.shape),
        "coarse_grid_shape_xyz": list(coarse_ag.shape),
        "source_ag_volume_fraction": float(source_ag.mean()),
        "coarse_ag_volume_fraction": float(coarse_ag.mean()),
        "mixed_cell_fraction": float(np.mean((coarse_ag > positive_tolerance) & (coarse_ag < 1.0 - positive_tolerance))),
        "coarse_conductance_rule": "neighboring cells conduct when both projected Ag volume fractions are positive",
        "reference_contact_rule": "positive exact area where binary source QSGS Ag occupies both sides of the shared coarse face",
        "periodicity": {"x": True, "y": True, "z": False},
        "active_coarse_ag_cells": int(active_ids.size),
        "coarse_solver_connectivity": component_summary(coarse_parent),
        "direct_source_ag_contact_connectivity": component_summary(contact_parent),
        "axis_results": axes,
        "limitations": [
            "A false-positive local face link does not by itself quantify the global heat-flow error.",
            "This compares geometry topology only; it does not modify the conservative volume-fraction solver.",
            "The QSGS binary mask is treated as piecewise constant and aligned with the saved source grid.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spacing", choices=("5A", "10A"), default="10A")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    with np.load(MODEL_DIR / "qsgs_morphology.npz", allow_pickle=False) as data:
        source_ag = ~data["skeleton"].astype(bool)
        size_A = data["size_A"].astype(float)
    with np.load(HERE / f"qsgs_phase_fractions_{args.spacing}.npz", allow_pickle=False) as data:
        ag_fraction = data["ag_fraction"].astype(float)
        spacing_A = data["spacing_A"].astype(float)
    result = audit_face_connectivity(source_ag, ag_fraction, size_A, spacing_A)
    result["spacing_label"] = args.spacing
    output = args.output or HERE / f"qsgs_ag_face_connectivity_{args.spacing}.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
