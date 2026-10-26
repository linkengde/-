#!/usr/bin/env python3
"""Compare voxel-resolved and projected Ag conductance in one QSGS window.

This is a geometry-only steady resistor-network diagnostic. It uses constant
intrinsic Ag conductivity, insulated non-drive faces, and unit temperature
difference on opposing faces of each tested axis. It does not model interface
resistance or claim experimental thermal conductivity.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import cg

from qsgs_face_contacts import exact_face_contact_fractions
from qsgs_component_coarsegrain import build_component_graph, component_effective_conductance
from qsgs_projection import _overlap_weights, project_phase_mask


HERE = Path(__file__).resolve().parent
MODEL_DIR = HERE.parent
SOURCE = MODEL_DIR / "qsgs_morphology.npz"
ROI_CASES = (
    {"name": "regression_1", "start": (36, 60, 18), "shape": (24, 24, 16), "motif": "prior regression window"},
    {"name": "regression_2", "start": (56, 12, 0), "shape": (24, 24, 16), "motif": "prior regression window"},
    {"name": "regression_3", "start": (64, 76, 0), "shape": (24, 24, 16), "motif": "prior regression window"},
    {"name": "low_fill_through", "start": (32, 32, 32), "shape": (24, 24, 16), "motif": "lower-Ag, through-connected path"},
    {"name": "high_fill_through", "start": (8, 8, 32), "shape": (24, 24, 16), "motif": "higher-Ag, through-connected network"},
    {"name": "low_fill_all_axis", "start": (36, 28, 32), "shape": (24, 24, 16), "motif": "low-volume-fraction network spanning x/y/z"},
    {"name": "high_fill_all_axis", "start": (28, 0, 28), "shape": (24, 24, 16), "motif": "high-volume-fraction network spanning x/y/z"},
    {"name": "isolated_islands_1", "start": (80, 56, 38), "shape": (12, 12, 8), "motif": "disconnected Ag islands; no axis spans"},
    {"name": "isolated_islands_2", "start": (0, 88, 34), "shape": (12, 12, 8), "motif": "disconnected Ag islands; no axis spans"},
)
COMPONENT_GROUPS = {
    "component_lumped_3x3x3": (3, 3, 3),
    "component_lumped_2x2x2": (2, 2, 2),
    "component_lumped_1x1x2": (1, 1, 2),
    "component_lumped_1x2x1": (1, 2, 1),
    "component_lumped_2x1x1": (2, 1, 1),
}
INTRINSIC_K = 1.0


def _boundary_face_fraction(source_ag: np.ndarray, size_A: np.ndarray, target_shape: tuple[int, ...], side: int) -> np.ndarray:
    spacing = np.asarray(size_A, dtype=float) / np.asarray(target_shape, dtype=float)
    wy = _overlap_weights(source_ag.shape[1], target_shape[1], float(size_A[1]))
    wz = _overlap_weights(source_ag.shape[2], target_shape[2], float(size_A[2]))
    return wy @ source_ag[0 if side == 0 else -1].astype(float) @ wz.T


def effective_ag_conductance(
    source_ag: np.ndarray,
    size_A: np.ndarray,
    ag_fraction: np.ndarray,
    spacing_A: np.ndarray,
    contacts: tuple[np.ndarray, np.ndarray, np.ndarray],
    transport_axis: int = 0,
) -> dict[str, float | int]:
    """Return the effective hot-to-cold conductance for a nonperiodic ROI."""
    if transport_axis not in (0, 1, 2):
        raise ValueError("transport_axis must be 0, 1, or 2")
    if transport_axis != 0:
        order = [transport_axis] + [axis for axis in range(3) if axis != transport_axis]
        source_ag = np.transpose(source_ag, order)
        ag_fraction = np.transpose(ag_fraction, order)
        size_A = np.asarray(size_A)[order]
        spacing_A = np.asarray(spacing_A)[order]
        contacts = tuple(np.transpose(contacts[old_axis], order) for old_axis in order)
    ag = np.asarray(ag_fraction, dtype=float)
    spacing = np.asarray(spacing_A, dtype=float)
    size = np.asarray(size_A, dtype=float)
    shape = ag.shape
    active_full = np.flatnonzero(ag.ravel() > 1.0e-12)
    full_to_active = np.full(ag.size, -1, dtype=np.int64)
    full_to_active[active_full] = np.arange(active_full.size)
    indices = np.arange(ag.size, dtype=np.int64).reshape(shape)
    volume = float(np.prod(spacing))
    edge_i: list[np.ndarray] = []
    edge_j: list[np.ndarray] = []
    edge_g: list[np.ndarray] = []
    for axis in range(3):
        low = [slice(None)] * 3
        high = [slice(None)] * 3
        low[axis], high[axis] = slice(0, -1), slice(1, None)
        low_t, high_t = tuple(low), tuple(high)
        contact = np.asarray(contacts[axis][low_t], dtype=float)
        ids_low = indices[low_t].ravel()
        ids_high = indices[high_t].ravel()
        pair_active = (ag.ravel()[ids_low] > 1.0e-12) & (ag.ravel()[ids_high] > 1.0e-12)
        g = INTRINSIC_K * contact.ravel() * volume / spacing[axis] ** 2
        valid = pair_active & (g > 0.0)
        edge_i.append(full_to_active[ids_low[valid]])
        edge_j.append(full_to_active[ids_high[valid]])
        edge_g.append(g[valid])
    i = np.concatenate(edge_i) if edge_i else np.empty(0, dtype=np.int64)
    j = np.concatenate(edge_j) if edge_j else np.empty(0, dtype=np.int64)
    g = np.concatenate(edge_g) if edge_g else np.empty(0, dtype=float)

    hot_fraction = _boundary_face_fraction(source_ag, size, shape, 0)
    cold_fraction = _boundary_face_fraction(source_ag, size, shape, 1)
    face_area = spacing[1] * spacing[2]
    hot_power = np.zeros(active_full.size, dtype=float)
    cold_power = np.zeros(active_full.size, dtype=float)
    hot_ids = indices[0, :, :].ravel()
    cold_ids = indices[-1, :, :].ravel()
    hot_g = 2.0 * INTRINSIC_K * hot_fraction.ravel() * face_area / spacing[0]
    cold_g = 2.0 * INTRINSIC_K * cold_fraction.ravel() * face_area / spacing[0]
    hot_valid = (hot_g > 0.0) & (ag.ravel()[hot_ids] > 1.0e-12)
    cold_valid = (cold_g > 0.0) & (ag.ravel()[cold_ids] > 1.0e-12)
    hot_active = full_to_active[hot_ids[hot_valid]]
    cold_active = full_to_active[cold_ids[cold_valid]]
    hot_power[hot_active] = hot_g[hot_valid]
    cold_power[cold_active] = cold_g[cold_valid]

    diagonal = np.bincount(i, weights=g, minlength=active_full.size)
    diagonal += np.bincount(j, weights=g, minlength=active_full.size)
    diagonal += hot_power + cold_power
    rows = np.concatenate((np.arange(active_full.size), i, j))
    cols = np.concatenate((np.arange(active_full.size), j, i))
    values = np.concatenate((diagonal, -g, -g))
    matrix = coo_matrix((values, (rows, cols)), shape=(active_full.size, active_full.size)).tocsr()
    rhs = hot_power  # T_hot=1, T_cold=0
    preconditioner = coo_matrix(
        (1.0 / np.maximum(diagonal, 1.0e-30), (np.arange(diagonal.size), np.arange(diagonal.size))),
        shape=matrix.shape,
    ).tocsr()
    temperature, info = cg(matrix, rhs, M=preconditioner, rtol=1.0e-12, atol=1.0e-14, maxiter=20_000)
    if info != 0:
        raise RuntimeError(f"QSGS conductance CG failed (info={info})")
    hot_current = float(np.dot(hot_power, 1.0 - temperature))
    cold_current = float(np.dot(cold_power, temperature))
    effective_k_fraction = hot_current * size[0] / float(np.prod(size[1:]))
    if abs(effective_k_fraction) < 1.0e-12:
        hot_current = 0.0
        cold_current = 0.0
        effective_k_fraction = 0.0
    return {
        "active_cells": int(active_full.size),
        "conducting_faces": int(g.size),
        "hot_face_ag_area_fraction": float(hot_fraction.mean()),
        "cold_face_ag_area_fraction": float(cold_fraction.mean()),
        "effective_conductance_over_k": hot_current,
        "effective_k_over_intrinsic_k": effective_k_fraction,
        "hot_cold_flux_relative_balance": abs(hot_current - cold_current) / max(abs(hot_current), 1.0e-30),
        "linear_residual_linf": float(np.max(np.abs(matrix @ temperature - rhs))),
    }


def run() -> dict:
    with np.load(SOURCE, allow_pickle=False) as data:
        full_skeleton = np.asarray(data["skeleton"], dtype=bool)
        full_size_A = np.asarray(data["size_A"], dtype=float)
        seed = int(data["seed"])
    results: dict[str, object] = {
        "status": "QSGS_LOCAL_GEOMETRY_TRANSPORT_DIAGNOSTIC",
        "source": SOURCE.name,
        "source_seed": seed,
        "source_shape_xyz": list(full_skeleton.shape),
        "roi_cases": [{"name": case["name"], "start": case["start"], "shape": case["shape"], "motif": case["motif"]} for case in ROI_CASES],
        "component_graph_groupings": {name: list(groups) for name, groups in COMPONENT_GROUPS.items()},
        "intrinsic_k_is_arbitrary_normalized_value": INTRINSIC_K,
        "boundary_conditions": "unit temperature difference on opposing faces of the tested x/y/z axis; remaining faces insulated; no periodic wrap",
        "interpretation": "geometry-only resistor-network comparison; not interface transport or a measured composite conductivity",
        "regions": {},
    }
    for region_index, case in enumerate(ROI_CASES, start=1):
        roi_start = case["start"]
        roi_shape = case["shape"]
        slices = tuple(
            slice(start, start + count) for start, count in zip(roi_start, roi_shape)
        )
        skeleton = full_skeleton[slices]
        source_ag = ~skeleton
        size_A = full_size_A * np.asarray(roi_shape) / np.asarray(full_skeleton.shape)
        voxel_spacing = size_A / np.asarray(roi_shape)
        voxel_contacts = exact_face_contact_fractions(
            source_ag, source_ag.astype(float), size_A, voxel_spacing
        )
        component_graphs = {
            name: build_component_graph(source_ag, size_A, coarse_voxels=groups, intrinsic_k=INTRINSIC_K)
            for name, groups in COMPONENT_GROUPS.items()
        }
        region = {
            "motif": case["motif"],
            "roi_start_voxels_xyz": list(roi_start),
            "roi_shape_voxels_xyz": list(roi_shape),
            "roi_size_A": size_A.tolist(),
            "source_voxel_spacing_A": voxel_spacing.tolist(),
            "roi_ag_volume_fraction": float(source_ag.mean()),
            "directions": {},
            "component_graphs": {
                name: {
                    "coarse_voxels": list(COMPONENT_GROUPS[name]),
                    "component_nodes": int(graph.voxel_count.size),
                    "source_ag_voxels": int(source_ag.sum()),
                    "preserved_ag_voxel_fraction": float(graph.voxel_count.sum() / source_ag.sum()),
                    "edge_count": int(graph.edge_conductance.size),
                    "node_capacity_voxel_counts_min": int(graph.voxel_count.min()),
                    "node_capacity_voxel_counts_median": float(np.median(graph.voxel_count)),
                    "node_capacity_voxel_counts_max": int(graph.voxel_count.max()),
                }
                for name, graph in component_graphs.items()
            },
        }
        projected: dict[str, tuple[np.ndarray, np.ndarray, tuple[np.ndarray, ...], list[int]]] = {}
        for label, spacing_A in (("projected_10A", 10.0), ("projected_5A", 5.0)):
            shape = tuple(int(v) for v in np.ceil(size_A / spacing_A))
            projection = project_phase_mask(skeleton, size_A, spacing_A, shape)
            contacts = exact_face_contact_fractions(
                source_ag, projection.ag_fraction, size_A, projection.spacing_A
            )
            projected[label] = (
                projection.ag_fraction,
                projection.spacing_A,
                contacts,
                list(shape),
            )
        for axis, axis_name in enumerate("xyz"):
            direction: dict[str, object] = {}
            direction["source_voxel_reference"] = effective_ag_conductance(
                source_ag, size_A, source_ag.astype(float), voxel_spacing, voxel_contacts, axis
            )
            for name, graph in component_graphs.items():
                direction[name] = component_effective_conductance(
                    graph, axis, intrinsic_k=INTRINSIC_K
                )
            for label, (fraction, spacing, contacts, shape) in projected.items():
                result = effective_ag_conductance(
                    source_ag, size_A, fraction, spacing, contacts, axis
                )
                result["grid_shape_xyz"] = shape
                result["actual_spacing_A_xyz"] = spacing.tolist()
                result["ag_volume_fraction"] = float(fraction.mean())
                direction[label] = result
            reference = direction["source_voxel_reference"]["effective_k_over_intrinsic_k"]
            for label in (*COMPONENT_GROUPS, "projected_10A", "projected_5A"):
                result = direction[label]
                result["relative_effective_k_error_vs_voxel_reference"] = (
                    None if reference == 0.0 else
                    (result["effective_k_over_intrinsic_k"] - reference) / reference
                )
            region["directions"][axis_name] = direction
        results["regions"][f"roi_{region_index}"] = region
    return results


def main() -> None:
    report = run()
    output = HERE / "nested_grid_diagnostics" / "qsgs_local_transport_benchmark.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
