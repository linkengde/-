#!/usr/bin/env python3
"""Compare voxel-resolved and projected Ag conductance in one QSGS window.

This is a geometry-only steady resistor-network diagnostic. It uses constant
intrinsic Ag conductivity, insulated lateral boundaries, and unit temperature
difference on opposing x faces. It does not model interface resistance or
claim experimental thermal conductivity.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import cg

from qsgs_face_contacts import exact_face_contact_fractions
from qsgs_projection import _overlap_weights, project_phase_mask


HERE = Path(__file__).resolve().parent
MODEL_DIR = HERE.parent
SOURCE = MODEL_DIR / "qsgs_morphology.npz"
ROI_SHAPE = (24, 24, 16)
ROI_STARTS = ((36, 60, 18), (56, 12, 0), (64, 76, 0))
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
) -> dict[str, float | int]:
    """Return the effective hot-to-cold conductance for a nonperiodic ROI."""
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
    temperature, info = cg(matrix, rhs, M=preconditioner, rtol=1.0e-10, atol=1.0e-12, maxiter=20_000)
    if info != 0:
        raise RuntimeError(f"QSGS conductance CG failed (info={info})")
    hot_current = float(np.dot(hot_power, 1.0 - temperature))
    cold_current = float(np.dot(cold_power, temperature))
    effective_k_fraction = hot_current * size[0] / float(np.prod(size[1:]))
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
        "roi_shape_voxels_xyz": list(ROI_SHAPE),
        "intrinsic_k_is_arbitrary_normalized_value": INTRINSIC_K,
        "boundary_conditions": "unit temperature difference on opposing x faces; all other faces insulated; no periodic wrap",
        "interpretation": "geometry-only resistor-network comparison; not interface transport or a measured composite conductivity",
        "regions": {},
    }
    for region_index, roi_start in enumerate(ROI_STARTS, start=1):
        slices = tuple(
            slice(start, start + count) for start, count in zip(roi_start, ROI_SHAPE)
        )
        skeleton = full_skeleton[slices]
        source_ag = ~skeleton
        size_A = full_size_A * np.asarray(ROI_SHAPE) / np.asarray(full_skeleton.shape)
        voxel_spacing = size_A / np.asarray(ROI_SHAPE)
        voxel_contacts = exact_face_contact_fractions(
            source_ag, source_ag.astype(float), size_A, voxel_spacing
        )
        region = {
            "roi_start_voxels_xyz": list(roi_start),
            "roi_size_A": size_A.tolist(),
            "source_voxel_spacing_A": voxel_spacing.tolist(),
            "roi_ag_volume_fraction": float(source_ag.mean()),
            "grids": {},
        }
        region["grids"]["source_voxel_reference"] = effective_ag_conductance(
            source_ag, size_A, source_ag.astype(float), voxel_spacing, voxel_contacts
        )
        for label, spacing_A in (("projected_10A", 10.0), ("projected_5A", 5.0)):
            shape = tuple(int(v) for v in np.ceil(size_A / spacing_A))
            projection = project_phase_mask(skeleton, size_A, spacing_A, shape)
            contacts = exact_face_contact_fractions(
                source_ag, projection.ag_fraction, size_A, projection.spacing_A
            )
            result = effective_ag_conductance(
                source_ag, size_A, projection.ag_fraction, projection.spacing_A, contacts
            )
            result["grid_shape_xyz"] = list(shape)
            result["actual_spacing_A_xyz"] = projection.spacing_A.tolist()
            result["ag_volume_fraction"] = float(projection.ag_fraction.mean())
            region["grids"][label] = result
        reference = region["grids"]["source_voxel_reference"]["effective_k_over_intrinsic_k"]
        for label in ("projected_10A", "projected_5A"):
            result = region["grids"][label]
            result["relative_effective_k_error_vs_voxel_reference"] = (
                result["effective_k_over_intrinsic_k"] - reference
            ) / reference
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
