"""Conservative projection of the stored QSGS voxel mask to thermal grids."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ProjectedMask:
    ag_fraction: np.ndarray
    skeleton_fraction: np.ndarray
    size_A: np.ndarray
    spacing_A: np.ndarray
    requested_spacing_A: float


def _overlap_weights(source_count: int, target_count: int, length_A: float) -> np.ndarray:
    """Return target-cell averages of piecewise-constant source voxels."""
    if source_count <= 0 or target_count <= 0 or not np.isfinite(length_A) or length_A <= 0:
        raise ValueError("Grid counts and physical length must be positive")
    source_edges = np.linspace(0.0, length_A, source_count + 1)
    target_edges = np.linspace(0.0, length_A, target_count + 1)
    source_lo = source_edges[:-1][None, :]
    source_hi = source_edges[1:][None, :]
    target_lo = target_edges[:-1][:, None]
    target_hi = target_edges[1:][:, None]
    overlap = np.maximum(0.0, np.minimum(target_hi, source_hi) - np.maximum(target_lo, source_lo))
    target_width = np.diff(target_edges)[:, None]
    weights = overlap / target_width
    if not np.allclose(weights.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-12):
        raise RuntimeError("Projection weights do not cover every target cell")
    return weights


def project_phase_mask(
    skeleton_mask: np.ndarray,
    size_A: np.ndarray,
    requested_spacing_A: float,
) -> ProjectedMask:
    """Map a binary skeleton voxel field to conservative Ag/ceramic fractions.

    The input axes are x, y, z, matching the QSGS generator. Each target cell
    spans the full material box along each axis; the returned grid therefore
    has uniform spacing no larger than ``requested_spacing_A``. This is a
    geometric volume projection only; it does not infer interface physics.
    """
    skeleton = np.asarray(skeleton_mask, dtype=bool)
    size = np.asarray(size_A, dtype=float)
    if skeleton.ndim != 3 or any(n <= 0 for n in skeleton.shape):
        raise ValueError("skeleton_mask must be a non-empty 3D array")
    if size.shape != (3,) or not np.all(np.isfinite(size)) or np.any(size <= 0.0):
        raise ValueError("size_A must contain three finite positive lengths")
    if not np.isfinite(requested_spacing_A) or requested_spacing_A <= 0.0:
        raise ValueError("requested_spacing_A must be finite and positive")

    source_shape = np.asarray(skeleton.shape, dtype=int)
    target_shape = np.ceil(size / requested_spacing_A).astype(int)
    spacing = size / target_shape
    wx, wy, wz = (
        _overlap_weights(int(source_shape[axis]), int(target_shape[axis]), float(size[axis]))
        for axis in range(3)
    )
    skeleton_values = skeleton.astype(float)
    projected = np.einsum("ai,bj,ck,ijk->abc", wx, wy, wz, skeleton_values, optimize=True)
    projected = np.clip(projected, 0.0, 1.0)
    ag_fraction = 1.0 - projected

    source_voxel_volume = float(np.prod(size / source_shape))
    target_cell_volume = float(np.prod(spacing))
    source_skeleton_volume = float(skeleton_values.sum() * source_voxel_volume)
    target_skeleton_volume = float(projected.sum() * target_cell_volume)
    source_ag_volume = float((1.0 - skeleton_values).sum() * source_voxel_volume)
    target_ag_volume = float(ag_fraction.sum() * target_cell_volume)
    tolerance = max(float(np.prod(size)) * 1.0e-12, 1.0e-9)
    if abs(source_skeleton_volume - target_skeleton_volume) > tolerance:
        raise RuntimeError("Skeleton volume was not conserved by grid projection")
    if abs(source_ag_volume - target_ag_volume) > tolerance:
        raise RuntimeError("Ag volume was not conserved by grid projection")
    return ProjectedMask(
        ag_fraction=ag_fraction,
        skeleton_fraction=projected,
        size_A=size.copy(),
        spacing_A=spacing,
        requested_spacing_A=float(requested_spacing_A),
    )
