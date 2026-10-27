"""Connectivity metrics for liquid Ag on the native QSGS voxel geometry."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import label
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def liquid_cluster_metrics(
    liquid_mask: np.ndarray,
    spacing_A: np.ndarray | tuple[float, float, float],
    periodic_axes: tuple[int, ...] = (0, 1),
) -> dict[str, int | float | list[float] | None]:
    """Measure six-face liquid clusters, including periodic x/y seam links.

    ``liquid_mask`` must identify liquid Ag voxels; scaffold voxels must be
    false. The returned centroid is the volume-weighted centroid of the
    largest cluster. Periodic coordinates use circular means so a cluster
    crossing a lateral seam remains near that seam.
    """
    mask = np.asarray(liquid_mask, dtype=bool)
    spacing = np.asarray(spacing_A, dtype=float)
    periodic = tuple(sorted(set(int(axis) for axis in periodic_axes)))
    if mask.ndim != 3 or any(n <= 0 for n in mask.shape):
        raise ValueError("liquid_mask must be a non-empty three-dimensional field")
    if spacing.shape != (3,) or np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0):
        raise ValueError("spacing_A must contain three finite positive lengths")
    if any(axis not in (0, 1, 2) for axis in periodic):
        raise ValueError("periodic_axes may contain only 0, 1, or 2")

    active_voxels = int(mask.sum())
    if active_voxels == 0:
        return {
            "connected_liquid_cluster_count": 0,
            "liquid_voxel_count": 0,
            "largest_liquid_cluster_voxels": 0,
            "largest_liquid_cluster_volume_A3": 0.0,
            "largest_liquid_cluster_centroid_A": None,
        }

    structure = np.zeros((3, 3, 3), dtype=bool)
    structure[1, 1, 1] = True
    for axis in range(3):
        low = [1, 1, 1]
        high = [1, 1, 1]
        low[axis] = 0
        high[axis] = 2
        structure[tuple(low)] = True
        structure[tuple(high)] = True
    labels, local_count = label(mask, structure=structure)

    seam_pairs = []
    for axis in periodic:
        if mask.shape[axis] <= 1:
            continue
        first = np.take(labels, 0, axis=axis).ravel()
        last = np.take(labels, -1, axis=axis).ravel()
        valid = (first > 0) & (last > 0) & (first != last)
        if np.any(valid):
            seam_pairs.append(np.column_stack((first[valid], last[valid])))

    # Each local image label is a graph node. Only seam-spanning labels need
    # edges; isolated nodes remain individual connected components.
    if seam_pairs:
        pairs = np.unique(np.concatenate(seam_pairs, axis=0), axis=0)
        rows = np.concatenate((pairs[:, 0], pairs[:, 1]))
        cols = np.concatenate((pairs[:, 1], pairs[:, 0]))
        adjacency = coo_matrix(
            (np.ones(rows.size, dtype=np.uint8), (rows, cols)),
            shape=(local_count + 1, local_count + 1),
        ).tocsr()
        _, component_by_label = connected_components(adjacency, directed=False)
    else:
        component_by_label = np.arange(local_count + 1, dtype=np.int32)

    local_sizes = np.bincount(labels.ravel(), minlength=local_count + 1)
    label_components = component_by_label[1:]
    cluster_sizes = np.bincount(
        label_components,
        weights=local_sizes[1:],
        minlength=local_count + 1,
    ).astype(np.int64)
    occupied_components = np.flatnonzero(cluster_sizes)
    largest_component = int(occupied_components[np.argmax(cluster_sizes[occupied_components])])
    largest_size = int(cluster_sizes[largest_component])

    # Aggregate coordinate moments from local image labels to periodic-aware
    # cluster IDs without allocating a full three-coordinate voxel array.
    label_sums: list[np.ndarray] = []
    for axis in range(3):
        indices = np.arange(mask.shape[axis], dtype=float)
        shape = [1, 1, 1]
        shape[axis] = mask.shape[axis]
        coordinate = (indices.reshape(shape) + 0.5) * spacing[axis]
        if axis in periodic:
            angle = 2.0 * np.pi * coordinate / (mask.shape[axis] * spacing[axis])
            sin_local = np.bincount(
                labels.ravel(), weights=np.broadcast_to(np.sin(angle), mask.shape).ravel(), minlength=local_count + 1
            )
            cos_local = np.bincount(
                labels.ravel(), weights=np.broadcast_to(np.cos(angle), mask.shape).ravel(), minlength=local_count + 1
            )
            sin_cluster = np.bincount(label_components, weights=sin_local[1:], minlength=local_count + 1)
            cos_cluster = np.bincount(label_components, weights=cos_local[1:], minlength=local_count + 1)
            length = mask.shape[axis] * spacing[axis]
            centroid = (np.arctan2(sin_cluster, cos_cluster) % (2.0 * np.pi)) * length / (2.0 * np.pi)
        else:
            local_sum = np.bincount(
                labels.ravel(),
                weights=np.broadcast_to(coordinate, mask.shape).ravel(),
                minlength=local_count + 1,
            )
            cluster_sum = np.bincount(label_components, weights=local_sum[1:], minlength=local_count + 1)
            centroid = np.divide(
                cluster_sum,
                cluster_sizes,
                out=np.zeros_like(cluster_sum),
                where=cluster_sizes > 0,
            )
        label_sums.append(centroid)

    centroid_A = [float(values[largest_component]) for values in label_sums]
    return {
        "connected_liquid_cluster_count": int(occupied_components.size),
        "liquid_voxel_count": active_voxels,
        "largest_liquid_cluster_voxels": largest_size,
        "largest_liquid_cluster_volume_A3": float(largest_size * np.prod(spacing)),
        "largest_liquid_cluster_centroid_A": centroid_A,
    }
