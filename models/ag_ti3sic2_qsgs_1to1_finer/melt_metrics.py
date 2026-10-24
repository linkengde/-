"""Offline/online metrics feeding the adaptive Ag melt-pool controller.

Liquid flags use an Ag CNA unknown label plus displacement over a fixed
window. The default thresholds are calibrated against the repository Ag EAM
solid/liquid smoke trajectories; they remain provisional for the composite.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class ClusterMetrics:
    member_count: int
    largest_cluster_atoms: int
    largest_cluster_mask: np.ndarray


def classify_liquid_ag(
    cna_code: np.ndarray,
    displacement_A: np.ndarray,
    cna_unknown_code: int = 5,
    minimum_displacement_A: float = 1.0,
) -> np.ndarray:
    """Require both disordered CNA and ≥threshold displacement over ~1 ps."""
    cna = np.asarray(cna_code)
    displacement = np.asarray(displacement_A, dtype=float)
    if cna.shape != displacement.shape:
        raise ValueError("CNA and displacement arrays must have the same shape")
    if not np.all(np.isfinite(displacement)) or np.any(displacement < 0.0):
        raise ValueError("Displacements must be finite and non-negative")
    if minimum_displacement_A <= 0.0:
        raise ValueError("minimum_displacement_A must be positive")
    return (cna == cna_unknown_code) & (displacement >= minimum_displacement_A)


def largest_connected_cluster(
    xyz_A: np.ndarray,
    members: np.ndarray,
    box_A: tuple[float, float, float],
    cutoff_A: float,
    periodic_xy: bool = True,
) -> ClusterMetrics:
    """Return the largest geometric cluster under a cutoff graph.

    z is always non-periodic. x/y are minimum-image periodic when requested.
    scipy's tree is built only from selected members, which keeps liquid-pool
    analysis tractable when the molten region is local.
    """
    xyz = np.asarray(xyz_A, dtype=float)
    mask = np.asarray(members, dtype=bool)
    box = np.asarray(box_A, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or mask.shape != (len(xyz),):
        raise ValueError("xyz_A and members have incompatible shapes")
    if not np.all(np.isfinite(xyz)) or box.shape != (3,) or np.any(box <= 0.0):
        raise ValueError("Coordinates and box dimensions must be finite and positive")
    if not np.isfinite(cutoff_A) or cutoff_A <= 0.0:
        raise ValueError("cutoff_A must be finite and positive")
    selected = np.flatnonzero(mask)
    full_cluster = np.zeros(len(xyz), dtype=bool)
    if selected.size == 0:
        return ClusterMetrics(0, 0, full_cluster)
    if selected.size == 1:
        full_cluster[selected[0]] = True
        return ClusterMetrics(1, 1, full_cluster)

    points = xyz[selected].copy()
    if periodic_xy:
        points[:, 0] %= box[0]
        points[:, 1] %= box[1]
        zmin = float(points[:, 2].min())
        zmax = float(points[:, 2].max())
        z_period = max(zmax - zmin + 2.0 * cutoff_A + 1.0e-6, 4.0 * cutoff_A)
        points[:, 2] = points[:, 2] - zmin + cutoff_A
        tree_box = np.array([box[0], box[1], z_period])
    else:
        shift = points.min(axis=0)
        points -= shift
        z_period = max(float(points[:, 2].max()) + 2.0 * cutoff_A + 1.0e-6, 4.0 * cutoff_A)
        points[:, 2] += cutoff_A
        tree_box = np.array([max(float(points[:, 0].max()) + 2.0 * cutoff_A, 4.0 * cutoff_A),
                             max(float(points[:, 1].max()) + 2.0 * cutoff_A, 4.0 * cutoff_A),
                             z_period])
        points[:, 0] += cutoff_A
        points[:, 1] += cutoff_A

    pairs = cKDTree(points, boxsize=tree_box).query_pairs(cutoff_A, output_type="ndarray")
    nmember = selected.size
    if pairs.size == 0:
        full_cluster[selected[0]] = True
        return ClusterMetrics(nmember, 1, full_cluster)
    edge = np.ones(pairs.shape[0] * 2, dtype=np.int8)
    row = np.concatenate((pairs[:, 0], pairs[:, 1]))
    col = np.concatenate((pairs[:, 1], pairs[:, 0]))
    graph = coo_matrix((edge, (row, col)), shape=(nmember, nmember)).tocsr()
    _, labels = connected_components(graph, directed=False, return_labels=True)
    counts = np.bincount(labels)
    largest_label = int(np.argmax(counts))
    local_mask = labels == largest_label
    full_cluster[selected[local_mask]] = True
    return ClusterMetrics(int(nmember), int(counts[largest_label]), full_cluster)


def surface_recession_from_height(
    initial_top_A: float,
    current_top_A: float,
) -> float:
    """Positive value means retreat downward from the initial top surface."""
    if not np.isfinite(initial_top_A) or not np.isfinite(current_top_A):
        raise ValueError("Surface heights must be finite")
    return float(initial_top_A - current_top_A)
