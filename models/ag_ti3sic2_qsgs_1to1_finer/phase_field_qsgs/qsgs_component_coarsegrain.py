"""Connectivity-preserving component graph for coarse QSGS heat transport.

Each coarse cell is split into its native-voxel 6-connected Ag components.
Only components that touch through a real Ag voxel face are connected across
coarse-cell boundaries. A component is isothermal in this first model; its
thermal capacity is proportional to its exact Ag voxel volume.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import label
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import cg


_SIX_NEIGHBOR = np.zeros((3, 3, 3), dtype=np.uint8)
for _axis in range(3):
    _low = [1, 1, 1]
    _high = [1, 1, 1]
    _low[_axis] = 0
    _high[_axis] = 2
    _SIX_NEIGHBOR[tuple(_low)] = 1
    _SIX_NEIGHBOR[tuple(_high)] = 1


@dataclass(frozen=True)
class ComponentGraph:
    component_id: np.ndarray  # Native voxel -> component node, -1 for scaffold.
    voxel_count: np.ndarray
    block_index: np.ndarray
    edge_i: np.ndarray
    edge_j: np.ndarray
    edge_axis: np.ndarray
    edge_conductance: np.ndarray
    boundary_area: tuple[tuple[np.ndarray, np.ndarray], ...]  # axis, low/high node areas
    boundary_half_width_A: tuple[tuple[float, float], ...]
    source_size_A: np.ndarray
    source_spacing_A: np.ndarray
    coarse_voxels: tuple[int, int, int]
    periodic_axes: tuple[int, ...]


def build_component_graph(
    source_ag: np.ndarray,
    size_A: np.ndarray,
    coarse_voxels: tuple[int, int, int] = (3, 3, 3),
    intrinsic_k: float = 1.0,
    periodic_axes: tuple[int, ...] = (),
) -> ComponentGraph:
    """Partition each coarse cell into Ag components and build real face links."""
    ag = np.asarray(source_ag, dtype=bool)
    size = np.asarray(size_A, dtype=float)
    groups = np.asarray(coarse_voxels, dtype=int)
    if ag.ndim != 3 or any(n <= 0 for n in ag.shape):
        raise ValueError("source_ag must be a non-empty 3D mask")
    if size.shape != (3,) or np.any(~np.isfinite(size)) or np.any(size <= 0.0):
        raise ValueError("size_A must contain three finite positive lengths")
    if groups.shape != (3,) or np.any(groups <= 0):
        raise ValueError("coarse_voxels must contain three positive group sizes")
    if not np.isfinite(intrinsic_k) or intrinsic_k <= 0.0:
        raise ValueError("intrinsic_k must be finite and positive")
    periodic = tuple(sorted(set(int(axis) for axis in periodic_axes)))
    if any(axis not in (0, 1, 2) for axis in periodic):
        raise ValueError("periodic_axes may contain only 0, 1, or 2")

    pitch = size / np.asarray(ag.shape, dtype=float)
    component_id = np.full(ag.shape, -1, dtype=np.int64)
    component_block: np.ndarray | list[tuple[int, int, int]] = []
    voxel_count: np.ndarray | list[int] = []
    block_counts = np.ceil(np.asarray(ag.shape) / groups).astype(int)
    block_widths = [
        np.minimum(groups[axis], ag.shape[axis] - np.arange(block_counts[axis]) * groups[axis])
        * pitch[axis]
        for axis in range(3)
    ]
    block_centers = [
        (np.arange(block_counts[axis]) * groups[axis] + block_widths[axis] / (2.0 * pitch[axis]))
        * pitch[axis]
        for axis in range(3)
    ]

    node_count = 0
    if tuple(groups) == (1, 1, 1):
        # Native-voxel reference: each occupied voxel is its own thermal node.
        component_id[ag] = np.arange(int(ag.sum()), dtype=np.int64)
        node_count = int(ag.sum())
        voxel_count = np.ones(node_count, dtype=np.int64)
        component_block = np.argwhere(ag).astype(np.int32, copy=False)
    elif tuple(groups) == (2, 1, 1):
        # The production coarse grouping pairs only adjacent x voxels. A pair
        # is either empty or one 6-connected Ag component, so vectorize rather
        # than calling ndimage.label millions of times on tiny blocks.
        nx, ny, nz = ag.shape
        slots = np.zeros((block_counts[0], ny, nz, 2), dtype=bool)
        slots[..., 0] = ag[0::2, :, :]
        if nx > 1:
            slots[: ag[1::2, :, :].shape[0], :, :, 1] = ag[1::2, :, :]
        block_active = slots.any(axis=-1)
        block_id = np.full(block_active.shape, -1, dtype=np.int64)
        block_id[block_active] = np.arange(int(block_active.sum()), dtype=np.int64)
        component_id[0::2, :, :] = np.where(slots[..., 0], block_id, -1)
        if nx > 1:
            component_id[1::2, :, :] = np.where(slots[: component_id[1::2].shape[0], ..., 1], block_id[: component_id[1::2].shape[0]], -1)
        counts = slots.sum(axis=-1)
        voxel_count = counts[block_active].astype(np.int64, copy=False)
        component_block = np.argwhere(block_active).astype(np.int32, copy=False)
        node_count = int(voxel_count.size)
    else:
        for bx in range(block_counts[0]):
            for by in range(block_counts[1]):
                for bz in range(block_counts[2]):
                    starts = np.asarray([bx, by, bz]) * groups
                    stops = np.minimum(starts + groups, np.asarray(ag.shape))
                    region = tuple(slice(int(starts[a]), int(stops[a])) for a in range(3))
                    local_labels, nlocal = label(ag[region], structure=_SIX_NEIGHBOR)
                    if nlocal == 0:
                        continue
                    active = local_labels > 0
                    component_id[region][active] = local_labels[active] + node_count - 1
                    counts = np.bincount(local_labels.ravel(), minlength=nlocal + 1)[1:]
                    voxel_count.extend(int(value) for value in counts)
                    component_block.extend([(bx, by, bz)] * nlocal)
                    node_count += nlocal

    if node_count == 0:
        raise ValueError("QSGS source region contains no Ag voxels")

    edge_i_parts: list[np.ndarray] = []
    edge_j_parts: list[np.ndarray] = []
    edge_axis_parts: list[np.ndarray] = []
    edge_g_parts: list[np.ndarray] = []
    blocks = np.asarray(component_block, dtype=int)
    for axis in range(3):
        low = [slice(None)] * 3
        high = [slice(None)] * 3
        low[axis], high[axis] = slice(0, -1), slice(1, None)
        low_t, high_t = tuple(low), tuple(high)
        ids_low = component_id[low_t].ravel()
        ids_high = component_id[high_t].ravel()
        valid = (ids_low >= 0) & (ids_high >= 0) & (ids_low != ids_high)
        voxel_face_area = float(np.prod(np.delete(pitch, axis)))
        if np.any(valid):
            pairs = np.column_stack((ids_low[valid], ids_high[valid]))
            # Native voxels have at most one face per node pair and direction.
            if tuple(groups) == (1, 1, 1):
                unique_pairs = pairs
                counts = np.ones(len(pairs), dtype=np.int64)
            else:
                # Multiple native faces can connect the same pair of component nodes.
                pairs.sort(axis=1)
                unique_pairs, counts = np.unique(pairs, axis=0, return_counts=True)
            contact_area = counts.astype(float) * voxel_face_area
            block_i = blocks[unique_pairs[:, 0], axis]
            block_j = blocks[unique_pairs[:, 1], axis]
            separation = np.abs(block_centers[axis][block_j] - block_centers[axis][block_i])
            if np.any(separation <= 0.0):
                raise RuntimeError("A component edge did not cross a coarse-cell boundary")
            conductance = intrinsic_k * contact_area / separation
            edge_i_parts.append(unique_pairs[:, 0])
            edge_j_parts.append(unique_pairs[:, 1])
            edge_axis_parts.append(np.full(unique_pairs.shape[0], axis, dtype=np.int8))
            edge_g_parts.append(conductance)

        if axis in periodic and ag.shape[axis] > 1:
            ids_first = np.take(component_id, 0, axis=axis).ravel()
            ids_last = np.take(component_id, -1, axis=axis).ravel()
            seam_valid = (ids_first >= 0) & (ids_last >= 0) & (ids_first != ids_last)
            if np.any(seam_valid):
                seam_pairs = np.column_stack((ids_first[seam_valid], ids_last[seam_valid]))
                if tuple(groups) == (1, 1, 1):
                    seam_unique = seam_pairs
                    seam_counts = np.ones(len(seam_pairs), dtype=np.int64)
                else:
                    seam_pairs.sort(axis=1)
                    seam_unique, seam_counts = np.unique(seam_pairs, axis=0, return_counts=True)
                seam_area = seam_counts.astype(float) * voxel_face_area
                seam_block_i = blocks[seam_unique[:, 0], axis]
                seam_block_j = blocks[seam_unique[:, 1], axis]
                seam_distance = (
                    block_centers[axis][0]
                    + size[axis]
                    - block_centers[axis][-1]
                )
                if seam_distance <= 0.0:
                    raise RuntimeError("Periodic component centres have invalid separation")
                edge_i_parts.append(seam_unique[:, 0])
                edge_j_parts.append(seam_unique[:, 1])
                edge_axis_parts.append(np.full(seam_unique.shape[0], axis, dtype=np.int8))
                edge_g_parts.append(intrinsic_k * seam_area / seam_distance)

    if edge_i_parts:
        edge_i = np.concatenate(edge_i_parts).astype(np.int64, copy=False)
        edge_j = np.concatenate(edge_j_parts).astype(np.int64, copy=False)
        edge_axis = np.concatenate(edge_axis_parts)
        edge_conductance = np.concatenate(edge_g_parts)
    else:
        edge_i = np.empty(0, dtype=np.int64)
        edge_j = np.empty(0, dtype=np.int64)
        edge_axis = np.empty(0, dtype=np.int8)
        edge_conductance = np.empty(0, dtype=float)

    boundary_area: list[tuple[np.ndarray, np.ndarray]] = []
    boundary_half_width: list[tuple[float, float]] = []
    for axis in range(3):
        sides: list[np.ndarray] = []
        for side in (0, -1):
            face_ids = np.take(component_id, side, axis=axis).ravel()
            face_ids = face_ids[face_ids >= 0]
            area_per_voxel = float(np.prod(np.delete(pitch, axis)))
            sides.append(np.bincount(face_ids, minlength=node_count).astype(float) * area_per_voxel)
        boundary_area.append((sides[0], sides[1]))
        half_widths = (
            float(block_widths[axis][0] / 2.0),
            float(block_widths[axis][-1] / 2.0),
        )
        boundary_half_width.append(half_widths)

    return ComponentGraph(
        component_id=component_id,
        voxel_count=np.asarray(voxel_count, dtype=np.int64),
        block_index=np.asarray(component_block, dtype=np.int32),
        edge_i=edge_i,
        edge_j=edge_j,
        edge_axis=edge_axis,
        edge_conductance=edge_conductance,
        boundary_area=tuple(boundary_area),
        boundary_half_width_A=tuple(boundary_half_width),
        source_size_A=size,
        source_spacing_A=pitch,
        coarse_voxels=tuple(int(v) for v in groups),
        periodic_axes=periodic,
    )


def component_effective_conductance(
    graph: ComponentGraph,
    axis: int,
    intrinsic_k: float = 1.0,
    rtol: float = 1.0e-12,
) -> dict[str, float | int]:
    """Solve steady conduction between opposite faces of a component graph."""
    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1, or 2")
    nnode = int(graph.voxel_count.size)
    edge_mask = graph.edge_axis >= 0
    i, j, g = graph.edge_i[edge_mask], graph.edge_j[edge_mask], graph.edge_conductance[edge_mask]
    hot_area, cold_area = graph.boundary_area[axis]
    half_hot, half_cold = graph.boundary_half_width_A[axis]
    hot_g = intrinsic_k * hot_area / half_hot
    cold_g = intrinsic_k * cold_area / half_cold

    adjacency = coo_matrix(
        (np.ones(2 * len(i)), (np.concatenate((i, j)), np.concatenate((j, i)))),
        shape=(nnode, nnode),
    ).tocsr()
    _, labels = connected_components(adjacency, directed=False, return_labels=True)
    boundary_nodes = np.flatnonzero((hot_g > 0.0) | (cold_g > 0.0))
    if boundary_nodes.size == 0:
        return {
            "component_nodes": nnode,
            "active_boundary_nodes": 0,
            "conducting_edges": int(len(g)),
            "effective_conductance_over_k": 0.0,
            "effective_k_over_intrinsic_k": 0.0,
            "hot_cold_flux_relative_balance": 0.0,
            "linear_residual_linf": 0.0,
        }
    anchored_components = np.unique(labels[boundary_nodes])
    keep = np.isin(labels, anchored_components)
    old_to_new = np.full(nnode, -1, dtype=np.int64)
    old_to_new[keep] = np.arange(np.count_nonzero(keep))
    edge_keep = keep[i] & keep[j]
    i_new, j_new, g_new = old_to_new[i[edge_keep]], old_to_new[j[edge_keep]], g[edge_keep]
    hot_new, cold_new = hot_g[keep], cold_g[keep]
    n = int(np.count_nonzero(keep))

    diagonal = np.bincount(i_new, weights=g_new, minlength=n).astype(float)
    diagonal += np.bincount(j_new, weights=g_new, minlength=n)
    diagonal += hot_new + cold_new
    rhs = hot_new  # Thot=1, Tcold=0.
    rows = np.concatenate((np.arange(n), i_new, j_new))
    cols = np.concatenate((np.arange(n), j_new, i_new))
    values = np.concatenate((diagonal, -g_new, -g_new))
    matrix = coo_matrix((values, (rows, cols)), shape=(n, n)).tocsr()
    if not np.any(rhs):
        temperature = np.zeros(n, dtype=float)
        info = 0
    else:
        preconditioner = coo_matrix(
            (1.0 / np.maximum(diagonal, 1.0e-30), (np.arange(n), np.arange(n))),
            shape=(n, n),
        ).tocsr()
        temperature, info = cg(
            matrix,
            rhs,
            M=preconditioner,
            rtol=rtol,
            atol=1.0e-14,
            maxiter=max(500, 20 * n),
        )
    if info != 0:
        raise RuntimeError(f"Component graph conductance CG failed (info={info})")
    hot_flux = float(np.dot(hot_new, 1.0 - temperature))
    cold_flux = float(np.dot(cold_new, temperature))
    cross_section = float(np.prod(np.delete(graph.source_size_A, axis)))
    effective_k = hot_flux * float(graph.source_size_A[axis]) / cross_section
    if abs(effective_k / intrinsic_k) < 1.0e-12:
        hot_flux = 0.0
        cold_flux = 0.0
        effective_k = 0.0
    return {
        "component_nodes": nnode,
        "active_boundary_nodes": n,
        "conducting_edges": int(g_new.size),
        "effective_conductance_over_k": hot_flux,
        "effective_k_over_intrinsic_k": effective_k / intrinsic_k,
        "hot_cold_flux_relative_balance": abs(hot_flux - cold_flux) / max(abs(hot_flux), 1.0e-30),
        "linear_residual_linf": float(np.max(np.abs(matrix @ temperature - rhs))) if n else 0.0,
    }
