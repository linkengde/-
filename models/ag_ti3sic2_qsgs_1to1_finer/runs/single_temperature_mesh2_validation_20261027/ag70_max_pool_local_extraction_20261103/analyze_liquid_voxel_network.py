#!/usr/bin/env python3
"""Topology-only graph analysis of the extracted liquid-Ag voxel cluster.

This computes geometric paths and solid-contact counts only. It deliberately
does not assign pressure, viscosity, conductance, velocity, or flow rate.
"""

from __future__ import annotations

import csv
import heapq
import json
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parent
MASKS = ROOT / "maximum_pool_local_windows.npz"
GEOMETRY_SHA256 = "5e509f16716e85464f604fcda2a53eaefd4f15e2d18d6c0bf24e0aa06c6656ce"
CHECKPOINT_SHA256 = "8b15a8b4766ad1d2fc8e5a6ec3d927e1fbe3728563ddc3f65546bd2a446406c6"


def shortest_spanning_path(mask: np.ndarray, spacing_nm: tuple[float, float, float], axis: int) -> dict:
    """Return shortest face-to-face liquid path through one connected mask."""
    coords = np.argwhere(mask)
    low = int(coords[:, axis].min())
    high = int(coords[:, axis].max())
    local_id = np.full(mask.shape, -1, dtype=np.int32)
    local_id[tuple(coords.T)] = np.arange(len(coords), dtype=np.int32)
    sources = np.flatnonzero(coords[:, axis] == low)
    targets = set(np.flatnonzero(coords[:, axis] == high).tolist())
    dist = np.full(len(coords), np.inf)
    parent = np.full(len(coords), -1, dtype=np.int32)
    heap = []
    for node in sources:
        dist[node] = 0.0
        heapq.heappush(heap, (0.0, int(node)))
    target = None
    steps = ((0, 1), (1, 1), (2, 1))
    while heap:
        d, node = heapq.heappop(heap)
        if d != dist[node]:
            continue
        if node in targets:
            target = node
            break
        x, y, z = coords[node]
        for ax, sign in steps:
            for delta in (-1, 1):
                nxt = [x, y, z]
                nxt[ax] += delta
                if nxt[ax] < 0 or nxt[ax] >= mask.shape[ax]:
                    continue
                other = int(local_id[tuple(nxt)])
                if other < 0:
                    continue
                nd = d + spacing_nm[ax]
                if nd < dist[other]:
                    dist[other] = nd
                    parent[other] = node
                    heapq.heappush(heap, (nd, other))
    if target is None:
        return {"spans": False, "shortest_path_length_nm": None, "straight_span_nm": None, "tortuosity": None, "path_edges_xyz": None}
    path = [target]
    while parent[path[-1]] >= 0:
        path.append(int(parent[path[-1]]))
    path.reverse()
    edge_counts = [0, 0, 0]
    for a, b in zip(path, path[1:]):
        edge_counts[int(np.flatnonzero(coords[a] != coords[b])[0])] += 1
    physical_length = float(sum(edge_counts[i] * spacing_nm[i] for i in range(3)))
    straight = float((high - low) * spacing_nm[axis])
    return {
        "spans": True,
        "shortest_path_length_nm": physical_length,
        "straight_span_nm": straight,
        "tortuosity": physical_length / straight if straight else 1.0,
        "path_edges_xyz": edge_counts,
    }


def paired_count(a: np.ndarray, b: np.ndarray) -> int:
    return int(np.count_nonzero(a & b))


def analyze_threshold(data: np.lib.npyio.NpzFile, threshold: float, spacing_nm: tuple[float, float, float]) -> dict:
    ag = data["pad8_ag_mask"].astype(bool)
    skeleton = data["pad8_skeleton_mask"].astype(bool)
    liquid_fraction = data["pad8_ag_liquid_fraction"]
    liquid = ag & (liquid_fraction >= threshold)
    structure = ndimage.generate_binary_structure(3, 1)
    labels, component_count = ndimage.label(liquid, structure=structure)
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    component_id = int(np.argmax(counts)) if component_count else 0
    largest = labels == component_id if component_id else np.zeros_like(liquid)
    coords = np.argwhere(largest)
    if coords.size == 0:
        raise ValueError(f"no liquid Ag voxels at threshold {threshold}")
    lo, hi = coords.min(axis=0), coords.max(axis=0)

    degree = np.zeros(largest.shape, dtype=np.uint8)
    edge_counts = []
    skeleton_contact_faces = []
    face_areas = [spacing_nm[1] * spacing_nm[2], spacing_nm[0] * spacing_nm[2], spacing_nm[0] * spacing_nm[1]]
    for axis in range(3):
        left = [slice(None)] * 3
        right = [slice(None)] * 3
        left[axis] = slice(0, -1)
        right[axis] = slice(1, None)
        left, right = tuple(left), tuple(right)
        internal_edges = largest[left] & largest[right]
        edge_count = int(np.count_nonzero(internal_edges))
        edge_counts.append(edge_count)
        degree[left] += internal_edges
        degree[right] += internal_edges
        contact = (largest[left] & skeleton[right]) | (skeleton[left] & largest[right])
        skeleton_contact_faces.append(int(np.count_nonzero(contact)))

    degree_hist = {str(k): int(np.count_nonzero(largest & (degree == k))) for k in range(7)}
    paths = [shortest_spanning_path(largest, spacing_nm, axis) for axis in range(3)]
    local_top_face = bool(np.any(largest[:, :, -1]))
    original_top_face = bool(np.any(largest[:, :, -1] & (data["pad8_source_z_index"][-1] == 49)))
    voxel_volume = float(np.prod(spacing_nm))
    face_area_nm2 = [float(n * area) for n, area in zip(skeleton_contact_faces, face_areas)]
    return {
        "liquid_fraction_threshold": threshold,
        "liquid_voxels_all_components": int(liquid.sum()),
        "six_face_components": int(component_count),
        "largest_component_voxels": int(coords.shape[0]),
        "largest_component_volume_nm3": float(coords.shape[0] * voxel_volume),
        "largest_component_bbox_cells_xyz": (hi - lo + 1).astype(int).tolist(),
        "largest_component_bbox_size_nm_xyz": ((hi - lo + 1) * np.asarray(spacing_nm)).astype(float).tolist(),
        "largest_component_internal_face_edges_xyz": edge_counts,
        "largest_component_mean_graph_degree": float(2 * sum(edge_counts) / coords.shape[0]),
        "largest_component_degree_histogram": degree_hist,
        "shortest_face_to_face_paths": {axis: paths[i] for i, axis in enumerate("xyz")},
        "skeleton_contact_face_counts_xyz": skeleton_contact_faces,
        "skeleton_contact_area_estimate_nm2_xyz": face_area_nm2,
        "skeleton_contact_area_estimate_nm2_total": float(sum(face_area_nm2)),
        "touches_local_window_z_top": local_top_face,
        "touches_original_sample_z_top": original_top_face,
    }


def main() -> None:
    with np.load(MASKS) as data:
        spacing_nm = (0.299635, 0.301857484616085, 0.317223)
        results = [analyze_threshold(data, q, spacing_nm) for q in (0.1, 0.5, 0.9, 0.99)]
        output = {
            "status": "TOPOLOGY_ONLY_LIQUID_AG_VOXEL_NETWORK_ANALYSIS",
            "method": "six-face voxel graph; Euclidean edge lengths for shortest paths; no pressure or flow solver",
            "source_geometry_sha256": GEOMETRY_SHA256,
            "source_checkpoint_sha256": CHECKPOINT_SHA256,
            "source_checkpoint_time_ps": 0.8495860511825916,
            "grid_spacing_nm_xyz": spacing_nm,
            "scope": "geometric connectivity, path tortuosity, and scaffold contacts only; no viscosity, conductance, velocity, flow rate, or crater prediction",
            "threshold_results": results,
        }
    (ROOT / "liquid_voxel_network_metrics.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    with (ROOT / "liquid_voxel_network_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ["threshold", "liquid_voxels", "components", "largest_voxels", "largest_volume_nm3", "mean_degree", "tortuosity_x", "tortuosity_y", "tortuosity_z", "scaffold_contact_area_nm2", "touches_original_z_top"]
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerow({
                "threshold": result["liquid_fraction_threshold"],
                "liquid_voxels": result["liquid_voxels_all_components"],
                "components": result["six_face_components"],
                "largest_voxels": result["largest_component_voxels"],
                "largest_volume_nm3": result["largest_component_volume_nm3"],
                "mean_degree": result["largest_component_mean_graph_degree"],
                "tortuosity_x": result["shortest_face_to_face_paths"]["x"]["tortuosity"],
                "tortuosity_y": result["shortest_face_to_face_paths"]["y"]["tortuosity"],
                "tortuosity_z": result["shortest_face_to_face_paths"]["z"]["tortuosity"],
                "scaffold_contact_area_nm2": result["skeleton_contact_area_estimate_nm2_total"],
                "touches_original_z_top": result["touches_original_sample_z_top"],
            })
    print(json.dumps({"output": "liquid_voxel_network_metrics.json", "threshold_results": results}, indent=2))


if __name__ == "__main__":
    main()
