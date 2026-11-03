#!/usr/bin/env python3
"""Read-only topology and throat-resolution audit of the existing Ag70 mask."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import ndimage


OUT = Path(__file__).resolve().parent
PROJECT = OUT.parents[2]
GEOMETRY = PROJECT / "runs/single_temperature_mesh2_validation_20261027/qsgs_morphology_ag70vol.npz"
REPORT = OUT / "ag70_geometry_audit.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def periodic_labels(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Six-face components, periodic x/y seams, open z faces."""
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

    roots = np.fromiter((root(i) for i in range(count + 1)), dtype=np.int64)
    labels = roots[labels]
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    sizes[0] = 0
    return labels, sizes


def spans_axis(mask: np.ndarray, axis: int) -> bool:
    labels, sizes = periodic_labels(mask)
    if not np.any(sizes):
        return False
    lo = np.take(labels, 0, axis=axis).ravel()
    hi = np.take(labels, -1, axis=axis).ravel()
    lo_ids = np.unique(lo[lo > 0])
    hi_ids = np.unique(hi[hi > 0])
    return bool(np.intersect1d(lo_ids, hi_ids, assume_unique=True).size)


def phase_connectivity(mask: np.ndarray) -> dict[str, object]:
    labels, sizes = periodic_labels(mask)
    positive = sizes[sizes > 0]
    largest_id = int(np.argmax(sizes)) if positive.size else 0
    seam_pairs = {
        "x": int(np.count_nonzero(np.take(mask, 0, axis=0) & np.take(mask, -1, axis=0))),
        "y": int(np.count_nonzero(np.take(mask, 0, axis=1) & np.take(mask, -1, axis=1))),
    }
    open_component_count = int(ndimage.label(mask)[1])
    return {
        "periodic_xy_face_connected_components": int(positive.size),
        "open_boundary_face_connected_components": open_component_count,
        "largest_component_cells": int(sizes[largest_id]) if positive.size else 0,
        "largest_component_fraction": float(sizes[largest_id] / np.count_nonzero(mask))
        if positive.size
        else 0.0,
        "percolates_x_periodically": spans_axis(mask, 0),
        "percolates_y_periodically": spans_axis(mask, 1),
        "spans_open_z_from_bottom_to_top": spans_axis(mask, 2),
        "direct_xy_seam_face_pairs": seam_pairs,
    }


def widest_percolating_throat(
    ag: np.ndarray,
    clearance_A: np.ndarray,
    spacing_A: np.ndarray,
    axis: int,
) -> dict[str, object]:
    """Threshold the Ag clearance field to estimate the widest path bottleneck.

    The critical EDT value is a voxel-center clearance to the nearest solid
    voxel center. The reported diameter subtracts one minimum cell width as a
    first-order voxel-boundary correction; it is an image-resolution estimate,
    not a molecular pore aperture measurement.
    """
    levels = np.unique(clearance_A[ag])
    low, high = 0, len(levels) - 1
    best: float | None = None
    while low <= high:
        mid = (low + high) // 2
        level = float(levels[mid])
        core = ag & (clearance_A >= level - 1.0e-12)
        if spans_axis(core, axis):
            best = level
            low = mid + 1
        else:
            high = mid - 1
    if best is None:
        return {
            "percolating_path_found": False,
            "critical_voxel_center_clearance_A": None,
            "estimated_bottleneck_diameter_A": None,
        }
    return {
        "percolating_path_found": True,
        "critical_voxel_center_clearance_A": best,
        "estimated_bottleneck_diameter_A": max(0.0, 2.0 * best - float(np.min(spacing_A))),
    }


def main() -> None:
    if not GEOMETRY.is_file():
        raise FileNotFoundError(GEOMETRY)
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        skeleton = np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
        seed = int(archive["seed"])

    ag = ~skeleton
    shape = np.asarray(ag.shape, dtype=int)
    spacing_A = size_A / shape
    cell_volume_A3 = float(np.prod(spacing_A))

    # Internal Ag/Ti3SiC2 interface: each x/y periodic face is counted once;
    # z is open and only internal material-material faces are included.
    interfaces = {
        "normal_x_faces": int(np.count_nonzero(ag != np.roll(ag, -1, axis=0))),
        "normal_y_faces": int(np.count_nonzero(ag != np.roll(ag, -1, axis=1))),
        "normal_z_faces": int(np.count_nonzero(ag[:, :, 1:] != ag[:, :, :-1])),
    }
    interface_area_A2 = (
        interfaces["normal_x_faces"] * spacing_A[1] * spacing_A[2]
        + interfaces["normal_y_faces"] * spacing_A[0] * spacing_A[2]
        + interfaces["normal_z_faces"] * spacing_A[0] * spacing_A[1]
    )

    # Tile in x/y before EDT so clearance near lateral borders respects the
    # periodic domain. z remains open; the solid skeleton defines the obstacles.
    tiled_ag = np.tile(ag, (3, 3, 1))
    tiled_clearance = ndimage.distance_transform_edt(tiled_ag, sampling=spacing_A)
    nx, ny, _ = ag.shape
    clearance_A = tiled_clearance[nx : 2 * nx, ny : 2 * ny, :]
    del tiled_ag, tiled_clearance
    minimum_clearance_A = float(np.min(clearance_A[ag]))
    minimum_clearance_cells = int(
        np.count_nonzero(np.isclose(clearance_A[ag], minimum_clearance_A, rtol=0.0, atol=1.0e-9))
    )

    throat = {
        "axis_order": "x, y, z",
        "minimum_interface_cell_clearance": {
            "minimum_center_clearance_to_scaffold_A": minimum_clearance_A,
            "one_cell_scale_width_estimate_A": max(
                0.0, 2.0 * minimum_clearance_A - float(np.min(spacing_A))
            ),
            "Ag_cells_at_minimum_clearance": minimum_clearance_cells,
            "note": "interface-adjacent cell clearance, not a measured channel throat; about one voxel and not a resolved continuum aperture",
        },
        "x_periodic": widest_percolating_throat(ag, clearance_A, spacing_A, 0),
        "y_periodic": widest_percolating_throat(ag, clearance_A, spacing_A, 1),
        "z_open_bottom_to_top": widest_percolating_throat(ag, clearance_A, spacing_A, 2),
        "method": (
            "largest EDT-clearance threshold retaining a six-face connected Ag path; "
            "approximate diameter = 2*critical center-clearance - minimum voxel spacing"
        ),
        "resolution_limit_A": float(np.min(spacing_A)),
    }

    report = {
        "status": "READ_ONLY_AG70_GEOMETRY_AUDIT_COMPLETE",
        "geometry_path": str(GEOMETRY.relative_to(PROJECT)),
        "geometry_sha256": sha256(GEOMETRY),
        "qsgs_seed": seed,
        "grid_xyz": shape.tolist(),
        "domain_size_nm_xyz": (size_A / 10.0).tolist(),
        "voxel_size_A_xyz": spacing_A.tolist(),
        "phase_volume_fraction": {
            "Ag": float(ag.mean()),
            "Ti3SiC2": float(skeleton.mean()),
        },
        "phase_cells": {
            "Ag": int(np.count_nonzero(ag)),
            "Ti3SiC2": int(np.count_nonzero(skeleton)),
        },
        "Ag_connectivity": phase_connectivity(ag),
        "Ag_percolating_throat_estimates": throat,
        "internal_Ag_Ti3SiC2_interface": {
            "phase_boundary_face_counts": interfaces,
            "area_nm2_voxel_face_estimate": float(interface_area_A2 / 100.0),
            "area_density_nm2_per_nm3": float(
                (interface_area_A2 / 100.0)
                / (np.prod(size_A / 10.0))
            ),
            "method": "anisotropic voxel-face area; lateral x/y periodic seams included, outer z surfaces excluded",
        },
        "continuum_resolution_assessment": {
            "smallest_voxel_nm": float(np.min(spacing_A) / 10.0),
            "grid_cells_across_1nm": float(1.0 / (np.min(spacing_A) / 10.0)),
            "interpretation": (
                "The approximately 0.30 nm cell spacing is atomic-scale. The mask does not provide a scale-separated, "
                "resolved pore-fluid mesh for quantitative continuum Navier-Stokes flow."
            ),
        },
        "limitations": [
            "The throat estimate is limited by voxelization and is not a validated hydrodynamic pore radius.",
            "Face-area interface is a stair-step voxel estimate; atomistic roughness and reconstruction are absent.",
            "Connectivity is topological only; it does not provide permeability, pressure drop, slip length, or flow velocity.",
            "The source geometry is a binary QSGS thermal mask, not a relaxed atomistic Ag/Ti3SiC2 interface.",
        ],
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
