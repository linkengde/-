#!/usr/bin/env python3
"""Audit counts, phase density, and coarse connectivity of the v8 copy."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent
TTM = HERE.parent
MODEL = TTM.parent
V6 = MODEL / "interface_reconstruction_v6_exploratory"
CANDIDATE = HERE / "candidate_xy202A_atoms_minus1pct"
sys.path.insert(0, str(TTM))
from diagnose_agw_ttm_mapping_v7 import read_dump, read_lammps_data_atoms  # noqa: E402

OUT = CANDIDATE / "candidate_structure_screen.json"
MATERIAL_Z_A = 200.0
NX = NY = NZ = 40
W_VOXEL_FRACTION = 181504 / (80**3)
MASS_AG = 107.8682
MASS_W = 183.84


def phase_connectivity(xyz: np.ndarray, type_mask: np.ndarray, box_xy: float) -> dict:
    use = type_mask & (xyz[:, 2] >= 0.0) & (xyz[:, 2] <= MATERIAL_Z_A)
    points = xyz[use]
    grid = np.zeros((NX, NY, NZ), dtype=bool)
    idx = np.column_stack((
        np.floor(np.mod(points[:, 0], box_xy) / box_xy * NX).astype(int).clip(0, NX - 1),
        np.floor(np.mod(points[:, 1], box_xy) / box_xy * NY).astype(int).clip(0, NY - 1),
        np.floor(points[:, 2] / MATERIAL_Z_A * NZ).astype(int).clip(0, NZ - 1),
    ))
    grid[tuple(idx.T)] = True
    structure = ndimage.generate_binary_structure(3, 1)
    labels, nlabels = ndimage.label(grid, structure=structure)
    parent = np.arange(nlabels + 1, dtype=np.int64)
    def root(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = int(parent[a])
        return a
    def union(a: int, b: int) -> None:
        ra, rb = root(a), root(b)
        if ra != rb:
            parent[rb] = ra
    for left, right in ((labels[0, :, :], labels[-1, :, :]), (labels[:, 0, :], labels[:, -1, :])):
        pairs = np.column_stack((left.ravel(), right.ravel()))
        pairs = pairs[(pairs[:, 0] > 0) & (pairs[:, 1] > 0)]
        for a, b in np.unique(pairs, axis=0):
            union(int(a), int(b))
    sizes: dict[int, int] = {}
    for label, size in zip(*np.unique(labels[labels > 0], return_counts=True)):
        r = root(int(label))
        sizes[r] = sizes.get(r, 0) + int(size)
    if not sizes:
        raise RuntimeError("Empty phase occupancy")
    largest_root = max(sizes, key=sizes.get)
    largest = labels.copy()
    for label in range(1, nlabels + 1):
        largest[labels == label] = root(label)
    voxels = np.argwhere(largest == largest_root)
    return {
        "grid": [NX, NY, NZ],
        "connectivity": "6-neighbor voxel occupancy with periodic x/y component merging; open z",
        "occupied_voxels": int(grid.sum()),
        "component_count_after_xy_periodic_merge": len(sizes),
        "largest_component_voxels": int(sizes[largest_root]),
        "largest_component_share_of_occupied_voxels": float(sizes[largest_root] / grid.sum()),
        "largest_component_touches_x_both_sides": bool(np.any(voxels[:, 0] == 0) and np.any(voxels[:, 0] == NX - 1)),
        "largest_component_touches_y_both_sides": bool(np.any(voxels[:, 1] == 0) and np.any(voxels[:, 1] == NY - 1)),
        "largest_component_spans_z": bool(np.any(voxels[:, 2] == 0) and np.any(voxels[:, 2] == NZ - 1)),
        "atoms_in_material_z_range": int(len(points)),
        "screen_only_not_atomic_neighbor_topology": True,
    }


def audit(label: str, ids, types, xyz, box_xy: float) -> dict:
    n_ag = int(np.count_nonzero(types == 2))
    n_w = int(np.count_nonzero(types == 1))
    mass_ratio = n_ag * MASS_AG / (n_w * MASS_W)
    material_volume = box_xy * box_xy * MATERIAL_Z_A
    geom_w_volume = material_volume * W_VOXEL_FRACTION
    geom_ag_volume = material_volume * (1.0 - W_VOXEL_FRACTION)
    return {
        "label": label,
        "box_xy_A": box_xy,
        "atoms_total": int(len(ids)),
        "atoms_Ag": n_ag,
        "atoms_W": n_w,
        "mass_ratio_Ag_over_W": float(mass_ratio),
        "mass_ratio_error_fraction": float(mass_ratio - 1.0),
        "within_1_to_1_plusminus_0p1_percent": bool(abs(mass_ratio - 1.0) <= 1.0e-3),
        "phase_number_density_using_fixed_v5_geometry_Ag_atoms_A3": {
            "Ag": float(n_ag / geom_ag_volume),
            "W": float(n_w / geom_w_volume),
        },
        "phase_topology_screen_5A": {
            "Ag": phase_connectivity(xyz, types == 2, box_xy),
            "W": phase_connectivity(xyz, types == 1, box_xy),
        },
    }


def main() -> None:
    ids0, types0, _ = read_lammps_data_atoms(MODEL / "interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz")
    ids1, types1, xyz1 = read_dump(V6 / "nve_fixedids_after_nvt_final.dump.gz")
    ids2, types2, xyz2 = read_dump(CANDIDATE / "candidate.restart.dump.gz")
    ids3, types3, xyz3 = read_dump(CANDIDATE / "cold_after_nve.dump.gz")
    if not (np.array_equal(ids0, ids1) and np.array_equal(types0, types1)):
        raise RuntimeError("Baseline template and restart ID/type arrays differ")
    type_by_id = dict(zip(ids0.tolist(), types0.tolist()))
    for label, check_ids, check_types in (("initial", ids2, types2), ("cold_end", ids3, types3)):
        if len(check_ids) != len(np.unique(check_ids)):
            raise RuntimeError(f"Candidate {label} dump contains duplicate atom IDs")
        try:
            expected = np.fromiter((type_by_id[int(i)] for i in check_ids), dtype=int, count=len(check_ids))
        except KeyError as exc:
            raise RuntimeError(f"Candidate {label} dump contains unknown atom ID {exc.args[0]}") from exc
        if not np.array_equal(check_types, expected):
            raise RuntimeError(f"Candidate {label} dump changed an atom type")
    result = {
        "title": "Ag-W v8 tuned-candidate count, density, and coarse connectivity screen",
        "baseline": audit("frozen v6 restart", ids1, types1, xyz1, 200.0),
        "candidate_initial": audit("1% XY dilation plus 1% species-matched atom removal, before cold check", ids2, types2, xyz2, 202.0),
        "candidate_after_cold_check": audit("same candidate after 0.5 ps NVT and 0.05 ps NVE", ids3, types3, xyz3, 202.0),
        "interpretation": [
            "The 5 A voxel connectivity calculation is a coarse geometry screen, not a bond-network or percolation proof.",
            "The phase-density estimate assumes the fixed v5 QSGS geometry volumes and therefore reflects the atom-count/box perturbation.",
            "A candidate mapping PASS does not establish cold equilibrium or static EFS acceptance.",
        ],
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "baseline_mass_ratio": result["baseline"]["mass_ratio_Ag_over_W"],
        "candidate_mass_ratio": result["candidate_after_cold_check"]["mass_ratio_Ag_over_W"],
        "candidate_phase_density": result["candidate_after_cold_check"]["phase_number_density_using_fixed_v5_geometry_Ag_atoms_A3"],
        "candidate_topology": result["candidate_after_cold_check"]["phase_topology_screen_5A"],
    }, indent=2))


if __name__ == "__main__":
    main()
