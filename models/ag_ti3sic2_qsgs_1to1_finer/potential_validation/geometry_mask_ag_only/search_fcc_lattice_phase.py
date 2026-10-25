#!/usr/bin/env python3
"""Bounded FCC occupancy registration against the QSGS geometric mask.

Only box-compatible in-plane cubic orientations are considered. Every
candidate is a rigid transformation of the existing Ag network, so it keeps
the Ag count, internal nearest-neighbor distances, and graph connectivity.
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt, map_coordinates

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent.parent
MORPHOLOGY = MODEL / "potential_validation/relaxed_backbone_candidate/qsgs_morphology.npz"


def load_points():
    with gzip.open(HERE / "Ag_only_id_map.csv.gz", "rt", newline="") as f:
        rows = csv.reader(f)
        next(rows)
        indices = np.asarray([[int(row[6]), int(row[7]), int(row[8])] for row in rows], dtype=np.int64)
    # Always search relative to the unshifted EAM-FCC remapping, so rerunning
    # preparation with the same result JSON is idempotent.
    return indices.astype(float) * (4.09 / 2.0)


def make_sdf(skeleton: np.ndarray, size: np.ndarray):
    grid = np.asarray(skeleton.shape, dtype=int)
    spacing = size / grid
    tiled = np.tile(skeleton, (3, 3, 1))
    pad_z = 12
    padded = np.pad(tiled, ((0, 0), (0, 0), (pad_z, pad_z)), constant_values=False)
    outside = distance_transform_edt(~padded, sampling=spacing)
    inside = distance_transform_edt(padded, sampling=spacing)
    half_voxel = 0.5 * float(np.min(spacing))
    signed = np.where(padded, -(inside - half_voxel), outside - half_voxel)
    signed = signed[grid[0]:2 * grid[0], grid[1]:2 * grid[1], pad_z:pad_z + grid[2]].astype(np.float32)
    return signed, spacing


def transform_points(points: np.ndarray, orientation: str, original_size: np.ndarray):
    xyz = points.copy()
    lx, ly = original_size[:2]
    if orientation == "rotate_z_90":
        xyz[:, 0], xyz[:, 1] = points[:, 1], np.mod(-points[:, 0], lx)
        size = np.array([ly, lx, original_size[2]])
    elif orientation == "rotate_z_180":
        xyz[:, 0], xyz[:, 1] = np.mod(-points[:, 0], lx), np.mod(-points[:, 1], ly)
        size = original_size.copy()
    elif orientation == "rotate_z_270":
        xyz[:, 0], xyz[:, 1] = np.mod(-points[:, 1], ly), points[:, 0]
        size = np.array([ly, lx, original_size[2]])
    elif orientation == "identity":
        size = original_size.copy()
    else:
        raise ValueError(orientation)
    return xyz, size


def main():
    points = load_points()
    archive = np.load(HERE / "qsgs_signed_distance_A.npz")
    original_size = archive["size_A"].astype(float)
    original_phi = archive["phi_A"].astype(np.float32)
    original_spacing = archive["spacing_A"].astype(float)
    morphology = np.load(MORPHOLOGY)
    skeleton = morphology["skeleton"].astype(bool)
    a0 = 4.09
    max_shift_component = a0 / 2.0
    wall_k = 5.0
    orientations = ("identity", "rotate_z_90", "rotate_z_180", "rotate_z_270")

    geometries = {"identity": (original_phi, original_spacing, original_size)}
    for orientation in orientations[1:]:
        _, box = transform_points(points[:1], orientation, original_size)
        if np.array_equal(box, original_size):
            geometries[orientation] = (original_phi, original_spacing, original_size)
        else:
            sdf, spacing = make_sdf(skeleton, box)
            geometries[orientation] = (sdf, spacing, box)

    def measure(orientation, dx, dy):
        xyz, size = transform_points(points, orientation, original_size)
        xyz[:, 0] = np.mod(xyz[:, 0] + dx, size[0])
        xyz[:, 1] = np.mod(xyz[:, 1] + dy, size[1])
        phi, spacing, _ = geometries[orientation]
        q = np.empty((3, len(xyz)), dtype=float)
        q[0] = xyz[:, 0] / spacing[0] - 0.5
        q[1] = xyz[:, 1] / spacing[1] - 0.5
        q[2] = np.clip(xyz[:, 2] / spacing[2] - 0.5, 0.0, phi.shape[2] - 1.0)
        signed = map_coordinates(phi, q, order=1, mode="grid-wrap", prefilter=False)
        signed[(xyz[:, 2] < 0.0) | (xyz[:, 2] >= size[2])] = 1.0e6
        depth = np.maximum(-signed, 0.0)
        mobile = xyz[:, 2] > 20.0
        return {
            "orientation": orientation,
            "box_size_A": [float(v) for v in size],
            "dx_A": float(dx),
            "dy_A": float(dy),
            "translation_norm_A": float(np.hypot(dx, dy)),
            "penetrating_atoms_all": int(np.count_nonzero(signed < 0.0)),
            "penetrating_atoms_mobile": int(np.count_nonzero((signed < 0.0) & mobile)),
            "all_penetration_integral_A2": float(np.dot(depth, depth)),
            "mobile_wall_energy_eV_at_k5": float(0.5 * wall_k * np.dot(depth[mobile], depth[mobile])),
            "max_penetration_A": float(depth.max(initial=0.0)),
        }

    def best(candidates):
        return min(candidates, key=lambda r: (
            r["penetrating_atoms_all"],
            r["all_penetration_integral_A2"],
            r["translation_norm_A"],
            0 if r["orientation"] == "identity" else 1,
        ))

    offsets = np.linspace(-max_shift_component, max_shift_component, 17)
    all_candidates, orientation_best = [], {}
    for orientation in orientations:
        coarse = [measure(orientation, dx, dy) for dx in offsets for dy in offsets]
        coarse_best = best(coarse)
        refine_halfwidth = max_shift_component / 16.0
        rx = np.linspace(max(-max_shift_component, coarse_best["dx_A"] - refine_halfwidth), min(max_shift_component, coarse_best["dx_A"] + refine_halfwidth), 9)
        ry = np.linspace(max(-max_shift_component, coarse_best["dy_A"] - refine_halfwidth), min(max_shift_component, coarse_best["dy_A"] + refine_halfwidth), 9)
        refined = [measure(orientation, dx, dy) for dx in rx for dy in ry]
        baseline = measure(orientation, 0.0, 0.0)
        all_candidates.extend([*coarse, *refined, baseline])
        orientation_best[orientation] = best([*coarse, *refined, baseline])

    baseline = measure("identity", 0.0, 0.0)
    selected = best(all_candidates)
    result = {
        "status": "BOUNDED_FCC_PHASE_AND_ORIENTATION_SEARCH",
        "model": "Ag-only EAM network; existing QSGS skeleton signed-distance mask",
        "equilibrium_fcc_a0_A": a0,
        "orientations_tested": list(orientations),
        "translation_search": {
            "axes": ["x", "y"],
            "max_abs_component_A": max_shift_component,
            "max_translation_norm_A": float(np.sqrt(2.0) * max_shift_component),
            "coarse_grid_count_per_orientation": len(offsets) ** 2,
            "refinement_grid_count_per_orientation": 81,
            "score_order": ["minimize all-Ag mask penetration count", "minimize integrated squared penetration", "minimize translation norm", "prefer identity on exact ties"],
            "rigid_transforms_preserve_ag_count_and_internal_neighbor_distances": True,
            "periodic_box_compatible_connectivity_preserved": True,
        },
        "atom_count": int(len(points)),
        "baseline": baseline,
        "best_per_orientation": orientation_best,
        "selected": selected,
        "selected_change_vs_baseline": {
            "penetration_count_delta": int(selected["penetrating_atoms_all"] - baseline["penetrating_atoms_all"]),
            "mobile_penetration_count_delta": int(selected["penetrating_atoms_mobile"] - baseline["penetrating_atoms_mobile"]),
            "penetration_integral_delta_A2": float(selected["all_penetration_integral_A2"] - baseline["all_penetration_integral_A2"]),
        },
        "candidate_ranking": sorted(all_candidates, key=lambda r: (
            r["penetrating_atoms_all"], r["all_penetration_integral_A2"], r["translation_norm_A"],
            0 if r["orientation"] == "identity" else 1,
        ))[:24],
        "limits": [
            "Candidates are rigid transforms of the existing Ag network; no atoms are added or deleted.",
            "The 90/270-degree candidates swap in-plane cell dimensions to remain periodic-compatible; the QSGS voxel field is resampled to the swapped cell size.",
            "Mask overlap is geometric and does not represent adhesion, wetting, interfacial heat transfer, or ceramic response.",
        ],
    }
    path = HERE / "fcc_lattice_phase_search.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
