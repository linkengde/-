#!/usr/bin/env python3
"""Assess v5 morphology mapping and bounded size/count perturbations on copies."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
TTM = HERE.parent
MODEL = TTM.parent
PACK = MODEL.parents[1]
REPO_ROOT = MODEL.parents[3]
sys.path.insert(0, str(TTM))
sys.path.insert(0, str(HERE))

from agw_ttm_reference import (  # noqa: E402
    AG_ATOMIC_VOLUME_A3,
    W_ATOMIC_VOLUME_A3,
    AtomGridMap,
    AgWMappingError,
    _raw_cic_volume_power,
    build_agw_grid_fields_joint,
    gaussian_ag_source_power,
)
from diagnose_agw_ttm_mapping_v7 import V6, local_capacity_flow, read_dump, read_lammps_data_atoms  # noqa: E402
from scan_agw_ttm_mapping_v8 import geometry_fraction_at_offset, regularize_v5_mask, sha256  # noqa: E402

OUT = HERE / "agw_ttm_v8_geometry_size_sensitivity.json"
DELETE_IDS_OUT = HERE / "v8_candidate_xy1pct_atoms1pct_remove_ids.csv"
SHAPE = (20, 20, 16)
BASE_SPACING = np.array([10.0, 10.0, 10.0])
TRANSFER_A = 20.0
OFFSET = (0.0, 2.5)


def source_for_box(xyz: np.ndarray, box_xy: float) -> np.ndarray:
    dxy = xyz[:, :2] - np.array([box_xy / 2.0, box_xy / 2.0])
    dxy = (dxy + box_xy / 2.0) % box_xy - box_xy / 2.0
    depth = np.maximum(0.0, 200.0 - xyz[:, 2])
    return (5.0 / 20.0) * np.exp(-np.sum(dxy**2, axis=1) / 30.0**2 - depth / 20.0)


def raw_fields(xyz_ag, xyz_w, q_ag, spacing):
    map_ag, vol_ag, src = _raw_cic_volume_power(xyz_ag, q_ag, SHAPE, tuple(spacing), AG_ATOMIC_VOLUME_A3)
    map_w, vol_w, _ = _raw_cic_volume_power(xyz_w, np.zeros(len(xyz_w)), SHAPE, tuple(spacing), W_ATOMIC_VOLUME_A3)
    return map_ag, vol_ag, src, map_w, vol_w


def geometry_guided_mapping(mask_w: np.ndarray, xyz_ag: np.ndarray, xyz_w: np.ndarray,
                            q_ag: np.ndarray, spacing: tuple[float, float, float],
                            offset_xy_A: tuple[float, float], sigma_cells: float = 0.0) -> dict:
    """Scale v5 morphology fractions to atom-volume totals and remap CIC support."""
    shift_x = int(round(offset_xy_A[0] / 2.5))
    shift_y = int(round(offset_xy_A[1] / 2.5))
    geom_w = geometry_fraction_at_offset(mask_w, shift_x, shift_y)
    if sigma_cells > 0.0:
        geom_w = gaussian_filter(geom_w, sigma=(sigma_cells, sigma_cells, sigma_cells),
                                 mode=("wrap", "wrap", "nearest"))
    geom_ag = 1.0 - geom_w
    cell_volume = float(np.prod(spacing))
    expected_ag = len(xyz_ag) * AG_ATOMIC_VOLUME_A3
    expected_w = len(xyz_w) * W_ATOMIC_VOLUME_A3
    scale_ag = expected_ag / float(geom_ag.sum() * cell_volume)
    scale_w = expected_w / float(geom_w.sum() * cell_volume)
    frac_ag = geom_ag * scale_ag
    frac_w = geom_w * scale_w
    if np.any(frac_ag + frac_w > 1.0 + 1e-12):
        raise RuntimeError("Scaled morphology fractions exceed cell capacity")

    map_ag, raw_ag, src_raw, map_w, raw_w = raw_fields(xyz_ag, xyz_w, q_ag, spacing)
    ncell = int(np.prod(SHAPE))
    cell_idx = np.indices(SHAPE).reshape(3, -1).T
    centres = (cell_idx + 0.5) * np.asarray(spacing)
    lengths = np.asarray(SHAPE) * np.asarray(spacing)

    def redirect(mapping: AtomGridMap, fraction: np.ndarray, phase: str):
        active = np.flatnonzero(fraction.ravel() > 0.0)
        if not len(active):
            raise RuntimeError(f"No active geometry cells for {phase}")
        active_xyz = centres[active]
        tiled_xyz, tiled_ids = [], []
        for sx in (-lengths[0], 0.0, lengths[0]):
            for sy in (-lengths[1], 0.0, lengths[1]):
                copy = active_xyz.copy()
                copy[:, 0] += sx
                copy[:, 1] += sy
                tiled_xyz.append(copy)
                tiled_ids.append(active)
        tree = cKDTree(np.concatenate(tiled_xyz, axis=0))
        original_ids = np.arange(ncell)
        invalid = np.flatnonzero(fraction.ravel() <= 0.0)
        if len(invalid):
            distances, nearest = tree.query(centres[invalid], k=1, workers=-1)
            nearest_active = np.concatenate(tiled_ids)[nearest]
            remap = original_ids.copy()
            remap[invalid] = nearest_active
            cell_ids = remap[mapping.cell_ids]
            moved = np.isin(mapping.cell_ids, invalid) & (mapping.weights > 0.0)
            moved_weights = mapping.weights[moved]
            moved_distances = distances[np.searchsorted(invalid, mapping.cell_ids[moved])]
            weighted_dist = float(np.dot(moved_weights, moved_distances) / moved_weights.sum()) if moved_weights.sum() else 0.0
            max_dist = float(moved_distances.max()) if len(moved_distances) else 0.0
            affected_atoms = int(np.count_nonzero(np.any(moved, axis=1)))
            moved_weight_fraction = float(moved_weights.sum() / mapping.weights.sum())
        else:
            cell_ids = mapping.cell_ids.copy()
            weighted_dist = max_dist = moved_weight_fraction = 0.0
            affected_atoms = 0
        moved_map = AtomGridMap(cell_ids, mapping.weights.copy(), mapping.grid_shape)
        support = np.zeros(ncell, dtype=bool)
        nonzero = moved_map.weights > 0.0
        support[moved_map.cell_ids[nonzero]] = True
        return moved_map, {
            "phase": phase,
            "geometry_active_cells": int(len(active)),
            "active_cells_without_atom_support": int(np.count_nonzero((fraction.ravel() > 0.0) & ~support)),
            "raw_CIC_shares_in_zero_geometry_cells": int(np.count_nonzero((mapping.weights > 0.0) & (fraction.ravel()[mapping.cell_ids] == 0.0))),
            "affected_atoms": affected_atoms,
            "moved_weight_fraction": moved_weight_fraction,
            "weighted_mean_relocation_A": weighted_dist,
            "maximum_relocation_A": max_dist,
            "atom_map_row_sum_max_error": float(np.max(np.abs(moved_map.weights.sum(axis=1) - 1.0))),
        }

    ag_map, ag_reloc = redirect(map_ag, frac_ag, "Ag")
    w_map, w_reloc = redirect(map_w, frac_w, "W")
    source = np.bincount(
        ag_map.cell_ids.ravel(),
        weights=np.repeat(AG_ATOMIC_VOLUME_A3 * q_ag, ag_map.cell_ids.shape[1]) * ag_map.weights.ravel(),
        minlength=ncell,
    ).reshape(SHAPE)
    bad_source = int(np.count_nonzero((source > 1e-14) & (frac_ag <= 0.0)))
    mapped_ag = np.bincount(ag_map.cell_ids.ravel(), weights=np.repeat(
        np.full(len(xyz_ag), AG_ATOMIC_VOLUME_A3), ag_map.cell_ids.shape[1]) * ag_map.weights.ravel(), minlength=ncell).reshape(SHAPE)
    mapped_w = np.bincount(w_map.cell_ids.ravel(), weights=np.repeat(
        np.full(len(xyz_w), W_ATOMIC_VOLUME_A3), w_map.cell_ids.shape[1]) * w_map.weights.ravel(), minlength=ncell).reshape(SHAPE)
    return {
        "status": "DIAGNOSTIC_ONLY",
        "method": "v5 regularized mask fractions scaled to phase atomic-volume totals; CIC atom map shares redirected to nearest active same-phase geometry cell",
        "morphology_smoothing_sigma_cells": float(sigma_cells),
        "fraction_scale_Ag": float(scale_ag),
        "fraction_scale_W": float(scale_w),
        "fraction_Ag_volume_A3": float(frac_ag.sum() * cell_volume),
        "fraction_W_volume_A3": float(frac_w.sum() * cell_volume),
        "fraction_sum_min": float((frac_ag + frac_w).min()),
        "fraction_sum_max": float((frac_ag + frac_w).max()),
        "phase_volume_relative_conservation_error": {
            "Ag": float((frac_ag.sum() * cell_volume - expected_ag) / expected_ag),
            "W": float((frac_w.sum() * cell_volume - expected_w) / expected_w),
        },
        "source_power_eV_ps": float(source.sum()),
        "source_expected_eV_ps": float(np.dot(np.full(len(xyz_ag), AG_ATOMIC_VOLUME_A3), q_ag)),
        "source_on_zero_Ag_fraction_cells": bad_source,
        "relocation": {"Ag": ag_reloc, "W": w_reloc},
        "atom_CIC_volume_vs_geometry_fraction_L1_A3": {
            "Ag": float(np.abs(mapped_ag - frac_ag * cell_volume).sum()),
            "W": float(np.abs(mapped_w - frac_w * cell_volume).sum()),
        },
        "atom_CIC_max_cell_fill_after_support_relocation": float(((mapped_ag + mapped_w) / cell_volume).max()),
        "mapping_is_not_a_validated_local_density_field": True,
    }


def deletion_score(mapping: AtomGridMap, raw_total: np.ndarray) -> np.ndarray:
    flat = raw_total.ravel()
    return np.sum(mapping.weights * flat[mapping.cell_ids], axis=1)


def main() -> None:
    manifest = json.loads((V6 / "exploratory_ttm_state_manifest.json").read_text(encoding="utf-8"))
    data_path = REPO_ROOT / manifest["data_file"]
    ids_t, types_t, xyz_t = read_lammps_data_atoms(data_path)
    ids_r, types_r, xyz_r = read_dump(V6 / "nve_fixedids_after_nvt_final.dump.gz")
    if not (np.array_equal(ids_t, ids_r) and np.array_equal(types_t, types_r)):
        raise RuntimeError("Template and frozen restart ID/type arrays do not match")
    mask_w = regularize_v5_mask(np.load(MODEL / "qsgs_morphology.npz")["W_backbone"].astype(bool))
    candidate = (xyz_t[:, 2] >= 40.0 - 1e-9) & (xyz_t[:, 2] <= 300.0 + 1e-9)
    mapped = candidate & (xyz_r[:, 2] >= 40.0 - 1e-9) & (xyz_r[:, 2] <= 200.0 + 1e-9)
    idx_ag = np.flatnonzero(mapped & (types_r == 2))
    idx_w = np.flatnonzero(mapped & (types_r == 1))
    offset = np.asarray(OFFSET)
    xy = xyz_r[:, :2].copy()
    xy[:, 0] = (xy[:, 0] - offset[0]) % 200.0
    xy[:, 1] = (xy[:, 1] - offset[1]) % 200.0
    xyz_ag = np.column_stack((xy[idx_ag], xyz_r[idx_ag, 2] - 40.0))
    xyz_w = np.column_stack((xy[idx_w], xyz_r[idx_w, 2] - 40.0))
    q_ag = gaussian_ag_source_power(xyz_r[idx_ag])
    geometry_candidates = [
        geometry_guided_mapping(mask_w, xyz_ag, xyz_w, q_ag, (10.0, 10.0, 10.0), OFFSET, sigma)
        for sigma in (0.0, 0.25, 0.5, 0.75, 1.0)
    ]

    full_ag = int(np.count_nonzero(types_t == 2))
    full_w = int(np.count_nonzero(types_t == 1))
    original_mass_ratio = full_ag * 107.8682 / (full_w * 183.84)
    size_scan = []
    for strain in (0.0, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05):
        factor = 1.0 + strain
        spacing = BASE_SPACING * np.array([factor, factor, 1.0])
        a = xyz_ag.copy(); w = xyz_w.copy()
        a[:, :2] *= factor; w[:, :2] *= factor
        scaled_abs_xy = xyz_r[idx_ag, :2] * factor
        q_scaled = source_for_box(np.column_stack((scaled_abs_xy, xyz_r[idx_ag, 2])), 200.0 * factor)
        map_a, vol_a, src, map_w, vol_w = raw_fields(a, w, q_scaled, spacing)
        fill = (vol_a + vol_w) / float(np.prod(spacing))
        flow = local_capacity_flow(vol_a + vol_w, tuple(spacing), TRANSFER_A)
        try:
            joint = build_agw_grid_fields_joint(a, w, q_scaled, SHAPE, tuple(spacing), TRANSFER_A)
            mapper = {"status": "PASS", "diagnostics": joint["diagnostics"]}
        except AgWMappingError as exc:
            mapper = {"status": "REJECTED", "reason": str(exc)}
        size_scan.append({
            "xy_linear_strain_fraction": strain,
            "box_xy_A": [200.0 * factor, 200.0 * factor],
            "grid_spacing_A": spacing.tolist(),
            "raw_combined_max_fill": float(fill.max()),
            "raw_overfilled_cells": int(np.count_nonzero(fill > 1.0 + 1e-10)),
            "raw_overflow_A3": float(np.maximum(vol_a + vol_w - np.prod(spacing), 0.0).sum()),
            "maxflow_20A": flow,
            "joint_mapper": mapper,
            "Ag_source_power_eV_ps": float(src.sum()),
            "source_power_relative_change": float(src.sum() / np.dot(np.full(len(q_ag), AG_ATOMIC_VOLUME_A3), q_ag) - 1.0),
            "mapping_only_sensitivity_copy": True,
            "affine_xy_dilation_requires_mechanical_relaxation_before_any_md": True,
        })

    # Atom-count candidates remove the highest local CIC-density-score atoms in
    # each phase, in a separate diagnostic copy. No atoms are deleted from the
    # frozen restart. Removal fractions are matched by species to preserve the
    # already-established approximately equal-mass inventory.
    raw_ag_map, vol_ag, src_base, raw_w_map, vol_w = raw_fields(xyz_ag, xyz_w, q_ag, BASE_SPACING)
    total_fill = (vol_ag + vol_w) / np.prod(BASE_SPACING)
    score_ag = deletion_score(raw_ag_map, total_fill)
    score_w = deletion_score(raw_w_map, total_fill)
    count_scan = []
    deletion_masks: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    for fraction in (0.001, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02):
        remove_ag = int(round(full_ag * fraction))
        remove_w = int(round(full_w * fraction))
        keep_ag = np.ones(len(xyz_ag), dtype=bool)
        keep_w = np.ones(len(xyz_w), dtype=bool)
        if remove_ag >= len(keep_ag) or remove_w >= len(keep_w):
            continue
        keep_ag[np.argsort(score_ag, kind="stable")[-remove_ag:]] = False
        keep_w[np.argsort(score_w, kind="stable")[-remove_w:]] = False
        deletion_masks[fraction] = (keep_ag, keep_w)
        a, w, qa = xyz_ag[keep_ag], xyz_w[keep_w], q_ag[keep_ag]
        ma, va, src, mw, vw = raw_fields(a, w, qa, BASE_SPACING)
        fill = (va + vw) / np.prod(BASE_SPACING)
        flow = local_capacity_flow(va + vw, tuple(BASE_SPACING), TRANSFER_A)
        try:
            joint = build_agw_grid_fields_joint(a, w, qa, SHAPE, tuple(BASE_SPACING), TRANSFER_A)
            mapper = {"status": "PASS", "diagnostics": joint["diagnostics"]}
        except AgWMappingError as exc:
            mapper = {"status": "REJECTED", "reason": str(exc)}
        new_ag = full_ag - remove_ag
        new_w = full_w - remove_w
        ratio = new_ag * 107.8682 / (new_w * 183.84)
        count_scan.append({
            "requested_fraction_per_species": fraction,
            "removed_atoms": {"Ag": remove_ag, "W": remove_w},
            "removed_fraction_of_total_atoms": float((remove_ag + remove_w) / (full_ag + full_w)),
            "remaining_total_atoms": int(new_ag + new_w),
            "mass_ratio_Ag_over_W": float(ratio),
            "mass_ratio_error_fraction": float(ratio - 1.0),
            "within_1_to_1_plusminus_0p1_percent": bool(abs(ratio - 1.0) <= 1e-3),
            "raw_combined_max_fill": float(fill.max()),
            "raw_overfilled_cells": int(np.count_nonzero(fill > 1.0 + 1e-10)),
            "raw_overflow_A3": float(np.maximum(va + vw - np.prod(BASE_SPACING), 0.0).sum()),
            "maxflow_20A": flow,
            "joint_mapper": mapper,
            "source_power_eV_ps": float(src.sum()),
            "source_power_relative_change": float(src.sum() / src_base.sum() - 1.0),
            "selection": "highest atom CIC contribution score in densest mapped cells, separate by phase",
            "mapping_only_sensitivity_copy": True,
            "physical_warning": "Vacancy creation is not accepted without atomistic relaxation, same-phase density checks, and connectivity checks.",
        })

    combined_scan = []
    for strain in (0.0, 0.005, 0.01, 0.015, 0.02, 0.03):
        factor = 1.0 + strain
        spacing = BASE_SPACING * np.array([factor, factor, 1.0])
        for fraction in (0.0, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02):
            if fraction == 0.0:
                ka, kw = np.ones(len(xyz_ag), dtype=bool), np.ones(len(xyz_w), dtype=bool)
            else:
                ka, kw = deletion_masks[fraction]
            a = xyz_ag[ka].copy(); w = xyz_w[kw].copy(); qa = q_ag[ka].copy()
            a[:, :2] *= factor; w[:, :2] *= factor
            scaled_abs_xy = xyz_r[idx_ag[ka], :2] * factor
            q_scaled = source_for_box(np.column_stack((scaled_abs_xy, xyz_r[idx_ag[ka], 2])), 200.0 * factor)
            _, va, src, _, vw = raw_fields(a, w, q_scaled, spacing)
            flow = local_capacity_flow(va + vw, tuple(spacing), TRANSFER_A)
            ratio = (full_ag - int(round(full_ag * fraction))) * 107.8682 / (
                (full_w - int(round(full_w * fraction))) * 183.84
            )
            try:
                joint = build_agw_grid_fields_joint(a, w, q_scaled, SHAPE, tuple(spacing), TRANSFER_A)
                mapper = {"status": "PASS", "diagnostics": joint["diagnostics"]}
            except AgWMappingError as exc:
                mapper = {"status": "REJECTED", "reason": str(exc)}
            combined_scan.append({
                "xy_linear_strain_fraction": strain,
                "box_xy_A": [200.0 * factor, 200.0 * factor],
                "grid_spacing_A": spacing.tolist(),
                "removed_fraction_per_species": fraction,
                "removed_atoms": {"Ag": int(round(full_ag * fraction)), "W": int(round(full_w * fraction))},
                "mass_ratio_Ag_over_W": float(ratio),
                "mass_ratio_within_1_to_1_plusminus_0p1_percent": bool(abs(ratio - 1.0) <= 1.0e-3),
                "maxflow_20A": flow,
                "joint_mapper": mapper,
                "source_power_eV_ps": float(src.sum()),
                "source_relative_change_from_frozen": float(src.sum() / src_base.sum() - 1.0),
                "mapping_only_sensitivity_copy": True,
                "post_relaxation_density_and_connectivity_not_yet_checked": True,
            })

    keep_ag_1pct, keep_w_1pct = deletion_masks[0.01]
    delete_ids = np.concatenate((ids_r[idx_ag[~keep_ag_1pct]], ids_r[idx_w[~keep_w_1pct]])).astype(np.int64)
    delete_types = np.concatenate((np.full(np.count_nonzero(~keep_ag_1pct), 2, dtype=np.int8),
                                   np.full(np.count_nonzero(~keep_w_1pct), 1, dtype=np.int8)))
    delete_rows = np.column_stack((delete_ids, delete_types))
    np.savetxt(DELETE_IDS_OUT, delete_rows, fmt="%d", delimiter=",", header="atom_id,type", comments="")

    result = {
        "title": "Ag-W TTM v8 v5-geometry and box/atom-count sensitivity",
        "frozen_input_sha256": {
            "template_data": sha256(data_path),
            "restart_dump": sha256(V6 / "nve_fixedids_after_nvt_final.dump.gz"),
        },
        "baseline": {
            "mapped_Ag_atoms": int(len(idx_ag)),
            "mapped_W_atoms": int(len(idx_w)),
            "all_atoms": int(len(types_t)),
            "mass_ratio_Ag_over_W": float(original_mass_ratio),
            "grid_spacing_A": BASE_SPACING.tolist(),
            "grid_shape": list(SHAPE),
            "phase_offset_xy_A": list(OFFSET),
        },
        "v5_geometry_guided_mapping_candidates": geometry_candidates,
        "xy_box_scale_scan": size_scan,
        "atom_count_scan": count_scan,
        "combined_box_and_count_scan": combined_scan,
        "candidate_id_list": {
            "file": DELETE_IDS_OUT.name,
            "sha256": sha256(DELETE_IDS_OUT),
            "candidate": "1.0% xy affine dilation plus removal of 1.0% of each full-system species from highest-CIC-density mapped atoms",
            "removed_Ag": int(np.count_nonzero(~keep_ag_1pct)),
            "removed_W": int(np.count_nonzero(~keep_w_1pct)),
        },
        "decision": "No perturbation is applied to the frozen structure. Candidate tuning is accepted only if local mapping and actual LAMMPS preflight pass, and post-relaxation density, mass ratio, and bicontinuity remain intact.",
        "limitations": [
            "V5 geometry fractions are a QSGS morphology reference, not a measured Voronoi atomic-density field.",
            "XY scaling is a uniform affine strain; every candidate needs mechanical relaxation before MD.",
            "Atom deletion candidates are diagnostic vacancy perturbations only and are not structurally validated.",
            "The failed EFS static audit and provisional Ag/W electronic interface conductance are unchanged.",
            "No TTM time step or laser source was applied.",
        ],
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "v5_geometry": [{"sigma_cells": g["morphology_smoothing_sigma_cells"], "fraction_sum_max": g["fraction_sum_max"], "relocation_Ag": g["relocation"]["Ag"], "relocation_W": g["relocation"]["W"], "atom_CIC_max_fill": g["atom_CIC_max_cell_fill_after_support_relocation"], "source_on_zero_Ag_cells": g["source_on_zero_Ag_fraction_cells"]} for g in geometry_candidates],
        "size_scan": [{"strain": x["xy_linear_strain_fraction"], "raw_max": x["raw_combined_max_fill"], "overflow": x["raw_overflow_A3"], "unresolved": x["maxflow_20A"].get("locally_unassignable_overflow_A3"), "mapper": x["joint_mapper"]["status"]} for x in size_scan],
        "count_scan": [{"fraction": x["requested_fraction_per_species"], "mass_ratio": x["mass_ratio_Ag_over_W"], "unresolved": x["maxflow_20A"].get("locally_unassignable_overflow_A3"), "mapper": x["joint_mapper"]["status"], "source_change": x["source_power_relative_change"]} for x in count_scan],
        "combined_pass_candidates": [{"strain": x["xy_linear_strain_fraction"], "removed_fraction": x["removed_fraction_per_species"], "mass_ratio": x["mass_ratio_Ag_over_W"], "unresolved": x["maxflow_20A"].get("locally_unassignable_overflow_A3"), "mapper": x["joint_mapper"]["status"]} for x in combined_scan if x["maxflow_20A"].get("status") == "PASS" and x["joint_mapper"]["status"] == "PASS"],
    }, indent=2))


if __name__ == "__main__":
    main()
