#!/usr/bin/env python3
"""Scan 10 A CIC mesh phases and compare atom fills with the v5 QSGS mask."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import gaussian_filter

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
    AgWMappingError,
    _raw_cic_volume_power,
    build_agw_grid_fields_joint,
    gaussian_ag_source_power,
)
from diagnose_agw_ttm_mapping_v7 import (  # noqa: E402
    V6,
    local_capacity_flow,
    read_dump,
    read_lammps_data_atoms,
)

OUTPUT = HERE / "agw_ttm_mapping_v8_report.json"
MASK_OUTPUT = HERE / "v5_regularized_qsgs_w_mask_2p5A.npz"
GRID_A = 10.0
TRANSFER_A = 20.0
OFFSET_VALUES_A = (0.0, 2.5, 5.0, 7.5)
GRID_SHAPE = (20, 20, 16)
SPACING = (10.0, 10.0, 10.0)
ORIGIN_BASE = np.array([0.0, 0.0, 40.0])
PHASE_BINS_A = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 60.0, 100.0, 150.0])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def regularize_v5_mask(source_mask: np.ndarray) -> np.ndarray:
    """Reproduce v5 sigma=2 thresholding with fixed source W voxel count."""
    source_mask = np.asarray(source_mask, dtype=bool)
    if source_mask.shape != (80, 80, 80):
        raise ValueError(f"Expected 80^3 QSGS mask, got {source_mask.shape}")
    occupancy = gaussian_filter(
        source_mask.astype(np.float32), sigma=2.0,
        mode=("wrap", "wrap", "nearest"),
    )
    n_w_voxels = int(source_mask.sum())
    take = np.argsort(occupancy.ravel(), kind="stable")[-n_w_voxels:]
    result = np.zeros(occupancy.size, dtype=bool)
    result[take] = True
    return result.reshape(source_mask.shape)


def geometry_fraction_at_offset(mask_w: np.ndarray, shift_x_vox: int, shift_y_vox: int) -> np.ndarray:
    """Aggregate 4x4x4 2.5 A morphology voxels into the 10 A mesh."""
    ix = (shift_x_vox + np.arange(80)) % 80
    iy = (shift_y_vox + np.arange(80)) % 80
    iz = np.arange(16, 80)  # z=40..200 A, no z phase shift
    fine = mask_w[np.ix_(ix, iy, iz)]
    return fine.reshape(20, 4, 20, 4, 16, 4).mean(axis=(1, 3, 5), dtype=np.float64)


def source_spatial_summary(power: np.ndarray, dx_A: float, dy_A: float) -> dict:
    nx, ny, nz = power.shape
    x = (np.arange(nx) + 0.5) * GRID_A + dx_A
    y = (np.arange(ny) + 0.5) * GRID_A + dy_A
    x = np.mod(x, 200.0)
    y = np.mod(y, 200.0)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    dxx = (xx - 100.0 + 100.0) % 200.0 - 100.0
    dyy = (yy - 100.0 + 100.0) % 200.0 - 100.0
    r = np.broadcast_to(np.sqrt(dxx**2 + dyy**2)[:, :, None], power.shape)
    depth = np.broadcast_to(
        (200.0 - (40.0 + (np.arange(nz) + 0.5) * GRID_A))[None, None, :], power.shape
    )
    total = float(power.sum())
    if total <= 0.0:
        raise RuntimeError("Ag source field has zero total power")
    radial_bins = np.digitize(r, PHASE_BINS_A[1:-1], right=False)
    radial_power = np.bincount(
        radial_bins.ravel(), weights=power.ravel(), minlength=len(PHASE_BINS_A) - 1
    )
    mean_r2 = float(np.sum(power * r**2) / total)
    mean_depth = float(np.sum(power * depth) / total)
    depth_variance = float(np.sum(power * (depth - mean_depth) ** 2) / total)
    return {
        "total_power_eV_ps": total,
        "radial_bin_edges_A": PHASE_BINS_A.tolist(),
        "radial_power_fraction": (radial_power / total).tolist(),
        "radial_rms_A": float(np.sqrt(mean_r2)),
        "depth_mean_A": mean_depth,
        "depth_rms_about_mean_A": float(np.sqrt(depth_variance)),
    }


def fraction_comparison(geometry_w: np.ndarray, raw_ag: np.ndarray, raw_w: np.ndarray) -> dict:
    geom_ag = 1.0 - geometry_w
    cell_volume = float(np.prod(SPACING))
    geometry_ag_volume = float(geom_ag.sum() * cell_volume)
    geometry_w_volume = float(geometry_w.sum() * cell_volume)
    atomic_ag = raw_ag * float(np.prod(SPACING))
    atomic_w = raw_w * float(np.prod(SPACING))

    def errors(reference: np.ndarray, observed: np.ndarray) -> dict:
        delta = np.abs(reference - observed)
        denom = float(np.sum(np.abs(reference)))
        corr = float(np.corrcoef(reference.ravel(), observed.ravel())[0, 1])
        return {
            "mean_absolute_fraction_difference": float(delta.mean()),
            "p95_absolute_fraction_difference": float(np.percentile(delta, 95)),
            "max_absolute_fraction_difference": float(delta.max()),
            "relative_L1_difference": float(np.abs(reference - observed).sum() / denom) if denom else None,
            "pearson_correlation": corr,
        }

    return {
        "geometry_phase_fraction_sum_min": float((geom_ag + geometry_w).min()),
        "geometry_phase_fraction_sum_max": float((geom_ag + geometry_w).max()),
        "geometry_Ag_volume_A3": geometry_ag_volume,
        "geometry_W_volume_A3": geometry_w_volume,
        "atom_CIC_Ag_volume_A3": float(atomic_ag.sum()),
        "atom_CIC_W_volume_A3": float(atomic_w.sum()),
        "atom_CIC_combined_fill_max": float((raw_ag + raw_w).max()),
        "geometry_vs_atom_CIC": {
            "Ag": errors(geom_ag, raw_ag),
            "W": errors(geometry_w, raw_w),
        },
        "geometry_W_vs_atom_derived_local_W_share": errors(
            geometry_w,
            np.divide(raw_w, raw_ag + raw_w, out=np.zeros_like(raw_w), where=(raw_ag + raw_w) > 0.0),
        ),
    }


def main() -> None:
    manifest = json.loads((V6 / "exploratory_ttm_state_manifest.json").read_text(encoding="utf-8"))
    data_path = REPO_ROOT / manifest["data_file"]
    restart_path = V6 / "nve_fixedids_after_nvt_final.dump.gz"
    ids_t, types_t, xyz_t = read_lammps_data_atoms(data_path)
    ids_r, types_r, xyz_r = read_dump(restart_path)
    if not (np.array_equal(ids_t, ids_r) and np.array_equal(types_t, types_r)):
        raise RuntimeError("Template and frozen restart ID/type arrays do not match")
    static_newton = (xyz_t[:, 2] >= 40.0 - 1e-9) & (xyz_t[:, 2] <= 300.0 + 1e-9)
    mapped = static_newton & (xyz_r[:, 2] >= 40.0 - 1e-9) & (xyz_r[:, 2] <= 200.0 + 1e-9)
    idx_ag = np.flatnonzero(mapped & (types_r == 2))
    idx_w = np.flatnonzero(mapped & (types_r == 1))
    xyz_ag_abs, xyz_w_abs = xyz_r[idx_ag], xyz_r[idx_w]
    q_ag = gaussian_ag_source_power(xyz_ag_abs)
    expected_source = float(np.dot(np.full(len(idx_ag), AG_ATOMIC_VOLUME_A3), q_ag))

    mask_archive = np.load(MODEL / "qsgs_morphology.npz")
    source_mask = mask_archive["W_backbone"].astype(bool)
    v5_mask = regularize_v5_mask(source_mask)
    source_mask_path = MODEL / "qsgs_morphology.npz"
    v5_report_path = MODEL / "interface_reconstruction_v5/trials/smoothed2_gap1p35_decay10_normal2/geometry_precondition_report.json"
    v5_report = json.loads(v5_report_path.read_text(encoding="utf-8"))
    changed_voxels = int(np.count_nonzero(v5_mask != source_mask))
    v5_expected = int(v5_report["source"]["morphology_voxels_reassigned"])
    if changed_voxels != v5_expected or int(v5_mask.sum()) != int(source_mask.sum()):
        raise RuntimeError("Recreated v5 regularized mask does not match the recorded voxel count/change count")
    np.savez_compressed(
        MASK_OUTPUT,
        W_backbone_v5_regularized=v5_mask.astype(np.uint8),
        source_qsgs_sha256=np.asarray(sha256(source_mask_path)),
        v5_mask_parameters=np.asarray(json.dumps(v5_report["parameters"], sort_keys=True)),
        morphology_voxels_reassigned=np.asarray(changed_voxels, dtype=np.int64),
    )

    baseline_profile = None
    phase_reports: list[dict] = []
    geometry_profiles: list[np.ndarray] = []
    any_mapper_pass = False
    for dx in OFFSET_VALUES_A:
        for dy in OFFSET_VALUES_A:
            rel_ag = xyz_ag_abs - ORIGIN_BASE
            rel_w = xyz_w_abs - ORIGIN_BASE
            rel_ag = rel_ag.copy(); rel_w = rel_w.copy()
            rel_ag[:, 0] -= dx; rel_ag[:, 1] -= dy
            rel_w[:, 0] -= dx; rel_w[:, 1] -= dy
            map_ag, volume_ag, source_raw = _raw_cic_volume_power(rel_ag, q_ag, GRID_SHAPE, SPACING, AG_ATOMIC_VOLUME_A3)
            map_w, volume_w, _ = _raw_cic_volume_power(rel_w, np.zeros(len(rel_w)), GRID_SHAPE, SPACING, W_ATOMIC_VOLUME_A3)
            cell_volume = float(np.prod(SPACING))
            frac_ag, frac_w = volume_ag / cell_volume, volume_w / cell_volume
            fill = frac_ag + frac_w
            flow = local_capacity_flow(volume_ag + volume_w, SPACING, TRANSFER_A)
            raw_profile = source_spatial_summary(source_raw, dx, dy)
            if baseline_profile is None:
                baseline_profile = np.asarray(raw_profile["radial_power_fraction"], dtype=float)
            raw_profile["radial_profile_relative_L1_vs_zero_offset"] = float(
                np.abs(np.asarray(raw_profile["radial_power_fraction"]) - baseline_profile).sum()
            )
            try:
                joint = build_agw_grid_fields_joint(
                    rel_ag, rel_w, q_ag, GRID_SHAPE, SPACING, TRANSFER_A
                )
                mapper_status = "PASS"
                any_mapper_pass = True
                joint_profile = source_spatial_summary(joint["source_power_ag_eV_ps"], dx, dy)
                joint_diagnostics = joint["diagnostics"]
            except AgWMappingError as exc:
                mapper_status = "REJECTED"
                joint_profile = None
                joint_diagnostics = {"reason": str(exc)}

            shift_x_vox = int(round(dx / 2.5))
            shift_y_vox = int(round(dy / 2.5))
            geom_w = geometry_fraction_at_offset(v5_mask, shift_x_vox, shift_y_vox)
            geometry_profiles.append(geom_w.astype(np.float32))
            geometry_comparison = fraction_comparison(geom_w, frac_ag, frac_w)
            phase_reports.append({
                "offset_xy_A": [dx, dy],
                "raw_atom_CIC": {
                    "Ag_max_fill": float(frac_ag.max()),
                    "W_max_fill": float(frac_w.max()),
                    "combined_max_fill": float(fill.max()),
                    "Ag_overfilled_cells": int(np.count_nonzero(frac_ag > 1.0 + 1e-10)),
                    "W_overfilled_cells": int(np.count_nonzero(frac_w > 1.0 + 1e-10)),
                    "combined_overfilled_cells": int(np.count_nonzero(fill > 1.0 + 1e-10)),
                    "overflow_volume_A3": float(np.maximum(volume_ag + volume_w - cell_volume, 0.0).sum()),
                    "Ag_volume_A3": float(volume_ag.sum()),
                    "W_volume_A3": float(volume_w.sum()),
                    "Ag_source_power_eV_ps": float(source_raw.sum()),
                    "expected_Ag_source_power_eV_ps": expected_source,
                },
                "local_capacity_flow_20A": flow,
                "joint_mapper": {"status": mapper_status, "diagnostics": joint_diagnostics},
                "Ag_source_spatial_distribution_raw": raw_profile,
                "Ag_source_spatial_distribution_joint": joint_profile,
                "v5_geometry_volume_fraction_comparison": geometry_comparison,
            })

    def phase_rank(row: dict) -> tuple[float, float, float, float]:
        flow = row["local_capacity_flow_20A"]
        return (
            float(flow.get("locally_unassignable_overflow_A3", float("inf"))),
            float(row["raw_atom_CIC"]["combined_overfilled_cells"]),
            float(row["raw_atom_CIC"]["combined_max_fill"]),
            float(row["Ag_source_spatial_distribution_raw"]["radial_profile_relative_L1_vs_zero_offset"]),
        )

    best = min(phase_reports, key=phase_rank)
    geometry_stack = np.stack(geometry_profiles, axis=0)
    np.savez_compressed(
        HERE / "v5_geometry_fractions_10A_all_offsets.npz",
        offsets_A=np.asarray([[x, y] for x in OFFSET_VALUES_A for y in OFFSET_VALUES_A], dtype=float),
        fraction_w=geometry_stack,
        fraction_ag=1.0 - geometry_stack,
    )
    result = {
        "title": "Ag-W TTM v8 mesh-phase and v5 geometry-fraction audit",
        "status": "V8_PHASE_WITH_REAL_ATOM_MAP_PASS" if any_mapper_pass else "V8_NO_REAL_ATOM_MAP_PASSES_20A_LIMIT",
        "frozen_inputs": {
            "template_data": manifest["data_file"],
            "template_data_sha256": manifest["data_sha256"],
            "restart_dump": "models/ag_w_qsgs_200A/interface_reconstruction_v6_exploratory/nve_fixedids_after_nvt_final.dump.gz",
            "restart_sha256": sha256(restart_path),
            "atom_count_unchanged": int(len(ids_r)),
            "mapped_atoms": int(mapped.sum()),
            "mapped_Ag_atoms": int(len(idx_ag)),
            "mapped_W_atoms": int(len(idx_w)),
            "grid_shape": list(GRID_SHAPE),
            "grid_spacing_A": list(SPACING),
            "transfer_limit_A": TRANSFER_A,
        },
        "phase_scan": {
            "offset_values_A": list(OFFSET_VALUES_A),
            "offset_count": len(phase_reports),
            "best_diagnostic_offset_xy_A": best["offset_xy_A"],
            "best_diagnostic_rank_key": list(phase_rank(best)),
            "all_offsets": phase_reports,
        },
        "v5_regularized_mask": {
            "source_mask_file": "models/ag_w_qsgs_200A/qsgs_morphology.npz",
            "source_mask_sha256": sha256(source_mask_path),
            "v5_geometry_report": "models/ag_w_qsgs_200A/interface_reconstruction_v5/trials/smoothed2_gap1p35_decay10_normal2/geometry_precondition_report.json",
            "v5_geometry_report_sha256": sha256(v5_report_path),
            "parameters": v5_report["parameters"],
            "source_W_voxels": int(source_mask.sum()),
            "regularized_W_voxels": int(v5_mask.sum()),
            "reassigned_voxels": changed_voxels,
            "changed_fraction": float(changed_voxels / v5_mask.size),
            "mask_archive": MASK_OUTPUT.name,
            "mask_array_sha256": hashlib.sha256(v5_mask.astype(np.uint8).tobytes(order="C")).hexdigest(),
            "coarse_fraction_archive": "v5_geometry_fractions_10A_all_offsets.npz",
            "geometry_fractions_are_diagnostic_only": True,
        },
        "box_atom_tuning": {
            "authorization": "User permits fine adjustment of box size and atom count.",
            "applied_in_this_audit": False,
            "reason": "First evaluate mesh phase and geometry-mask agreement on the frozen coordinates; any changed-box or changed-count candidate will be a separate copy with mass ratio and connectivity rechecked.",
        },
        "decision": "Geometry fractions are diagnostic and do not replace atom-derived fractions. Use only a real LAMMPS mapping-only preflight after selecting a phase that satisfies capacity/conservation; otherwise proceed to separate box/count sensitivity copies.",
        "limitations": [
            "The 20 A transfer radius is an exploratory locality screen, not a calibrated physical length.",
            "The v5 mask is a reconstructed geometry with exactly preserved W voxel count; its phase fractions are not verified as the frozen restart's local atomic-density field.",
            "No TTM step or laser source is executed by this scan.",
            "The v5 EFS static audit remains failed and Ag/W interface conductance remains exploratory.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "best_diagnostic_offset_xy_A": result["phase_scan"]["best_diagnostic_offset_xy_A"],
        "offset_summaries": [{
            "offset_xy_A": row["offset_xy_A"],
            "raw_max_fill": row["raw_atom_CIC"]["combined_max_fill"],
            "overfilled_cells": row["raw_atom_CIC"]["combined_overfilled_cells"],
            "overflow_A3": row["raw_atom_CIC"]["overflow_volume_A3"],
            "locally_unassignable_A3": row["local_capacity_flow_20A"].get("locally_unassignable_overflow_A3"),
            "mapper": row["joint_mapper"]["status"],
            "source_profile_L1": row["Ag_source_spatial_distribution_raw"]["radial_profile_relative_L1_vs_zero_offset"],
            "geom_W_MAE": row["v5_geometry_volume_fraction_comparison"]["geometry_vs_atom_CIC"]["W"]["mean_absolute_fraction_difference"],
        } for row in phase_reports],
        "v5_mask_reassigned_voxels": changed_voxels,
    }, indent=2))


if __name__ == "__main__":
    main()
