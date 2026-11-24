#!/usr/bin/env python3
"""V11 padded-z mapping scan; no atom edits, MD integration, TTM step, or source."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
REPO = Path(__file__).resolve().parents[5]
TTM = MODEL / "ttm"
sys.path.insert(0, str(TTM))

from agw_ttm_reference import (  # noqa: E402
    AG_ATOMIC_VOLUME_A3,
    W_ATOMIC_VOLUME_A3,
    _raw_cic_volume_power,
    gaussian_ag_source_power,
    joint_capacity_maxflow,
)
from diagnose_agw_ttm_mapping_v7 import read_dump, read_lammps_data_atoms  # noqa: E402

DATA = MODEL / "interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz"
RESTART = MODEL / "interface_reconstruction_v9_hybrid/v9_nve2_end.dump.gz"
MANIFEST = MODEL / "interface_reconstruction_v9_hybrid/v9_exploratory_state_manifest.json"
OUT = HERE / "v11_padded_z_mapping_scan.json"

BOX_XY_A = 200.0
Z_MIN_A = 35.0
Z_MAX_A = 205.0
GRID_A = 5.0
TRANSFER_A = 20.0
PHASES = ((0.0, 0.0), (0.0, 2.5), (2.5, 0.0), (2.5, 2.5))
RADIAL_BINS_A = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 60.0, 100.0, 150.0])


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_profile(power: np.ndarray, dx: float, dy: float) -> dict:
    nx, ny, nz = power.shape
    xx = (dx + (np.arange(nx) + 0.5) * GRID_A) % BOX_XY_A
    yy = (dy + (np.arange(ny) + 0.5) * GRID_A) % BOX_XY_A
    x, y = np.meshgrid(xx, yy, indexing="ij")
    x = (x - BOX_XY_A / 2.0 + BOX_XY_A / 2.0) % BOX_XY_A - BOX_XY_A / 2.0
    y = (y - BOX_XY_A / 2.0 + BOX_XY_A / 2.0) % BOX_XY_A - BOX_XY_A / 2.0
    radius = np.broadcast_to(np.sqrt(x*x + y*y)[:, :, None], power.shape)
    total = float(power.sum())
    bins = np.digitize(radius.ravel(), RADIAL_BINS_A[1:-1])
    radial = np.bincount(bins, weights=power.ravel(), minlength=len(RADIAL_BINS_A)-1)
    zc = Z_MIN_A + (np.arange(nz) + 0.5) * GRID_A
    depth = np.maximum(0.0, 200.0 - zc)[None, None, :]
    return {
        "total_eV_ps": total,
        "radial_bins_A": RADIAL_BINS_A.tolist(),
        "radial_power_fraction": (radial / total).tolist() if total else None,
        "radial_RMS_A": float(np.sqrt(np.sum(power * radius**2) / total)) if total else None,
        "mean_depth_from_nominal_z200_A": float(np.sum(power * depth) / total) if total else None,
    }


def transfer_summary(flow: dict, raw_total: np.ndarray) -> dict:
    transfers = flow.get("transfer_fractions", {})
    distances, amounts = [], []
    nx, ny, nz = raw_total.shape
    for donor, targets in transfers.items():
        ix, iy, iz = int(donor)//(ny*nz), (int(donor)//nz)%ny, int(donor)%nz
        for target, fraction in targets:
            tx, ty, tz = int(target)//(ny*nz), (int(target)//nz)%ny, int(target)%nz
            dx = min((tx-ix)%nx, (ix-tx)%nx) * GRID_A
            dy = min((ty-iy)%ny, (iy-ty)%ny) * GRID_A
            dz = (tz-iz) * GRID_A
            distances.append(float(np.sqrt(dx*dx + dy*dy + dz*dz)))
            amounts.append(float(raw_total.ravel()[int(donor)] * fraction))
    if not distances:
        return {"edges": 0, "transferred_volume_A3": 0.0, "mean_distance_A": 0.0, "max_distance_A": 0.0}
    return {
        "edges": len(distances),
        "transferred_volume_A3": float(np.sum(amounts)),
        "mean_distance_A": float(np.average(distances, weights=amounts)),
        "max_distance_A": float(max(distances)),
    }


def main() -> None:
    for path in (DATA, RESTART, MANIFEST):
        if not path.is_file():
            raise FileNotFoundError(path)

    data_ids, data_types, data_xyz = read_lammps_data_atoms(DATA)
    ids, types, xyz = read_dump(RESTART)
    if not (np.array_equal(data_ids, ids) and np.array_equal(data_types, types)):
        raise RuntimeError("V9 restart does not preserve the frozen template atom IDs/types")

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if sha256(RESTART).lower() != str(manifest["restart_dump_sha256"]).lower():
        raise RuntimeError("V9 restart SHA-256 does not match its frozen manifest")

    # Preserve the exact static Newton-group membership used by V7–V10.
    static_newton = (data_xyz[:, 2] >= 40.0 - 1e-9) & (data_xyz[:, 2] <= 300.0 + 1e-9)
    inside = static_newton & (xyz[:, 2] >= Z_MIN_A - 1e-9) & (xyz[:, 2] <= Z_MAX_A + 1e-9)
    outside = static_newton & ~inside
    outside_by_side = {
        "below_z_min": int(np.count_nonzero(static_newton & (xyz[:, 2] < Z_MIN_A - 1e-9))),
        "above_z_max": int(np.count_nonzero(static_newton & (xyz[:, 2] > Z_MAX_A + 1e-9))),
    }
    idx_ag = np.flatnonzero(inside & (types == 2))
    idx_w = np.flatnonzero(inside & (types == 1))
    if not len(idx_ag) or not len(idx_w):
        raise RuntimeError("Padded mapping domain must contain both Ag and W")

    shape = (int(round(BOX_XY_A / GRID_A)), int(round(BOX_XY_A / GRID_A)),
             int(round((Z_MAX_A - Z_MIN_A) / GRID_A)))
    if shape != (40, 40, 34):
        raise RuntimeError(f"Unexpected V11 padded mesh shape: {shape}")
    spacing = (GRID_A, GRID_A, GRID_A)
    q_ag = gaussian_ag_source_power(
        xyz[idx_ag],
        center_xy_A=(BOX_XY_A/2.0, BOX_XY_A/2.0),
        surface_z_A=200.0,
        allow_above_surface=True,
        periodic_box_xy_A=BOX_XY_A,
    )
    exact_source = float(AG_ATOMIC_VOLUME_A3 * q_ag.sum())
    phase_rows = []

    for dx, dy in PHASES:
        origin = np.array([dx, dy, Z_MIN_A])
        _, v_ag, source = _raw_cic_volume_power(
            xyz[idx_ag] - origin, q_ag, shape, spacing, AG_ATOMIC_VOLUME_A3
        )
        _, v_w, _ = _raw_cic_volume_power(
            xyz[idx_w] - origin, np.zeros(len(idx_w)), shape, spacing, W_ATOMIC_VOLUME_A3
        )
        raw_total = v_ag + v_w
        flow = joint_capacity_maxflow(raw_total, spacing, TRANSFER_A)
        over = np.maximum(raw_total - GRID_A**3, 0.0)
        layer_over = over.sum(axis=(0, 1))
        top = []
        for flat in np.argsort(over.ravel())[-30:][::-1]:
            value = float(over.ravel()[flat])
            if value <= 0.0:
                break
            i, j, k = np.unravel_index(int(flat), over.shape)
            top.append({
                "ijk": [int(i), int(j), int(k)],
                "center_A": [dx+(i+0.5)*GRID_A, dy+(j+0.5)*GRID_A, Z_MIN_A+(k+0.5)*GRID_A],
                "raw_fill": float(raw_total[i, j, k] / GRID_A**3),
                "overflow_A3": value,
            })
        phase_rows.append({
            "offset_xy_A": [dx, dy],
            "raw_max_fill": float(raw_total.max() / GRID_A**3),
            "raw_overflow_A3": float(over.sum()),
            "overfilled_cells": int(np.count_nonzero(over > 1e-12)),
            "overflow_by_z_layer_A3": layer_over.tolist(),
            "top_30_overfilled_cells": top,
            "Ag_volume_conservation_error_A3": float(v_ag.sum() - len(idx_ag)*AG_ATOMIC_VOLUME_A3),
            "W_volume_conservation_error_A3": float(v_w.sum() - len(idx_w)*W_ATOMIC_VOLUME_A3),
            "source_grid_power_eV_ps": float(source.sum()),
            "source_atom_power_eV_ps": exact_source,
            "source_conservation_error_eV_ps": float(source.sum() - exact_source),
            "source_profile": source_profile(source, dx, dy),
            "capacity_flow": {k:v for k,v in flow.items() if k != "transfer_fractions"},
            "transfer_summary": transfer_summary(flow, raw_total),
        })
        print(json.dumps({
            "phase_xy_A": [dx, dy],
            "flow_status": flow["status"],
            "unassigned_A3": flow["locally_unassignable_overflow_A3"],
            "z_outside_atoms": int(outside.sum()),
        }), flush=True)

    best = min(phase_rows, key=lambda row: row["capacity_flow"]["locally_unassignable_overflow_A3"])
    all_flow_pass = all(row["capacity_flow"]["status"] == "PASS" for row in phase_rows)
    domain_pass = int(outside.sum()) == 0
    result = {
        "title": "Ag-W V11 padded-z mapping phase scan",
        "status": "PASS_MAPPING_SCAN" if all_flow_pass and domain_pass else "FAIL_MAPPING_SCAN",
        "purpose": "Conservative mapping-only diagnostic; no atoms modified, no MD integration, no TTM step, no laser source.",
        "input": {
            "data": str(DATA.relative_to(REPO)),
            "data_sha256": sha256(DATA),
            "restart": str(RESTART.relative_to(MODEL)),
            "restart_sha256": sha256(RESTART),
            "atoms_total": int(len(ids)),
            "Ag_total": int(np.count_nonzero(types == 2)),
            "W_total": int(np.count_nonzero(types == 1)),
            "Ag_W_mass_ratio": float(np.count_nonzero(types == 2)*107.8682 /
                                     (np.count_nonzero(types == 1)*183.84)),
            "box_A": [200.0, 200.0, 300.0],
            "coordinates_or_atoms_modified": False,
        },
        "mapping": {
            "z_bounds_A": [Z_MIN_A, Z_MAX_A],
            "grid_shape": list(shape),
            "spacing_A": list(spacing),
            "transfer_limit_A": TRANSFER_A,
            "phase_xy_A": [list(p) for p in PHASES],
            "static_newton_atoms": int(static_newton.sum()),
            "mapped_atoms": int(inside.sum()),
            "excluded_static_newton_atoms": int(outside.sum()),
            "excluded_by_side": outside_by_side,
            "domain_coverage_pass": domain_pass,
            "guard_cells": {"below_old_z40_boundary": 1, "above_old_z200_boundary": 1},
            "atomic_volumes_unchanged": True,
        },
        "phase_results": phase_rows,
        "best_phase": {
            "offset_xy_A": best["offset_xy_A"],
            "unassigned_A3": best["capacity_flow"]["locally_unassignable_overflow_A3"],
            "status": best["capacity_flow"]["status"],
        },
        "all_four_phases_pass": all_flow_pass,
        "source_surface": {
            "nominal_z_A": 200.0,
            "allow_above_surface_clamp": True,
            "Ag_atoms_above_nominal_surface": int(np.count_nonzero(xyz[idx_ag, 2] > 200.0 + 1e-9)),
            "note": "Exploratory flat-plane approximation; atoms above z=200 have source depth clamped to zero. This is not a moving-free-surface model.",
        },
        "ttm_steps_executed": 0,
        "laser_source_applied": False,
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "best_phase": result["best_phase"],
        "domain_coverage_pass": domain_pass,
        "output": str(OUT),
    }, indent=2))


if __name__ == "__main__":
    main()
