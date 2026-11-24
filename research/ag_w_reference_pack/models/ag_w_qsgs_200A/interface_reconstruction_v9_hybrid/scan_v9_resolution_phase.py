#!/usr/bin/env python3
"""Full-atom V9 mesh-resolution and 2.5 A x/y phase scan, no structure edits."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np

MODEL = Path(__file__).resolve().parents[1]
TTM = MODEL / "ttm"
OUTDIR = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(TTM))
sys.path.insert(0, str(MODEL / "ttm/v8_grid_geometry_mapping"))

from agw_ttm_reference import (  # noqa: E402
    AG_ATOMIC_VOLUME_A3,
    W_ATOMIC_VOLUME_A3,
    _raw_cic_volume_power,
    gaussian_ag_source_power,
)
from diagnose_agw_ttm_mapping_v7 import (  # noqa: E402
    V6,
    local_capacity_flow,
    read_dump,
    read_lammps_data_atoms,
)

OUT = OUTDIR / "v9_resolution_phase_scan.json"
TRANSFER_A = 20.0
SOURCE_BINS_A = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 60.0, 100.0, 150.0])
RESTART = V6 / "nve_fixedids_after_nvt_final.dump.gz"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def profile(values: np.ndarray, dx: float, dy: float, spacing: float, nz: int, exact_atoms=None) -> dict:
    nx, ny = values.shape[:2]
    centers_x = (dx + (np.arange(nx) + 0.5) * spacing) % 200.0
    centers_y = (dy + (np.arange(ny) + 0.5) * spacing) % 200.0
    xx, yy = np.meshgrid(centers_x, centers_y, indexing="ij")
    dxx = (xx - 100.0 + 100.0) % 200.0 - 100.0
    dyy = (yy - 100.0 + 100.0) % 200.0 - 100.0
    r = np.sqrt(dxx**2 + dyy**2)[:, :, None]
    depth = (200.0 - (40.0 + (np.arange(nz) + 0.5) * spacing))[None, None, :]
    total = float(values.sum())
    radial = np.digitize(np.broadcast_to(r, values.shape).ravel(), SOURCE_BINS_A[1:-1])
    radial_power = np.bincount(radial, weights=values.ravel(), minlength=len(SOURCE_BINS_A) - 1)
    return {
        "total_eV_ps": total,
        "radial_bins_A": SOURCE_BINS_A.tolist(),
        "radial_power_fraction": (radial_power / total).tolist() if total else None,
        "radial_RMS_A": float(np.sqrt(np.sum(values * r**2) / total)) if total else None,
        "depth_mean_A": float(np.sum(values * depth) / total) if total else None,
    }


def main() -> None:
    manifest = json.loads((V6 / "exploratory_ttm_state_manifest.json").read_text(encoding="utf-8"))
    data = REPO / manifest["data_file"]
    ids0, types0, xyz0 = read_lammps_data_atoms(data)
    ids, types, xyz = read_dump(RESTART)
    if not (np.array_equal(ids0, ids) and np.array_equal(types0, types)):
        raise RuntimeError("V7 full-atom template and frozen restart IDs/types differ")
    static_newton = (xyz0[:, 2] >= 40.0 - 1e-9) & (xyz0[:, 2] <= 300.0 + 1e-9)
    inside = static_newton & (xyz[:, 2] >= 40.0 - 1e-9) & (xyz[:, 2] <= 200.0 + 1e-9)
    idx_ag = np.flatnonzero(inside & (types == 2))
    idx_w = np.flatnonzero(inside & (types == 1))
    q_ag = gaussian_ag_source_power(xyz[idx_ag])
    exact_source = float(np.dot(np.full(len(idx_ag), AG_ATOMIC_VOLUME_A3), q_ag))
    phases = []

    for grid in (5.0, 10.0, 20.0):
        nxy = int(round(200.0 / grid))
        nz = int(round(160.0 / grid))
        shape = (nxy, nxy, nz)
        spacing = (grid, grid, grid)
        offsets = np.arange(0.0, grid, 2.5)
        for dx in offsets:
            for dy in offsets:
                rel_ag = xyz[idx_ag] - np.array([dx, dy, 40.0])
                rel_w = xyz[idx_w] - np.array([dx, dy, 40.0])
                map_ag, vol_ag, source_ag = _raw_cic_volume_power(
                    rel_ag, q_ag, shape, spacing, AG_ATOMIC_VOLUME_A3
                )
                map_w, vol_w, _ = _raw_cic_volume_power(
                    rel_w, np.zeros(len(idx_w)), shape, spacing, W_ATOMIC_VOLUME_A3
                )
                cell_volume = grid**3
                combined = vol_ag + vol_w
                raw_profile = profile(source_ag, dx, dy, grid, nz)
                phases.append({
                    "grid_spacing_A": grid,
                    "shape": shape,
                    "offset_xy_A": [float(dx), float(dy)],
                    "raw_fill": {
                        "Ag_max": float((vol_ag / cell_volume).max()),
                        "W_max": float((vol_w / cell_volume).max()),
                        "combined_max": float((combined / cell_volume).max()),
                        "overfilled_cells": int(np.count_nonzero(combined > cell_volume + 1e-9)),
                        "overflow_A3": float(np.maximum(combined - cell_volume, 0.0).sum()),
                        "Ag_volume_A3": float(vol_ag.sum()),
                        "W_volume_A3": float(vol_w.sum()),
                    },
                    "source_raw": raw_profile,
                    "source_power_conservation_error_eV_ps": float(source_ag.sum() - exact_source),
                })

    # Same direct maximum-flow diagnostic is applied to every 10/20 A phase
    # and to the best raw-fill 5 A phase. This keeps the expensive 5 A graph
    # bounded while still comparing the physically relevant candidate.
    flow_candidates = []
    for grid in (5.0, 10.0, 20.0):
        rows = [r for r in phases if r["grid_spacing_A"] == grid]
        if grid == 5.0:
            selected = sorted(rows, key=lambda r: (
                r["raw_fill"]["overflow_A3"],
                r["raw_fill"]["overfilled_cells"],
                r["raw_fill"]["combined_max"],
            ))[:1]
        else:
            selected = rows
        for row in selected:
            dx, dy = row["offset_xy_A"]
            grid = float(row["grid_spacing_A"])
            shape = tuple(int(v) for v in row["shape"])
            rel_ag = xyz[idx_ag] - np.array([dx, dy, 40.0])
            rel_w = xyz[idx_w] - np.array([dx, dy, 40.0])
            map_ag, vol_ag, source_ag = _raw_cic_volume_power(
                rel_ag, q_ag, shape, (grid,)*3, AG_ATOMIC_VOLUME_A3
            )
            map_w, vol_w, _ = _raw_cic_volume_power(
                rel_w, np.zeros(len(idx_w)), shape, (grid,)*3, W_ATOMIC_VOLUME_A3
            )
            raw_total = vol_ag + vol_w
            flow = local_capacity_flow(raw_total, (grid,)*3, TRANSFER_A)
            row["direct_local_capacity_flow_20A"] = flow
            row["flow_candidate"] = True
            flow_candidates.append({
                "grid_spacing_A": grid,
                "shape": shape,
                "offset_xy_A": [dx, dy],
                "flow": flow,
            })

    for grid in (5.0, 10.0, 20.0):
        candidates = [r for r in phases if r["grid_spacing_A"] == grid]
        for row in candidates:
            row.setdefault("flow_candidate", False)
            row.setdefault("direct_local_capacity_flow_20A", None)

    result = {
        "title": "Ag-W V9 full-atom grid resolution and phase scan",
        "frozen_structure": {
            "template_data": manifest["data_file"],
            "template_sha256": manifest["data_sha256"],
            "restart_dump": str(RESTART.relative_to(MODEL)),
            "restart_sha256": sha256(RESTART),
            "atoms_total": int(len(ids)),
            "Ag_total": int(np.count_nonzero(types == 2)),
            "W_total": int(np.count_nonzero(types == 1)),
            "mapped_Ag": int(len(idx_ag)),
            "mapped_W": int(len(idx_w)),
            "Ag_W_mass_ratio": float(np.count_nonzero(types == 2) * 107.8682 / (np.count_nonzero(types == 1) * 183.84)),
            "box_A": [200.0, 200.0, 300.0],
            "source_reference_power_S1_eV_ps": exact_source,
            "atom_count_or_coordinates_modified": False,
        },
        "scan_definition": {
            "grid_sizes_A": [5.0, 10.0, 20.0],
            "phase_step_A": 2.5,
            "offsets_per_axis": {"5": [0.0, 2.5], "10": [0.0, 2.5, 5.0, 7.5], "20": list(np.arange(0, 20, 2.5))},
            "xy_periodic_z_open": True,
            "max_transfer_distance_A": TRANSFER_A,
            "phase_count": len(phases),
            "flow_policy": "Direct donor-to-original-free-recipient maximum flow; no cascaded redistribution.",
        },
        "flow_candidates": flow_candidates,
        "phases": phases,
        "interpretation": {
            "20A_resolution_warning": "At 20 A, the 30 A Gaussian radius spans 1.5 cells and 20 A skin depth spans one cell; use only if the spatial source profile passes sensitivity review.",
            "5A_flow_scan_scope": "Only the best raw-fill 5 A phase was sent to max flow due to graph size; all phase raw fills are retained.",
            "no_TTM_or_source_run": True,
        },
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "SCAN_COMPLETE",
        "phase_count": len(phases),
        "flow_candidates": flow_candidates,
        "best_raw_by_grid": [{
            "grid": g,
            "best_offset": min((r for r in phases if r["grid_spacing_A"] == g), key=lambda r: r["raw_fill"]["overflow_A3"])["offset_xy_A"],
            "raw_max": min((r for r in phases if r["grid_spacing_A"] == g), key=lambda r: r["raw_fill"]["overflow_A3"])["raw_fill"]["combined_max"],
            "overflow": min((r for r in phases if r["grid_spacing_A"] == g), key=lambda r: r["raw_fill"]["overflow_A3"])["raw_fill"]["overflow_A3"],
        } for g in (5.0, 10.0, 20.0)],
    }, indent=2))


if __name__ == "__main__":
    main()
