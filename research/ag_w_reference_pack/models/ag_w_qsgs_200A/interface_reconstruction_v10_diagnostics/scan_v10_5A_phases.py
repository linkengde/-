#!/usr/bin/env python3
"""Run the V9 production direct-capacity solver on all four V10 5 A phases."""
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

V9 = MODEL / "interface_reconstruction_v9_hybrid"
V6 = MODEL / "interface_reconstruction_v6_exploratory"
DATA = MODEL / "interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz"
RESTART = V9 / "v9_nve2_end.dump.gz"
OUT = HERE / "v10_5A_phase_scan.json"
BIN_EDGES = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 60.0, 100.0, 150.0])


def source_profile(power: np.ndarray, dx: float, dy: float, h: float) -> dict:
    nx, ny, nz = power.shape
    xx = (dx + (np.arange(nx) + 0.5) * h) % 200.0
    yy = (dy + (np.arange(ny) + 0.5) * h) % 200.0
    x, y = np.meshgrid(xx, yy, indexing="ij")
    x = (x - 100.0 + 100.0) % 200.0 - 100.0
    y = (y - 100.0 + 100.0) % 200.0 - 100.0
    r = np.broadcast_to(np.sqrt(x*x + y*y)[:, :, None], power.shape)
    total = float(power.sum())
    bins = np.digitize(r.ravel(), BIN_EDGES[1:-1])
    radial = np.bincount(bins, weights=power.ravel(), minlength=len(BIN_EDGES)-1)
    z_centers = 40.0 + (np.arange(nz) + 0.5) * h
    depth = (200.0 - z_centers)[None, None, :]
    return {
        "power_eV_ps": total,
        "radial_bins_A": BIN_EDGES.tolist(),
        "radial_power_fraction": (radial / total).tolist() if total else None,
        "radial_RMS_A": float(np.sqrt(np.sum(power * r*r) / total)) if total else None,
        "mean_depth_A": float(np.sum(power * depth) / total) if total else None,
    }


def transfer_distance_summary(flow: dict, raw_total: np.ndarray, h: float) -> dict:
    transfers = flow.get("transfer_fractions", {})
    if not transfers:
        return {"actual_transfer_edges": 0, "actual_max_distance_A": None, "actual_mean_distance_A": None}
    ny, nz = raw_total.shape[1:]
    nx = raw_total.shape[0]
    distances, amounts = [], []
    for donor, targets in transfers.items():
        ix, iy, iz = int(donor) // (ny*nz), (int(donor)//nz) % ny, int(donor)%nz
        for target, frac in targets:
            tx, ty, tz = int(target)//(ny*nz), (int(target)//nz)%ny, int(target)%nz
            dx = min((tx-ix)%nx, (ix-tx)%nx) * h
            dy = min((ty-iy)%ny, (iy-ty)%ny) * h
            dz = (tz-iz) * h
            distances.append(float(np.sqrt(dx*dx+dy*dy+dz*dz)))
            amounts.append(float(raw_total.ravel()[int(donor)] * frac))
    if not distances:
        return {"actual_transfer_edges": 0, "actual_max_distance_A": None, "actual_mean_distance_A": None}
    return {
        "actual_transfer_edges": len(distances),
        "actual_max_distance_A": max(distances),
        "actual_mean_distance_A": float(np.average(distances, weights=amounts)),
        "actual_transferred_volume_A3": float(np.sum(amounts)),
    }


def main() -> None:
    ids_data, types_data, xyz_data = read_lammps_data_atoms(DATA)
    ids, types, xyz = read_dump(RESTART)
    if not (np.array_equal(ids_data, ids) and np.array_equal(types_data, types)):
        raise RuntimeError("V9 final restart does not preserve template IDs/types")
    static_newton = (xyz_data[:, 2] >= 40.0-1e-9) & (xyz_data[:, 2] <= 300.0+1e-9)
    inside = static_newton & (xyz[:, 2] >= 40.0-1e-9) & (xyz[:, 2] <= 200.0+1e-9)
    idx_ag = np.flatnonzero(inside & (types == 2))
    idx_w = np.flatnonzero(inside & (types == 1))
    q_ag = gaussian_ag_source_power(xyz[idx_ag])
    exact_source = float(np.dot(np.full(len(idx_ag), AG_ATOMIC_VOLUME_A3), q_ag))
    rows = []

    for dx, dy in ((0.0, 0.0), (0.0, 2.5), (2.5, 0.0), (2.5, 2.5)):
        origin = np.array([dx, dy, 40.0])
        _, v_ag, source = _raw_cic_volume_power(
            xyz[idx_ag]-origin, q_ag, (40, 40, 32), (5.0,)*3, AG_ATOMIC_VOLUME_A3
        )
        _, v_w, _ = _raw_cic_volume_power(
            xyz[idx_w]-origin, np.zeros(len(idx_w)), (40, 40, 32), (5.0,)*3, W_ATOMIC_VOLUME_A3
        )
        total = v_ag + v_w
        flow = joint_capacity_maxflow(total, (5.0,)*3, 20.0)
        distances = transfer_distance_summary(flow, total, 5.0)
        flow_without_arrays = {k:v for k,v in flow.items() if k != "transfer_fractions"}
        fill = (total / 125.0)
        overfill = np.maximum(total-125.0, 0.0)
        rows.append({
            "offset_xy_A": [dx, dy],
            "raw_fill_max_Ag": float((v_ag/125.0).max()),
            "raw_fill_max_W": float((v_w/125.0).max()),
            "raw_fill_max_combined": float(fill.max()),
            "raw_overfilled_cells": int(np.count_nonzero(fill > 1.0+1e-12)),
            "raw_overflow_demand_A3": float(overfill.sum()),
            "global_free_capacity_A3": float(np.maximum(125.0-total,0.0).sum()),
            "Ag_volume_conservation_error_A3": float(v_ag.sum()-len(idx_ag)*AG_ATOMIC_VOLUME_A3),
            "W_volume_conservation_error_A3": float(v_w.sum()-len(idx_w)*W_ATOMIC_VOLUME_A3),
            "source_power_conservation_error_eV_ps": float(source.sum()-exact_source),
            "source_spatial_profile": source_profile(source, dx, dy, 5.0),
            "production_capacity_flow": flow_without_arrays,
            "transfer_distance_diagnostics": distances,
        })
        print(json.dumps({"phase":[dx,dy],"status":flow["status"],"unassigned_A3":flow["locally_unassignable_overflow_A3"]}), flush=True)

    best = min(rows, key=lambda x: x["production_capacity_flow"]["locally_unassignable_overflow_A3"])
    v9_manifest = json.loads((V9/"v9_exploratory_state_manifest.json").read_text())
    result = {
        "title": "Ag-W V10 complete 5 A phase capacity-flow scan",
        "frozen_v9_commit": "579c855f25f8286a3cdb535ea7336be95962f994",
        "input": {
            "template_data": str(DATA.relative_to(REPO)),
            "template_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
            "restart": str(RESTART.relative_to(MODEL)),
            "restart_sha256": hashlib.sha256(RESTART.read_bytes()).hexdigest(),
            "atoms_total": int(len(ids)),
            "Ag_total": int(np.count_nonzero(types==2)),
            "W_total": int(np.count_nonzero(types==1)),
            "mapped_Ag": int(len(idx_ag)),
            "mapped_W": int(len(idx_w)),
            "mass_ratio_full_Ag_over_W": float(np.count_nonzero(types==2)*107.8682/(np.count_nonzero(types==1)*183.84)),
            "box_A": [200.0,200.0,300.0],
            "box_or_atom_changes": False,
            "atomic_volumes_unchanged": True,
            "transfer_limit_A": 20.0,
            "mapping_temperature_state_manifest_sha256": hashlib.sha256((V9/"v9_exploratory_state_manifest.json").read_bytes()).hexdigest(),
        },
        "grid_A": 5.0,
        "phase_step_A": 2.5,
        "phases_tested": len(rows),
        "algorithm": "joint_capacity_maxflow shared with V9 production mapper; direct donors to initially free cells; xy periodic, z open; no cascade",
        "phase_results": rows,
        "best_by_unassigned_overflow": {
            "offset_xy_A": best["offset_xy_A"],
            "unassigned_A3": best["production_capacity_flow"]["locally_unassignable_overflow_A3"],
            "overflow_demand_A3": best["production_capacity_flow"]["overflow_demand_A3"],
            "capacity_status": best["production_capacity_flow"]["status"],
        },
        "all_phases_pass": all(x["production_capacity_flow"]["status"]=="PASS" for x in rows),
        "no_ttm_or_laser_step": True,
    }
    OUT.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"output":str(OUT),"best":result["best_by_unassigned_overflow"],"all_pass":result["all_phases_pass"]},indent=2))


if __name__=="__main__": main()
