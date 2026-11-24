#!/usr/bin/env python3
"""Recompute selected full-atom phases with V9's production capacity solver."""
from __future__ import annotations

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
from diagnose_agw_ttm_mapping_v7 import V6, read_dump, read_lammps_data_atoms  # noqa: E402

SCAN = HERE / "v9_resolution_phase_scan.json"
OUT = HERE / "v9_shared_capacity_flow_verification.json"
RESTART = V6 / "nve_fixedids_after_nvt_final.dump.gz"
TRANSFER_A = 20.0


def main() -> None:
    scan = json.loads(SCAN.read_text(encoding="utf-8"))
    manifest = json.loads((V6 / "exploratory_ttm_state_manifest.json").read_text(encoding="utf-8"))
    data = REPO / manifest["data_file"]
    ids0, types0, xyz0 = read_lammps_data_atoms(data)
    ids, types, xyz = read_dump(RESTART)
    if not (np.array_equal(ids0, ids) and np.array_equal(types0, types)):
        raise RuntimeError("Frozen template and restart atom IDs/types differ")
    mask = ((xyz0[:, 2] >= 40.0 - 1.0e-9) & (xyz0[:, 2] <= 300.0 + 1.0e-9)
            & (xyz[:, 2] >= 40.0 - 1.0e-9) & (xyz[:, 2] <= 200.0 + 1.0e-9))
    i_ag = np.flatnonzero(mask & (types == 2))
    i_w = np.flatnonzero(mask & (types == 1))
    q_ag = gaussian_ag_source_power(xyz[i_ag])
    rows_out = []

    for h in (5.0, 10.0, 20.0):
        rows = [r for r in scan["phases"] if r["grid_spacing_A"] == h]
        flowed = [r for r in rows if r.get("direct_local_capacity_flow_20A")]
        if not flowed:
            raise RuntimeError(f"No prior direct-flow candidate found for {h:g} A")
        selected = min(
            flowed,
            key=lambda r: r["direct_local_capacity_flow_20A"]["locally_unassignable_overflow_A3"],
        )
        dx, dy = selected["offset_xy_A"]
        nxy, nz = int(round(200.0 / h)), int(round(160.0 / h))
        shape = (nxy, nxy, nz)
        origin = np.array([dx, dy, 40.0])
        _, v_ag, s_ag = _raw_cic_volume_power(
            xyz[i_ag] - origin, q_ag, shape, (h, h, h), AG_ATOMIC_VOLUME_A3
        )
        _, v_w, _ = _raw_cic_volume_power(
            xyz[i_w] - origin, np.zeros(len(i_w)), shape, (h, h, h), W_ATOMIC_VOLUME_A3
        )
        raw_total = v_ag + v_w
        flow = joint_capacity_maxflow(raw_total, (h, h, h), TRANSFER_A)
        flow.pop("transfer_fractions", None)
        rows_out.append({
            "grid_A": h,
            "shape": list(shape),
            "offset_xy_A": [dx, dy],
            "phase_selection": "best prior direct-flow unassigned-volume candidate at this resolution",
            "prior_diagnostic_resolution_A3": selected["direct_local_capacity_flow_20A"]["capacity_resolution_A3"],
            "prior_unassigned_A3": selected["direct_local_capacity_flow_20A"]["locally_unassignable_overflow_A3"],
            "raw_combined_max_fill": float(raw_total.max() / h**3),
            "raw_overflow_A3": float(np.maximum(raw_total - h**3, 0.0).sum()),
            "Ag_volume_conservation_error_A3": float(v_ag.sum() - len(i_ag) * AG_ATOMIC_VOLUME_A3),
            "W_volume_conservation_error_A3": float(v_w.sum() - len(i_w) * W_ATOMIC_VOLUME_A3),
            "Ag_source_power_eV_ps": float(s_ag.sum()),
            "shared_production_solver": flow,
        })

    result = {
        "title": "Ag-W V9 selected phases checked with the shared production max-flow solver",
        "input": {
            "template_data": manifest["data_file"],
            "template_sha256": manifest["data_sha256"],
            "restart_dump": str(RESTART.relative_to(MODEL)),
            "restart_sha256": manifest["restart_dump_sha256"],
            "atoms_total": int(len(ids)),
            "mapped_Ag": int(len(i_ag)),
            "mapped_W": int(len(i_w)),
            "coordinates_or_atoms_modified": False,
            "transfer_limit_A": TRANSFER_A,
        },
        "results": rows_out,
        "decision": "FAIL if no mapping resolution has a complete direct assignment; 20 A source resolution also requires spatial-profile review.",
        "ttm_or_laser_run": False,
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUT),
        "results": [{
            "grid_A": x["grid_A"], "offset_xy_A": x["offset_xy_A"],
            "status": x["shared_production_solver"]["status"],
            "unassignable_A3": x["shared_production_solver"]["locally_unassignable_overflow_A3"],
            "resolution_A3": x["shared_production_solver"]["capacity_resolution_A3"],
        } for x in rows_out],
    }, indent=2))


if __name__ == "__main__":
    main()
