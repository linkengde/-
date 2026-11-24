#!/usr/bin/env python3
"""Summarize the second V9 cold NVT/NVE segment."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(HERE))
from analyze_v9_cold_baseline import read_dump_loose, read_thermo  # noqa: E402
sys.path.insert(0, str(MODEL / "ttm"))
from diagnose_agw_ttm_mapping_v7 import read_dump, read_lammps_data_atoms  # noqa: E402

LOG = HERE / "v9_cold_followup.log"
START = HERE / "v9_nve_end.dump.gz"
NVE_END = HERE / "v9_nve2_end.dump.gz"
TEMPLATE = MODEL / "interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz"
OUT = HERE / "v9_cold_followup_metrics.json"


def main() -> None:
    header, rows = read_thermo(LOG)
    ix = {name: header.index(name) for name in (
        "Step", "Atoms", "PotEng", "TotEng", "c_T_mobile", "c_T_W", "c_T_Ag",
        "c_Fmax_mobile", "v_Frms_mobile", "c_z_max",
    )}
    by_step = {int(row[ix["Step"]]): row for row in rows}
    initial, nvt, nve = by_step[0], by_step[2000], by_step[2200]
    ids0, types0, xyz0 = read_dump_loose(START)
    idsnve, typesnve, xyznve = read_dump_loose(NVE_END)
    data_ids, data_types, data_xyz = read_lammps_data_atoms(TEMPLATE)
    if not (np.array_equal(ids0, idsnve) and np.array_equal(types0, typesnve)
            and np.array_equal(data_ids, ids0) and np.array_equal(data_types, types0)):
        raise RuntimeError("Atom IDs/types changed during the V9 cold follow-up")
    fixed_ids = data_ids[data_xyz[:, 2] < 20.0]
    fixed = np.isin(idsnve, fixed_ids)
    if int(fixed.sum()) != 46117:
        raise RuntimeError("Original fixed-ID set changed")
    delta = xyznve - xyz0
    delta[:, :2] = (delta[:, :2] + 100.0) % 200.0 - 100.0
    disp = np.linalg.norm(delta, axis=1)
    sample025 = by_step[1000]
    sample030 = by_step[1200]
    result = {
        "title": "Ag-W V9 second cold NVT/NVE continuation",
        "settings": {"NVT_ps": 0.5, "NVE_ps": 0.05, "dt_fs": 0.25, "fixed_bottom_atoms": int(fixed.sum()), "TTM_or_source": False},
        "atom_integrity": {
            "start_end_atoms": [int(len(ids0)), int(len(idsnve))],
            "ids_types_unchanged": True,
            "box_A": [200.0, 200.0, 300.0],
            "fixed_id_max_displacement_A": float(np.max(np.linalg.norm(xyznve[fixed] - xyz0[fixed], axis=1))),
            "mobile_max_displacement_A": float(disp[~fixed].max()),
            "mobile_atoms_over_1A": int(np.count_nonzero(disp[~fixed] > 1.0)),
            "z_max_start_A": float(xyz0[:, 2].max()),
            "z_max_NVT_end_A": float(nvt[ix["c_z_max"]]),
            "z_max_NVE_end_A": float(xyznve[:, 2].max()),
        },
        "NVT": {
            "mobile_temperature_start_K": float(initial[ix["c_T_mobile"]]),
            "mobile_temperature_end_K": float(nvt[ix["c_T_mobile"]]),
            "Ag_temperature_end_K": float(nvt[ix["c_T_Ag"]]),
            "W_temperature_end_K": float(nvt[ix["c_T_W"]]),
            "Ag_W_temperature_difference_K": float(nvt[ix["c_T_W"]] - nvt[ix["c_T_Ag"]]),
            "mobile_temperature_change_last_0p1ps_K": float(nvt[ix["c_T_mobile"]] - by_step[1600][ix["c_T_mobile"]]),
            "mobile_temperature_change_last_0p2ps_K": float(nvt[ix["c_T_mobile"]] - by_step[1200][ix["c_T_mobile"]]),
            "potential_energy_change_eV": float(nvt[ix["PotEng"]] - initial[ix["PotEng"]]),
            "early_temperature_points_K": [float(sample025[ix["c_T_mobile"]]), float(sample030[ix["c_T_mobile"]])],
        },
        "NVE": {
            "total_energy_start_eV": float(nvt[ix["TotEng"]]),
            "total_energy_end_eV": float(nve[ix["TotEng"]]),
            "total_energy_drift_eV": float(nve[ix["TotEng"]] - nvt[ix["TotEng"]]),
            "mobile_temperature_end_K": float(nve[ix["c_T_mobile"]]),
            "Ag_temperature_end_K": float(nve[ix["c_T_Ag"]]),
            "W_temperature_end_K": float(nve[ix["c_T_W"]]),
            "maximum_mobile_force_eV_A": float(nve[ix["c_Fmax_mobile"]]),
            "rms_mobile_force_eV_A": float(nve[ix["v_Frms_mobile"]]),
        },
        "decision": {
            "cold_baseline_gate": "FAIL_PHASE_TEMPERATURE_SPLIT_AND_CONTINUING_POTENTIAL_RELAXATION",
            "ttm_pair_allowed": False,
            "reason": "The global mobile temperature is near 300 K, but Ag/W temperatures remain separated and potential energy continues to decrease during NVT.",
        },
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
