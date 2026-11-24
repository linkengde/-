#!/usr/bin/env python3
"""Summarize V9 NVT/NVE state, atom retention, fixed IDs and displacements."""
from __future__ import annotations

import gzip
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
REPO = Path(__file__).resolve().parents[5]
TTM = MODEL / "ttm"
sys.path.insert(0, str(TTM))
from diagnose_agw_ttm_mapping_v7 import read_dump, read_lammps_data_atoms  # noqa: E402

LOG = HERE / "v9_cold_nvt_nve.log"
START_DUMP = MODEL / "interface_reconstruction_v6_exploratory/nve_fixedids_after_nvt_final.dump.gz"
NVE_END = HERE / "v9_nve_end.dump.gz"
TEMPLATE = MODEL / "interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz"
OUT = HERE / "v9_cold_baseline_metrics.json"


def read_thermo(path: Path) -> tuple[list[str], list[list[float]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    header = None
    rows: list[list[float]] = []
    for line in lines:
        fields = line.split()
        if fields and fields[0] == "Step" and "c_T_mobile" in fields:
            header = fields
            continue
        if header is None or not fields:
            continue
        try:
            values = [float(x) for x in fields]
        except ValueError:
            continue
        if len(values) == len(header):
            rows.append(values)
    if not rows or header is None:
        raise RuntimeError(f"No thermostat thermo rows found in {path}")
    # Run boundaries repeat the same step; retain the last row at each step.
    by_step: dict[int, list[float]] = {}
    for row in rows:
        by_step[int(row[0])] = row
    return header, [by_step[k] for k in sorted(by_step)]


def read_dump_loose(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if path.suffix == ".gz":
        return read_dump(path)
    with path.open("rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("ITEM: ATOMS"):
                columns = line.split()[2:]
                break
        else:
            raise ValueError(f"No ATOMS section in {path}")
        table = np.loadtxt(stream)
    ix = {name: columns.index(name) for name in ("id", "type", "x", "y", "z")}
    order = np.argsort(table[:, ix["id"]].astype(np.int64))
    return (table[order, ix["id"]].astype(np.int64),
            table[order, ix["type"]].astype(np.int32),
            table[order][:, [ix["x"], ix["y"], ix["z"]]].astype(float))


def main() -> None:
    header, rows = read_thermo(LOG)
    idx = {name: header.index(name) for name in (
        "Step", "Time", "Atoms", "PotEng", "KinEng", "TotEng", "c_T_mobile",
        "c_T_W", "c_T_Ag", "c_Fmax_mobile", "v_Frms_mobile", "c_z_min", "c_z_max",
    )}
    table = np.asarray(rows, dtype=float)
    nvt = table[(table[:, idx["Step"]] >= 0) & (table[:, idx["Step"]] <= 2000)]
    at = {int(r[idx["Step"]]): r for r in table}
    if not {1200, 1600, 2000, 2200}.issubset(at):
        raise RuntimeError(f"Missing expected NVT/NVE thermo steps: {sorted(at)}")
    ids0, types0, xyz0 = read_dump(START_DUMP)
    ids_nve, types_nve, xyz_nve = read_dump_loose(NVE_END)
    data_ids, data_types, data_xyz = read_lammps_data_atoms(TEMPLATE)
    if not (np.array_equal(ids0, ids_nve) and np.array_equal(types0, types_nve)
            and np.array_equal(data_ids, ids0) and np.array_equal(data_types, types0)):
        raise RuntimeError("V9 run changed atom IDs or element types")
    fixed_ids = data_ids[data_xyz[:, 2] < 20.0]
    fixed = np.isin(ids_nve, fixed_ids)
    if int(fixed.sum()) != 46117:
        raise RuntimeError(f"Fixed-ID count changed: {int(fixed.sum())}")
    fixed_max_move = float(np.max(np.linalg.norm(xyz_nve[fixed] - xyz0[fixed], axis=1)))
    delta = xyz_nve - xyz0
    delta[:, 0:2] = (delta[:, 0:2] + 100.0) % 200.0 - 100.0
    displacement = np.linalg.norm(delta, axis=1)
    end = at[2200]
    initial = at[0]
    nvt_end = at[2000]
    v6_ratio = float(np.count_nonzero(types_nve == 2) * 107.8682 / (np.count_nonzero(types_nve == 1) * 183.84))
    nvt_series = nvt[:, idx["c_T_mobile"]]
    result = {
        "title": "Ag-W V9 cold NVT/NVE continuation",
        "status": "ANALYZED",
        "settings": {
            "starting_restart": str(START_DUMP.relative_to(MODEL)),
            "starting_step": 4200,
            "fixed_group_source": "original data IDs whose initial z < 20 A",
            "fixed_atoms": int(fixed.sum()),
            "mobile_atoms": int(len(ids_nve) - fixed.sum()),
            "NVT_duration_ps": 0.5,
            "NVE_duration_ps": 0.05,
            "timestep_fs": 0.25,
            "TTM_or_source": False,
        },
        "atom_integrity": {
            "atoms_start_end": [int(len(ids0)), int(len(ids_nve))],
            "ids_and_types_unchanged": True,
            "Ag_W_mass_ratio": v6_ratio,
            "original_box_A": [200.0, 200.0, 300.0],
            "fixed_id_max_displacement_A": fixed_max_move,
            "mobile_max_displacement_A": float(displacement[~fixed].max()),
            "mobile_atoms_displaced_over_1A": int(np.count_nonzero(displacement[~fixed] > 1.0)),
            "mobile_atoms_displaced_over_2A": int(np.count_nonzero(displacement[~fixed] > 2.0)),
            "z_max_start_A": float(xyz0[:, 2].max()),
            "z_max_NVT_end_A": float(nvt_end[idx["c_z_max"]]),
            "z_max_NVE_end_A": float(xyz_nve[:, 2].max()),
        },
        "NVT": {
            "mobile_temperature_start_K": float(initial[idx["c_T_mobile"]]),
            "mobile_temperature_end_K": float(nvt_end[idx["c_T_mobile"]]),
            "Ag_temperature_end_K": float(nvt_end[idx["c_T_Ag"]]),
            "W_temperature_end_K": float(nvt_end[idx["c_T_W"]]),
            "mobile_temperature_change_last_0p1ps_K": float(nvt_end[idx["c_T_mobile"]] - at[1600][idx["c_T_mobile"]]),
            "mobile_temperature_change_last_0p2ps_K": float(nvt_end[idx["c_T_mobile"]] - at[1200][idx["c_T_mobile"]]),
            "mobile_temperature_range_K": [float(nvt_series.min()), float(nvt_series.max())],
            "potential_energy_change_eV": float(nvt_end[idx["PotEng"]] - initial[idx["PotEng"]]),
        },
        "NVE": {
            "start_total_energy_eV": float(nvt_end[idx["TotEng"]]),
            "end_total_energy_eV": float(end[idx["TotEng"]]),
            "total_energy_drift_eV": float(end[idx["TotEng"]] - nvt_end[idx["TotEng"]]),
            "mobile_temperature_start_K": float(nvt_end[idx["c_T_mobile"]]),
            "mobile_temperature_end_K": float(end[idx["c_T_mobile"]]),
            "maximum_mobile_force_eV_A": float(end[idx["c_Fmax_mobile"]]),
            "rms_mobile_force_eV_A": float(end[idx["v_Frms_mobile"]]),
        },
        "interpretation": {
            "temperature_target_range_285_to_315K": bool(285.0 <= nvt_end[idx["c_T_mobile"]] <= 315.0),
            "temperature_plateau": "review NVT tail deltas; a 300 K endpoint alone is not a balance pass",
            "cold_baseline_gate": "PENDING_REVIEW",
            "ttm_pair_authorized_by_this_result": False,
        },
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
