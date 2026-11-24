#!/usr/bin/env python3
"""Summarize the candidate's no-source NVT/NVE run and exact energy drift."""
from __future__ import annotations

import gzip
import csv
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from diagnose_agw_ttm_mapping_v7 import read_dump, read_lammps_data_atoms  # noqa: E402

OUT = HERE / "cold_baseline_metrics.json"
MASS = {1: 183.84, 2: 107.8682}
BOX_XY = 202.0
FIXED_Z_CEILING = 20.0


def parse_thermo(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def parse_energy(path: Path) -> dict[int, dict[str, float]]:
    result = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split()
            if len(parts) != 4:
                continue
            try:
                step = int(parts[0])
                pe, ke, etot = map(float, parts[1:])
            except ValueError:
                continue
            result[step] = {"potential_eV": pe, "kinetic_eV": ke, "total_eV": etot}
    return result


def displacement_report(start_path: Path, end_path: Path, fixed_ids: np.ndarray) -> dict:
    ids0, types0, xyz0 = read_dump(start_path)
    ids1, types1, xyz1 = read_dump(end_path)
    if not np.array_equal(ids0, ids1) or not np.array_equal(types0, types1):
        raise RuntimeError("Cold-check dumps changed the atom ID/type set")
    delta = xyz1 - xyz0
    delta[:, :2] = (delta[:, :2] + BOX_XY / 2.0) % BOX_XY - BOX_XY / 2.0
    disp = np.linalg.norm(delta, axis=1)
    fixed = np.isin(ids0, fixed_ids)
    mobile = ~fixed
    out = {
        "all_max_A": float(disp.max()),
        "mobile_max_A": float(disp[mobile].max()),
        "mobile_rms_A": float(np.sqrt(np.mean(disp[mobile] ** 2))),
        "mobile_atoms_over_1A": int(np.count_nonzero(disp[mobile] > 1.0)),
        "fixed_max_A": float(disp[fixed].max()),
        "by_phase": {},
    }
    for typ, element in ((2, "Ag"), (1, "W")):
        phase = (types0 == typ) & mobile
        out["by_phase"][element] = {
            "mobile_max_A": float(disp[phase].max()),
            "mobile_rms_A": float(np.sqrt(np.mean(disp[phase] ** 2))),
            "mobile_atoms_over_1A": int(np.count_nonzero(disp[phase] > 1.0)),
        }
    return out


def main() -> None:
    rows = parse_thermo(HERE / "cold_baseline_thermo.csv")
    by_step = {int(row["step"]): row for row in rows}
    start = by_step[4200]
    nvt_end = by_step[6200]
    nve_end = by_step[6400]
    precise = parse_energy(HERE / "energy_precision_full.lammps.log")
    e_start, e_nvt, e_nve = precise[4200], precise[6200], precise[6400]
    nve_drift = e_nve["total_eV"] - e_nvt["total_eV"]
    data_ids, _, data_xyz = read_lammps_data_atoms(HERE / "candidate.data.gz")
    fixed_ids = data_ids[(data_xyz[:, 2] >= 0.0) & (data_xyz[:, 2] < FIXED_Z_CEILING)]
    displacements = displacement_report(HERE / "cold_start.dump.gz", HERE / "cold_after_nve.dump.gz", fixed_ids)
    report = {
        "title": "Ag-W v8 tuned-copy no-source cold baseline check",
        "status": "COLD_BASELINE_FAIL_INTERPHASE_TEMPERATURE_MISMATCH",
        "input_candidate": "candidate_xy202A_atoms_minus1pct/candidate.data.gz plus candidate.restart.dump.gz",
        "settings": {
            "boundary": "p p f",
            "box_A": [202.0, 202.0, 300.0],
            "frozen_bottom_IDs": "selected once from candidate.data.gz at 0<=z<20 A",
            "fixed_atom_count": 46117,
            "mobile_atom_count": 431658,
            "temperature_target_K": 300.0,
            "NVT_duration_ps": 0.5,
            "NVT_damping_ps": 0.1,
            "NVE_duration_ps": 0.05,
            "timestep_fs": 0.25,
            "laser_source": "none",
            "TTM_steps": 0,
        },
        "atom_count_start_nvt_end_nve_end": [477775, 477775, 477775],
        "thermo": {
            "NVT_start_step4200": start,
            "NVT_end_step6200": nvt_end,
            "NVE_end_step6400": nve_end,
            "sampled_rows": rows,
            "sampling_note": "Captured from the LAMMPS execution output while it was running.",
        },
        "precise_energy_eV": {
            "NVT_start": e_start,
            "NVT_end_NVE_start": e_nvt,
            "NVE_end": e_nve,
            "thermostat_energy_added_during_NVT": float(e_nvt["total_eV"] - e_start["total_eV"]),
            "NVE_total_energy_drift_eV": float(nve_drift),
            "NVE_total_energy_drift_eV_per_atom": float(nve_drift / 477775),
        },
        "displacement_from_NVT_start_to_NVE_end": displacements,
        "decision": {
            "atom_count_preserved_during_run": True,
            "energy_stability": "NVE drift is small at 0.123 eV total, about 2.58e-7 eV/atom over 0.05 ps.",
            "phase_temperature_difference_K_at_NVE_end": float(nve_end["T_W_K"] - nve_end["T_Ag_K"]),
            "temperature_gate_pass": False,
            "TTM_or_laser_authorized_by_this_result": False,
        },
        "artifacts": {
            "lammps_input": "in.cold_baseline_0p5ps",
            "sampled_thermo": "cold_baseline_thermo.csv",
            "energy_precision_input": "in.energy_precision_check",
            "energy_precision_log": "energy_precision_full.lammps.log",
            "snapshots": ["cold_start.dump.gz", "cold_after_nvt.dump.gz", "cold_after_nve.dump.gz"],
            "full_cold_run_log_retained": False,
            "log_note": "The default LAMMPS log was replaced by the subsequent zero-step precision energy check; endpoint energies are retained at full precision and sampled thermo rows were preserved from execution output.",
        },
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "NVE_energy_drift_eV": nve_drift,
        "NVE_phase_temperatures_K": {"Ag": nve_end["T_Ag_K"], "W": nve_end["T_W_K"]},
        "NVE_phase_temperature_difference_K": report["decision"]["phase_temperature_difference_K_at_NVE_end"],
        "displacement": displacements,
    }, indent=2))


if __name__ == "__main__":
    main()
