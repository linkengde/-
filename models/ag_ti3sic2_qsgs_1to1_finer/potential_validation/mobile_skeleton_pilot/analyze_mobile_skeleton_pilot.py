#!/usr/bin/env python3
"""Summarize the short mobile-skeleton NVE pilot without claiming equilibration."""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1]
BOX = np.array([301.9555789947, 304.2005152739, 259.3915614456])
ELEMENT = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}


def load_dump(name: str) -> np.ndarray:
    path = HERE / name
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return np.loadtxt(handle, skiprows=9)


def parse_thermo() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    log = (HERE / "log.mobile_skeleton_nve_0p025ps").read_text(errors="replace")
    for line in log.splitlines():
        match = re.match(r"^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$", line)
        if match:
            vals = [float(v) for v in match.groups()]
            rows.append(dict(zip(("step", "time_ps", "T_mobile_K", "T_skeleton_K", "T_Ag_K", "PE_eV", "KE_eV", "E_total_eV", "Fmax_mobile_eV_A", "Fmax_skeleton_eV_A", "Fmax_Ag_eV_A"), vals)))
    if len(rows) < 2:
        raise RuntimeError(f"Expected initial and final thermo records; got {len(rows)}")
    return rows


def main() -> None:
    start = load_dump("start.dump.gz")
    end = load_dump("final.dump.gz")
    if start.shape != end.shape or not np.array_equal(start[:, :2], end[:, :2]):
        raise RuntimeError("Start/end atom IDs or types differ")
    typ = start[:, 1].astype(int)
    mobile = start[:, 4] >= 20.0
    dx = end[:, 2:5] - start[:, 2:5]
    for axis in (0, 1):
        dx[:, axis] -= BOX[axis] * np.rint(dx[:, axis] / BOX[axis])
    disp = np.linalg.norm(dx, axis=1)
    force = np.linalg.norm(end[:, 8:11], axis=1)
    thermo = parse_thermo()
    t0, tf = thermo[0], thermo[-1]
    element_results = {}
    for type_id, name in ELEMENT.items():
        mask = (typ == type_id) & mobile
        if not np.any(mask):
            continue
        element_results[name] = {
            "mobile_atoms": int(mask.sum()),
            "displacement_A_p50_p95_p99_max": [float(v) for v in np.percentile(disp[mask], [50, 95, 99, 100])],
            "fraction_displacement_gt_0p5A": float(np.mean(disp[mask] > 0.5)),
            "fraction_displacement_gt_1A": float(np.mean(disp[mask] > 1.0)),
            "end_force_eV_A_p50_p95_p99_max": [float(v) for v in np.percentile(force[mask], [50, 95, 99, 100])],
        }
    report = {
        "status": "EXPLORATORY_MOBILE_SKELETON_NVE_SCREEN",
        "decision": "COLD_BASELINE_NOT_ESTABLISHED; DO_NOT_START_HEAT_TRANSFER_OR_TTM_FROM_THIS_ENDPOINT",
        "simulation_duration_ps": float(tf["time_ps"] - t0["time_ps"]),
        "steps": int(tf["step"] - t0["step"]),
        "timestep_fs": 0.25,
        "boundary": "p p f",
        "fixed_region_A": "all atoms with initial z < 20",
        "mobile_counts": {ELEMENT[t]: int(np.sum((typ == t) & mobile)) for t in ELEMENT},
        "thermo_initial": t0,
        "thermo_final": tf,
        "total_energy_drift_eV": float(tf["E_total_eV"] - t0["E_total_eV"]),
        "total_energy_drift_fraction_of_abs_initial_eV": float(abs(tf["E_total_eV"] - t0["E_total_eV"]) / abs(t0["E_total_eV"])),
        "potential_energy_change_eV": float(tf["PE_eV"] - t0["PE_eV"]),
        "kinetic_energy_change_eV": float(tf["KE_eV"] - t0["KE_eV"]),
        "mobile_temperature_change_K": float(tf["T_mobile_K"] - t0["T_mobile_K"]),
        "maximum_mobile_force_change_eV_A": float(tf["Fmax_mobile_eV_A"] - t0["Fmax_mobile_eV_A"]),
        "mobile_atoms_displaced_over_1A": int(np.sum(disp[mobile] > 1.0)),
        "mobile_species_metrics": element_results,
        "interface_model": "Ag-skeleton WCA geometric wall from Model B; exploratory, repulsive-only, not calibrated",
        "decision_limits": [
            "No-source short mobility test only; no temperature gradient, TTM, heating, or melting test.",
            "The fixed bottom anchors atoms mechanically and is not a thermal reservoir.",
            "Short duration cannot establish solid-state stability or heat-sink performance.",
            "Nonzero PE-to-KE relaxation and high residual force must be interpreted as initial-state relaxation, not thermal transport.",
        ],
    }
    (HERE / "mobile_skeleton_pilot_metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
