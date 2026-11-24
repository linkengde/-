#!/usr/bin/env python3
"""Analyze the 0.5 ps NVT plus 0.05 ps NVE cold-baseline follow-up."""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
BOX = np.array([200.0, 200.0, 1000.0])
ELEMENT = {1: "W", 2: "Ag"}


def read_dump(name: str) -> np.ndarray:
    path = HERE / name
    with gzip.open(path, "rt") as f:
        return np.loadtxt(f, skiprows=9)


def distance_metrics(pos: np.ndarray, typ: np.ndarray) -> dict[str, float]:
    pos = pos.copy()
    pos[:, :2] %= 200.0
    w = pos[typ == 1]
    ag = pos[typ == 2]
    wt = cKDTree(w, boxsize=BOX)
    at = cKDTree(ag, boxsize=BOX)
    return {
        "W_W_min_A": float(wt.query(w, k=2, workers=-1)[0][:, 1].min()),
        "Ag_Ag_min_A": float(at.query(ag, k=2, workers=-1)[0][:, 1].min()),
        "Ag_W_min_A": float(wt.query(ag, workers=-1)[0].min()),
    }


def thermo_rows() -> list[dict[str, float]]:
    names = ("step", "time_ps", "atoms", "PE_eV", "KE_eV", "Etot_eV", "T_mobile_K", "T_W_K", "T_Ag_K", "Fmax_mobile_eV_A", "Frms_mobile_eV_A", "zmin_A", "zmax_A")
    rows = []
    for line in (HERE / "log.nvt_0p5ps_nve_0p05ps").read_text(errors="replace").splitlines():
        fields = line.split()
        if len(fields) != len(names) or not fields[0].isdigit():
            continue
        try:
            rows.append(dict(zip(names, map(float, fields))))
        except ValueError:
            continue
    if len(rows) < 10:
        raise RuntimeError(f"Expected NVT/NVE thermo history, found {len(rows)} rows")
    return rows


def motion(start: np.ndarray, end: np.ndarray) -> dict[str, object]:
    typ = start[:, 1].astype(int)
    p0, p1 = start[:, 2:5].copy(), end[:, 2:5].copy()
    d = p1 - p0
    for ax in (0, 1):
        d[:, ax] -= 200.0 * np.rint(d[:, ax] / 200.0)
    disp = np.linalg.norm(d, axis=1)
    mobile = p0[:, 2] > 20.0
    out: dict[str, object] = {"mobile_atoms": int(mobile.sum()), "mobile_max_displacement_A": float(disp[mobile].max()), "mobile_atoms_over_1A": int(np.sum(disp[mobile] > 1.0))}
    for tid, name in ELEMENT.items():
        mask = mobile & (typ == tid)
        out[name] = {
            "mobile_atoms": int(mask.sum()),
            "displacement_p50_p95_p99_max_A": [float(v) for v in np.percentile(disp[mask], [50, 95, 99, 100])],
        }
    return out


def main() -> None:
    start = read_dump("nvt_start.dump.gz")
    nvt_end = read_dump("nvt_end.dump.gz")
    nve_end = read_dump("nve_after_nvt_final.dump.gz")
    if not (np.array_equal(start[:, :2], nvt_end[:, :2]) and np.array_equal(start[:, :2], nve_end[:, :2])):
        raise RuntimeError("Atom count/IDs/types differ among snapshots")
    typ = start[:, 1].astype(int)
    rows = thermo_rows()
    nvt_start_row = next(r for r in rows if r["step"] == 0)
    nvt_end_row = next(r for r in rows if r["step"] == 2000)
    nve_rows = [r for r in rows if r["step"] >= 2000]
    nve_start_row, nve_final_row = nve_rows[0], nve_rows[-1]
    nvt_rows = [r for r in rows if 0 <= r["step"] <= 2000]
    report = {
        "title": "Ag-W v6 cold-baseline NVT/NVE follow-up",
        "status": "NVT_TARGET_NOT_REACHED; SHORT_NVE_CHECK_COMPLETED",
        "decision": "LEVEL_2_COLD_BASELINE_NOT_ESTABLISHED; DO_NOT_START_TTM_OR_LASER",
        "settings": {"boundary": "p p f", "fixed_region": "atoms selected by z<=20 A at NVT start; group retained for both segments", "NVT_temperature_K": 300.0, "NVT_damping_ps": 0.1, "NVT_duration_ps": 0.5, "NVE_duration_ps": 0.05, "timestep_fs": 0.25},
        "fixed_and_mobile_counts": {"fixed_bottom": 46232, "mobile": 436369},
        "atom_counts_start_nvt_end_nvt_end_nve": [int(len(start)), int(len(nvt_end)), int(len(nve_end))],
        "NVT_start": nvt_start_row,
        "NVT_end": nvt_end_row,
        "NVT_mobile_temperature_min_max_K": [float(min(r["T_mobile_K"] for r in nvt_rows)), float(max(r["T_mobile_K"] for r in nvt_rows))],
        "NVT_total_energy_change_eV": float(nvt_end_row["Etot_eV"] - nvt_start_row["Etot_eV"]),
        "NVT_PE_change_eV": float(nvt_end_row["PE_eV"] - nvt_start_row["PE_eV"]),
        "NVT_KE_change_eV": float(nvt_end_row["KE_eV"] - nvt_start_row["KE_eV"]),
        "NVE_start": nve_start_row,
        "NVE_end": nve_final_row,
        "NVE_total_energy_drift_eV": float(nve_final_row["Etot_eV"] - nve_start_row["Etot_eV"]),
        "NVE_mobile_temperature_change_K": float(nve_final_row["T_mobile_K"] - nve_start_row["T_mobile_K"]),
        "distance_minima_A": {
            "NVT_start": distance_metrics(start[:, 2:5], typ),
            "NVT_end": distance_metrics(nvt_end[:, 2:5], typ),
            "NVE_end": distance_metrics(nve_end[:, 2:5], typ),
        },
        "displacement_from_NVT_start": {"NVT_end": motion(start, nvt_end), "NVE_end": motion(start, nve_end)},
        "limitations": [
            "The thermostat did not bring the mobile system close to 300 K within 0.5 ps; its final mobile temperature was about 197 K.",
            "Most net energy supplied during NVT went into the potential-energy increase, consistent with continuing structural relaxation.",
            "The subsequent NVE segment checks energy integration near 197 K, not a 300 K equilibrium baseline.",
            "No laser source, TTM, melting, or cooling run was performed.",
        ],
    }
    (HERE / "nvt_followup_metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
