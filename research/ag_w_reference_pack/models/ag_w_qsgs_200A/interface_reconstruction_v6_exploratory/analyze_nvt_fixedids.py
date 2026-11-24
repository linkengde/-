#!/usr/bin/env python3
"""Analyze the corrected original-fixed-ID NVT/NVE continuation."""
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
    with gzip.open(HERE / name, "rt") as f:
        return np.loadtxt(f, skiprows=9)


def thermo_rows() -> list[dict[str, float]]:
    names = ("step", "time_ps", "atoms", "PE_eV", "KE_eV", "Etot_eV", "T_mobile_K", "T_W_K", "T_Ag_K", "Fmax_mobile_eV_A", "Frms_mobile_eV_A", "zmin_A", "zmax_A")
    rows = []
    for line in (HERE / "log.nvt_followup_1ps_fixed_ids").read_text(errors="replace").splitlines():
        fields = line.split()
        if len(fields) != len(names) or not fields[0].isdigit():
            continue
        try:
            rows.append(dict(zip(names, map(float, fields))))
        except ValueError:
            continue
    if len(rows) < 8:
        raise RuntimeError(f"Expected thermo history, found {len(rows)} rows")
    return rows


def distance_metrics(pos: np.ndarray, typ: np.ndarray) -> dict[str, float]:
    p = pos.copy()
    p[:, :2] %= 200.0
    w, ag = p[typ == 1], p[typ == 2]
    wt, at = cKDTree(w, boxsize=BOX), cKDTree(ag, boxsize=BOX)
    return {
        "W_W_min_A": float(wt.query(w, k=2, workers=-1)[0][:, 1].min()),
        "Ag_Ag_min_A": float(at.query(ag, k=2, workers=-1)[0][:, 1].min()),
        "Ag_W_min_A": float(wt.query(ag, workers=-1)[0].min()),
    }


def motion(start: np.ndarray, end: np.ndarray, mobile_mask: np.ndarray) -> dict[str, object]:
    typ = start[:, 1].astype(int)
    d = end[:, 2:5] - start[:, 2:5]
    for ax in (0, 1):
        d[:, ax] -= 200.0 * np.rint(d[:, ax] / 200.0)
    disp = np.linalg.norm(d, axis=1)
    out: dict[str, object] = {
        "mobile_atoms": int(mobile_mask.sum()),
        "mobile_atoms_over_1A": int(np.sum(disp[mobile_mask] > 1.0)),
        "mobile_max_displacement_A": float(disp[mobile_mask].max()),
    }
    for tid, name in ELEMENT.items():
        mask = mobile_mask & (typ == tid)
        out[name] = {
            "mobile_atoms": int(mask.sum()),
            "displacement_p50_p95_p99_max_A": [float(v) for v in np.percentile(disp[mask], [50, 95, 99, 100])],
            "atoms_over_1A": int(np.sum(disp[mask] > 1.0)),
        }
    return out


def main() -> None:
    # The original fixed IDs are reconstructed from the first v6 NVE start.
    original = read_dump("nve_start.dump.gz")
    start = read_dump("nvt_fixedids_start.dump.gz")
    nvt_end = read_dump("nvt_fixedids_end.dump.gz")
    nve_end = read_dump("nve_fixedids_after_nvt_final.dump.gz")
    for a in (nvt_end, nve_end):
        if not np.array_equal(start[:, :2], a[:, :2]):
            raise RuntimeError("Atom IDs/types changed across continuation")
    fixed_ids = set(original[original[:, 4] <= 20.0, 0].astype(int))
    mobile_mask = ~np.isin(start[:, 0].astype(int), np.fromiter(fixed_ids, dtype=int))
    typ = start[:, 1].astype(int)
    rows = thermo_rows()
    nvt_rows = [r for r in rows if 0 <= r["step"] <= 4000]
    nve_start_row = [r for r in rows if r["step"] == 4000][-1]
    nve_end_row = rows[-1]
    nvt_start_row = next(r for r in nvt_rows if r["step"] == 0)
    nvt_end_row = [r for r in nvt_rows if r["step"] == 4000][-1]
    row_by_step = {int(r["step"]): r for r in nvt_rows}
    report = {
        "title": "Ag-W v6 corrected-fixed-ID 1 ps NVT plus 0.05 ps NVE",
        "status": "NVT_TARGET_REACHED_NEAR_1PS; NVE_CHECK_COMPLETED",
        "decision": "EXPLORATORY_COLD_BASELINE_CANDIDATE; REVIEW_MOTION_AND_TREND_BEFORE_TTM",
        "settings": {"boundary": "p p f", "fixed_ID_source": "original v6 NVE-start atoms with z<=20 A", "fixed_bottom_count": len(fixed_ids), "mobile_count": int(mobile_mask.sum()), "NVT_K": 300.0, "Tdamp_ps": 0.1, "NVT_duration_ps": 1.0, "NVE_duration_ps": 0.05, "dt_fs": 0.25},
        "restored_previously_frozen_atoms": 115,
        "atom_counts_start_nvt_end_nve_end": [int(len(start)), int(len(nvt_end)), int(len(nve_end))],
        "NVT_start": nvt_start_row,
        "NVT_end": nvt_end_row,
        "NVT_mobile_temperature_min_max_K": [float(min(r["T_mobile_K"] for r in nvt_rows)), float(max(r["T_mobile_K"] for r in nvt_rows))],
        "NVT_temperature_change_last_0p1ps_K": float(row_by_step[4000]["T_mobile_K"] - row_by_step[3600]["T_mobile_K"]),
        "NVT_temperature_change_last_0p2ps_K": float(row_by_step[4000]["T_mobile_K"] - row_by_step[3200]["T_mobile_K"]),
        "NVT_PE_change_eV": float(nvt_end_row["PE_eV"] - nvt_start_row["PE_eV"]),
        "NVT_KE_change_eV": float(nvt_end_row["KE_eV"] - nvt_start_row["KE_eV"]),
        "NVT_Etot_change_eV": float(nvt_end_row["Etot_eV"] - nvt_start_row["Etot_eV"]),
        "NVE_start": nve_start_row,
        "NVE_end": nve_end_row,
        "NVE_Etot_drift_eV": float(nve_end_row["Etot_eV"] - nve_start_row["Etot_eV"]),
        "NVE_mobile_temperature_change_K": float(nve_end_row["T_mobile_K"] - nve_start_row["T_mobile_K"]),
        "minimum_pair_distances_A": {
            "NVT_start": distance_metrics(start[:, 2:5], typ),
            "NVT_end": distance_metrics(nvt_end[:, 2:5], typ),
            "NVE_end": distance_metrics(nve_end[:, 2:5], typ),
        },
        "displacement_from_NVT_start": {"NVT_end": motion(start, nvt_end, mobile_mask), "NVE_end": motion(start, nve_end, mobile_mask)},
        "z_bounds_A": {"NVT_start": [float(start[:, 4].min()), float(start[:, 4].max())], "NVT_end": [float(nvt_end[:, 4].min()), float(nvt_end[:, 4].max())], "NVE_end": [float(nve_end[:, 4].min()), float(nve_end[:, 4].max())]},
        "invalid_zero_interaction_attempt": "Excluded: first corrected-ID input omitted pair style/coeff; zero pair energy and forces. See retained invalid input/log.",
        "limitations": [
            "The NVT temperature was still increasing during its final 0.1 ps; the NVE check is the better indicator of post-thermostat stability.",
            "This remains an exploratory result using the v5 volume-regularized QSGS candidate and Dai EFS; it is not a validated quantitative Ag-W model.",
            "No TTM or laser source was used in this continuation.",
        ],
    }
    (HERE / "nvt_fixedids_metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
