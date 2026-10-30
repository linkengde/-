#!/usr/bin/env python3
"""Create a grid-scale interface-resistance regime sweep from candidate k(T)."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RUN_DIR = HERE.parent.parent
PROJECT_DIR = RUN_DIR.parents[1]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
TABLE = PHASE_DIR / "thermophysical_candidates.csv"
OUTPUT = HERE / "interface_conductance_resistance_sweep.csv"
G_VALUES = (0.0, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, float("inf"))


def load_k(material: str, phase: str) -> tuple[np.ndarray, np.ndarray]:
    with TABLE.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row["material"] == material and row["phase"] == phase]
    return (
        np.asarray([float(row["temperature_K"]) for row in rows]),
        np.asarray([float(row["thermal_conductivity_W_m_K"]) for row in rows]),
    )


def main() -> None:
    with np.load(GEOMETRY, allow_pickle=False) as data:
        spacing_m = np.asarray(data["size_A"], dtype=float) / np.asarray(data["skeleton"].shape, dtype=float) * 1e-10
    ag_t, ag_solid_k = load_k("Ag", "solid")
    sc_t, sc_k = load_k("Ti3SiC2", "solid")
    tm = float(ag_t[-1])
    states = (
        ("300K_candidates", 300.0, float(np.interp(300.0, ag_t, ag_solid_k)), float(np.interp(300.0, sc_t, sc_k))),
        ("near_Ag_Tm", tm, float(np.interp(tm, ag_t, ag_solid_k)), float(np.interp(tm, sc_t, sc_k))),
        ("near_CaseB_scaffold_peak", 1431.5269152051005, 175.0, float(np.interp(1431.5269152051005, sc_t, sc_k))),
    )
    rows = []
    for state, temperature, k_ag, k_sc in states:
        for axis, dx in enumerate(spacing_m):
            r_bulk = float(dx / (2.0 * k_ag) + dx / (2.0 * k_sc))
            g_equal_bulk = 1.0 / r_bulk
            for g in G_VALUES:
                if g == 0.0:
                    r_interface = float("inf")
                    ratio = float("inf")
                    flux_ratio = 0.0
                elif np.isinf(g):
                    r_interface = 0.0
                    ratio = 0.0
                    flux_ratio = 1.0
                else:
                    r_interface = 1.0 / g
                    ratio = r_interface / r_bulk
                    flux_ratio = r_bulk / (r_bulk + r_interface)
                rows.append({
                    "state": state,
                    "temperature_K_for_k_anchor": temperature,
                    "face_normal_axis": "xyz"[axis],
                    "voxel_pitch_nm": dx * 1e9,
                    "candidate_k_Ag_W_m_K": k_ag,
                    "candidate_k_Ti3SiC2_W_m_K": k_sc,
                    "bulk_half_cell_areal_resistance_m2_K_W": r_bulk,
                    "G_equal_to_bulk_resistance_W_m2_K": g_equal_bulk,
                    "interface_conductance_W_m2_K": "0" if g == 0 else ("inf" if np.isinf(g) else f"{g:.8g}"),
                    "interface_areal_resistance_m2_K_W": "inf" if np.isinf(r_interface) else f"{r_interface:.8g}",
                    "R_interface_over_R_bulk": "inf" if np.isinf(ratio) else f"{ratio:.8g}",
                    "face_flux_over_perfect_contact": flux_ratio,
                    "interpretation": "numerical regime sweep only; no measured Ag/Ti3SiC2 interface G_K is selected",
                })
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {OUTPUT} ({len(rows)} rows)")
    print("grid-scale G_K where interface resistance equals two bulk half-cell resistances:")
    for state, _, k_ag, k_sc in states:
        values = [float(dx / (2.0 * k_ag) + dx / (2.0 * k_sc)) for dx in spacing_m]
        print(state, [1.0 / r for r in values])


if __name__ == "__main__":
    main()
