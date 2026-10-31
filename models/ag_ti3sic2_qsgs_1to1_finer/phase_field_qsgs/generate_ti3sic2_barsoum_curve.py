#!/usr/bin/env python3
"""Generate the Ti3SiC2 correlation-derived property grid.

These are formula-evaluated values, not digitized raw experimental points.
Run from this directory:
    python generate_ti3sic2_barsoum_curve.py
"""
from __future__ import annotations

import csv
from pathlib import Path

MOLAR_MASS_KG_MOL = 0.195708
DENSITY_KG_M3 = 4520.0
OUTPUT = Path(__file__).with_name("ti3sic2_barsoum_correlation_curve_300_2000K.csv")


def cp_status(t: float) -> str:
    if abs(t - 300.0) < 1e-9:
        return "room_temperature_correlation_value_not_digitized_raw_point"
    if t < 373.15:
        return "correlation_evaluation_below_Barsoum_DSC_range"
    if t <= 1273.15:
        return "correlation_evaluation_within_Barsoum_DSC_range_100_to_1000C_not_raw_digitized_data"
    if t <= 1473.15:
        return "Barsoum_reported_extrapolation_to_about_1200C"
    if t <= 1873.15:
        return "Qiu_secondary_review_stated_range_only_not_independent_measurement"
    return "extrapolation_beyond_cited_Cp_ranges"


def main() -> None:
    temperatures = {float(t) for t in range(300, 2001, 25)}
    temperatures.update({373.15, 1273.15, 1473.15, 1873.15, 2000.0})
    fields = [
        "temperature_K", "temperature_C", "Cp_molar_J_mol_K",
        "Cp_mass_J_kg_K", "thermal_conductivity_W_m_K", "density_kg_m3",
        "Cp_data_status", "k_data_status", "density_data_status",
    ]
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for t in sorted(temperatures):
            cp_molar = 164.4 - 16419.0 / t
            k = 38.6 - 0.0045 * t
            writer.writerow({
                "temperature_K": f"{t:.2f}",
                "temperature_C": f"{t - 273.15:.2f}",
                "Cp_molar_J_mol_K": f"{cp_molar:.8f}",
                "Cp_mass_J_kg_K": f"{cp_molar / MOLAR_MASS_KG_MOL:.8f}",
                "thermal_conductivity_W_m_K": f"{k:.8f}",
                "density_kg_m3": f"{DENSITY_KG_M3:.2f}",
                "Cp_data_status": cp_status(t),
                "k_data_status": (
                    "Barsoum_reported_fit_range_25_to_1200C"
                    if t <= 1473.15
                    else "extrapolation_beyond_Barsoum_reported_k_fit_range"
                ),
                "density_data_status": "fixed_reference_density_assumption_not_temperature_dependent",
            })


if __name__ == "__main__":
    main()
