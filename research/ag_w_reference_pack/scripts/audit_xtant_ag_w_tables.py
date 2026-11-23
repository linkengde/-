#!/usr/bin/env python3
"""Inventory and compare the user-supplied XTANT-3 Ag/W TTM source tables.

This is a provenance/schema audit, not a fit or an MD/TTM run. It preserves
the distinct Ce datasets included with the dynamical tables and K tables.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


PACK = Path(__file__).resolve().parents[1]
DATA = PACK / "ttm_data" / "XTANT3_TTM_parameters_v3"
LIN_W = PACK / "ttm_data" / "W_TTM_Lin2008_digitized_300_20000K.csv"
OUTPUT = PACK / "ttm_data" / "XTANT3_AG_W_SOURCE_AUDIT.json"

FILES = {
    "Ce_Ag": ("Ce_metals/Ce_Ag.txt", 2, "J/(m^3 K)"),
    "Ce_W": ("Ce_metals/Ce_W.txt", 2, "J/(m^3 K)"),
    "G_Ag": ("G_e-ph_metals/G_e-ph_Ag.txt", 2, "W/(m^3 K)"),
    "G_W": ("G_e-ph_metals/G_e-ph_W.txt", 2, "W/(m^3 K)"),
    "K_Ag": ("K_metals/K_Ag.dat", 6, "W/(m K), eV, J/(m^3 K)"),
    "K_W": ("K_metals/K_W.dat", 6, "W/(m K), eV, J/(m^3 K)"),
}
KEY_TEMPERATURES = [300, 500, 600, 1000, 3000, 5000, 11000, 15000, 19000, 20000, 23000, 50000]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([float(value) for value in line.split()])
    return rows


def at_temperature(rows: list[list[float]], temperature: int) -> list[float] | None:
    for row in rows:
        if int(round(row[0])) == temperature:
            return row[1:]
    return None


def summarize_table(label: str, filename: str, expected_columns: int, units: str) -> dict:
    path = DATA / filename
    rows = read_table(path)
    temperatures = [row[0] for row in rows]
    column_counts = sorted({len(row) for row in rows})
    positive_numeric_values = all(value > 0 for row in rows for value in row[1:])
    finite_numeric_values = all(value == value and abs(value) != float("inf") for row in rows for value in row)
    steps = sorted({round(b - a, 10) for a, b in zip(temperatures, temperatures[1:])})
    return {
        "filename": filename,
        "sha256": sha256(path),
        "row_count": len(rows),
        "columns_per_row": column_counts,
        "expected_columns": expected_columns,
        "column_count_matches_expected": column_counts == [expected_columns],
        "temperature_min_K": temperatures[0],
        "temperature_max_K": temperatures[-1],
        "temperature_steps_K": steps,
        "strictly_increasing_temperature": all(b > a for a, b in zip(temperatures, temperatures[1:])),
        "all_values_finite": finite_numeric_values,
        "all_parameter_values_positive": positive_numeric_values,
        "units": units,
        "key_values": {
            str(t): at_temperature(rows, t)
            for t in KEY_TEMPERATURES
            if at_temperature(rows, t) is not None
        },
        "_rows": rows,
    }


def read_lin_w() -> dict[int, dict[str, float]]:
    result: dict[int, dict[str, float]] = {}
    with LIN_W.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            result[int(float(row["Te_K"]))] = {
                key: float(value) for key, value in row.items() if key != "Te_K"
            }
    return result


def main() -> None:
    tables = {
        label: summarize_table(label, filename, columns, units)
        for label, (filename, columns, units) in FILES.items()
    }
    lin_w = read_lin_w()
    comparisons = {"W_Ce_XTANT_dynamic_vs_Lin_digitized": {}, "W_G_XTANT_vs_Lin_digitized": {},
                   "W_Ce_XTANT_dynamic_vs_K_internal_Ce": {}, "Ag_Ce_XTANT_dynamic_vs_K_internal_Ce": {}}

    for temperature in KEY_TEMPERATURES:
        xtant_ce = at_temperature(tables["Ce_W"]["_rows"], temperature)
        xtant_g = at_temperature(tables["G_W"]["_rows"], temperature)
        k_w = at_temperature(tables["K_W"]["_rows"], temperature)
        k_ag = at_temperature(tables["K_Ag"]["_rows"], temperature)
        ag_ce = at_temperature(tables["Ce_Ag"]["_rows"], temperature)

        if temperature in lin_w and xtant_ce:
            lin = lin_w[temperature]["Ce_Lin2008_J_m-3_K-1"]
            comparisons["W_Ce_XTANT_dynamic_vs_Lin_digitized"][str(temperature)] = {
                "XTANT_Ce_J_m-3_K-1": xtant_ce[0], "Lin_digitized_Ce_J_m-3_K-1": lin,
                "XTANT_over_Lin": xtant_ce[0] / lin,
            }
        if temperature in lin_w and xtant_g:
            lin = lin_w[temperature]["G_Lin2008_digitized_W_m-3_K-1"]
            comparisons["W_G_XTANT_vs_Lin_digitized"][str(temperature)] = {
                "XTANT_G_W_m-3_K-1": xtant_g[0], "Lin_digitized_G_W_m-3_K-1": lin,
                "XTANT_over_Lin": xtant_g[0] / lin,
                "XTANT_over_Duffy_experimental_G0": xtant_g[0] / 1.4e17,
            }
        if xtant_ce and k_w:
            comparisons["W_Ce_XTANT_dynamic_vs_K_internal_Ce"][str(temperature)] = {
                "Ce_dynamic_J_m-3_K-1": xtant_ce[0], "Ce_K_table_J_m-3_K-1": k_w[4],
                "dynamic_over_K_table": xtant_ce[0] / k_w[4],
            }
        if ag_ce and k_ag:
            comparisons["Ag_Ce_XTANT_dynamic_vs_K_internal_Ce"][str(temperature)] = {
                "Ce_dynamic_J_m-3_K-1": ag_ce[0], "Ce_K_table_J_m-3_K-1": k_ag[4],
                "dynamic_over_K_table": ag_ce[0] / k_ag[4],
            }

    # The reported total is consistent with the reciprocal combination of the
    # electron-phonon and electron-electron scattering-limited conductivities.
    k_combination = {}
    for species in ("Ag", "W"):
        for temperature in (300, 1000, 5000, 11000, 19000):
            row = at_temperature(tables[f"K_{species}"]["_rows"], temperature)
            if row:
                k_total, k_ep, k_ee = row[:3]
                harmonic = 1.0 / (1.0 / k_ep + 1.0 / k_ee)
                k_combination[f"{species}_{temperature}K"] = {
                    "reported_k_total_W_m-1_K-1": k_total,
                    "harmonic_combination_k_ep_k_ee_W_m-1_K-1": harmonic,
                    "relative_difference": abs(harmonic - k_total) / k_total,
                }

    for table in tables.values():
        del table["_rows"]

    result = {
        "title": "Schema, provenance, and cross-source audit of user-supplied XTANT-3 Ag/W tables",
        "source_metadata": "ttm_data/XTANT3_TTM_parameters_v3/README.md",
        "user_supplied_source_readme": "ttm_data/XTANT3_TTM_parameters_v3/SOURCE_README_user_upload.md",
        "license_file": "ttm_data/XTANT3_TTM_parameters_v3/LICENSE_user_upload",
        "lin_comparison_source": "ttm_data/W_TTM_Lin2008_digitized_300_20000K.csv",
        "tables": tables,
        "cross_source_comparisons": comparisons,
        "K_total_check": k_combination,
        "interpretation": [
            "Ce_Ag/W are dynamical-simulation values; Ce columns in K_Ag/W are from ideal crystals with a 7x7x7 k-point mesh and are a distinct dataset.",
            "K tables extend to 50000 K, but standalone Ce/G tables do not; do not extrapolate Ce/G merely because K has higher coverage.",
            "The supplied G files contain only Te and G columns, so they provide no explicit lattice-temperature dependence.",
            "The W XTANT G curve is substantially below the Lin 2008 digitized curve and the Duffy low-temperature anchor; preserve as an alternate sensitivity model.",
            "No Ag/W-specific electronic interfacial conductance or validated mixed-phase closure is supplied by these files.",
        ],
        "scope": "Provenance/schema audit only; no parameter blending, fitting, or TTM trajectory is performed.",
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(PACK)}")


if __name__ == "__main__":
    main()
