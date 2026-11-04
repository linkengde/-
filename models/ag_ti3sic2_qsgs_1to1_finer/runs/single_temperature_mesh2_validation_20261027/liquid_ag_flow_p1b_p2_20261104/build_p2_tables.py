#!/usr/bin/env python3
"""Rebuild auditable liquid-Ag property tables from cited fits and points."""

from __future__ import annotations

import csv
import math
from pathlib import Path


OUT = Path(__file__).resolve().parent

VISCOSITY_POINTS = [
    (1273, 3.520, 0.372),
    (1323, 2.904, 0.225),
    (1373, 2.812, 0.323),
    (1423, 2.649, 0.165),
    (1473, 2.540, 0.110),
    (1523, 2.498, 0.229),
    (1573, 2.399, 0.176),
    (1623, 2.332, 0.195),
    (1673, 2.214, 0.635),
]


def mu_mpas(temperature_k: float) -> float:
    return math.exp(-0.4631 + 2089.0 / temperature_k)


def write_csv(filename: str, header: list[str], rows: list[list[object]]) -> None:
    with (OUT / filename).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def build_viscosity_tables() -> dict[str, float]:
    raw_rows = []
    residuals = []
    log_residuals = []
    for temperature, measured, reported_aux in VISCOSITY_POINTS:
        fitted = mu_mpas(temperature)
        residual_pct = 100.0 * (fitted - measured) / measured
        raw_rows.append([
            "Dubberstein et al. (2016), Table 4",
            "10.1007/s10765-016-2104-7",
            temperature,
            f"{measured:.6f}",
            f"{reported_aux:.6f}",
            f"{fitted:.9f}",
            f"{residual_pct:.6f}",
            "primary-paper Table 4; user reports NIST JSON's nine Ag values match",
            "1273-1673",
            "externally checked in chat as reported by user; raw JSON not yet archived",
            "JSON reportedly has per-point expanded-uncertainty fields; values, units, and coverage not yet locally extracted",
            "mPa s; Table 4 auxiliary eta_sigma column retained verbatim; definition not stated in inspected text",
        ])
        residuals.append(residual_pct)
        log_residuals.append(math.log(fitted / measured))

    write_csv(
        "ag_viscosity_experimental.csv",
        [
            "source", "doi", "temperature_K", "measured_dynamic_viscosity_mPa_s",
            "reported_eta_sigma_mPa_s_undefined", "fit_prediction_mPa_s",
            "fit_minus_measured_percent", "data_provenance", "reported_temperature_range_K",
            "nist_json_nine_value_match_status", "nist_expanded_uncertainty_status", "notes",
        ],
        raw_rows,
    )

    grid = list(range(1273, 1674, 25))
    grid = sorted(set(grid + [temperature for temperature, _, _ in VISCOSITY_POINTS]))
    write_csv(
        "ag_viscosity_fit.csv",
        ["temperature_K", "mu_fit_mPa_s", "mu_fit_Pa_s", "fit_status", "equation"],
        [[
            temperature,
            f"{mu_mpas(temperature):.9f}",
            f"{mu_mpas(temperature) * 1e-3:.12g}",
            "fitted correlation within reported 1273-1673 K experiment range",
            "ln(mu/(mPa s)) = -0.4631 + 2089/T",
        ] for temperature in grid],
    )

    return {
        "mape_pct": sum(abs(value) for value in residuals) / len(residuals),
        "max_abs_error_pct": max(abs(value) for value in residuals),
        "mean_signed_error_pct": sum(residuals) / len(residuals),
        "log_rmse": math.sqrt(sum(value * value for value in log_residuals) / len(log_residuals)),
    }


def build_density_table() -> list[dict[str, float]]:
    rows: list[list[object]] = []
    comparisons: list[dict[str, float]] = []
    abdullaev_points = [1234.93, 1300.0, 1500.0, 1573.0, 1600.0, 1631.0]
    abdullaev_reported_errors = {1234.93: 0.25, 1300.0: 0.26, 1600.0: 0.34}
    for temperature in abdullaev_points:
        density = 9327.1 - 0.9956 * (temperature - 1234.93)
        error = abdullaev_reported_errors.get(temperature)
        rows.append([
            "Abdullaev, Khairulin & Stankus (2020), Eq. 2",
            "10.1088/1742-6596/1677/1/012161",
            f"{temperature:.2f}",
            f"{density:.6f}",
            "9327.1 - 0.9956*(T - 1234.93)",
            "1234.93-1631",
            "primary experimental-source correlation",
            "" if error is None else f"{error:.2f}",
            "kg m^-3; reported error coverage convention unspecified; blank means no pointwise error quoted here",
        ])

    alternate_points = [1235.0, 1300.0, 1500.0, 1573.0]
    for temperature in alternate_points:
        density = 9577.0 - 1.41 * (temperature - 1235.0)
        abdullaev = 9327.1 - 0.9956 * (temperature - 1234.93)
        difference_pct = 100.0 * (density - abdullaev) / abdullaev
        comparisons.append({"temperature_K": temperature, "difference_pct": difference_pct})
        rows.append([
            "Dubberstein & Heller (2014), as cited in Dubberstein et al. (2016) Table 3",
            "no DOI independently verified in this audit",
            f"{temperature:.2f}",
            f"{density:.6f}",
            "9577 - 1.41*(T - 1235)",
            "1235-1573",
            "secondary citation in inspected primary paper; candidate cross-check, not selected anchor",
            "",
            f"kg m^-3; differs from Abdullaev fit by {difference_pct:.4f}% at this T; source paper reference 21",
        ])

    write_csv(
        "ag_density_reference.csv",
        [
            "source", "doi_or_reference_status", "temperature_K", "rho_fit_kg_m3",
            "equation", "reported_temperature_range_K", "evidence_status",
            "reported_error_percent", "notes",
        ],
        rows,
    )
    return comparisons


def build_surface_tension_table() -> None:
    rows: list[list[object]] = []
    red_temps = list(range(1219, 1735, 25)) + [1734]
    red_temps = sorted(set(red_temps))
    for temperature in red_temps:
        sigma = 984.0 - 0.298 * (temperature - 1234.0)
        rows.append([
            "Morohoshi et al. (2008)", "10.2320/jinstmet.72.708",
            "Ar/He + 5% H2; H2O=2.6 ppm; apparent pO2 varies with temperature",
            temperature, f"{sigma:.6f}", "-0.298", "1219-1734",
            "measured oscillating-drop linear fit; use only within range",
        ])

    oxygen_temps = list(range(1350, 1522, 25)) + [1521]
    oxygen_temps = sorted(set(oxygen_temps))
    for temperature in oxygen_temps:
        sigma = 698.0 + 0.546 * (temperature - 1234.0)
        rows.append([
            "Morohoshi et al. (2008)", "10.2320/jinstmet.72.708",
            "pO2=100 Pa",
            temperature, f"{sigma:.6f}", "+0.546", "1350-1521",
            "measured oscillating-drop linear fit; use only within range",
        ])

    rows.extend([
        ["Hibiya & Ozawa (2013)", "10.1002/crat.201200514", "pO2=1e-10 Pa",
         "", "", "-0.245", "temperature range not transcribed in current audit",
         "reported slope-only contextual cross-check; no curve generated"],
        ["Hibiya & Ozawa (2013)", "10.1002/crat.201200514", "pO2=1 Pa",
         "", "", "+0.609", "temperature range not transcribed in current audit",
         "reported slope-only contextual cross-check; no curve generated"],
    ])

    write_csv(
        "ag_surface_tension_candidates.csv",
        [
            "source", "doi", "atmosphere_condition", "temperature_K",
            "surface_tension_mN_m", "d_sigma_dT_mN_m_K", "validity_range_K", "evidence_status",
        ],
        rows,
    )


def main() -> None:
    viscosity = build_viscosity_tables()
    density_comparisons = build_density_table()
    build_surface_tension_table()
    assert len(VISCOSITY_POINTS) == 9
    assert abs(mu_mpas(1273) - 3.247565) < 1e-5
    assert all(math.isfinite(mu_mpas(t)) and mu_mpas(t) > 0 for t in range(1273, 1674, 10))
    print(f"viscosity_summary={viscosity}")
    print("density_candidate_differences_pct=" + repr(density_comparisons))
    print("wrote auditable CSVs to " + str(OUT))


if __name__ == "__main__":
    main()
