#!/usr/bin/env python3
"""Build P1P Ag properties from the user-supplied structured source pack.

No network access is used. Saturation pressure and the derivative-based latent
heat are generated from equations transcribed in the supplied README; the five
alpha_ev knots are read from its Table 8 CSV. The resulting alpha values away
from those knots remain blank: the C kernel interpolates between finite knots
and refuses alpha extrapolation unless an explicit sensitivity fallback is
passed by the caller.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source_data" / "user_source_pack"
OUT = ROOT / "data"
R = 8.31446261815324  # J mol-1 K-1
M_AG = 0.1078682  # kg mol-1
N_A = 6.02214076e23  # mol-1


def read_alpha_knots(path: Path) -> dict[float, float]:
    knots: dict[float, float] = {}
    with path.open(newline="", encoding="utf-8-sig") as fp:
        for row in csv.DictReader(fp):
            t = float(row["T_K"])
            alpha = float(row["alpha_ev_dimensionless"])
            if not 0.0 <= alpha <= 1.0:
                raise ValueError(f"alpha out of range at {t} K")
            knots[t] = alpha
    if len(knots) != 5 or min(knots) != 1373 or max(knots) != 1773:
        raise ValueError("expected the five supplied alpha_ev anchors, 1373-1773 K")
    return knots


def alpha_for_temperature(t: float, knots: dict[float, float]) -> float | None:
    points = sorted(knots.items())
    if t < points[0][0] or t > points[-1][0]:
        return None
    for ti, ai in points:
        if t == ti:
            return ai
    for (t0, a0), (t1, a1) in zip(points, points[1:]):
        if t0 < t < t1:
            return a0 + (a1 - a0) * (t - t0) / (t1 - t0)
    raise AssertionError(t)


def wang_psat(t: float) -> float:
    return 10.0 ** (-14400.0 / t - 0.85 * math.log10(t) + 13.825)


def nist_psat(t: float) -> float | None:
    if not 1823.0 <= t <= 2425.0:
        return None
    log10_p_bar = 1.95303 - 2505.533 / (t - 1194.947)
    return 1.0e5 * 10.0**log10_p_bar


def lv_mass(t: float) -> float:
    # Clausius-Clapeyron derivative of the Wang fit; ideal-vapor estimate.
    return R * (14400.0 * math.log(10.0) - 0.85 * t) / M_AG


def main() -> None:
    alpha_path = SOURCE / "Ag_evaporation_P1P_P1Q_ag_evaporation_coefficient_table8.csv"
    knots = read_alpha_knots(alpha_path)
    OUT.mkdir(parents=True, exist_ok=True)
    temperatures = set(float(t) for t in range(1234, 2435, 25))
    temperatures.update((1234.0, 1235.0, 1373.0, 1473.0, 1573.0, 1673.0,
                         1773.0, 1823.0, 2000.0, 2200.0, 2300.0, 2425.0, 2436.0))
    temps = sorted(temperatures)

    main_csv = OUT / "Ag_evaporation_properties.csv"
    solver_csv = OUT / "Ag_evaporation_solver_table.csv"
    with main_csv.open("w", newline="", encoding="utf-8") as fp, \
         solver_csv.open("w", newline="", encoding="utf-8") as sfp:
        writer = csv.writer(fp)
        solver = csv.writer(sfp)
        writer.writerow(("temperature_K", "p_sat_Wang2024_Pa", "p_sat_NIST_Antoine_Pa",
                         "Lv_derived_J_kg", "alpha_ev", "alpha_ev_status",
                         "p_sat_status", "Lv_status"))
        solver.writerow(("temperature_K", "p_sat_Pa", "L_v_J_kg", "alpha_ev"))
        for t in temps:
            alpha = alpha_for_temperature(t, knots)
            alpha_status = "" if alpha is None else (
                "Table8 anchor" if t in knots else
                "linear interpolation between Table8 experimental anchors")
            p_wang = wang_psat(t)
            p_nist = nist_psat(t)
            lv = lv_mass(t)
            p_status = "Wang2024 liquid-Ag fit (1234-2436 K)"
            lv_status = "derived by Clausius-Clapeyron; ideal-vapor approximation; not calorimetry"
            writer.writerow((f"{t:.0f}", f"{p_wang:.12g}",
                             "" if p_nist is None else f"{p_nist:.12g}",
                             f"{lv:.12g}", "" if alpha is None else f"{alpha:.9g}",
                             alpha_status, p_status, lv_status))
            solver.writerow((f"{t:.0f}", f"{p_wang:.12g}", f"{lv:.12g}",
                             "" if t not in knots else f"{knots[t]:.9g}"))

    audit_path = OUT / "Ag_evaporation_property_audit.csv"
    with audit_path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp)
        writer.writerow(("temperature_K", "Wang_p_sat_Pa", "NIST_p_sat_Pa",
                         "NIST_minus_Wang_percent", "Lv_derived_kJ_mol", "Lv_derived_MJ_kg",
                         "alpha_ev", "alpha_status"))
        for t in (1373, 1473, 1573, 1673, 1773, 1823, 2000, 2200, 2300, 2425):
            alpha = alpha_for_temperature(float(t), knots)
            pn = nist_psat(float(t))
            writer.writerow((t, f"{wang_psat(t):.9g}", "" if pn is None else f"{pn:.9g}",
                             "" if pn is None else f"{100.*(pn/wang_psat(t)-1.):.9g}",
                             f"{lv_mass(t)*M_AG/1000.:.9g}", f"{lv_mass(t)/1.e6:.9g}",
                             "" if alpha is None else f"{alpha:.9g}",
                             "outside Table8 alpha range" if alpha is None else
                             ("Table8 anchor" if float(t) in knots else "interpolated")))

    # Engineering diagnostics: the high-temperature alpha values are brackets
    # borrowed from the measured interval, not extrapolated measurements.
    # p_Ag,v=0 is a free-evaporation upper-bound scenario only.
    flux_path = OUT / "Ag_evaporation_flux_scale_diagnostic.csv"
    r_ag = R / M_AG
    with flux_path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp)
        writer.writerow(("temperature_K", "alpha_scenario", "p_Ag_v_Pa",
                         "HK_flux_kg_m2_s", "rho_assumed_kg_m3",
                         "equivalent_recession_100ps_nm", "expected_Ag_atoms_100nm2_100ps",
                         "expected_Ag_atoms_100nm2_5ps", "scenario_status"))
        area = 100.e-18
        for alpha in (0.0834, 0.1428, 0.2715):
            t = 2300.0
            flux = alpha * wang_psat(t) / math.sqrt(2.0 * math.pi * r_ag * t)
            recession_nm = flux * 100.0e-12 / 9000.0 * 1.0e9
            mass_100ps = flux*area*100.e-12
            atoms_100ps = mass_100ps/M_AG*N_A
            writer.writerow((t, alpha, 0.0, f"{flux:.12g}", 9000.0,
                             f"{recession_nm:.12g}", f"{atoms_100ps:.12g}",
                             f"{atoms_100ps*0.05:.12g}",
                             "free-evaporation upper bound; alpha bracket borrowed from Table8 1373-1773 K, not validated at 2300 K"))
        writer.writerow((2300, "", f"{wang_psat(2300):.12g}", 0.0, 9000, 0.0,
                         0.0, 0.0,
                         "zero net flux when p_Ag,v equals p_sat; illustrative boundary condition"))

    # Convert the user-supplied Table 5 rates to SI and independently check
    # the Table 8 theoretical maximum rates against the supplied p_sat fit.
    measured_in = SOURCE.parent / "user_implementation_pack" / \
        "P1P_Ag_evaporation_implementation_pack_measured_evaporation_rates_table5.csv"
    if measured_in.exists():
        measured_out = OUT / "Ag_measured_evaporation_rates_SI.csv"
        with measured_in.open(newline="", encoding="utf-8-sig") as src, \
             measured_out.open("w", newline="", encoding="utf-8") as dst:
            reader = csv.DictReader(src)
            writer = csv.writer(dst)
            writer.writerow(("temperature_K", "system_pressure_Pa", "rate_kg_m2_s",
                             "interpretation"))
            for row in reader:
                writer.writerow((row["T_K"], row["system_pressure_Pa"],
                                 row["rate_kg_m2_s"],
                                 "system pressure is not local Ag vapor partial pressure"))

    alpha_in = alpha_path
    theory_path = OUT / "Ag_table8_theoretical_rate_check.csv"
    with alpha_in.open(newline="", encoding="utf-8-sig") as src, \
         theory_path.open("w", newline="", encoding="utf-8") as dst:
        writer = csv.writer(dst)
        writer.writerow(("temperature_K", "published_alpha_ev", "omega_theory_recomputed_g_cm2_min",
                         "omega_theory_Table8_g_cm2_min", "relative_difference_percent",
                         "status_note"))
        for row in csv.DictReader(src):
            t = float(row["T_K"])
            alpha = float(row["alpha_ev_dimensionless"])
            theoretical_flux = wang_psat(t)/math.sqrt(2.*math.pi*r_ag*t)
            # 1 kg m^-2 s^-1 = 6 g cm^-2 min^-1.
            omega_recalc = theoretical_flux*6.
            omega_table = float(row["omega_max_theoretical_g_cm-2_min-1"])
            difference = 100.*(omega_recalc/omega_table-1.)
            writer.writerow((f"{t:.0f}",f"{alpha:.7g}",f"{omega_recalc:.9g}",
                             f"{omega_table:.9g}",f"{difference:.7g}",
                             "formula reconstruction is ~1.2% above rounded Table8; retain published alpha"))

    # Cross-check the supplied tabulations/equations so source transcription
    # errors are caught before the CSV is used by the C solver.
    reference = {1373: 4.6773247, 1473: 22.701798, 1573: 89.804625,
                 1673: 300.42506, 1773: 874.52084, 2000: 6593.6289,
                 2200: 27446.45, 2300: 50893.436}
    for t, expected in reference.items():
        if abs(wang_psat(t) / expected - 1.0) > 2.0e-7:
            raise ValueError(f"supplied Wang pressure check failed at {t} K")
    if abs(lv_mass(2300) / 2.405062e6 - 1.0) > 1.0e-6:
        raise ValueError("derived 2300 K latent-heat check failed")
    print(f"generated {len(temps)} rows; source alpha anchors={len(knots)}; "
          f"alpha validity={min(knots):.0f}-{max(knots):.0f} K; "
          f"property fit={temps[0]:.0f}-{temps[-1]:.0f} K")


if __name__ == "__main__":
    main()
