#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""按随附原始论文数据生成可审计的 Ag/Ti₃SiC₂ 材料CSV表。"""
from __future__ import annotations

import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
MOLAR_MASS_AG = 0.1078682
MOLAR_MASS_TI = 0.195708
TM_AG = 1235.08
RHO_AG_REF = 10499.4
RHO_TI_REF = 4500.0


def cp_ag_solid_molar(T_K: float) -> float:
    """Pankratz p.32: cal/mol/K to J/mol/K; T is kelvin."""
    return 4.184 * (5.272 + 1.836e-3 * T_K + 22200.0 / (T_K * T_K))


def h_ag_solid_mass(T_K: float) -> float:
    """Ag solid sensible enthalpy from 300 K, J/kg."""
    tref = 300.0
    integral_cal_mol = (5.272 * (T_K - tref)
                        + 0.000918 * (T_K * T_K - tref * tref)
                        + 22200.0 * (1.0 / tref - 1.0 / T_K))
    return integral_cal_mol * 4.184 / MOLAR_MASS_AG


def rho_ag_solid(T_K: float) -> float:
    d = T_K - 293.15
    return 10503.8 - 0.6453 * d - 6.94e-5 * d * d


def write_csv(name: str, columns: list[str], rows: list[list]):
    path = HERE / name
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


def main():
    solid_k = [
        (300.0, 4.29), (400.0, 4.25), (500.0, 4.19), (600.0, 4.12),
        (700.0, 4.04), (800.0, 3.96), (900.0, 3.88), (1000.0, 3.79),
        (1100.0, 3.70), (1200.0, 3.61), (1235.08, 3.58),
    ]
    solid_rows = []
    for t, k_raw in solid_k:
        cp_mol = cp_ag_solid_molar(t)
        solid_rows.append([
            f"{t:.2f}", f"{cp_mol / 4.184:.8f}", f"{cp_mol:.8f}",
            f"{cp_mol / MOLAR_MASS_AG:.8f}",
            f"{h_ag_solid_mass(t) * MOLAR_MASS_AG / 4.184:.8f}",
            f"{h_ag_solid_mass(t):.8f}", f"{k_raw:.5f}", f"{k_raw * 100.0:.5f}",
            f"{rho_ag_solid(t):.5f}", "recommended_table_value",
            "Ho et al. 1972; Pankratz 1982; Abdullaev et al. 2020",
        ])
    write_csv("Ag_solid_properties.csv", [
        "temperature_K", "Cp_molar_cal_mol_K", "Cp_molar_J_mol_K", "Cp_mass_J_kg_K",
        "sensible_enthalpy_from_300K_cal_mol", "sensible_enthalpy_from_300K_J_kg", "Ho_table_k_W_cm_K",
        "thermal_conductivity_W_m_K", "density_kg_m3", "data_status", "source_ids",
    ], solid_rows)

    liquid_k = [(1235.08, 1.75), (1273.2, 1.77), (1300.0, 1.78),
                (1400.0, 1.83), (1473.2, 1.86), (1500.0, 1.87),
                (1600.0, 1.90), (1700.0, 1.93), (1800.0, 1.95),
                (1900.0, 1.96), (2000.0, 1.97)]
    liquid_rows = []
    for t, k_raw in liquid_k:
        rho = 9327.1 - 0.9956 * (t - 1234.93)
        liquid_rows.append([
            f"{t:.2f}", "8.00000000", "33.47200000", f"{33.472 / MOLAR_MASS_AG:.8f}",
            f"{k_raw:.5f}", f"{k_raw * 100.0:.5f}", f"{rho:.5f}",
            "provisional_estimate_uncertainty_about_20_percent",
            "Pankratz 1982; Ho et al. 1972; Abdullaev et al. 2020",
        ])
    write_csv("Ag_liquid_properties.csv", [
        "temperature_K", "Cp_molar_cal_mol_K", "Cp_molar_J_mol_K", "Cp_mass_J_kg_K",
        "Ho_table_k_W_cm_K", "thermal_conductivity_W_m_K", "density_kg_m3",
        "data_status", "source_ids",
    ], liquid_rows)

    thermodynamics = [
        ["A_Pankratz_main", "solid", "298.15-1235.08", "5.272+0.001836*T+22200/T^2",
         "cal/mol/K", "4.184*(5.272+0.001836*T+22200/T^2)", "J/mol/K",
         "2.700", "kcal/mol", "11296.8", "J/mol", "solid formula integrated from 300 K; latent once at Tm",
         "compiled_reference", "Pankratz 1982 Bulletin 672, printed p.32, PDF p.40"],
        ["A_Pankratz_main", "liquid", "1235.08-2000", "8.000", "cal/mol/K", "8.000*4.184", "J/mol/K",
         "2.700", "kcal/mol", "11296.8", "J/mol",
         "table uses constant Cp; values above 1600 K extrapolated in source",
         "compiled_reference_with_extrapolation", "Pankratz 1982 Bulletin 672, printed p.32, PDF p.40"],
        ["B_Feber_Hultgren_solid", "liquid", "1243-1588", "8.111", "cal/mol/K", "8.111*4.184/0.1078682", "J/kg/K", "2590", "cal/mol", "2590*4.184/0.1078682", "J/kg", "liquid enthalpy fit; ΔH uses Hultgren solid at Tm",
         "drop_calorimetry", "Feber et al. 1969, printed pp.171-173, PDF pp.3 and 5, Eq.(1), Tables 1-2"],
        ["B_Feber_author_solid_extrapolation", "liquid", "1243-1588", "8.111", "cal/mol/K", "8.111*4.184/0.1078682", "J/kg/K", "2710", "cal/mol", "2710*4.184/0.1078682", "J/kg", "separate ΔH from authors' solid fit; do not average",
         "drop_calorimetry_mixed_solid_fit", "Feber et al. 1969, printed pp.171-173, PDF pp.3 and 5, Eq.(1), Table 2"],
        ["C_Cagran_2006", "solid_near_melting", "1000-1234", "165.6", "J/kg/K", "165.6", "J/kg/K", "101.59", "kJ/kg", "101.59", "kJ/kg",
         "Eq.(5) solid fit H=17.2367+0.1656T kJ/kg; Table 2 Hs differs",
         "pulse_heating_fit", "Cagran et al. 2006, printed pp.107-109, supplied PDF pp.4-6, Eqs.(5)-(6), Table 2"],
        ["C_Cagran_2006", "liquid", "1234-2000", "259.8", "J/kg/K", "259.8", "J/kg/K", "101.59", "kJ/kg", "101.59", "kJ/kg",
         "Eq.(6) H=2.5807+0.2598T kJ/kg; Table 2 values Hs=231.49, Hl=323.17",
         "pulse_heating_fit_with_internal_discrepancy", "Cagran et al. 2006, printed pp.107-109, supplied PDF pp.4-6, Eqs.(5)-(6), Table 2"],
    ]
    write_csv("Ag_thermodynamic_alternatives.csv", [
        "material_set", "phase", "source_temperature_range_K", "source_Cp_or_model_expression",
        "source_Cp_units", "SI_Cp_or_model_expression", "SI_Cp_units",
        "source_latent_heat_value_or_expression", "source_latent_heat_units",
        "SI_latent_heat_value_or_expression", "SI_latent_heat_units",
        "interpretation", "evidence_class", "source_location",
    ], thermodynamics)

    matula = [
        (300, 1.629, "solid_total_recommended"), (350, 1.932, "solid_total_recommended"),
        (400, 2.241, "solid_total_recommended"), (500, 2.875, "solid_total_recommended"),
        (600, 3.531, "solid_total_recommended"), (700, 4.209, "solid_total_recommended"),
        (800, 4.912, "solid_total_recommended"), (900, 5.638, "solid_total_recommended"),
        (1000, 6.396, "solid_total_recommended"), (1100, 7.215, "solid_total_recommended"),
        (1200, 8.089, "solid_total_recommended"), (1235.08, 8.415, "solid_total_recommended"),
        (1235.08, 17.30, "liquid_total_recommended"), (1300, 17.85, "liquid_total_recommended"),
        (1400, 18.69, "liquid_total_recommended"), (1500, 19.54, "liquid_total_recommended"),
        (1600, 20.38, "liquid_total_recommended"), (1673, 21.00, "liquid_total_recommended"),
    ]
    resist_rows = [[t, f"{rho_1e8:.5f}", f"{rho_1e8 * 1e-8:.10e}", "Matula 1979 Table 11",
                    state, "recommended; total resistivity; corrected for thermal expansion"]
                   for t, rho_1e8, state in matula]
    resist_rows += [
        [1100, "", "", "Cagran Eq.(7)", "solid_initial_geometry_fit",
         "rho_microohm_m=-0.0029415+7.2481e-5*T; T in K"],
        [1234, "", "", "Cagran Eq.(8)", "liquid_initial_geometry_fit",
         "rho_microohm_m=0.0554+8.3336e-5*T; T in K"],
    ]
    write_csv("Ag_electrical_resistivity.csv", [
        "temperature_K", "resistivity_value_1e-8_ohm_m", "resistivity_ohm_m",
        "source_or_equation", "state_or_fit", "status_and_note",
    ], resist_rows)

    # Abdullaev density checkpoints plus directly reported liquid viscosity and
    # surface-tension coefficients. Blank sigma is deliberate: Fig. 2 is not a table.
    fluid_rows = []
    for t in [300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1234.93]:
        fluid_rows.append(["density", "solid", t, f"{rho_ag_solid(t):.5f}", "kg/m^3", f"{rho_ag_solid(t):.5f}", "", "",
                           "equation (1), measured density fit", "Abdullaev et al. 2020, Table 1, Eq.(1)"])
    for t in [1234.93, 1300, 1400, 1500, 1600, 1631]:
        rho = 9327.1 - 0.9956 * (t - 1234.93)
        fluid_rows.append(["density", "liquid", t, f"{rho:.5f}", "kg/m^3", f"{rho:.5f}", "", "",
                           "equation (2), measured density fit", "Abdullaev et al. 2020, Table 1, Eq.(2)"])
    for t, eta, eta_sd in [(1273, 3.520, 0.372), (1323, 2.904, 0.225),
                           (1373, 2.812, 0.323), (1423, 2.649, 0.165),
                           (1473, 2.540, 0.110), (1523, 2.498, 0.229),
                           (1573, 2.399, 0.176), (1623, 2.332, 0.195),
                           (1673, 2.214, 0.635)]:
        fluid_rows.append(["dynamic_viscosity", "liquid", t, eta, "mPa s", f"{eta*1e-3:.8f}", "", eta_sd,
                           "measured; Table 4, 99.9 mass% Ag", "Dubberstein et al. 2016, Table 4"])
    fluid_rows.append(["dynamic_viscosity_fit", "liquid", "",
                       "ln(eta_mPa_s)=-0.4631+2089/T", "equation",
                       "ln(eta_Pa_s)=-7.370855+2089/T", "", "",
                       "empirical fit over the reported liquid-Ag measurement interval",
                       "Dubberstein et al. 2016, fit accompanying Table 4"])
    for oxygen, dsdT, note in [(1e-10, -0.245e-3, "low oxygen, linear clean-surface branch"),
                               (1.0, 0.609e-3, "oxygen adsorption; positive slope below sigma maximum")]:
        fluid_rows.append(["surface_tension_temperature_coefficient", "liquid", "",
                           dsdT, "N/(m K)", dsdT, oxygen, "", note,
                           "Hibiya and Ozawa 2013, printed p.211, Fig.2 discussion"])
    fluid_rows.append(["contact_angle", "liquid_Ag_on_Ti3SiC2", 1243.15, 72.0, "degree", f"{(72.0*3.141592653589793/180.0):.8f}", "vacuum", "",
                       "sessile drop at 970 C; wettability does not determine interface thermal resistance",
                       "Zhang et al. 2018, printed p.811, Fig.1 discussion"])
    write_csv("Ag_fluid_properties.csv", [
        "property", "phase_or_pair", "temperature_K", "value", "original_units", "value_SI",
        "oxygen_partial_pressure_Pa_or_atmosphere", "standard_deviation_or_reported_error",
        "data_status_and_note", "source_location",
    ], fluid_rows)

    ti_rows = []
    for t in [300.0, 373.15, 1273.15, 1473.15, 1873.15, 2300.0]:
        cp_mol = 164.4 - 16419.0 / t
        k = 38.6 - 0.0045 * t
        status = "primary_fit_or_measurement_range" if t <= 1273.15 else (
            "reported_fit_extrapolation" if t <= 1473.15 else "extrapolation_do_not_use_for_prediction")
        ti_rows.append([f"{t:.2f}", f"{cp_mol:.8f}", f"{cp_mol/MOLAR_MASS_TI:.8f}",
                        f"{k:.8f}", f"{RHO_TI_REF:.5f}", status,
                        "Barsoum et al. 1999, printed pp.435-436, Figs.5-6; density 4.5 g/cm^3"])
    write_csv("Ti3SiC2_solid_properties.csv", [
        "temperature_K", "Cp_molar_J_mol_K", "Cp_mass_J_kg_K", "thermal_conductivity_W_m_K",
        "reference_density_kg_m3", "temperature_status", "source_location",
    ], ti_rows)

    payload = {
        "schema_version": 2,
        "database_scope": "single-temperature enthalpy FVM properties; Ag flow properties stored for a later stage",
        "Ag": {
            "molar_mass_kg_mol": MOLAR_MASS_AG,
            "melting_temperature_K": TM_AG,
            "latent_heat_fusion_J_mol": 11296.8,
            "latent_heat_fusion_J_kg": 11296.8 / MOLAR_MASS_AG,
            "solid_cp_formula_J_mol_K": "4.184*(5.272+0.001836*T+22200/T**2)",
            "solid_cp_range_K": [298.15, 1235.08],
            "liquid_cp_J_mol_K": 33.472,
            "liquid_cp_range_K": [1235.08, 2000.0],
            "solid_reference_density_fixed_kg_m3": RHO_AG_REF,
            "solid_density_equation_not_used_in_first_solver": "10503.8-0.6453*(T-293.15)-6.94e-5*(T-293.15)^2",
            "liquid_density_equation_not_used_in_first_solver": "9327.1-0.9956*(T-1234.93)",
            "density_source_transition_change_percent": {"value": 5.16, "uncertainty": 0.14},
            "thermal_conductivity_policy": "separate solid/liquid interpolation tables; no interpolation across Tm; linear phase-fraction blend only inside 1 K numerical interval",
            "liquid_conductivity_status": "Ho recommended/provisional estimate; no direct molten-Ag thermal-conductivity experiment in supplied Ho copy; estimated uncertainty about 20 percent",
        },
        "Ti3SiC2": {
            "molar_mass_kg_mol": MOLAR_MASS_TI,
            "reference_density_kg_m3": RHO_TI_REF,
            "cp_formula_J_mol_K": "164.4-16419/T",
            "conductivity_formula_W_m_K": "38.6-0.0045*T",
            "primary_solver_upper_temperature_K": 1473.15,
            "density_policy": "fixed 4500 kg/m^3; no thermal expansion or structural change",
            "expansion_note": "CSV value alpha_volume=8.90e-6 K^-1 is not the tensor-derived volumetric coefficient; 2 alpha_a+alpha_c is about 26.95e-6 K^-1",
        },
        "solver_policy": {
            "reference_temperature_K": 300.0,
            "Ag_mushy_interval_width_baseline_K": 1.0,
            "fixed_density_in_enthalpy_solver": True,
            "interface_resistance_baseline_m2_K_W": 0.0,
            "interface_resistance_sensitivity_m2_K_W": [1e-10, 1e-9, 1e-8],
        },
        "source_notes": [
            "Pankratz compiled heat of fusion is 11296.8 J/mol; Feber and Cagran experimental values are preserved as alternatives, not averaged.",
            "Stage 4 database and checksums remain unchanged; v2 is an independent candidate.",
            "Pankratz data above 1600 K are marked as source extrapolation; Ti3SiC2 solver stops at 1473.15 K.",
        ],
    }
    (HERE / "material_parameters.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote v2 property tables to {HERE}")


if __name__ == "__main__":
    main()
