#!/usr/bin/env python3
"""Audit the selected Ag–W TTM Baseline v1 and run one W energy-ledger pilot."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


PACK = Path(__file__).resolve().parents[1]
DATA = PACK / "ttm_data"
BASE = DATA / "AgW_TTM_baseline_v1"
W_LIN_CSV = DATA / "W_TTM_Lin2008_digitized_300_20000K.csv"
AG_LIN_CSV = BASE / "Ag_TTM_Lin2008_digitized_300_20000K.csv"
AG_K_TABLE = DATA / "XTANT3_TTM_parameters_v3/K_metals/K_Ag.dat"
OUT_TRAJECTORY = BASE / "W_Lin_baseline_energy_pilot.csv"
OUT_AG_KE = BASE / "Ag_ke_Lin_Drude_room_anchor.csv"
AGW_GRID_SMOKE = PACK / "models/ag_w_qsgs_200A/ttm/agw_ttm_grid_smoke_metrics.json"
OUT_JSON = BASE / "baseline_v1_validation.json"
OUT_MD = BASE / "BASELINE_V1_VALIDATION_20261124.md"

T0 = 500.0
Q_INPUT = 68_650.0 * 25.0
PULSE_FWHM_PS = 0.100
PULSE_CENTER_PS = 0.200
DT_PS = 0.001
T_END_PS = 20.0
WRITE_EVERY = 10
M_W = 0.18384
G0_DUFFY_EXP = 1.4e17
G0_DUFFY_SIGMA = 0.3e17
G0_DUFFY_THEORY = 1.65e17
KE0_W = 0.75 * 174.0
KE0_AG = 429.0  # User-supplied 300 K bulk total; Ag is electronically dominated.
AG_W_INTERFACE_H = 1.0e9  # W m^-2 K^-1, provisional method analogue only


def read_csv_table(path: Path, ce_col: str, g_col: str) -> tuple[list[float], list[float], list[float]]:
    te, ce, g = [], [], []
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            te.append(float(row["Te_K"]))
            ce.append(float(row[ce_col]))
            g.append(float(row[g_col]))
    return te, ce, g


def read_ag_k() -> tuple[list[float], list[float], list[float], list[float]]:
    te, k_total, k_eph, k_ee = [], [], [], []
    with AG_K_TABLE.open(encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            columns = line.split()
            te.append(float(columns[0]))
            k_total.append(float(columns[1]))
            k_eph.append(float(columns[2]))
            k_ee.append(float(columns[3]))
    return te, k_total, k_eph, k_ee


def interp(xs: list[float], ys: list[float], x: float) -> float:
    if x < xs[0] - 1e-9 or x > xs[-1] + 1e-9:
        raise ValueError(f"{x:g} K lies outside the source grid {xs[0]:g}–{xs[-1]:g} K")
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    f = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] + f * (ys[hi] - ys[lo])


def integrate_table(xs: list[float], ys: list[float], start: float, end: float) -> float:
    if abs(end - start) < 1e-14:
        return 0.0
    sign = 1.0 if end > start else -1.0
    a, b = (start, end) if end > start else (end, start)
    points = [a] + [x for x in xs if a < x < b] + [b]
    return sign * sum(0.5 * (interp(xs, ys, l) + interp(xs, ys, r)) * (r - l) for l, r in zip(points, points[1:]))


def cp_w_molar(t: float) -> float:
    # Tolias (2017), recommended solid-W White–Minges fit, valid through 3080 K.
    return 21.868372 + 8.068661e-3 * t - 3.756196e-6 * t**2 + 1.075862e-9 * t**3 + 1.406637e4 / t**2


def rho_w(t: float) -> float:
    d = t - 293.15
    return (19.25 - 2.66207e-4 * d - 3.0595e-9 * d**2 - 9.5185e-12 * d**3) * 1000.0


def cl_w(t: float) -> float:
    return cp_w_molar(t) * rho_w(t) / M_W


def integrate_cl(start: float, end: float) -> float:
    if abs(end - start) < 1e-14:
        return 0.0
    sign = 1.0 if end > start else -1.0
    a, b = (start, end) if end > start else (end, start)
    n = 100
    h = (b - a) / n
    values = [cl_w(a + i * h) for i in range(n + 1)]
    result = values[0] + values[-1] + 4 * sum(values[1:-1:2]) + 2 * sum(values[2:-1:2])
    return sign * result * h / 3.0


def pulse_sigma() -> float:
    return PULSE_FWHM_PS / (2.0 * math.sqrt(2.0 * math.log(2.0)))


def pulse_rate(t: float) -> float:
    sigma = pulse_sigma()
    lo = math.erf((0.0 - PULSE_CENTER_PS) / (math.sqrt(2.0) * sigma))
    hi = math.erf((T_END_PS - PULSE_CENTER_PS) / (math.sqrt(2.0) * sigma))
    normalization = 0.5 * (hi - lo)
    return Q_INPUT / (sigma * math.sqrt(2.0 * math.pi) * normalization) * math.exp(-0.5 * ((t - PULSE_CENTER_PS) / sigma) ** 2)


def pulse_energy(t: float) -> float:
    sigma = pulse_sigma()
    lo = math.erf((0.0 - PULSE_CENTER_PS) / (math.sqrt(2.0) * sigma))
    hi = math.erf((T_END_PS - PULSE_CENTER_PS) / (math.sqrt(2.0) * sigma))
    value = math.erf((t - PULSE_CENTER_PS) / (math.sqrt(2.0) * sigma))
    return Q_INPUT * min(1.0, max(0.0, (value - lo) / (hi - lo)))


def ke_w(te: float, tl: float, w_t: list[float], w_ce: list[float]) -> float:
    ce_te = interp(w_t, w_ce, te)
    ce_300 = interp(w_t, w_ce, 300.0)
    return KE0_W * ce_te / ce_300 * (300.0 / tl) ** 1.5


def ke_ag(te: float, tl: float, ag_t: list[float], ag_ce: list[float], a_over_b: float) -> float:
    """Lin-style Drude closure, normalized to the 300 K Ag room-temperature anchor."""
    ce_te = interp(ag_t, ag_ce, te)
    ce_300 = interp(ag_t, ag_ce, 300.0)
    scattering_ratio = (300.0 + a_over_b * 300.0**2) / (tl + a_over_b * te**2)
    return KE0_AG * (ce_te / ce_300) * scattering_ratio


def rhs(t: float, te: float, tl: float, t_grid: list[float], ce_grid: list[float], g_grid: list[float]) -> tuple[float, float]:
    ce = interp(t_grid, ce_grid, te)
    g = interp(t_grid, g_grid, te)
    exchange_j_m3_ps = g * 1.0e-12 * (te - tl)
    return (pulse_rate(t) - exchange_j_m3_ps) / ce, exchange_j_m3_ps / cl_w(tl)


def rk4(t: float, te: float, tl: float, dt: float, t_grid: list[float], ce_grid: list[float], g_grid: list[float]) -> tuple[float, float]:
    a1, b1 = rhs(t, te, tl, t_grid, ce_grid, g_grid)
    a2, b2 = rhs(t + dt/2, te + dt*a1/2, tl + dt*b1/2, t_grid, ce_grid, g_grid)
    a3, b3 = rhs(t + dt/2, te + dt*a2/2, tl + dt*b2/2, t_grid, ce_grid, g_grid)
    a4, b4 = rhs(t + dt, te + dt*a3, tl + dt*b3, t_grid, ce_grid, g_grid)
    return te + dt*(a1 + 2*a2 + 2*a3 + a4)/6, tl + dt*(b1 + 2*b2 + 2*b3 + b4)/6


def run_energy_pilot(t_grid: list[float], ce_grid: list[float], g_grid: list[float]) -> dict[str, float]:
    te = tl = T0
    steps = round(T_END_PS / DT_PS)
    rows = []
    max_residual = 0.0
    max_te = te
    t50 = None
    t90 = None
    for i in range(steps + 1):
        t = i * DT_PS
        ue = integrate_table(t_grid, ce_grid, T0, te)
        ul = integrate_cl(T0, tl)
        ein = pulse_energy(t)
        residual = (ue + ul - ein) / Q_INPUT
        max_residual = max(max_residual, abs(residual))
        max_te = max(max_te, te)
        fraction_l = ul / Q_INPUT
        if t >= PULSE_CENTER_PS and t50 is None and fraction_l >= 0.5:
            t50 = t
        if t >= PULSE_CENTER_PS and t90 is None and fraction_l >= 0.9:
            t90 = t
        if i % WRITE_EVERY == 0:
            rows.append((t, te, tl, ein, ue, ul, residual))
        if i < steps:
            te, tl = rk4(t, te, tl, DT_PS, t_grid, ce_grid, g_grid)

    with OUT_TRAJECTORY.open("w", encoding="utf-8", newline="") as stream:
        fields = ["time_ps", "Te_K", "Tl_K", "input_energy_J_m3", "electron_energy_J_m3", "lattice_energy_J_m3", "energy_ledger_residual_fraction"]
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)
    return {
        "pulse_energy_density_J_m3": Q_INPUT,
        "pulse_FWHM_ps": PULSE_FWHM_PS,
        "time_step_ps": DT_PS,
        "end_time_ps": T_END_PS,
        "max_Te_K": max_te,
        "final_Te_K": te,
        "final_Tl_K": tl,
        "final_lattice_energy_J_m3": integrate_cl(T0, tl),
        "final_lattice_energy_fraction": integrate_cl(T0, tl) / Q_INPUT,
        "time_to_50pct_pulse_energy_in_lattice_ps": t50,
        "time_to_90pct_pulse_energy_in_lattice_ps": t90,
        "max_abs_energy_ledger_residual_fraction": max_residual,
        "energy_ledger_pass_threshold_fraction": 1e-8,
        "energy_ledger_status": "PASS" if max_residual < 1e-8 else "FAIL",
    }


def main() -> None:
    w_t, w_ce, w_g = read_csv_table(W_LIN_CSV, "Ce_Lin2008_J_m-3_K-1", "G_Lin2008_digitized_W_m-3_K-1")
    ag_t, ag_ce, ag_g = read_csv_table(AG_LIN_CSV, "Ce_Lin2008_digitized_J_m-3_K-1", "G_Lin2008_digitized_W_m-3_K-1")
    ag_k_t, ag_k_total, ag_k_eph, ag_k_ee = read_ag_k()
    row_300 = (ag_k_t.index(300.0) if 300.0 in ag_k_t else min(range(len(ag_k_t)), key=lambda i: abs(ag_k_t[i] - 300.0)))
    if abs(ag_k_t[row_300] - 300.0) > 1e-9:
        raise ValueError("The Ag XTANT K table has no exact 300 K anchor row")
    # Lin et al. write tau_e^-1 = A*Te^2 + B*Tl. The XTANT K file lists
    # separate e-ph and e-e conductivities; under the channel model, their
    # inverse ratio at Te=Tl=300 K gives A/B = (k_eph/k_ee)/300.
    ag_a_over_b = ag_k_eph[row_300] / ag_k_ee[row_300] / 300.0

    g_w_300 = interp(w_t, w_g, 300.0)
    exp_low, exp_high = G0_DUFFY_EXP - G0_DUFFY_SIGMA, G0_DUFFY_EXP + G0_DUFFY_SIGMA
    g0_check = {
        "Te_K": 300.0,
        "Lin_digitized_G_W_m3_K": g_w_300,
        "Lin_300K_value_has_low_temperature_extrapolation": True,
        "Daraszewicz_experimental_G0_W_m3_K": G0_DUFFY_EXP,
        "Daraszewicz_experimental_1sigma_interval_W_m3_K": [exp_low, exp_high],
        "Daraszewicz_ab_initio_G0_W_m3_K": G0_DUFFY_THEORY,
        "Lin_over_experiment_nominal": g_w_300 / G0_DUFFY_EXP,
        "Lin_relative_difference_from_experiment_nominal": (g_w_300 - G0_DUFFY_EXP) / G0_DUFFY_EXP,
        "Lin_exceeds_experiment_upper_1sigma_by_fraction": (g_w_300 - exp_high) / exp_high,
        "Lin_relative_difference_from_ab_initio": (g_w_300 - G0_DUFFY_THEORY) / G0_DUFFY_THEORY,
        "within_experiment_1sigma": exp_low <= g_w_300 <= exp_high,
        "interpretation": "Lin's digitized low-Te value does not reproduce the Daraszewicz experimental G0 within its quoted 1-sigma range. The Lin value is itself extrapolated below 600 K; retain the discrepancy as a model/source limitation rather than rescaling the Lin curve.",
    }

    ke_checks = {
        "formula": "ke_W(Te,Tl)=130.5*[Ce_W_Lin(Te)/Ce_W_Lin(300)]*(300/Tl)^1.5 W m^-1 K^-1",
        "ke0_W_m_K": KE0_W,
        "ke_W_at_Te_Tl_300K_W_m_K": ke_w(300.0, 300.0, w_t, w_ce),
        "ke_W_at_Te_10000_Tl_300K_W_m_K": ke_w(10000.0, 300.0, w_t, w_ce),
        "ke_W_at_Te_10000_Tl_1000K_W_m_K": ke_w(10000.0, 1000.0, w_t, w_ce),
        "reference_fit_range_K": [300, 2500],
        "high_temperature_use": "extrapolative closure; high-Te nonequilibrium validity is not independently established",
        "Ag_ke_source": "Lin et al. (2008), Sec. II C Drude relation, room-anchored with Ag Ce from Fig. 2(c)",
        "Ag_ke_formula": "429 * [Ce_Ag_Lin(Te)/Ce_Ag_Lin(300)] * [(300 + a*300^2)/(Tl + a*Te^2)] W m^-1 K^-1, with a=A/B",
        "Ag_ke0_W_m_K": KE0_AG,
        "Ag_ke_anchor_provenance": "429 W m^-1 K^-1 at 300 K as supplied by the user; total bulk conductivity is used as an electronically-dominated proxy",
        "Ag_A_over_B_K_minus1": ag_a_over_b,
        "Ag_A_over_B_method": "(XTANT K_Ag 300 K kappa_e_ph / kappa_e_e) / 300 K, assuming the listed channel conductivities follow Lin's inverse-relaxation-time decomposition",
        "Ag_XTANT_300K_channel_values_W_m_K": {"kappa_tot": ag_k_total[row_300], "kappa_e_ph": ag_k_eph[row_300], "kappa_e_e": ag_k_ee[row_300]},
        "Ag_ke_grid_K": [ag_k_t[0], ag_k_t[-1]],
        "Ag_XTANT_kappa_tot_at_300K_W_m_K": interp(ag_k_t, ag_k_total, 300.0),
        "Ag_ke_at_Te_Tl_300K_W_m_K": ke_ag(300.0, 300.0, ag_t, ag_ce, ag_a_over_b),
        "Ag_ke_at_Te_10000_Tl_300K_W_m_K": ke_ag(10000.0, 300.0, ag_t, ag_ce, ag_a_over_b),
        "Ag_ke_at_Te_Tl_10000K_W_m_K": ke_ag(10000.0, 10000.0, ag_t, ag_ce, ag_a_over_b),
        "Ag_ke_status": "exploratory Lin-Drude closure; room anchored; A/B is inferred from the XTANT 300 K channel ratio and is not independently validated",
    }

    with OUT_AG_KE.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["Te_K", "Tl_K", "Ce_Ag_Lin_J_m3_K", "ke_Ag_W_m_K"])
        for te in ag_t:
            writer.writerow([f"{te:.9g}", "300", f"{interp(ag_t, ag_ce, te):.12g}", f"{ke_ag(te, 300.0, ag_t, ag_ce, ag_a_over_b):.12g}"])
            writer.writerow([f"{te:.9g}", f"{te:.9g}", f"{interp(ag_t, ag_ce, te):.12g}", f"{ke_ag(te, te, ag_t, ag_ce, ag_a_over_b):.12g}"])

    energy_pilot = run_energy_pilot(w_t, w_ce, w_g)
    grid_smoke = json.loads(AGW_GRID_SMOKE.read_text(encoding="utf-8")) if AGW_GRID_SMOKE.is_file() else None
    summary = {
        "baseline": "AgW_TTM_baseline_v1",
        "selection_status": "selected exploratory parameter baseline; W low-temperature experimental G0 consistency check is not passed",
        "parameter_branches": {
            "Ag_Ce": {"source": "Lin et al. (2008), Fig. 2(c), digitized candidate", "table": str(AG_LIN_CSV.relative_to(PACK)), "range_K": [ag_t[0], ag_t[-1]]},
            "Ag_G": {"source": "Lin et al. (2008), Fig. 2(d), digitized candidate", "table": str(AG_LIN_CSV.relative_to(PACK)), "range_K": [ag_t[0], ag_t[-1]]},
            "Ag_ke": {"source": "Lin et al. (2008), Sec. II C Drude form, normalized to the user-provided 300 K Ag room-temperature bulk conductivity anchor", "formula": ke_checks["Ag_ke_formula"], "a_over_b_K_minus1": ag_a_over_b, "Ag_ke_at_Te_Tl_300K_W_m_K": ke_checks["Ag_ke_at_Te_Tl_300K_W_m_K"], "Ag_A_over_B_K_minus1": ag_a_over_b, "table": str(OUT_AG_KE.relative_to(PACK)), "range_K": [ag_t[0], ag_t[-1]], "temperature_dependent_variables": ["Te", "Tl"], "status": ke_checks["Ag_ke_status"]},
            "W_Ce": {"source": "Lin et al. (2008), Fig. 7(c), digitized candidate", "table": str(W_LIN_CSV.relative_to(PACK)), "range_K": [w_t[0], w_t[-1]]},
            "W_G": {"source": "Lin et al. (2008), Fig. 7(d), digitized candidate", "table": str(W_LIN_CSV.relative_to(PACK)), "range_K": [w_t[0], w_t[-1]]},
            "W_ke": {"source": "Daraszewicz et al. (2014), Drude model equation (6)", **ke_checks},
            "W_G0_reference": {"source": "Daraszewicz et al. (2014), measured and ab-initio near-equilibrium anchors", **g0_check},
        },
        "timescale_at_reference_states": {
            "formula_tau_eph": "Ce/G",
            "formula_tau_eq": "Ce*Cl/[G*(Ce+Cl)] at constant coefficients",
            "W_300K": {"Te_K": 300.0, "Cl_J_m3_K": cl_w(300.0), "tau_eph_ps": interp(w_t, w_ce, 300.0)/g_w_300*1e12, "tau_eq_ps": (interp(w_t, w_ce, 300.0)*cl_w(300.0)/(g_w_300*(interp(w_t, w_ce, 300.0)+cl_w(300.0))))*1e12},
            "W_500K": {"Te_K": 500.0, "Cl_J_m3_K": cl_w(500.0), "tau_eph_ps": interp(w_t, w_ce, 500.0)/interp(w_t, w_g, 500.0)*1e12, "tau_eq_ps": (interp(w_t, w_ce, 500.0)*cl_w(500.0)/(interp(w_t, w_g, 500.0)*(interp(w_t, w_ce, 500.0)+cl_w(500.0))))*1e12},
            "Cl_source": "Tolias 2017 recommended solid-W White-Minges Cp and density fits; Cp is used as a narrow-range Cv proxy",
        },
        "single_branch_W_energy_pilot": energy_pilot,
        "synthetic_Ag_W_grid_smoke": grid_smoke,
        "Ag_W_solver_artifacts": {
            "reference_module": "models/ag_w_qsgs_200A/ttm/agw_ttm_reference.py",
            "gated_LAMMPS_driver": "models/ag_w_qsgs_200A/ttm/agw_ttm_md_pilot.py",
            "grid_smoke_script": "models/ag_w_qsgs_200A/ttm/run_agw_ttm_grid_smoke.py",
            "grid_smoke_status": grid_smoke.get("status") if grid_smoke else "NOT_RUN",
            "atomistic_TTM_MD_status": "NOT_RUN_STATIC_STRUCTURE_GATE_NOT_PASSED",
        },
        "provisional_Ag_W_electronic_interface": {
            "conductance_W_m2_K": AG_W_INTERFACE_H,
            "model": "finite face-based electron conductance; equal-and-opposite interphase energy transfer",
            "source_basis": "Li et al. (2024) Au/Pt and Au/Al method analogues; selected value is not a direct Ag/W measurement",
            "spatial_rule": "unlike-phase fractions across shared grid-cell faces; sub-cell interface area omitted",
            "status": "numerically implemented and smoke-tested; material value remains exploratory",
        },
        "output_trajectory": OUT_TRAJECTORY.name,
        "output_Ag_ke_table": OUT_AG_KE.name,
        "limitations": [
            "The W Lin G curve below 600 K is extrapolated and differs from Daraszewicz's experimental G0; no empirical rescaling was applied.",
            "The Ag Lin curves are digitized from plots; the low-temperature G plateau and short 5000–6607 K bridge are documented approximations.",
            "Ag ke is a room-anchored Lin Drude closure; the user-provided 429 W m^-1 K^-1 anchor approximates total as electronic-dominated, and the A/B scattering ratio is inferred from XTANT's room-temperature channel split.",
            "Ag ke high-Te nonequilibrium behavior is not independently validated; the raw XTANT kappa_tot table remains archived and is not used as the room-temperature amplitude.",
            "Ag/W electron-interface conductance is set to a provisional 1 GW m^-2 K^-1 from non-Ag/W analogues; its material value and omitted sub-cell area are not validated.",
            "The Ag/W solver is an exploratory Python reference bridge, not a compiled LAMMPS TTM fix; the gated atomistic driver has not run because no static structure passed.",
            "No Ag lattice-conductivity curve or jointly validated Ag/W interface conductance is supplied.",
            "The 0D energy pilot uses one normalized temporal Gaussian and tests energy accounting only; it is not the spatial laser calculation and does not test melting.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_report(summary)
    print(json.dumps({"G0_check": g0_check, "ke_checks": ke_checks, "energy_pilot": energy_pilot}, indent=2))
    print(f"Wrote {OUT_JSON.relative_to(PACK)} and {OUT_MD.relative_to(PACK)}")


def write_report(summary: dict) -> None:
    g = summary["parameter_branches"]["W_G0_reference"]
    k = summary["parameter_branches"]["W_ke"]
    e = summary["single_branch_W_energy_pilot"]
    t300 = summary["timescale_at_reference_states"]["W_300K"]
    t500 = summary["timescale_at_reference_states"]["W_500K"]
    ag_k = summary["parameter_branches"]["Ag_ke"]
    report = f"""# Ag–W TTM Baseline v1: validation record

Date: 2026-11-24
Generated by `scripts/validate_agw_ttm_baseline_v1.py`.

## Selected parameter set

- **Ag `Ce(Te), G(Te)`:** candidate curves digitized from Lin et al. (2008), Fig. 2(c,d), archived with an explicit digitization report.
- **Ag `ke(Te,Tl)`:** room-anchored Lin et al. (2008) Drude closure, `{ag_k['formula']}`. The 300 K value is `{ag_k['Ag_ke_at_Te_Tl_300K_W_m_K']:.1f} W m⁻¹ K⁻¹`. The room anchor uses the user's cited Ag bulk total `429 W m⁻¹ K⁻¹` as an electronic-dominant proxy. `A/B={ag_k['Ag_A_over_B_K_minus1']:.4g} K⁻¹` is inferred from the XTANT 300 K channel ratio; this is a provisional shape closure, not an Ag-specific measured `ke(Te,Tl)` surface.
- **W `Ce(Te), G(Te)`:** Lin et al. (2008), Fig. 7(c,d), candidate digitized tables.
- **W `ke(Te,Tl)`:** Daraszewicz et al. (2014), equation (6): `{k['formula']}`. It gives `{k['ke_W_at_Te_Tl_300K_W_m_K']:.1f} W m⁻¹ K⁻¹` at 300 K. The published exponent was fitted to equilibrium conductivity behavior over 300–2500 K; use beyond that interval is an extrapolation.
- **W low-temperature `G0` reference:** Daraszewicz's measured `1.4(3)×10¹⁷ W m⁻³ K⁻¹` and ab-initio `1.65×10¹⁷ W m⁻³ K⁻¹` remain validation anchors, not replacements for the Lin `G(Te)` curve.

## W low-temperature coupling check

At 300 K, the Lin digitization gives `G={g['Lin_digitized_G_W_m3_K']:.4g} W m⁻³ K⁻¹`. This is `{100*g['Lin_relative_difference_from_experiment_nominal']:.1f}%` above the experimental central value and `{100*g['Lin_exceeds_experiment_upper_1sigma_by_fraction']:.1f}%` above its upper 1σ bound. It is also `{100*g['Lin_relative_difference_from_ab_initio']:.1f}%` above the paper's ab-initio anchor.

**The Lin curve does not reproduce the Daraszewicz experimental `G0` within the quoted 1σ range.** Lin's 300 K value is itself in the documented below-600 K extrapolation region. This baseline therefore retains a visible source discrepancy; the Lin curve was not rescaled to force agreement.

## Ag/W electronic interface assumption

The first Ag–W grid implementation uses a face-based Ag/W electron conductance of `{AG_W_INTERFACE_H/1e9:.1f} GW m⁻² K⁻¹`. This provisional value lies between the Au/Pt and Au/Al method analogues reported by Li et al. (2024); it is not an Ag/W measurement. Transfers are equal and opposite between the phase electron grids. The mesh couples unlike-phase fractions across shared cell faces and omits sub-cell interface area, so this numerical closure is not a validated physical interface model.

## Single-branch W energy-accounting pilot

The one-branch pilot uses Lin W `Ce/G`, a 100 fs normalized temporal Gaussian with the same volumetric energy in the electron subsystem, and Tolias's solid-W lattice heat-capacity fit. It is a 0D energy-ledger check; spatial conduction, `ke`, melting, and the Ag/W interface do not enter this test.

- Initial state: `Te=Tl=500 K`; input energy density `{e['pulse_energy_density_J_m3']:.5g} J/m³`.
- Peak electron temperature: `{e['max_Te_K']:.3f} K`.
- Time for half of the input energy to appear in the lattice: `{e['time_to_50pct_pulse_energy_in_lattice_ps']:.3f} ps`; 90%: `{e['time_to_90pct_pulse_energy_in_lattice_ps']:.3f} ps`.
- Final lattice energy at 20 ps: `{e['final_lattice_energy_J_m3']:.5g} J/m³` (`{100*e['final_lattice_energy_fraction']:.3f}%` of input; the remaining energy stays in the electron reservoir at equilibrium).
- Maximum absolute total-energy ledger residual: `{e['max_abs_energy_ledger_residual_fraction']:.3g}` of input; threshold `{e['energy_ledger_pass_threshold_fraction']:.1g}`; **{e['energy_ledger_status']}**.

The two initial-state timescales are both reported. At 300 K, `Ce/G={t300['tau_eph_ps']:.4g} ps` and `τeq={t300['tau_eq_ps']:.4g} ps`; at 500 K, `Ce/G={t500['tau_eph_ps']:.4g} ps` and `τeq={t500['tau_eq_ps']:.4g} ps`. Here `τeq=Ce Cl/[G(Ce+Cl)]` is the constant-coefficient two-reservoir temperature-difference time; `Cl` comes from Tolias's solid-W `Cp(T)` and density fits, with `Cp` used as a narrow-range `Cv` proxy.

## Status and limits

`AgW_TTM_baseline_v1` is selected as the repository's exploratory material-parameter baseline. The W energy accounting passes. The independent W low-temperature experimental `G0` consistency check does not pass for the Lin digitized value; this mismatch remains attached to the baseline. Ag `G` is digitized and graph-based at low temperature. The corrected Ag conductivity closure now has a reasonable 300 K anchor, while its high-Te nonequilibrium shape remains exploratory and explicitly dependent on the inferred `A/B` ratio.

The Ag/W electronic interface is numerically closed by the provisional face rule above, but its conductance and sub-cell area approximation are not physically validated. The reused-core Ag/W grid smoke check status is **{summary['Ag_W_solver_artifacts']['grid_smoke_status']}**; it exercises a synthetic grid, not atoms or LAMMPS. The gated Ag/W LAMMPS coupling entry point is `models/ag_w_qsgs_200A/ttm/agw_ttm_md_pilot.py`; it remains unrun until a statically accepted Ag–W initial configuration and explicit pulse envelope are available.

## Files

- Parameter/source definition: `baseline_parameters.json`
- W Lin inputs: `../W_TTM_Lin2008_digitized_300_20000K.csv`
- Ag Lin digitization: `Ag_TTM_Lin2008_digitized_300_20000K.csv` and `Ag_TTM_Lin2008_digitization_report.json`
- Ag Drude conductivity grid: `Ag_ke_Lin_Drude_room_anchor.csv`
- Energy-pilot trajectory: `W_Lin_baseline_energy_pilot.csv`
- Machine-readable result: `baseline_v1_validation.json`
- Ag/W solver and gated pilot: `../../models/ag_w_qsgs_200A/ttm/agw_ttm_reference.py`, `../../models/ag_w_qsgs_200A/ttm/agw_ttm_md_pilot.py`
- Reproduction script: `../../scripts/validate_agw_ttm_baseline_v1.py`
"""
    OUT_MD.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
