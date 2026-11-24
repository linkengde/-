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


def read_csv_table(path: Path, ce_col: str, g_col: str) -> tuple[list[float], list[float], list[float]]:
    te, ce, g = [], [], []
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            te.append(float(row["Te_K"]))
            ce.append(float(row[ce_col]))
            g.append(float(row[g_col]))
    return te, ce, g


def read_ag_k() -> tuple[list[float], list[float]]:
    te, ke = [], []
    with AG_K_TABLE.open(encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            columns = line.split()
            te.append(float(columns[0]))
            ke.append(float(columns[1]))  # upstream kappa_tot, assembled by Matthiessen's rule
    return te, ke


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
    ag_k_t, ag_ke = read_ag_k()

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
        "Ag_ke_source": "XTANT-3 K_Ag.dat kappa_tot column; Te-only, not a ke(Te,Tl) surface",
        "Ag_ke_grid_K": [ag_k_t[0], ag_k_t[-1]],
        "Ag_ke_at_300K_W_m_K": interp(ag_k_t, ag_ke, 300.0),
    }

    energy_pilot = run_energy_pilot(w_t, w_ce, w_g)
    summary = {
        "baseline": "AgW_TTM_baseline_v1",
        "selection_status": "selected exploratory parameter baseline; W low-temperature experimental G0 consistency check is not passed",
        "parameter_branches": {
            "Ag_Ce": {"source": "Lin et al. (2008), Fig. 2(c), digitized candidate", "table": str(AG_LIN_CSV.relative_to(PACK)), "range_K": [ag_t[0], ag_t[-1]]},
            "Ag_G": {"source": "Lin et al. (2008), Fig. 2(d), digitized candidate", "table": str(AG_LIN_CSV.relative_to(PACK)), "range_K": [ag_t[0], ag_t[-1]]},
            "Ag_ke": {"source": "XTANT-3 supplied Ag K table kappa_tot column", "table": str(AG_K_TABLE.relative_to(PACK)), "range_K": [ag_k_t[0], ag_k_t[-1]], "temperature_dependent_variables": ["Te"]},
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
        "output_trajectory": OUT_TRAJECTORY.name,
        "limitations": [
            "The W Lin G curve below 600 K is extrapolated and differs from Daraszewicz's experimental G0; no empirical rescaling was applied.",
            "The Ag Lin curves are digitized from plots; the low-temperature G plateau and short 5000–6607 K bridge are documented approximations.",
            "Ag ke is Te-only and XTANT warns that room-temperature values may not be reliable.",
            "No Ag/W electron-interface conductance, mixed-cell electronic transport rule, or jointly validated Ag/W TTM-MD closure is supplied.",
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
    report = f"""# Ag–W TTM Baseline v1: validation record

Date: 2026-11-24
Generated by `scripts/validate_agw_ttm_baseline_v1.py`.

## Selected parameter set

- **Ag `Ce(Te), G(Te)`:** candidate curves digitized from Lin et al. (2008), Fig. 2(c,d), archived with an explicit digitization report.
- **Ag `ke(Te)`:** XTANT-3 Ag `K` table's `kappa_tot` column. This is Te-only; the source cautions that room-temperature values may not be reliable.
- **W `Ce(Te), G(Te)`:** Lin et al. (2008), Fig. 7(c,d), candidate digitized tables.
- **W `ke(Te,Tl)`:** Daraszewicz et al. (2014), equation (6): `{k['formula']}`. It gives `{k['ke_W_at_Te_Tl_300K_W_m_K']:.1f} W m⁻¹ K⁻¹` at 300 K. The published exponent was fitted to equilibrium conductivity behavior over 300–2500 K; use beyond that interval is an extrapolation.
- **W low-temperature `G0` reference:** Daraszewicz's measured `1.4(3)×10¹⁷ W m⁻³ K⁻¹` and ab-initio `1.65×10¹⁷ W m⁻³ K⁻¹` remain validation anchors, not replacements for the Lin `G(Te)` curve.

## W low-temperature coupling check

At 300 K, the Lin digitization gives `G={g['Lin_digitized_G_W_m3_K']:.4g} W m⁻³ K⁻¹`. This is `{100*g['Lin_relative_difference_from_experiment_nominal']:.1f}%` above the experimental central value and `{100*g['Lin_exceeds_experiment_upper_1sigma_by_fraction']:.1f}%` above its upper 1σ bound. It is also `{100*g['Lin_relative_difference_from_ab_initio']:.1f}%` above the paper's ab-initio anchor.

**The Lin curve does not reproduce the Daraszewicz experimental `G0` within the quoted 1σ range.** Lin's 300 K value is itself in the documented below-600 K extrapolation region. This baseline therefore retains a visible source discrepancy; the Lin curve was not rescaled to force agreement.

## Single-branch W energy-accounting pilot

The one-branch pilot uses Lin W `Ce/G`, a 100 fs normalized temporal Gaussian with the same volumetric energy in the electron subsystem, and Tolias's solid-W lattice heat-capacity fit. It is a 0D energy-ledger check; spatial conduction, `ke`, melting, and the Ag/W interface do not enter this test.

- Initial state: `Te=Tl=500 K`; input energy density `{e['pulse_energy_density_J_m3']:.5g} J/m³`.
- Peak electron temperature: `{e['max_Te_K']:.3f} K`.
- Time for half of the input energy to appear in the lattice: `{e['time_to_50pct_pulse_energy_in_lattice_ps']:.3f} ps`; 90%: `{e['time_to_90pct_pulse_energy_in_lattice_ps']:.3f} ps`.
- Final lattice energy at 20 ps: `{e['final_lattice_energy_J_m3']:.5g} J/m³` (`{100*e['final_lattice_energy_fraction']:.3f}%` of input; the remaining energy stays in the electron reservoir at equilibrium).
- Maximum absolute total-energy ledger residual: `{e['max_abs_energy_ledger_residual_fraction']:.3g}` of input; threshold `{e['energy_ledger_pass_threshold_fraction']:.1g}`; **{e['energy_ledger_status']}**.

The two initial-state timescales are both reported. At 300 K, `Ce/G={t300['tau_eph_ps']:.4g} ps` and `τeq={t300['tau_eq_ps']:.4g} ps`; at 500 K, `Ce/G={t500['tau_eph_ps']:.4g} ps` and `τeq={t500['tau_eq_ps']:.4g} ps`. Here `τeq=Ce Cl/[G(Ce+Cl)]` is the constant-coefficient two-reservoir temperature-difference time; `Cl` comes from Tolias's solid-W `Cp(T)` and density fits, with `Cp` used as a narrow-range `Cv` proxy.

## Status and limits

`AgW_TTM_baseline_v1` is selected as the repository's exploratory material-parameter baseline. Energy accounting passes. The independent W low-temperature experimental `G0` consistency check does not pass for the Lin digitized value; this mismatch remains attached to the baseline. Ag `G` is digitized and graph-based at low temperature. The selected Ag and W conductivity closures also have different source and temperature dependence, which is recorded in `baseline_parameters.json`.

This baseline does not close Ag/W electronic interface conductance, cross-phase electronic transport, or a mixed finite-volume-cell rule. The existing free-surface geometry still needs a compatible TTM solver and a structurally accepted Ag–W initial configuration before TTM-MD.

## Files

- Parameter/source definition: `baseline_parameters.json`
- W Lin inputs: `../W_TTM_Lin2008_digitized_300_20000K.csv`
- Ag Lin digitization: `Ag_TTM_Lin2008_digitized_300_20000K.csv` and `Ag_TTM_Lin2008_digitization_report.json`
- Energy-pilot trajectory: `W_Lin_baseline_energy_pilot.csv`
- Machine-readable result: `baseline_v1_validation.json`
- Reproduction script: `../../scripts/validate_agw_ttm_baseline_v1.py`
"""
    OUT_MD.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
