#!/usr/bin/env python3
"""Build a source-separated W Ce/G and electron-phonon time-scale comparison.

The outputs retain Lin, XTANT-3, and Daraszewicz inputs as distinct branches.
Daraszewicz supplies low-Te G0 anchors, not a full G(Te) curve; illustrative
anchor time scales use Lin's Ce values and are restricted to selected points.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path


PACK = Path(__file__).resolve().parents[1]
TTM = PACK / "ttm_data"
LIN_CSV = TTM / "W_TTM_Lin2008_digitized_300_20000K.csv"
XTANT_CE = TTM / "XTANT3_TTM_parameters_v3/Ce_metals/Ce_W.txt"
XTANT_G = TTM / "XTANT3_TTM_parameters_v3/G_e-ph_metals/G_e-ph_W.txt"
CSV_OUT = TTM / "W_eph_timescale_comparison.csv"
JSON_OUT = TTM / "W_eph_timescale_summary.json"
MD_OUT = TTM / "W_EPH_TIMESCALE_AUDIT_20261124.md"

DUFFY_DISPLAY_TEMPERATURES_K = (300, 500, 600, 1000)
DUFFY_G_EXPERIMENT = 1.4e17
DUFFY_G_EXPERIMENT_UNCERTAINTY = 0.3e17
DUFFY_G_AB_INITIO = 1.65e17


def load_two_column(path: Path) -> dict[int, float]:
    result: dict[int, float] = {}
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            te, value = (float(item) for item in line.split()[:2])
            result[int(round(te))] = value
    return result


def tau_ps(ce_j_m3_k: float, g_w_m3_k: float) -> float:
    # J/W = seconds; multiply by 1e12 for picoseconds.
    return ce_j_m3_k / g_w_m3_k * 1.0e12


def load_lin() -> dict[int, dict[str, float]]:
    rows: dict[int, dict[str, float]] = {}
    with LIN_CSV.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            te = int(float(row["Te_K"]))
            rows[te] = {
                "ce": float(row["Ce_Lin2008_J_m-3_K-1"]),
                "g": float(row["G_Lin2008_digitized_W_m-3_K-1"]),
            }
    return rows


def fmt(value: float | None, precision: int = 3) -> str:
    if value is None:
        return "—"
    return f"{value:.{precision}g}"


def main() -> None:
    lin = load_lin()
    xtant_ce = load_two_column(XTANT_CE)
    xtant_g = load_two_column(XTANT_G)
    temperatures = sorted(lin)

    csv_fields = [
        "Te_K",
        "Ce_Lin_digitized_J_m-3_K-1",
        "G_Lin_digitized_W_m-3_K-1",
        "tau_Lin_Ce_over_G_ps",
        "Ce_XTANT_dynamic_J_m-3_K-1",
        "G_XTANT_dynamic_W_m-3_K-1",
        "XTANT_over_Lin_Ce",
        "XTANT_over_Lin_G",
        "tau_XTANT_Ce_over_G_ps",
        "XTANT_over_Lin_tau",
        "G_Daraszewicz_experiment_G0_W_m-3_K-1",
        "tau_Lin_Ce_over_Daraszewicz_experiment_G0_ps",
        "tau_Duffy_experiment_G0_uncertainty_low_ps",
        "tau_Duffy_experiment_G0_uncertainty_high_ps",
        "G_Daraszewicz_ab_initio_G0_W_m-3_K-1",
        "tau_Lin_Ce_over_Daraszewicz_ab_initio_G0_ps",
    ]
    output_rows: list[dict[str, float | int | str]] = []
    for te in temperatures:
        ce_lin = lin[te]["ce"]
        g_lin = lin[te]["g"]
        ce_x = xtant_ce.get(te)
        g_x = xtant_g.get(te)
        row: dict[str, float | int | str] = {
            "Te_K": te,
            "Ce_Lin_digitized_J_m-3_K-1": ce_lin,
            "G_Lin_digitized_W_m-3_K-1": g_lin,
            "tau_Lin_Ce_over_G_ps": tau_ps(ce_lin, g_lin),
            "Ce_XTANT_dynamic_J_m-3_K-1": ce_x if ce_x is not None else "",
            "G_XTANT_dynamic_W_m-3_K-1": g_x if g_x is not None else "",
            "XTANT_over_Lin_Ce": ce_x / ce_lin if ce_x is not None else "",
            "XTANT_over_Lin_G": g_x / g_lin if g_x is not None else "",
            "tau_XTANT_Ce_over_G_ps": tau_ps(ce_x, g_x) if ce_x is not None and g_x is not None else "",
            "XTANT_over_Lin_tau": (tau_ps(ce_x, g_x) / tau_ps(ce_lin, g_lin)) if ce_x is not None and g_x is not None else "",
        }
        if te in DUFFY_DISPLAY_TEMPERATURES_K:
            g_low = DUFFY_G_EXPERIMENT - DUFFY_G_EXPERIMENT_UNCERTAINTY
            g_high = DUFFY_G_EXPERIMENT + DUFFY_G_EXPERIMENT_UNCERTAINTY
            row.update(
                {
                    "G_Daraszewicz_experiment_G0_W_m-3_K-1": DUFFY_G_EXPERIMENT,
                    "tau_Lin_Ce_over_Daraszewicz_experiment_G0_ps": tau_ps(ce_lin, DUFFY_G_EXPERIMENT),
                    "tau_Duffy_experiment_G0_uncertainty_low_ps": tau_ps(ce_lin, g_high),
                    "tau_Duffy_experiment_G0_uncertainty_high_ps": tau_ps(ce_lin, g_low),
                    "G_Daraszewicz_ab_initio_G0_W_m-3_K-1": DUFFY_G_AB_INITIO,
                    "tau_Lin_Ce_over_Daraszewicz_ab_initio_G0_ps": tau_ps(ce_lin, DUFFY_G_AB_INITIO),
                }
            )
        else:
            row.update({field: "" for field in csv_fields[10:]})
        output_rows.append(row)

    with CSV_OUT.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=csv_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    sample_temperatures = (300, 500, 600, 1000, 3000, 5000, 11000, 15000, 19000, 20000)
    by_te = {int(row["Te_K"]): row for row in output_rows}
    comparisons = []
    for te in sample_temperatures:
        row = by_te[te]
        xtant_tau = row["tau_XTANT_Ce_over_G_ps"]
        lin_tau = row["tau_Lin_Ce_over_G_ps"]
        comparisons.append(
            {
                "Te_K": te,
                "tau_Lin_ps": lin_tau,
                "tau_XTANT_ps": xtant_tau if xtant_tau != "" else None,
                "XTANT_over_Lin_tau": xtant_tau / lin_tau if xtant_tau != "" else None,
                "tau_Duffy_experiment_G0_using_Lin_Ce_ps": row.get("tau_Lin_Ce_over_Daraszewicz_experiment_G0_ps") or None,
                "tau_Duffy_experiment_G0_uncertainty_ps": [
                    row.get("tau_Duffy_experiment_G0_uncertainty_low_ps"),
                    row.get("tau_Duffy_experiment_G0_uncertainty_high_ps"),
                ] if te in DUFFY_DISPLAY_TEMPERATURES_K else None,
                "tau_Duffy_ab_initio_G0_using_Lin_Ce_ps": row.get("tau_Lin_Ce_over_Daraszewicz_ab_initio_G0_ps") or None,
            }
        )

    lin_tau_values = [float(row["tau_Lin_Ce_over_G_ps"]) for row in output_rows]
    xtant_tau_values = [float(row["tau_XTANT_Ce_over_G_ps"]) for row in output_rows if row["tau_XTANT_Ce_over_G_ps"] != ""]
    common_rows = [row for row in output_rows if row["tau_XTANT_Ce_over_G_ps"] != ""]
    tau_ratios = [(float(row["XTANT_over_Lin_tau"]), int(row["Te_K"])) for row in common_rows]
    ce_ratios = [float(row["XTANT_over_Lin_Ce"]) for row in common_rows]
    g_ratios = [float(row["XTANT_over_Lin_G"]) for row in common_rows]
    report = {
        "title": "Source-separated tungsten electron-phonon time-scale comparison",
        "formula": "tau_e-ph = Ce(Te)/G(Te); for Ce in J m^-3 K^-1 and G in W m^-3 K^-1, tau is seconds; output values are converted to ps.",
        "interpretation_limit": "Ce/G is a characteristic electron energy-transfer time for a simple fixed-lattice estimate. A two-reservoir relaxation eigen-time also depends on lattice heat capacity; transport and nonlinear temperature evolution are omitted.",
        "branches": {
            "Lin_2008": "Both Ce and G are candidate curves digitized from the paper; G below 600 K is extrapolated and digitization uncertainty is about 2-5% before model/source uncertainty.",
            "XTANT3": "Ce and G are the supplied dynamical tables, retained as a separate branch.",
            "Daraszewicz_2014": "The paper obtains Ce using the Lin et al. method and reports a near-equilibrium experimental G0=1.4(3)e17 and ab-initio G0=1.65e17 W m^-3 K^-1. The displayed tau values combine Lin-digitized Ce with those G0 anchors only at selected comparison points; they are not a high-Te G(Te) model.",
        },
        "duffy_anchor_display_temperatures_K": list(DUFFY_DISPLAY_TEMPERATURES_K),
        "source_files": {
            "Lin_Ce_G": "ttm_data/W_TTM_Lin2008_digitized_300_20000K.csv",
            "XTANT_Ce": "ttm_data/XTANT3_TTM_parameters_v3/Ce_metals/Ce_W.txt",
            "XTANT_G": "ttm_data/XTANT3_TTM_parameters_v3/G_e-ph_metals/G_e-ph_W.txt",
            "Daraszewicz": "papers/Daraszewicz_2014_Tungsten_electron_phonon_coupling.pdf",
        },
        "temperature_ranges_K": {
            "Lin_digitized_Ce_G": [min(lin), max(lin)],
            "XTANT_dynamic_Ce_G": [max(min(xtant_ce), min(xtant_g)), min(max(xtant_ce), max(xtant_g))],
        },
        "tau_ranges_ps_on_native_grids": {
            "Lin_2008_digitized": [min(lin_tau_values), max(lin_tau_values)],
            "XTANT3_dynamic": [min(xtant_tau_values), max(xtant_tau_values)],
        },
        "pointwise_XTANT_over_Lin_ratios_on_common_grid": {
            "temperature_range_K": [min(int(row["Te_K"]) for row in common_rows), max(int(row["Te_K"]) for row in common_rows)],
            "Ce_ratio_min_max": [min(ce_ratios), max(ce_ratios)],
            "G_ratio_min_max": [min(g_ratios), max(g_ratios)],
            "tau_ratio_min_max": [min(tau_ratios)[0], max(tau_ratios)[0]],
            "tau_ratio_min_at_K": min(tau_ratios)[1],
            "tau_ratio_max_at_K": max(tau_ratios)[1],
            "ratio_definition": "(Ce_XTANT/G_XTANT)/(Ce_Lin/G_Lin); values are pointwise at matching Te, no interpolation",
        },
        "sample_comparisons": comparisons,
        "outputs": {
            "csv": CSV_OUT.name,
            "json": JSON_OUT.name,
            "markdown": MD_OUT.name,
        },
    }
    JSON_OUT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    table = [
        "| Te (K) | Lin τ (ps) | XTANT τ (ps) | XTANT/Lin τ | Duffy exp G₀ + Lin Ce τ (ps) | Duffy ab initio G₀ + Lin Ce τ (ps) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in comparisons:
        exp_tau = row["tau_Duffy_experiment_G0_using_Lin_Ce_ps"]
        exp_range = row["tau_Duffy_experiment_G0_uncertainty_ps"]
        exp_text = "—" if exp_tau is None else f"{exp_tau:.3g} ({exp_range[0]:.3g}–{exp_range[1]:.3g})"
        theory = row["tau_Duffy_ab_initio_G0_using_Lin_Ce_ps"]
        theory_text = "—" if theory is None else f"{theory:.3g}"
        ratio = "—" if row["XTANT_over_Lin_tau"] is None else f"{row['XTANT_over_Lin_tau']:.3g}×"
        table.append(
            f"| {row['Te_K']} | {row['tau_Lin_ps']:.3g} | {fmt(row['tau_XTANT_ps'])} | {ratio} | {exp_text} | {theory_text} |"
        )

    md = f"""# W electron–phonon time-scale source comparison

Date: 2026-11-24
Generated by `scripts/compare_w_eph_timescales.py`.

## Definition and source branches

The comparison metric is `tau_e-ph ≈ Ce(Te)/G(Te)`. Since `Ce` is in J m⁻³ K⁻¹ and `G` is in W m⁻³ K⁻¹, the ratio is seconds; the files report it in picoseconds. It is a simple electron energy-transfer time with the lattice treated as a fixed reservoir. A true two-reservoir relaxation time also depends on lattice heat capacity; diffusion, interfaces, and nonlinear coupled evolution are not represented by this ratio.

- **Lin et al. (2008):** both `Ce` and `G` are digitized candidate curves. They are kept from the same source; `G` below 600 K contains an explicit extrapolation, and graph reading carries an estimated 2–5% uncertainty.
- **XTANT-3:** the supplied standalone dynamical `Ce_W.txt` and `G_e-ph_W.txt` tables form a separate native branch.
- **Daraszewicz et al. (2014):** reports near-equilibrium `G₀=1.4(3)×10¹⁷ W m⁻³ K⁻¹` from experiment and `1.65×10¹⁷ W m⁻³ K⁻¹` from ab initio. Their methods cite Lin et al. for `Ce`; the comparison below therefore uses Lin-digitized `Ce` with each `G₀` anchor. These are low-temperature scale comparisons, not extensions of Duffy `G₀` into a high-`Te` curve. The experimental uncertainty range is propagated through `tau=Ce/G`.

## Selected temperatures

{chr(10).join(table)}

The Duffy anchor columns are shown only at 300, 500, 600, and 1,000 K to keep their near-equilibrium scope visible. A dash means that the branch has no value at that temperature or the anchor comparison is not tabulated there.

## What the comparison shows

- Across the supplied XTANT native grid (500–19,000 K), `Ce/G` ranges from {min(xtant_tau_values):.3g} to {max(xtant_tau_values):.3g} ps. The Lin digitized native grid (300–20,000 K) ranges from {min(lin_tau_values):.3g} to {max(lin_tau_values):.3g} ps. These extrema occur at different temperatures; they are not a single pointwise uncertainty interval.
- On their exact common 500–19,000 K grid, the pointwise XTANT/Lin `Ce/G` ratio ranges from {min(tau_ratios)[0]:.3g}× at {min(tau_ratios)[1]} K to {max(tau_ratios)[0]:.3g}× at {max(tau_ratios)[1]} K. The CSV also exposes the separate `Ce` and `G` ratios, so the time-scale difference is not mistaken for a coupling-only comparison.
- At 500 K, XTANT's time scale is {float(by_te[500]['tau_XTANT_Ce_over_G_ps']) / float(by_te[500]['tau_Lin_Ce_over_G_ps']):.2f} times Lin's; at 11,000 K it is {float(by_te[11000]['tau_XTANT_Ce_over_G_ps']) / float(by_te[11000]['tau_Lin_Ce_over_G_ps']):.2f} times Lin's. The two source branches therefore predict substantially different energy-transfer rates.
- Daraszewicz's experimental `G₀` anchor yields `tau=Ce_Lin/G₀` at 300 K of {by_te[300]['tau_Lin_Ce_over_Daraszewicz_experiment_G0_ps']:.3g} ps (reported-G uncertainty range {by_te[300]['tau_Duffy_experiment_G0_uncertainty_low_ps']:.3g}–{by_te[300]['tau_Duffy_experiment_G0_uncertainty_high_ps']:.3g} ps). The ab-initio anchor gives {by_te[300]['tau_Lin_Ce_over_Daraszewicz_ab_initio_G0_ps']:.3g} ps at the same `Ce`.

Do not average the source branches. Keep their names with each later TTM run and use them for sensitivity cases. No Ag/W trajectory is started by this analysis; the full-size QSGS geometry still fails the static-start audit.
"""
    MD_OUT.write_text(md, encoding="utf-8")
    print(f"Wrote {CSV_OUT.relative_to(PACK)}, {JSON_OUT.relative_to(PACK)}, and {MD_OUT.relative_to(PACK)}")


if __name__ == "__main__":
    main()
