#!/usr/bin/env python3
"""Digitize candidate W TTM curves from Lin, Zhigilei & Celli (2008), Fig. 7.

This is a plot digitizer, not a reproduction of the authors' electronic-structure
calculation. The output is exploratory and retains alternative G0 anchors.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "papers/Lin_Zhigilei_Celli_2008_Electron_phonon_coupling_and_heat_capacity.pdf"
OUT = ROOT / "ttm_data"
OUT.mkdir(exist_ok=True)

T_MIN_K = 300
T_MAX_K = 20_000
T_STEP_K = 100

# Pixel calibration for the 180-dpi render of PDF page 13 (printed page 13).
X0, X1 = 861.0, 1309.0
CE_Y0, CE_Y35 = 552.0, 176.0
G_Y2, G_Y12 = 1083.0, 710.0
PIXEL_THRESHOLD = 120

J_M3_TO_EV_A3 = 6.241509074e-12
W_M3_TO_EV_PS_A3 = 6.241509074e-24
W_M_TO_EV_PS_A = 6.241509074e-4


def render_figure() -> np.ndarray:
    with tempfile.TemporaryDirectory(prefix="lin2008_wfig_") as temp:
        prefix = Path(temp) / "lin"
        subprocess.run(
            ["pdftoppm", "-f", "13", "-l", "13", "-png", "-r", "180", str(PDF), str(prefix)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        page = next(Path(temp).glob("lin-*.png"))
        return np.asarray(Image.open(page).convert("L"))


def digitize_ce(image: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x_values, y_values = [], []
    # Ignore the first 15 px where the y-axis tick marks reach into the data panel.
    # The right plot border obscures the last two columns; stop at x=1307.
    for x in range(int(X0) + 15, int(X1) - 1):
        ys = np.flatnonzero(image[176:550, x] < PIXEL_THRESHOLD) + 176
        # The legend is printed inside the plot near its top. Excluding y<200
        # removes it while retaining the high-Te Ce curve (near y≈205 px).
        ys = ys[(ys > 200) & (ys < 550)]
        if len(ys):
            x_values.append(x)
            # Lin's DOS curve lies above the dashed Sommerfeld reference.
            y_values.append(float(ys.min()))
    return np.asarray(x_values), np.asarray(y_values)


def g_candidates(image: np.ndarray, x: int) -> np.ndarray:
    ys = np.flatnonzero(image[710:1082, x] < PIXEL_THRESHOLD) + 710
    return ys[(ys > 709) & (ys < 1081)]


def digitize_g(image: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    seed_x = 920
    seed_candidates = g_candidates(image, seed_x)
    if not len(seed_candidates):
        raise RuntimeError("Could not find the Fig. 7(d) G curve seed")
    seed_y = float(seed_candidates[np.argmin(abs(seed_candidates - 1055))])
    path: dict[int, float] = {seed_x: seed_y}

    for direction in (1, -1):
        prev_x = seed_x
        slope = -1.0 if direction == 1 else 1.0
        # The right plot border obscures the curve after x=1307.
        stop = int(X1) - 1 if direction == 1 else int(X0) + 4
        for x in range(seed_x + direction, stop, direction):
            ys = g_candidates(image, x)
            if not len(ys):
                path[x] = path[prev_x] + slope * direction
                prev_x = x
                continue
            prediction = path[prev_x] + slope * direction
            distances = abs(ys - prediction)
            if float(distances.min()) > 25.0:
                # A plotted reference line or annotation is not the smooth G curve.
                y = float(prediction)
                path[x] = y
                prev_x = x
                continue
            y = float(ys[np.argmin(distances)])
            local_slope = (y - path[prev_x]) / direction
            slope = max(-3.0, min(3.0, 0.75 * local_slope + 0.25 * slope))
            path[x] = y
            prev_x = x

    xs = np.asarray(sorted(path), dtype=float)
    ys = np.asarray([path[int(x)] for x in xs], dtype=float)
    return xs, ys


def main() -> None:
    image = render_figure()
    if image.shape[1] < 1310:
        raise RuntimeError(f"Unexpected rendered page dimensions {image.shape}")

    ce_x, ce_y = digitize_ce(image)
    g_x, g_y = digitize_g(image)
    temperature = np.arange(T_MIN_K, T_MAX_K + T_STEP_K, T_STEP_K, dtype=float)
    x_grid = X0 + (X1 - X0) * temperature / 20_000.0

    ce_y_grid = np.interp(x_grid, ce_x, ce_y)
    ce_j_m3_k = (CE_Y0 - ce_y_grid) / (CE_Y0 - CE_Y35) * 35.0e5
    # The plot's leftmost ticks obscure the curve; use the paper's stated
    # Sommerfeld coefficient for the low-Te linear region.
    ce_j_m3_k[temperature <= 600.0] = 137.3 * temperature[temperature <= 600.0]

    g_y_grid = np.interp(x_grid, g_x, g_y)
    g_w_m3_k = 2.0e17 + (G_Y2 - g_y_grid) / (G_Y2 - G_Y12) * 10.0e17
    # Below the visible graph's lower-temperature edge, extend the near-zero-Te
    # Lin curve from its plotted ~1.9e17 value to the first digitized pixels.
    low = temperature < 600.0
    g_600 = float(np.interp(X0 + (X1 - X0) * 600.0 / 20_000.0, g_x, g_y))
    g_600_val = 2.0e17 + (G_Y2 - g_600) / (G_Y2 - G_Y12) * 10.0e17
    g_w_m3_k[low] = np.interp(temperature[low], [0.0, 600.0], [1.9e17, g_600_val])

    # Duffy et al.'s Drude-style expression, shown at Te=Tl for this table.
    # In a TTM solve use the equation with independent Te and Tl, not this curve.
    ce_300 = float(np.interp(300.0, temperature, ce_j_m3_k))
    ke0_w_m_k = 0.75 * 174.0
    ke_equilibrium_w_m_k = ke0_w_m_k * (ce_j_m3_k / ce_300) * (300.0 / temperature) ** 1.5

    ce_ev_a3_k = ce_j_m3_k * J_M3_TO_EV_A3
    g_ev_ps_a3_k = g_w_m3_k * W_M3_TO_EV_PS_A3
    ke_ev_ps_a_k = ke_equilibrium_w_m_k * W_M_TO_EV_PS_A
    ue_j_m3 = np.concatenate(([0.0], np.cumsum(0.5 * (ce_j_m3_k[1:] + ce_j_m3_k[:-1]) * np.diff(temperature))))
    ue_ev_a3 = ue_j_m3 * J_M3_TO_EV_A3

    csv_path = OUT / "W_TTM_Lin2008_digitized_300_20000K.csv"
    header = (
        "Te_K,Ce_Lin2008_J_m-3_K-1,Ce_eV_A-3_K-1,"
        "G_Lin2008_digitized_W_m-3_K-1,G_eV_ps-1_A-3_K-1,"
        "G_Duffy2014_experiment_constant_W_m-3_K-1,G_Duffy2014_ab_initio_constant_W_m-3_K-1,"
        "ke_DuffyDrude_Te_eq_Tl_W_m-1_K-1,ke_eV_ps-1_A-1_K-1,"
        "Ue_delta_from_300K_eV_A-3"
    )
    values = np.column_stack(
        (
            temperature,
            ce_j_m3_k,
            ce_ev_a3_k,
            g_w_m3_k,
            g_ev_ps_a3_k,
            np.full_like(temperature, 1.4e17),
            np.full_like(temperature, 1.65e17),
            ke_equilibrium_w_m_k,
            ke_ev_ps_a_k,
            ue_ev_a3,
        )
    )
    np.savetxt(csv_path, values, delimiter=",", header=header, comments="", fmt="%.10g")

    report = {
        "title": "Candidate W TTM data digitized from Lin et al. (2008), Fig. 7(c,d)",
        "source_pdf": str(PDF.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(PDF.read_bytes()).hexdigest(),
        "source_citation": "Z. Lin, L. V. Zhigilei, and V. Celli, Phys. Rev. B 77, 075133 (2008), DOI 10.1103/PhysRevB.77.075133",
        "digitization": {
            "pdf_page": 13,
            "render_dpi": 180,
            "pixel_threshold_grayscale": PIXEL_THRESHOLD,
            "temperature_range_K": [T_MIN_K, T_MAX_K],
            "step_K": T_STEP_K,
            "estimated_plot_reading_uncertainty": "approximately 2-5%; curve digitization and source-model uncertainty are not included in this pixel estimate",
            "caveat": "Digitized curves are exploratory candidates, not original raw author data. G below 600 K is a stated low-Te extrapolation from the lower plotted range.",
        },
        "cross_checks": {
            "Ce_gamma_W_J_m-3_K-2": 137.3,
            "Ce_at_300K_J_m-3_K-1": float(ce_j_m3_k[0]),
            "Ce_at_11000K_over_gammaTe": float(np.interp(11000.0, temperature, ce_j_m3_k) / (137.3 * 11000.0)),
            "G_at_300K_Lin_digitized_W_m-3_K-1": float(g_w_m3_k[0]),
            "G_at_300K_Duffy_experiment_W_m-3_K-1": 1.4e17,
            "G_at_300K_Duffy_ab_initio_W_m-3_K-1": 1.65e17,
            "G_at_20000K_Lin_digitized_W_m-3_K-1": float(g_w_m3_k[-1]),
            "ke0_Duffy_W_m-1_K-1": ke0_w_m_k,
        },
        "ke_closure": "Duffy et al. (2014): ke(Te,Tl)=130.5*(Ce(Te)/Ce(300))*(300/Tl)^1.5 W m^-1 K^-1. The CSV ke column evaluates this only at Te=Tl; high-Te nonequilibrium validity is not independently verified.",
        "G_alternatives": "CSV includes the Lin 2008 digitized G curve plus separate constant Duffy 2014 experimental and ab initio anchors. These are alternate sensitivity cases, not values to average.",
        "output_csv": csv_path.name,
    }
    report_path = OUT / "W_TTM_Lin2008_digitization_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
