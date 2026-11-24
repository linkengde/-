#!/usr/bin/env python3
"""Digitize candidate Ag Ce(Te) and G(Te) curves from Lin et al. (2008), Fig. 2.

This extracts plotted data and explicitly labels the low-temperature G plateau
and short bridge as graph-based approximations; it does not recreate the authors'
electronic-structure calculation or claim original tabulated values.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image


PACK = Path(__file__).resolve().parents[1]
PDF = PACK / "papers/Lin_Zhigilei_Celli_2008_Electron_phonon_coupling_and_heat_capacity.pdf"
OUT = PACK / "ttm_data/AgW_TTM_baseline_v1"
OUT.mkdir(parents=True, exist_ok=True)
CSV_OUT = OUT / "Ag_TTM_Lin2008_digitized_300_20000K.csv"
JSON_OUT = OUT / "Ag_TTM_Lin2008_digitization_report.json"
PLOT_OUT = OUT / "Ag_TTM_Lin2008_digitization_preview.png"

PAGE = 7
RENDER_DPI = 180
T_MIN_K, T_MAX_K, T_STEP_K = 300, 20_000, 100

# Pixel coordinates calibrated to the panel axes on rendered PDF page 7.
X0, X1 = 852.0, 1300.0
CE_Y0, CE_Y30 = 557.0, 196.0  # Ce axis: 0 to 30 x 10^5 J m^-3 K^-1
G_Y0, G_Y3 = 1086.0, 677.0  # G axis: 0 to 3 x 10^17 W m^-3 K^-1
GRAY_THRESHOLD = 120
AG_GAMMA_EXP_J_M3_K2 = 63.3
G_LOW_PLATEAU_W_M3_K = 0.25e17
G_LOW_PLATEAU_END_K = 5000.0
G_CURVE_SEED_X = 1000
G_CURVE_SEED_Y = 1045.0


def render_page() -> np.ndarray:
    with tempfile.TemporaryDirectory(prefix="lin2008_ag_fig2_") as temp:
        prefix = Path(temp) / "lin"
        subprocess.run(
            ["pdftoppm", "-f", str(PAGE), "-l", str(PAGE), "-png", "-r", str(RENDER_DPI), str(PDF), str(prefix)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        page = next(Path(temp).glob("lin-*.png"))
        return np.asarray(Image.open(page).convert("L"))


def temperature_from_x(x: float) -> float:
    return (x - X0) / (X1 - X0) * T_MAX_K


def digitize_ce(image: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    seed_x, seed_y = 1000, 474.0
    path: dict[int, float] = {seed_x: seed_y}
    for direction in (1, -1):
        previous_x = seed_x
        slope = -1.4
        stop = int(X1) - 1 if direction == 1 else int(X0) + 15
        for x in range(seed_x + direction, stop, direction):
            candidates = np.flatnonzero(image[148:556, x] < GRAY_THRESHOLD) + 148
            candidates = candidates[(candidates > 148) & (candidates < 556)]
            prediction = path[previous_x] + slope * direction
            if not len(candidates):
                path[x] = prediction
                previous_x = x
                continue
            distances = np.abs(candidates - prediction)
            if float(distances.min()) > 32.0:
                # Skip legend strokes and stray marks rather than jumping curves.
                path[x] = prediction
                previous_x = x
                continue
            y = float(candidates[np.argmin(distances)])
            local_slope = (y - path[previous_x]) / direction
            slope = max(-3.5, min(0.3, 0.75 * local_slope + 0.25 * slope))
            path[x] = y
            previous_x = x
    xs = np.asarray(sorted(path), dtype=float)
    ys = np.asarray([path[int(x)] for x in xs], dtype=float)
    return xs, ys


def g_candidates(image: np.ndarray, x: int) -> np.ndarray:
    ys = np.flatnonzero(image[700:1080, x] < GRAY_THRESHOLD) + 700
    # Exclude the x/y axes and keep the plotted curve region.
    return ys[(ys > 700) & (ys < 1080)]


def digitize_g(image: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    path: dict[int, float] = {G_CURVE_SEED_X: G_CURVE_SEED_Y}
    for direction in (1, -1):
        previous_x = G_CURVE_SEED_X
        slope = -1.3 if direction == 1 else -1.3
        stop = int(X1) - 2 if direction == 1 else int(X0) + 18
        for x in range(G_CURVE_SEED_X + direction, stop, direction):
            candidates = g_candidates(image, x)
            prediction = path[previous_x] + slope * direction
            if not len(candidates):
                path[x] = prediction
                previous_x = x
                continue
            distances = np.abs(candidates - prediction)
            if float(distances.min()) > 28.0:
                path[x] = prediction
                previous_x = x
                continue
            y = float(candidates[np.argmin(distances)])
            local_slope = (y - path[previous_x]) / direction
            slope = max(-3.5, min(1.0, 0.75 * local_slope + 0.25 * slope))
            path[x] = y
            previous_x = x
    xs = np.asarray(sorted(path), dtype=float)
    ys = np.asarray([path[int(x)] for x in xs], dtype=float)
    return xs, ys


def main() -> None:
    image = render_page()
    if image.shape[1] < 1301:
        raise RuntimeError(f"Unexpected rendered PDF page shape {image.shape}")

    ce_x, ce_y = digitize_ce(image)
    g_x, g_y = digitize_g(image)
    temperatures = np.arange(T_MIN_K, T_MAX_K + T_STEP_K, T_STEP_K, dtype=float)
    x_grid = X0 + (X1 - X0) * temperatures / T_MAX_K
    ce_y_grid = np.interp(x_grid, ce_x, ce_y)
    ce = (CE_Y0 - ce_y_grid) / (CE_Y0 - CE_Y30) * 30.0e5
    # The paper reports the experimental low-Te coefficient gamma_exp=63.3.
    # The text says the plotted DOS result follows the low-Te linear form below
    # several thousand kelvin; use the reported gamma through 5000 K.
    ce[temperatures <= 5000.0] = AG_GAMMA_EXP_J_M3_K2 * temperatures[temperatures <= 5000.0]

    seed_temperature = temperature_from_x(G_CURVE_SEED_X)
    seed_g = (G_Y0 - G_CURVE_SEED_Y) / (G_Y0 - G_Y3) * 3.0e17
    g_x_temperature = np.asarray([temperature_from_x(x) for x in g_x])
    g_plot = (G_Y0 - g_y) / (G_Y0 - G_Y3) * 3.0e17
    g = np.interp(temperatures, g_x_temperature, g_plot)
    low = temperatures <= G_LOW_PLATEAU_END_K
    bridge = (temperatures > G_LOW_PLATEAU_END_K) & (temperatures < seed_temperature)
    g[low] = G_LOW_PLATEAU_W_M3_K
    g[bridge] = np.interp(temperatures[bridge], [G_LOW_PLATEAU_END_K, seed_temperature], [G_LOW_PLATEAU_W_M3_K, seed_g])

    values = np.column_stack((temperatures, ce, g))
    np.savetxt(
        CSV_OUT,
        values,
        delimiter=",",
        header="Te_K,Ce_Lin2008_digitized_J_m-3_K-1,G_Lin2008_digitized_W_m-3_K-1",
        comments="",
        fmt="%.10g",
    )

    report = {
        "title": "Candidate silver TTM Ce and G digitized from Lin et al. (2008), Fig. 2(c,d)",
        "source_pdf": str(PDF.relative_to(PACK)),
        "source_sha256": hashlib.sha256(PDF.read_bytes()).hexdigest(),
        "citation": "Z. Lin, L. V. Zhigilei, and V. Celli, Phys. Rev. B 77, 075133 (2008), DOI 10.1103/PhysRevB.77.075133",
        "digitization": {
            "pdf_page": PAGE,
            "render_dpi": RENDER_DPI,
            "temperature_range_K": [T_MIN_K, T_MAX_K],
            "step_K": T_STEP_K,
            "pixel_gray_threshold": GRAY_THRESHOLD,
            "calibration": {
                "temperature_axis_pixels": [X0, X1],
                "Ce_zero_and_30e5_axis_pixels": [CE_Y0, CE_Y30],
                "G_zero_and_3e17_axis_pixels": [G_Y0, G_Y3],
            },
            "estimated_graph_reading_uncertainty": "approximately 5% for plotted portions; digitization and model uncertainty are not included",
            "Ce_low_temperature_rule": "For Te<=5000 K, use the paper's quoted experimental gamma_exp=63.3 J m^-3 K^-2 multiplied by Te. The text says the plotted DOS result follows the low-Te linear form below several thousand kelvin; above 5000 K use the digitized Fig. 2(c) solid curve.",
            "G_low_temperature_rule": "Fig. 2(d) solid curve is approximately on a 0.25e17 W m^-3 K^-1 low-Te plateau. Hold this plotted plateau through 5000 K, then bridge linearly to the curve-tracing seed at 6607 K because the low-T trace overlaps experimental reference lines. Above the seed use traced Fig. 2(d).",
            "caveat": "These are extracted candidate curves, not original raw author tables. In particular the low-temperature G plateau and 5000–6607 K bridge are explicit graph-based approximations.",
        },
        "cross_checks": {
            "Ce_gamma_exp_J_m3_K2": AG_GAMMA_EXP_J_M3_K2,
            "Ce_at_300K_J_m3_K": float(ce[0]),
            "Ce_at_5000K_J_m3_K": float(np.interp(5000.0, temperatures, ce)),
            "Ce_at_10000K_J_m3_K": float(np.interp(10000.0, temperatures, ce)),
            "G_at_300K_W_m3_K": float(g[0]),
            "G_at_5000K_W_m3_K": float(np.interp(5000.0, temperatures, g)),
            "G_at_10000K_W_m3_K": float(np.interp(10000.0, temperatures, g)),
            "G_at_20000K_W_m3_K": float(g[-1]),
        },
        "output_csv": CSV_OUT.name,
        "preview_png": PLOT_OUT.name,
    }
    JSON_OUT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_preview(temperatures, ce, g)
    print(json.dumps(report["cross_checks"], indent=2))
    print(f"Wrote {CSV_OUT.relative_to(PACK)}, {JSON_OUT.relative_to(PACK)}, and {PLOT_OUT.relative_to(PACK)}")


def write_preview(temperatures: np.ndarray, ce: np.ndarray, g: np.ndarray) -> None:
    import os

    config = Path("/tmp/ag_w_mplconfig")
    cache = Path("/tmp/ag_w_xdg_cache")
    config.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(config))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.0), sharex=True, constrained_layout=True)
    axes[0].plot(temperatures, ce / 1e5, color="#1768ac")
    axes[0].set_ylabel(r"$C_e$ ($10^5$ J m$^{-3}$ K$^{-1}$)")
    axes[0].set_title("Ag Lin 2008 digitized TTM candidate")
    axes[1].plot(temperatures, g / 1e17, color="#b23a48")
    axes[1].set_ylabel(r"$G$ ($10^{17}$ W m$^{-3}$ K$^{-1}$)")
    axes[1].set_xlabel("Electron temperature (K)")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(PLOT_OUT, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
