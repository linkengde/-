"""Plot QSGS geometry and temperature/liquid slices from saved 3x3 checkpoints."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize


RUN_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUN_DIR.parents[1]
GEOMETRY = PROJECT_DIR / "qsgs_morphology.npz"
HEAT_CHECKPOINT = RUN_DIR / "heat/checkpoints/caseA_tile3_heat_0p6ps.npz"
COOL_CHECKPOINT = RUN_DIR / "cooling_from_0p6/checkpoints/caseA_tile3_cool_from_0p6_1p15ps.npz"
OUTPUT = RUN_DIR / "tile3_snapshot_heating_0p6_vs_cooling_1p15ps.png"


def load_checkpoint(path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    with np.load(path, allow_pickle=False) as archive:
        temperature = np.asarray(archive["temperature_K"])
        liquid = np.asarray(archive["ag_liquid_fraction"])
        metadata = json.loads(str(archive["metadata_json"].item()))
    return temperature, liquid, metadata


def masked_image(ax, values: np.ndarray, ag: np.ndarray, extent: tuple[float, float, float, float],
                 norm: Normalize, title: str, xlabel: str, ylabel: str):
    cmap = plt.get_cmap("inferno").copy()
    cmap.set_bad("#d8dde3")
    shown = np.ma.masked_where(~ag, values)
    im = ax.imshow(
        shown.T, origin="lower", extent=extent, aspect="equal", cmap=cmap,
        norm=norm, interpolation="nearest",
    )
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return im


def main() -> None:
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        skeleton_base = np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
    ag_base = ~skeleton_base
    ag = np.tile(ag_base, (3, 3, 1))
    spacing_A = size_A / np.asarray(ag_base.shape)
    temperature_hot, liquid_hot, hot_metadata = load_checkpoint(HEAT_CHECKPOINT)
    temperature_cool, liquid_cool, cool_metadata = load_checkpoint(COOL_CHECKPOINT)
    central = (slice(100, 200), slice(100, 200), slice(None))
    central_ag = ag[central]
    hot_roi = temperature_hot[central]
    liquid_roi = liquid_hot[central]
    cool_roi = temperature_cool[central]
    liquid_cool_roi = liquid_cool[central]

    # Native-cell coordinates for the central repeated geometry tile.
    lx_A, ly_A, lz_A = size_A
    source_local_x = float(hot_metadata["identity"]["source_center_xy_A"][0]) - lx_A
    source_local_y = float(hot_metadata["identity"]["source_center_xy_A"][1]) - ly_A
    norm = Normalize(vmin=300.0, vmax=2400.0)
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), constrained_layout=True)
    fig.suptitle(
        "Ag–Ti₃SiC₂ QSGS | 3×3 periodic domain | Case A",
        fontsize=16,
    )

    # Whole-domain top layer: central-cell box and the single heat-source center.
    im = masked_image(
        axes[0, 0], temperature_hot[:, :, -1], ag[:, :, -1],
        (0.0, 3.0 * lx_A / 10.0, 0.0, 3.0 * ly_A / 10.0), norm,
        f"Whole domain, top layer · heating t={hot_metadata['time_ps']:.2f} ps",
        "x (nm)", "y (nm)",
    )
    axes[0, 0].add_patch(plt.Rectangle(
        (lx_A / 10.0, ly_A / 10.0), lx_A / 10.0, ly_A / 10.0,
        linewidth=1.5, edgecolor="#41d6ff", facecolor="none", linestyle="--",
    ))
    axes[0, 0].plot(
        float(hot_metadata["identity"]["source_center_xy_A"][0]) / 10.0,
        float(hot_metadata["identity"]["source_center_xy_A"][1]) / 10.0,
        marker="+", color="cyan", markersize=10, markeredgewidth=2,
    )
    axes[0, 0].text(
        0.02, 0.97, "dashed box: central source cell", transform=axes[0, 0].transAxes,
        va="top", color="white", fontsize=9,
        bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none"},
    )

    hot_surface = hot_roi[:, :, -1]
    im = masked_image(
        axes[0, 1], hot_surface, central_ag[:, :, -1],
        (0.0, lx_A / 10.0, 0.0, ly_A / 10.0), norm,
        "Central cell, top layer · fₗ ≥ 0.5 contour",
        "local x (nm)", "local y (nm)",
    )
    axes[0, 1].plot(source_local_x / 10.0, source_local_y / 10.0,
                    marker="+", color="cyan", markersize=9, markeredgewidth=1.7)
    if np.any((liquid_roi[:, :, -1] >= 0.5) & central_ag[:, :, -1]):
        axes[0, 1].contour(
            np.linspace(0.5 * spacing_A[0] / 10.0, lx_A / 10.0 - 0.5 * spacing_A[0] / 10.0, 100),
            np.linspace(0.5 * spacing_A[1] / 10.0, ly_A / 10.0 - 0.5 * spacing_A[1] / 10.0, 100),
            ((liquid_roi[:, :, -1] >= 0.5) & central_ag[:, :, -1]).T.astype(float),
            levels=[0.5], colors=["#4dfcff"], linewidths=0.85,
        )
    axes[0, 1].text(
        0.02, 0.97, "cyan outline: liquid threshold", transform=axes[0, 1].transAxes,
        va="top", color="white", fontsize=9,
        bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none"},
    )

    # Cross-section through the center of the local Gaussian source.
    j_center = int(np.clip(round(source_local_y / spacing_A[1] - 0.5), 0, 99))
    im = masked_image(
        axes[1, 0], hot_roi[:, j_center, :], central_ag[:, j_center, :],
        (0.0, lx_A / 10.0, 0.0, lz_A / 10.0), norm,
        f"Central x–z section at y={source_local_y / 10.0:.2f} nm · heating",
        "local x (nm)", "z from bottom (nm)",
    )
    axes[1, 0].axvline(source_local_x / 10.0, color="#4dfcff", linewidth=0.7, linestyle="--")
    hot_cross_liquid = (liquid_roi[:, j_center, :] >= 0.5) & central_ag[:, j_center, :]
    if np.any(hot_cross_liquid):
        axes[1, 0].contour(
            np.linspace(0.5 * spacing_A[0] / 10.0, lx_A / 10.0 - 0.5 * spacing_A[0] / 10.0, 100),
            np.linspace(0.5 * spacing_A[2] / 10.0, lz_A / 10.0 - 0.5 * spacing_A[2] / 10.0, 50),
            hot_cross_liquid.T.astype(float), levels=[0.5], colors=["#4dfcff"], linewidths=0.85,
        )

    im = masked_image(
        axes[1, 1], cool_roi[:, :, -1], central_ag[:, :, -1],
        (0.0, lx_A / 10.0, 0.0, ly_A / 10.0), norm,
        f"Central cell, top layer · after heat-off t={cool_metadata['time_ps']:.2f} ps",
        "local x (nm)", "local y (nm)",
    )
    liquid_count = int(np.count_nonzero((liquid_cool_roi >= 0.5) & central_ag))
    axes[1, 1].text(
        0.02, 0.97, f"fₗ ≥ 0.5 Ag voxels: {liquid_count}", transform=axes[1, 1].transAxes,
        va="top", color="white", fontsize=9,
        bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none"},
    )

    colorbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.86, pad=0.025)
    colorbar.set_label("Ag temperature (K)")
    fig.text(
        0.5, 0.005,
        "Gray: Ti₃SiC₂ insulating mask · blue outline: fₗ ≥ 0.5 · Case A has zero scaffold heat capacity/conductivity",
        ha="center", fontsize=9,
    )
    fig.savefig(OUTPUT, dpi=180, bbox_inches="tight")
    print(OUTPUT)


if __name__ == "__main__":
    main()
