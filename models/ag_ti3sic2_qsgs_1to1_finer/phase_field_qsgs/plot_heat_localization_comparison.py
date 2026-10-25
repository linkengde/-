"""Plot 0.1 ps surface-plane temperature rises for localization pilots."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qsgs-matplotlib-cache")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
import numpy as np

HERE = Path(__file__).resolve().parent

CASES = (
    ("base_volume_fraction", "10 Å · original box · volume-fraction links"),
    ("base_exact_contacts", "10 Å · original box · QSGS face contacts"),
    ("3x3_volume_fraction", "10 Å · 3×3 repeat box · volume-fraction links"),
    ("3x3_exact_contacts", "10 Å · 3×3 repeat box · QSGS face contacts"),
    ("5A_exact_contacts", "5 Å · original box · QSGS face contacts"),
)


def _panel(ax, path: Path, field_name: str, cmap: str, vmin: float, vmax: float, title: str) -> None:
    with np.load(path, allow_pickle=False) as data:
        temp = np.asarray(data[f"t_0.1ps_{field_name}_temperature_K"], dtype=float)
        ag = np.asarray(data["ag_fraction"], dtype=float)
        spacing = np.asarray(data["spacing_A"], dtype=float)
        size = np.asarray(data["size_A"], dtype=float)
    rise = np.maximum(temp[:, :, -1] - 300.0, 0.0).T
    mask = (ag[:, :, -1].T < 1.0e-12) | (rise < vmin)
    image = ax.imshow(
        np.ma.array(rise, mask=mask),
        origin="lower",
        extent=(0.0, size[0] / 10.0, 0.0, size[1] / 10.0),
        interpolation="nearest",
        cmap=cmap,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        aspect="equal",
    )
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    if size[0] > 400.0:
        tile_x, tile_y = size[0] / 30.0, size[1] / 30.0
        cx, cy = size[0] / 20.0, size[1] / 20.0
        ax.add_patch(
            Rectangle(
                (cx - tile_x / 2.0, cy - tile_y / 2.0),
                tile_x,
                tile_y,
                fill=False,
                edgecolor="white",
                linewidth=0.8,
                linestyle="--",
            )
        )
    return image


def main() -> None:
    fig, axes = plt.subplots(2, 5, figsize=(19.0, 7.0), constrained_layout=True)
    electron_image = lattice_image = None
    for column, (stem, title) in enumerate(CASES):
        path = HERE / f"qsgs_localization_{stem}_fields.npz"
        electron_image = _panel(
            axes[0, column], path, "electron", "magma", 1.0, 6000.0, title
        )
        lattice_image = _panel(
            axes[1, column], path, "lattice", "viridis", 0.01, 5.0, title
        )
    axes[0, 0].set_ylabel("Electron ΔT (K)\ny (nm)")
    axes[1, 0].set_ylabel("Lattice ΔT (K)\ny (nm)")
    fig.suptitle(
        "QSGS local-heating screen at 0.1 ps (top z plane; dashed box = central 30 nm tile)",
        fontsize=13,
    )
    fig.colorbar(electron_image, ax=axes[0, :], label="Electron temperature rise (K)", shrink=0.82)
    fig.colorbar(lattice_image, ax=axes[1, :], label="Lattice temperature rise (K)", shrink=0.82)
    output = HERE / "qsgs_heat_localization_screen_0p1ps.png"
    fig.savefig(output, dpi=220)
    print(output)


if __name__ == "__main__":
    main()
