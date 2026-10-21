#!/usr/bin/env python3
"""Render a visual comparison of the bicontinuous reference and stopped slab run."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from matplotlib.colors import ListedColormap
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
AG = "#D7A83B"
TSC = "#168B91"


def colors_for(atoms):
    return np.where(np.asarray(atoms.arrays["phase"]).astype(str) == "Ag", AG, TSC)


def main():
    bicont_path = ROOT / "exploratory_melt_cool_20261018/input/NEW_MODEL_OVITO.extxyz"
    stopped_path = ROOT / "mace_v17_auto_melt_cool_20261019_continuation7/HEATING_CHECKPOINT.extxyz"
    bicont = read(bicont_path, format="extxyz")
    stopped = read(stopped_path, format="extxyz")

    fig = plt.figure(figsize=(16, 8.4), facecolor="white", layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=(1.25, 1.0))

    # Full 3D scatter of the geometry-only 39,486-atom bicontinuous reference.
    ax3d = fig.add_subplot(gs[0, 0], projection="3d")
    p = bicont.positions
    ph = np.asarray(bicont.arrays["phase"]).astype(str)
    # The open-x vacuum is omitted from the view so the material morphology is legible.
    keep = p[:, 0] >= 50.0
    mask_ag = keep & (ph == "Ag")
    mask_tsc = keep & (ph != "Ag")
    for mask, color in ((mask_tsc, TSC), (mask_ag, AG)):
        ax3d.scatter(p[mask, 0], p[mask, 1], p[mask, 2], s=1.0, c=color,
                     alpha=0.52, linewidths=0, depthshade=False, rasterized=True)
    ax3d.set_xlim(50, 110)
    ax3d.set_ylim(0, 98.1)
    ax3d.set_zlim(0, 100.9)
    ax3d.set_box_aspect((60, 98.1, 100.9))
    ax3d.view_init(elev=22, azim=-55)
    ax3d.set_xlabel("X (Å; open surface)", labelpad=8)
    ax3d.set_ylabel("Y (Å; periodic)", labelpad=8)
    ax3d.set_zlabel("Z (Å; periodic)", labelpad=8)
    ax3d.set_title("Bicontinuous Scheme A reference\n39,486 atoms; geometry-only source", pad=14, weight="bold")

    # Side projection of the actual 1,084-atom structure used by continuation 7.
    ax2d = fig.add_subplot(gs[0, 1])
    p = stopped.positions
    ph = np.asarray(stopped.arrays["phase"]).astype(str)
    for mask, color, label in ((ph != "Ag", TSC, "Ti₃SiC₂"), (ph == "Ag", AG, "Ag")):
        ax2d.scatter(p[mask, 0], p[mask, 2], s=15, c=color, alpha=0.78,
                     linewidths=0, rasterized=True, label=label)
    ax2d.set_aspect("equal", adjustable="datalim")
    ax2d.set_xlabel("X (Å; in-plane)")
    ax2d.set_ylabel("Z (Å; surface-normal)")
    ax2d.set_title("Structure actually heated, now stopped\n1,084 atoms; layered Ag/Ti₃SiC₂", pad=14, weight="bold")
    ax2d.legend(loc="upper right", frameon=True, markerscale=1.2)
    ax2d.grid(alpha=0.18)

    fig.suptitle("Ag–Ti₃SiC₂ structures: required bicontinuous geometry vs. stopped heating input",
                 fontsize=16, weight="bold")
    fig.text(0.5, -0.008,
             "Left: true two-phase connected networks (Ag gold, Ti₃SiC₂ teal); displayed material starts at X=50 Å, vacuum omitted.  "
             "Right: last saved frame of the separate layered run; it is not a bicontinuous model.",
             ha="center", va="top", fontsize=9.5)
    fig.savefig(OUT / "MODEL_GEOMETRY_COMPARISON.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # A comparable phase-map rendering: periodic YZ slices through the actual
    # 39,486-atom file, assigning each pixel the phase of its nearest atom in a
    # narrow X slab. This reveals the periodic morphology more clearly than a
    # single perspective projection of every atom.
    p = bicont.positions
    phase = np.asarray(bicont.arrays["phase"]).astype(str)
    ly, lz = np.asarray(bicont.cell.lengths())[1:3]
    ny = nz = 300
    y = np.linspace(0.0, ly, ny, endpoint=False)
    z = np.linspace(0.0, lz, nz, endpoint=False)
    yy, zz = np.meshgrid(y, z, indexing="xy")
    query = np.column_stack((yy.ravel(), zz.ravel()))
    cmap = ListedColormap(("#d5d7dc", "#2864d2"))
    fig, axes = plt.subplots(3, 3, figsize=(10, 10), dpi=180)
    for ax, x0 in zip(axes.ravel(), (60, 65, 70, 75, 80, 85, 90, 95, 100)):
        selected = np.abs(p[:, 0] - x0) <= 1.2
        yz = p[selected, 1:3].copy()
        yz[:, 0] %= ly
        yz[:, 1] %= lz
        is_ag = (phase[selected] == "Ag").astype(np.uint8)
        tree = cKDTree(yz, boxsize=(ly, lz))
        _, nearest = tree.query(query, k=1)
        phase_map = is_ag[nearest].reshape(nz, ny)
        ax.imshow(phase_map, origin="lower", extent=(0, ly, 0, lz),
                  interpolation="nearest", cmap=cmap, vmin=0, vmax=1,
                  aspect="equal")
        ax.set_title(f"X = {x0} Å")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("39,486-atom model: periodic YZ phase slices\nBlue = Ag; gray = Ti₃SiC₂",
                 fontsize=14, weight="bold")
    fig.text(0.5, 0.015,
             "Each panel samples atoms within X = listed value ±1.2 Å; pixels take the nearest sampled atom's phase. "
             "This shows the actual model's periodic morphology, not a continuous reconstructed phase field.",
             ha="center", va="bottom", fontsize=8.5)
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(OUT / "BICONNECTIVE_YZ_SLICES.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
