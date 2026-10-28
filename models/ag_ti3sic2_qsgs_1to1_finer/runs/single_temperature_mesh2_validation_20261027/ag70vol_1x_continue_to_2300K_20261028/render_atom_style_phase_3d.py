#!/usr/bin/env python3
"""Render the saved Ag70 phase field as one sphere glyph per native voxel."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy import ndimage


BASE = Path(__file__).resolve().parent
RUN_DIR = BASE.parent
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
CHECKPOINT = BASE / "checkpoints" / "caseA_ag70_1x_heat_18p636ps.npz"
OUTPUT = BASE / "ag70_atom_style_phase_spheres_18p636ps.png"

COLORS = {
    "scaffold": "#263f51",
    "solid_ag": "#83b7d5",
    "partial_ag": "#ffd166",
    "liquid_ag": "#f28e2b",
}


def centres(mask: np.ndarray, spacing_nm: np.ndarray) -> np.ndarray:
    coords = np.argwhere(mask).astype(float) + 0.5
    return coords * spacing_nm[None, :]


def add_sample_outline(ax, size_nm: np.ndarray) -> None:
    corners = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    ], dtype=float) * size_nm
    edges = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
             (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7))
    for i, j in edges:
        ax.plot(*corners[[i, j]].T, color="#71808a", linewidth=0.7, alpha=0.55)


def main() -> None:
    if not GEOMETRY.exists() or not CHECKPOINT.exists():
        raise FileNotFoundError(f"Required geometry/checkpoint not found: {CHECKPOINT}")
    with np.load(GEOMETRY, allow_pickle=False) as data:
        scaffold = np.asarray(data["skeleton"], dtype=bool)
        size_nm = np.asarray(data["size_A"], dtype=float) / 10.0
    with np.load(CHECKPOINT, allow_pickle=False) as data:
        temperature = np.asarray(data["temperature_K"], dtype=float)
        liquid_fraction = np.asarray(data["ag_liquid_fraction"], dtype=float)
        metadata = json.loads(str(data["metadata_json"].item()))

    ag = ~scaffold
    spacing_nm = size_nm / np.asarray(scaffold.shape, dtype=float)
    liquid = ag & (liquid_fraction >= 0.5)
    partial = ag & (liquid_fraction > 1e-12) & (liquid_fraction < 0.5)
    solid_ag = ag & ~liquid & ~partial

    labels, _ = ndimage.label(liquid, ndimage.generate_binary_structure(3, 1))
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    pool_id = int(np.argmax(counts))
    pool = labels == pool_id
    pool_voxels = int(counts[pool_id])
    pool_volume_nm3 = pool_voxels * float(np.prod(spacing_nm))
    coords = np.where(pool)
    lo = np.array([index.min() for index in coords])
    hi = np.array([index.max() for index in coords])
    pool_bbox_nm = (hi - lo + 1) * spacing_nm

    # Remove one display octant to reveal interior spheres; no source cells or
    # phase values are changed. Every other voxel remains represented.
    visible = np.ones(scaffold.shape, dtype=bool)
    nx, ny, nz = scaffold.shape
    visible[nx // 2 :, ny // 2 :, nz // 2 :] = False

    fig = plt.figure(figsize=(15.0, 11.0))
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_proj_type("ortho")
    phase_specs = (
        (COLORS["scaffold"], "Ti₃SiC₂ scaffold"),
        (COLORS["solid_ag"], "unmelted Ag"),
        (COLORS["partial_ag"], "partially liquid Ag"),
        (COLORS["liquid_ag"], "liquid Ag"),
    )
    points = centres(visible, spacing_nm)
    phase_colors = np.array([matplotlib.colors.to_rgba(color) for color, _ in phase_specs])
    point_colors = phase_colors[np.where(scaffold, 0, np.where(partial, 2, np.where(liquid, 3, 1)))[visible]]
    ax.scatter(
        points[:, 0], points[:, 1], points[:, 2],
        s=48, marker="o", c=point_colors,
        depthshade=True, linewidths=0, rasterized=True,
    )

    add_sample_outline(ax, size_nm)
    ax.set_xlim(0, size_nm[0])
    ax.set_ylim(0, size_nm[1])
    ax.set_zlim(0, size_nm[2])
    ax.set_box_aspect(size_nm)
    ax.set_xlabel("x (nm)", labelpad=8)
    ax.set_ylabel("y (nm)", labelpad=8)
    ax.set_zlabel("z (nm)", labelpad=8)
    ax.view_init(elev=29, azim=-56)
    ax.grid(True, alpha=0.2)
    ax.legend(
        handles=[Line2D([0], [0], marker="o", linestyle="", markersize=9,
                        markerfacecolor=color, markeredgecolor="#26343d", label=label)
                 for color, label in phase_specs],
        loc="upper left", bbox_to_anchor=(0.02, 0.98), framealpha=0.96, fontsize=10,
    )
    time_ps = float(metadata["time_ps"])
    ax.text2D(
        0.02, 0.03,
        f"One sphere is drawn at every visible native-voxel centre ({np.count_nonzero(visible):,} of {visible.size:,}).\n"
        f"Largest connected liquid cluster: {pool_volume_nm3:.0f} nm³; "
        f"native-grid size {pool_bbox_nm[0]:.1f} × {pool_bbox_nm[1]:.1f} × {pool_bbox_nm[2]:.1f} nm.",
        transform=ax.transAxes, fontsize=9, ha="left", va="bottom",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "alpha": 0.88, "edgecolor": "#c8d0d6"},
    )
    fig.suptitle(
        f"Sphere-style 3D phase map | Ag–Ti₃SiC₂ 70:30 vol% | t={time_ps:.3f} ps\n"
        f"Peak temperature {temperature.max():.1f} K | liquid cluster={pool_voxels:,} native voxels",
        fontsize=17, y=0.96,
    )
    fig.text(
        0.5, 0.012,
        "Spheres are visualization glyphs at continuum voxel centres, not atomistic coordinates; "
        "the model does not resolve liquid flow or splashing.",
        ha="center", fontsize=9,
    )
    fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.06)
    fig.savefig(OUTPUT, dpi=200, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"saved={OUTPUT}")
    print(f"time_ps={time_ps:.6f} peak_K={temperature.max():.4f}")
    print(f"pool_voxels={pool_voxels} volume_nm3={pool_volume_nm3:.6f} bbox_nm={pool_bbox_nm.tolist()}")
    print(f"rendered_spheres={np.count_nonzero(visible)} of {visible.size}")


if __name__ == "__main__":
    main()
