#!/usr/bin/env python3
"""Render the 18.636 ps phase map as a 3D voxel-context figure."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy import ndimage


BASE = Path(__file__).resolve().parent
RUN_DIR = BASE.parent
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
CHECKPOINT = BASE / "checkpoints" / "caseA_ag70_1x_heat_18p636ps.npz"
OUTPUT = BASE / "ag70_phase_map_3d_18p636ps.png"

COLORS = {
    "scaffold": "#263f51",
    "solid_ag": "#83b7d5",
    "partial_ag": "#ffd166",
    "liquid_ag": "#f28e2b",
}


def exposed_faces(mask: np.ndarray, spacing_nm: np.ndarray) -> tuple[list[np.ndarray], list[str]]:
    """Return the exposed voxel faces for a mask, including material interfaces."""
    polygons: list[np.ndarray] = []
    colors: list[str] = []
    axis_colors = {
        (0, -1): "#db741d", (0, 1): "#ef8628",
        (1, -1): "#e47c20", (1, 1): "#f39232",
        (2, -1): "#e98224", (2, 1): "#ffab4c",
    }
    for axis in range(3):
        other = [d for d in range(3) if d != axis]
        for side in (-1, 1):
            neighbor = np.roll(mask, -side, axis=axis)
            edge = [slice(None)] * 3
            edge[axis] = -1 if side == 1 else 0
            neighbor[tuple(edge)] = False
            cells = np.argwhere(mask & ~neighbor)
            if not len(cells):
                continue
            quads = np.repeat(cells[:, None, :].astype(float), 4, axis=1)
            quads[:, :, axis] += 1.0 if side == 1 else 0.0
            for j, bits in enumerate(((0, 0), (1, 0), (1, 1), (0, 1))):
                quads[:, j, other[0]] += bits[0]
                quads[:, j, other[1]] += bits[1]
            quads *= spacing_nm[None, None, :]
            polygons.extend(quads)
            colors.extend([axis_colors[(axis, side)]] * len(quads))
    return polygons, colors


def draw_pool(ax, pool: np.ndarray, spacing_nm: np.ndarray) -> None:
    faces, shades = exposed_faces(pool, spacing_nm)
    ax.add_collection3d(Poly3DCollection(
        faces,
        facecolors=shades,
        edgecolors=(0.42, 0.20, 0.07, 0.28),
        linewidths=0.035,
        alpha=0.48,
        zsort="average",
        rasterized=True,
    ))


def sample_centres(mask: np.ndarray, spacing_nm: np.ndarray, maximum: int, seed: int) -> np.ndarray:
    flat = np.flatnonzero(mask)
    if not len(flat):
        return np.empty((0, 3), dtype=float)
    rng = np.random.default_rng(seed)
    if len(flat) > maximum:
        flat = rng.choice(flat, size=maximum, replace=False)
    coords = np.column_stack(np.unravel_index(flat, mask.shape)).astype(float) + 0.5
    return coords * spacing_nm[None, :]


def add_sample_outline(ax, size_nm: np.ndarray) -> None:
    corners = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    ], dtype=float) * size_nm
    edges = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
             (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7))
    for i, j in edges:
        ax.plot(*corners[[i, j]].T, color="#536574", linewidth=0.75, alpha=0.78)


def add_slice_plane(ax, phase: np.ndarray, z_index: int, spacing_nm: np.ndarray) -> None:
    """Add the exact full-resolution XY phase slice as a plane inside the 3D domain."""
    nx, ny = phase.shape[:2]
    z = (z_index + 0.5) * spacing_nm[2] + 0.012  # tiny offset avoids z-fighting only
    x, y = np.meshgrid(np.arange(nx + 1) * spacing_nm[0],
                       np.arange(ny + 1) * spacing_nm[1], indexing="ij")
    quads = np.empty((nx * ny, 4, 3), dtype=float)
    quads[:, 0, 0] = x[:-1, :-1].ravel()
    quads[:, 0, 1] = y[:-1, :-1].ravel()
    quads[:, 1, 0] = x[1:, :-1].ravel()
    quads[:, 1, 1] = y[1:, :-1].ravel()
    quads[:, 2, 0] = x[1:, 1:].ravel()
    quads[:, 2, 1] = y[1:, 1:].ravel()
    quads[:, 3, 0] = x[:-1, 1:].ravel()
    quads[:, 3, 1] = y[:-1, 1:].ravel()
    quads[:, :, 2] = z
    palette = np.array([
        matplotlib.colors.to_rgba(COLORS["scaffold"], 0.98),
        matplotlib.colors.to_rgba(COLORS["solid_ag"], 0.98),
        matplotlib.colors.to_rgba(COLORS["partial_ag"], 0.98),
        matplotlib.colors.to_rgba(COLORS["liquid_ag"], 0.98),
    ])
    facecolors = palette[phase[:, :, z_index].ravel()]
    ax.add_collection3d(Poly3DCollection(
        quads, facecolors=facecolors, edgecolors="none", zsort="average", rasterized=True
    ))


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
    pool_volume = pool_voxels * float(np.prod(spacing_nm))
    coords = np.where(pool)
    lo = np.array([c.min() for c in coords])
    hi = np.array([c.max() for c in coords])
    bbox = (hi - lo + 1) * spacing_nm
    time_ps = float(metadata["time_ps"])

    # Discrete native-grid phase values, with the same four labels as the source slice.
    phase = np.zeros(scaffold.shape, dtype=np.uint8)
    phase[solid_ag] = 1
    phase[partial] = 2
    phase[liquid] = 3

    fig = plt.figure(figsize=(15.0, 11.0))
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_proj_type("ortho")

    # Sparse centres give whole-volume context; the orange pool boundary is full resolution.
    for mask, color, label, maximum, seed, size, alpha in (
        (scaffold, COLORS["scaffold"], "Ti₃SiC₂ scaffold", 10500, 1, 5.5, 0.38),
        (solid_ag, COLORS["solid_ag"], "unmelted Ag", 12500, 2, 5.0, 0.34),
        (partial, COLORS["partial_ag"], "partially liquid Ag", 2700, 3, 8.0, 0.72),
        (pool, COLORS["liquid_ag"], "liquid Ag cluster", 4500, 4, 4.7, 0.44),
    ):
        points = sample_centres(mask, spacing_nm, maximum, seed)
        if len(points):
            ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=size, c=color,
                       alpha=alpha, depthshade=False, linewidths=0, rasterized=True,
                       label=label)
    draw_pool(ax, pool, spacing_nm)

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
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98), framealpha=0.96, fontsize=10)
    ax.text2D(
        0.02, 0.03,
        f"Full-resolution connected-pool surface: {pool_volume:.0f} nm³; "
        f"bounding size {bbox[0]:.1f} × {bbox[1]:.1f} × {bbox[2]:.1f} nm.\n"
        "Colored background phase centres are sampled for legibility; pool geometry uses native voxels.",
        transform=ax.transAxes, fontsize=9, ha="left", va="bottom",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "alpha": 0.88, "edgecolor": "#c8d0d6"},
    )
    fig.suptitle(
        f"3D Ag melt-pool shape with phase context | Ag–Ti₃SiC₂ 70:30 vol% | t={time_ps:.3f} ps\n"
        f"Peak temperature {temperature.max():.1f} K | pool={pool_voxels:,} native voxels",
        fontsize=17, y=0.96,
    )
    fig.text(
        0.5, 0.012,
        "Thermal enthalpy visualization only: it does not model liquid flow, free-surface motion, or splashing.",
        ha="center", fontsize=9,
    )
    fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.06)
    fig.savefig(OUTPUT, dpi=230, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"saved={OUTPUT}")
    print(f"t_ps={time_ps:.6f} peak_K={temperature.max():.4f} pool_voxels={pool_voxels}")
    print(f"pool_volume_nm3={pool_volume:.6f} bbox_nm={bbox.tolist()}")


if __name__ == "__main__":
    main()
