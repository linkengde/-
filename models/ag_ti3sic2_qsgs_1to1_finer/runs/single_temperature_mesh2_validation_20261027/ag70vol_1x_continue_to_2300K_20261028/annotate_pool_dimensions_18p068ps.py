#!/usr/bin/env python3
"""Re-render the saved 18.068 ps Case A pool with explicit 3D dimension arrows."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

import render_current_3d as render


BASE = Path(__file__).resolve().parent
GEOMETRY = BASE.parent / "qsgs_morphology_ag70vol.npz"
CHECKPOINT = BASE / "checkpoints" / "caseA_ag70_1x_heat_18p068ps.npz"
OUTPUT = BASE / "ag70_pool_shape_3d_18p068ps_dimensions.png"


def draw_box(ax, bounds: tuple[tuple[float, float], ...], color: str, linestyle: str = "--") -> None:
    (x0, x1), (y0, y1), (z0, z1) = bounds
    corners = np.array([
        [x0,y0,z0], [x1,y0,z0], [x1,y1,z0], [x0,y1,z0],
        [x0,y0,z1], [x1,y0,z1], [x1,y1,z1], [x0,y1,z1],
    ])
    edges = ((0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7))
    for i, j in edges:
        ax.plot(*corners[[i,j]].T, color=color, linewidth=1.2, linestyle=linestyle, alpha=0.9)


def add_dimension(ax, start: np.ndarray, end: np.ndarray, label: str, color: str) -> None:
    vector = end - start
    length = float(np.linalg.norm(vector))
    unit = vector / length
    head = min(0.8, length * 0.06)
    ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], color=color, linewidth=2.0)
    for point, direction in ((start, unit), (end, -unit)):
        ax.quiver(
            point[0], point[1], point[2],
            direction[0] * head, direction[1] * head, direction[2] * head,
            color=color, linewidth=1.8, arrow_length_ratio=0.85, normalize=False,
        )
    midpoint = (start + end) / 2.0
    ax.text(*midpoint, label, color=color, fontsize=10, fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": color, "alpha": 0.92, "pad": 2.5})


def main() -> None:
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        scaffold = np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
    with np.load(CHECKPOINT, allow_pickle=False) as archive:
        liquid_fraction = np.asarray(archive["ag_liquid_fraction"], dtype=float)
        metadata = json.loads(str(archive["metadata_json"].item()))

    ag = ~scaffold
    spacing_nm = (size_A / np.asarray(ag.shape, dtype=float)) / 10.0
    pool = ag & (liquid_fraction >= 0.5)
    labels, count = ndimage.label(pool, ndimage.generate_binary_structure(3, 1))
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    largest_id = int(np.argmax(counts))
    pool = labels == largest_id
    coords = np.where(pool)
    lo = np.asarray([axis.min() for axis in coords])
    hi = np.asarray([axis.max() for axis in coords])
    bounds = np.column_stack((lo * spacing_nm, (hi + 1) * spacing_nm))
    bbox_nm = (hi - lo + 1) * spacing_nm
    volume_nm3 = float(counts[largest_id]) * float(np.prod(spacing_nm))
    equivalent_diameter_nm = (6.0 * volume_nm3 / math.pi) ** (1.0 / 3.0)

    roi_lo = np.maximum(lo - 4, 0)
    roi_hi = np.minimum(hi + 5, np.asarray(pool.shape))
    roi = tuple(slice(int(a), int(b)) for a, b in zip(roi_lo, roi_hi, strict=True))
    fig = plt.figure(figsize=(12.5, 10.0))
    ax = fig.add_subplot(111, projection="3d")
    render.add_pool_surface(ax, pool[roi], spacing_nm, roi_lo, alpha=0.52)
    draw_box(ax, tuple((float(a), float(b)) for a, b in bounds), "#334e62")

    x0, x1 = bounds[0]
    y0, y1 = bounds[1]
    z0, z1 = bounds[2]
    pad = 1.5
    add_dimension(ax, np.array([x0, y0 - pad, z1 + pad]), np.array([x1, y0 - pad, z1 + pad]),
                  f"X ≈ {bbox_nm[0]:.2f} nm", "#b22222")
    add_dimension(ax, np.array([x1 + pad, y0, z1 + pad]), np.array([x1 + pad, y1, z1 + pad]),
                  f"Y ≈ {bbox_nm[1]:.2f} nm", "#176b3a")
    add_dimension(ax, np.array([x1 + pad, y0 - pad, z0]), np.array([x1 + pad, y0 - pad, z1]),
                  f"Z ≈ {bbox_nm[2]:.2f} nm", "#4b3f91")

    limits = (
        (float(x0 - 2.5), float(x1 + 5.0)),
        (float(y0 - 5.0), float(y1 + 3.0)),
        (float(z0 - 2.0), float(z1 + 5.0)),
    )
    render.set_nm_axes(ax, limits, bbox_nm + np.array([2.0, 2.0, 2.0]))
    ax.view_init(elev=27, azim=-54)
    time_ps = float(metadata["time_ps"])
    ax.set_title(
        f"Ag70 Case A · 1× · t={time_ps:.3f} ps · largest connected liquid-rich cluster (fₗ ≥ 0.5)",
        pad=18,
        fontsize=14,
    )
    fig.text(
        0.5, 0.905,
        f"Measured outer bounding box: X × Y × Z = {bbox_nm[0]:.2f} × {bbox_nm[1]:.2f} × {bbox_nm[2]:.2f} nm",
        ha="center", va="center", fontsize=12, fontweight="bold",
        bbox={"facecolor": "white", "edgecolor": "#334e62", "alpha": 0.96, "pad": 5},
    )
    fig.text(
        0.5, 0.035,
        f"Voxel cluster volume = {volume_nm3:.1f} nm³  |  equivalent-sphere diameter D_eq = {equivalent_diameter_nm:.2f} nm "
        "(not X or Y width)\n"
        "The irregular, percolating cluster does not fill a compact cuboid.",
        ha="center", va="bottom", fontsize=11,
        bbox={"facecolor": "white", "edgecolor": "#888888", "alpha": 0.95, "pad": 6},
    )
    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.15)
    fig.savefig(OUTPUT, dpi=220, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"saved={OUTPUT}")
    print(f"time_ps={time_ps:.6f} volume_nm3={volume_nm3:.6f} bbox_nm={bbox_nm.tolist()} equivalent_sphere_diameter_nm={equivalent_diameter_nm:.6f} components={count}")


if __name__ == "__main__":
    main()
