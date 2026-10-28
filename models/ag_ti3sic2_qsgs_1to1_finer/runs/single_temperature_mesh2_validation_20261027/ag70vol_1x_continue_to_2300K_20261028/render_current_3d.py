#!/usr/bin/env python3
"""Render the saved Ag70 Case A pool as voxel-boundary surfaces in 3D."""

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
OUTPUT_COMPOSITE = BASE / "ag70_full_model_and_pool_shape_3d_final.png"
OUTPUT_POOL = BASE / "ag70_pool_shape_3d_final.png"
OUTPUT_SLICE = BASE / "ag70_phase_slice_through_pool_final.png"


def boundary_quads(mask: np.ndarray, spacing_nm: np.ndarray, origin_index: np.ndarray) -> tuple[list[np.ndarray], list[str]]:
    """Build only exposed voxel faces; keeps a large pool fast to render."""
    polygons: list[np.ndarray] = []
    colors: list[str] = []
    axis_colors = {
        (0, -1): "#c7600e", (0, 1): "#e47a1c",
        (1, -1): "#d96e12", (1, 1): "#f28e2b",
        (2, -1): "#e98724", (2, 1): "#ffae45",
    }
    for axis in range(3):
        other = [d for d in range(3) if d != axis]
        for side in (-1, 1):
            neighbor = np.roll(mask, -side, axis=axis)
            edge = [slice(None)] * 3
            edge[axis] = -1 if side == 1 else 0
            neighbor[tuple(edge)] = False
            cells = np.argwhere(mask & ~neighbor)
            if cells.size == 0:
                continue
            base = cells + origin_index
            quads = np.repeat(base[:, None, :].astype(float), 4, axis=1)
            quads[:, :, axis] += 1.0 if side == 1 else 0.0
            for j, bits in enumerate(((0, 0), (1, 0), (1, 1), (0, 1))):
                quads[:, j, other[0]] += bits[0]
                quads[:, j, other[1]] += bits[1]
            quads *= spacing_nm[None, None, :]
            polygons.extend(quads)
            colors.extend([axis_colors[(axis, side)]] * len(quads))
    return polygons, colors


def add_pool_surface(
    ax, mask: np.ndarray, spacing_nm: np.ndarray, origin_index: np.ndarray, alpha: float = 0.98
) -> None:
    polygons, colors = boundary_quads(mask, spacing_nm, origin_index)
    surface = Poly3DCollection(
        polygons,
        facecolors=colors,
        edgecolors=(0.34, 0.16, 0.06, 0.42),
        linewidths=0.10,
        alpha=alpha,
        zsort="average",
        rasterized=True,
    )
    ax.add_collection3d(surface)


def add_sample_box(ax, size_nm: np.ndarray, color: str = "#596b7a", width: float = 0.8) -> None:
    x0, y0, z0 = 0.0, 0.0, 0.0
    x1, y1, z1 = size_nm
    corners = np.array([
        [x0,y0,z0], [x1,y0,z0], [x1,y1,z0], [x0,y1,z0],
        [x0,y0,z1], [x1,y0,z1], [x1,y1,z1], [x0,y1,z1],
    ])
    edges = ((0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7))
    for i, j in edges:
        ax.plot(*corners[[i, j]].T, color=color, linewidth=width)


def cloud_points(ax, mask: np.ndarray, spacing_nm: np.ndarray, color: str, stride: int, size: float, alpha: float) -> None:
    sampled = mask[::stride, ::stride, ::stride]
    points = np.argwhere(sampled)
    if points.size == 0:
        return
    points = points.astype(float) * stride + stride / 2.0
    shape = np.array(mask.shape, dtype=float)
    points = np.minimum(points, shape[None, :] - 0.5)
    points *= spacing_nm[None, :]
    ax.scatter(
        points[:, 0], points[:, 1], points[:, 2],
        s=size, c=color, alpha=alpha, depthshade=False, linewidths=0, rasterized=True,
    )


def set_nm_axes(ax, limits: tuple[tuple[float, float], ...], aspect: np.ndarray) -> None:
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_zlim(*limits[2])
    ax.set_box_aspect(aspect)
    ax.set_proj_type("ortho")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_zlabel("z (nm)")
    ax.grid(True, alpha=0.35)


def main() -> None:
    if not GEOMETRY.exists() or not CHECKPOINT.exists():
        raise FileNotFoundError(f"Geometry or checkpoint missing: {CHECKPOINT}")
    with np.load(GEOMETRY, allow_pickle=False) as archive:
        scaffold = np.asarray(archive["skeleton"], dtype=bool)
        size_A = np.asarray(archive["size_A"], dtype=float)
    with np.load(CHECKPOINT, allow_pickle=False) as archive:
        temperature = np.asarray(archive["temperature_K"], dtype=float)
        liquid_fraction = np.asarray(archive["ag_liquid_fraction"], dtype=float)
        metadata = json.loads(str(archive["metadata_json"].item()))

    ag = ~scaffold
    spacing_A = size_A / np.asarray(ag.shape, dtype=float)
    spacing_nm = spacing_A / 10.0
    liquid = ag & (liquid_fraction >= 0.5)
    labels, count = ndimage.label(liquid, ndimage.generate_binary_structure(3, 1))
    counts = np.bincount(labels.ravel())
    if counts.size:
        counts[0] = 0
    label_id = int(np.argmax(counts))
    pool = labels == label_id
    voxels = int(counts[label_id])
    volume_nm3 = voxels * float(np.prod(spacing_nm))
    coords = np.where(pool)
    lo = np.array([c.min() for c in coords])
    hi = np.array([c.max() for c in coords])
    bbox_nm = (hi - lo + 1) * spacing_nm
    centroid_nm = (np.array([c.mean() for c in coords]) + 0.5) * spacing_nm
    size_nm = size_A / 10.0
    time_ps = float(metadata["time_ps"])
    peak_K = float(temperature.max())
    ag_liquid_fraction = float(np.sum(ag * liquid_fraction) / np.sum(ag))
    roi_lo = np.maximum(lo - 5, 0)
    roi_hi = np.minimum(hi + 6, np.array(pool.shape))
    roi = tuple(slice(int(a), int(b)) for a, b in zip(roi_lo, roi_hi, strict=True))
    partial = ag[roi] & (liquid_fraction[roi] > 1.0e-12) & (liquid_fraction[roi] < 0.5)

    fig = plt.figure(figsize=(17.2, 8.8))
    fig.subplots_adjust(left=0.035, right=0.985, top=0.84, bottom=0.10, wspace=0.03)
    ax_full = fig.add_subplot(1, 2, 1, projection="3d")
    solid_ag = ag & (liquid_fraction < 0.5)
    partial_all = ag & (liquid_fraction > 1.0e-12) & (liquid_fraction < 0.5)
    cloud_points(ax_full, solid_ag, spacing_nm, "#83b7d5", stride=3, size=2.2, alpha=0.18)
    cloud_points(ax_full, scaffold, spacing_nm, "#263f51", stride=3, size=2.4, alpha=0.34)
    add_pool_surface(ax_full, pool, spacing_nm, np.zeros(3, dtype=int), alpha=0.48)
    partial_points = np.argwhere(partial_all)
    if partial_points.size:
        partial_points = partial_points[::max(1, len(partial_points) // 12000 + 1)].astype(float) + 0.5
        partial_points *= spacing_nm[None, :]
        ax_full.scatter(partial_points[:, 0], partial_points[:, 1], partial_points[:, 2],
                        s=4.0, c="#ffd166", alpha=0.62, depthshade=False,
                        linewidths=0, rasterized=True)
    add_sample_box(ax_full, size_nm)
    ax_full.scatter(*centroid_nm, color="#b43f16", s=34, depthshade=False, label="pool centroid")
    source_nm = np.array([size_nm[0] / 2, size_nm[1] / 2, size_nm[2]])
    ax_full.scatter(*source_nm, color="#2674ae", marker="*", s=110, depthshade=False, label="Gaussian source")
    set_nm_axes(ax_full, tuple((0.0, float(v)) for v in size_nm), size_nm)
    ax_full.set_title("Whole dual-phase QSGS sample", pad=10)
    ax_full.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#83b7d5", markersize=7, label="Ag matrix (sampled points)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#263f51", markersize=7, label="Ti3SiC2 scaffold (sampled points)"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#f28e2b", markersize=7, label="connected liquid Ag surface"),
        Line2D([0], [0], marker="*", color="w", markerfacecolor="#2674ae", markersize=10, label="Gaussian source"),
    ], loc="upper left", fontsize=8, framealpha=0.9)

    ax_zoom = fig.add_subplot(1, 2, 2, projection="3d")
    ag_roi = ag[roi] & ~pool[roi] & ~partial
    scaffold_roi = scaffold[roi]
    roi_spacing = spacing_nm
    roi_origin_nm = roi_lo * spacing_nm
    for mask, color, size, alpha in (
        (ag_roi, "#83b7d5", 5.0, 0.24),
        (scaffold_roi, "#263f51", 5.4, 0.34),
    ):
        pts = np.argwhere(mask)
        if pts.size:
            # Keep the close-up legible while preserving the original voxel coordinates.
            pts = pts[::max(1, len(pts) // 9000 + 1)].astype(float) + 0.5
            pts = (pts + roi_lo[None, :]) * roi_spacing[None, :]
            ax_zoom.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=size, c=color, alpha=alpha,
                            depthshade=False, linewidths=0, rasterized=True)
    pcoords = np.argwhere(partial)
    if pcoords.size:
        pcoords = pcoords[::max(1, len(pcoords) // 3500 + 1)].astype(float) + 0.5
        pcoords = (pcoords + roi_lo[None, :]) * roi_spacing[None, :]
        ax_zoom.scatter(pcoords[:, 0], pcoords[:, 1], pcoords[:, 2], s=9, c="#ffd166",
                        alpha=0.72, edgecolors="#8d6814", linewidths=0.2, depthshade=False, rasterized=True)
    add_pool_surface(ax_zoom, pool[roi], spacing_nm, roi_lo, alpha=0.38)
    lower_nm, upper_nm = roi_lo * spacing_nm, roi_hi * spacing_nm
    local_limits = tuple((float(lower_nm[i]), float(upper_nm[i])) for i in range(3))
    local_aspect = (roi_hi - roi_lo) * spacing_nm
    set_nm_axes(ax_zoom, local_limits, local_aspect)
    ax_zoom.view_init(elev=27, azim=-54)
    ax_zoom.set_title("Largest connected liquid Ag pool (f >= 0.5)", pad=10)

    ax_full.view_init(elev=24, azim=-53)
    fig.suptitle(
        f"Ag-Ti3SiC2 70:30 vol% | 1x Case A | t={time_ps:.3f} ps | peak T={peak_K:.1f} K\n"
        f"Pool: {voxels:,} voxels, {volume_nm3:.1f} nm^3, "
        f"bbox={bbox_nm[0]:.1f} x {bbox_nm[1]:.1f} x {bbox_nm[2]:.1f} nm | "
        f"Ag liquid fraction={100*ag_liquid_fraction:.2f}%",
        fontsize=15,
    )
    fig.text(
        0.5, 0.035,
        "Orange is the connected liquid cluster at f >= 0.5; grey/blue points show solid phase context. "
        "This is a thermal enthalpy result, not a fluid or spatter simulation.",
        ha="center", fontsize=9,
    )
    OUTPUT_COMPOSITE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_COMPOSITE, dpi=190, facecolor="white", bbox_inches="tight")
    plt.close(fig)

    fig2 = plt.figure(figsize=(9.5, 8.5))
    ax = fig2.add_subplot(1, 1, 1, projection="3d")
    pool_crop = pool[roi]
    pool_crop_size_nm = (roi_hi - roi_lo) * spacing_nm
    add_pool_surface(ax, pool_crop, spacing_nm, np.zeros(3, dtype=int))
    limits = tuple((0.0, float(v)) for v in pool_crop_size_nm)
    set_nm_axes(ax, limits, pool_crop_size_nm)
    ax.view_init(elev=27, azim=-54)
    ax.set_title(
        f"3D liquid-pool shape at {time_ps:.3f} ps\n"
        f"{volume_nm3:.1f} nm^3 | {bbox_nm[0]:.1f} x {bbox_nm[1]:.1f} x {bbox_nm[2]:.1f} nm",
        pad=14,
    )
    fig2.tight_layout()
    fig2.savefig(OUTPUT_POOL, dpi=220, facecolor="white", bbox_inches="tight")
    plt.close(fig2)

    # This slice shows the phases hidden behind the 3D pool surface.
    z_index = int(np.clip(round(coords[2].mean()), 0, ag.shape[2] - 1))
    phase = np.zeros(ag.shape[:2], dtype=np.uint8)
    phase[ag[:, :, z_index]] = 1
    phase[partial_all[:, :, z_index]] = 2
    phase[pool[:, :, z_index]] = 3
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    fig3, ax_slice = plt.subplots(figsize=(9.2, 8.2))
    cmap = ListedColormap(["#263f51", "#83b7d5", "#ffd166", "#f28e2b"])
    extent_xy = (0.0, float(size_nm[0]), 0.0, float(size_nm[1]))
    ax_slice.imshow(phase.T, origin="lower", extent=extent_xy, interpolation="nearest",
                    aspect="equal", cmap=cmap, vmin=0, vmax=3)
    ax_slice.contour(pool[:, :, z_index].T.astype(float), levels=[0.5], colors=["#9b3d08"],
                     linewidths=0.8, origin="lower", extent=extent_xy)
    ax_slice.set_xlabel("x (nm)")
    ax_slice.set_ylabel("y (nm)")
    ax_slice.set_title(
        f"Phase map through pool center: z={(z_index + 0.5) * spacing_nm[2]:.2f} nm, "
        f"t={time_ps:.3f} ps"
    )
    ax_slice.legend(handles=[
        Patch(color="#263f51", label="Ti3SiC2 scaffold"),
        Patch(color="#83b7d5", label="unmelted Ag"),
        Patch(color="#ffd166", label="partial-liquid Ag"),
        Patch(color="#f28e2b", label="liquid Ag (f >= 0.5)"),
    ], loc="upper right", framealpha=0.92)
    ax_slice.grid(False)
    fig3.tight_layout()
    fig3.savefig(OUTPUT_SLICE, dpi=220, facecolor="white", bbox_inches="tight")
    plt.close(fig3)

    print(f"checkpoint={CHECKPOINT}")
    print(f"saved={OUTPUT_COMPOSITE}")
    print(f"saved={OUTPUT_POOL}")
    print(f"saved={OUTPUT_SLICE}")
    print(f"time_ps={time_ps:.6f} peak_K={peak_K:.4f} pool_voxels={voxels} volume_nm3={volume_nm3:.6f}")
    print(f"bbox_nm={bbox_nm.tolist()} ag_liquid_fraction={ag_liquid_fraction:.9g} components={count}")


if __name__ == "__main__":
    main()
