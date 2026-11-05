#!/usr/bin/env python3
"""Render a dense, atom-like 3D view of a P1H continuum-grid snapshot."""
from __future__ import annotations

import argparse
import gzip
import struct
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D

HEADER = struct.Struct("<8s6I5d")
MAIN_MAGIC = b"P1HSNAP1"
MAIN_FIELD_COUNT = 20


def read_frame(path: Path):
    raw = gzip.open(path, "rb").read() if path.suffix == ".gz" else path.read_bytes()
    if len(raw) < HEADER.size:
        raise ValueError(f"short P1H frame: {path}")
    magic, version, nx, ny, nz, nfields, flags, time_ps, dx_m, x0, y0, z0 = HEADER.unpack_from(raw)
    if (magic, version, nfields) != (MAIN_MAGIC, 1, MAIN_FIELD_COUNT):
        raise ValueError(f"unexpected P1H frame header in {path}")
    expected = HEADER.size + 4 * nfields * nx * ny * nz
    if len(raw) != expected:
        raise ValueError(f"partial frame: expected {expected} bytes, got {len(raw)}")
    fields = np.frombuffer(raw, dtype="<f4", offset=HEADER.size).reshape(nfields, nx, ny, nz)
    if not np.isfinite(fields).all():
        raise ValueError("P1H frame contains NaN or Inf")
    return fields, (nx, ny, nz), flags, time_ps, dx_m, (x0, y0, z0)


def add_box_edges(ax, bounds, color="#aab3bb", linewidth=0.9):
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = bounds
    for y in (ymin, ymax):
        for z in (zmin, zmax):
            ax.plot([xmin, xmax], [y, y], [z, z], color=color, lw=linewidth)
    for x in (xmin, xmax):
        for z in (zmin, zmax):
            ax.plot([x, x], [ymin, ymax], [z, z], color=color, lw=linewidth)
    for x in (xmin, xmax):
        for y in (ymin, ymax):
            ax.plot([x, x], [y, y], [zmin, zmax], color=color, lw=linewidth)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--main", required=True, type=Path, help="Compressed/raw P1H thermal and material frame")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--pressure-mpa", type=float, default=0.)
    parser.add_argument("--ratio", default="70:30")
    parser.add_argument("--largest-cluster-cells", type=int, default=0)
    parser.add_argument("--largest-cluster-nm3", type=float, default=float("nan"))
    parser.add_argument("--bbox-nm", type=float, nargs=3, metavar=("X", "Y", "Z"), default=(float("nan"),)*3)
    args = parser.parse_args()

    # Noto Sans CJK is available in the cloud image; register it for Chinese labels.
    font_candidates = ["/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
                      "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc"]
    for font in font_candidates:
        if Path(font).exists():
            font_manager.fontManager.addfont(font)
            plt.rcParams["font.family"] = font_manager.FontProperties(fname=font).get_name()
            break
    plt.rcParams["axes.unicode_minus"] = False

    fields, (nx, ny, nz), _, time_ps, dx_m, origin_m = read_frame(args.main)
    liquid_fraction = fields[4]
    temperature = fields[5]
    material_id = np.rint(fields[19]).astype(np.uint8)
    ag_cells = material_id == 1
    skeleton = material_id == 2
    solid_ag = ag_cells & (liquid_fraction <= 0.01)
    partial_ag = ag_cells & (liquid_fraction > 0.01) & (liquid_fraction < 0.5)
    liquid_ag = ag_cells & (liquid_fraction >= 0.5)
    sample = ag_cells | skeleton

    # Render the sample's outer and pore-facing shell densely, as in an
    # atom/voxel model, while retaining the entire liquid cluster so the pool
    # remains visible through the translucent cell display.
    padded = np.pad(sample, 1, constant_values=False)
    interior = (padded[2:,1:-1,1:-1] & padded[:-2,1:-1,1:-1] &
                padded[1:-1,2:,1:-1] & padded[1:-1,:-2,1:-1] &
                padded[1:-1,1:-1,2:] & padded[1:-1,1:-1,:-2])
    surface_shell = sample & ~interior

    dx_nm = dx_m * 1e9
    origin_nm = np.asarray(origin_m, dtype=np.float64) * 1e9
    x = origin_nm[0] + (np.arange(nx) + 0.5) * dx_nm
    y = origin_nm[1] + (np.arange(ny) + 0.5) * dx_nm
    z = origin_nm[2] + (np.arange(nz) + 0.5) * dx_nm
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")

    occupied_z = np.flatnonzero(np.any(sample, axis=(0, 1)))
    if not occupied_z.size:
        raise ValueError("snapshot has no sample cells")
    bounds = ((origin_nm[0], origin_nm[0] + nx*dx_nm),
              (origin_nm[1], origin_nm[1] + ny*dx_nm),
              (origin_nm[2] + occupied_z[0]*dx_nm,
               origin_nm[2] + (occupied_z[-1]+1)*dx_nm))

    colors = {"skeleton":"#233b4a", "solid_ag":"#83b6d2",
              "partial_ag":"#ffed00", "liquid_ag":"#f28a22"}
    categories = [
        # Solid Ag is deliberately translucent so the subsurface liquid pool
        # can be seen through the porous matrix, as in the supplied reference.
        ("solid_ag", solid_ag & surface_shell, 20, 0.30, "#527e96"),
        ("skeleton", skeleton & surface_shell, 21, 0.96, "#162936"),
        ("partial_ag", partial_ag & surface_shell, 23, 1.00, "#cab900"),
        ("liquid_ag", liquid_ag & surface_shell, 24, 0.72, "#ad5112"),
    ]

    # Plot all occupied fixed-grid cell centers as shaded spheres. This is a
    # visualization of continuum cells, not an atomistic reconstruction.
    fig = plt.figure(figsize=(16, 10), facecolor="white")
    ax = fig.add_axes([0.10, 0.035, 0.86, 0.80], projection="3d")
    # The translucent internal liquid is drawn first; surface materials then
    # retain the shell colors and depth cues.
    inner_liquid = liquid_ag & ~surface_shell
    if np.any(inner_liquid):
        ax.scatter(xx[inner_liquid], yy[inner_liquid], zz[inner_liquid],
                   s=14, c=colors["liquid_ag"], alpha=0.18,
                   edgecolors="none", depthshade=False, rasterized=True)
    for name, mask, size, alpha, edge in categories:
        if np.any(mask):
            ax.scatter(xx[mask], yy[mask], zz[mask], s=size, c=colors[name],
                       alpha=alpha, edgecolors=edge, linewidths=0.12,
                       depthshade=True, rasterized=True)

    (xmin, xmax), (ymin, ymax), (zmin, zmax) = bounds
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin, zmax)
    ax.set_box_aspect((xmax-xmin, ymax-ymin, zmax-zmin))
    # Use a higher camera angle so the exposed top and the pool read clearly,
    # matching the reference render while preserving the physical aspect ratio.
    ax.view_init(elev=36, azim=-58)
    ax.set_proj_type("persp", focal_length=0.95)
    ax.set_axis_off()
    add_box_edges(ax, bounds)

    sc_peak = float(np.max(temperature[skeleton])) if np.any(skeleton) else float("nan")
    legend = [
        Line2D([0],[0], marker="s", color="none", label="Ti3SiC2 骨架", markerfacecolor=colors["skeleton"], markersize=13),
        Line2D([0],[0], marker="s", color="none", label="固态 Ag（f_l≈0）", markerfacecolor=colors["solid_ag"], markersize=13),
        Line2D([0],[0], marker="s", color="none", label="部分液化 Ag（0<f_l<0.5）", markerfacecolor=colors["partial_ag"], markersize=13),
        Line2D([0],[0], marker="s", color="none", label="液相占主导 Ag（f_l>=0.5）", markerfacecolor=colors["liquid_ag"], alpha=0.7, markersize=13),
    ]
    fig.legend(handles=legend, loc="upper left", bbox_to_anchor=(0.015, 0.985),
               frameon=False, ncol=1, fontsize=14, labelspacing=0.6, handletextpad=0.5)
    fig.text(0.58, 0.957,
             f"Ag–Ti3SiC2 三维相图 | {args.ratio} 体积比 | t = {time_ps:.6f} ps",
             ha="center", va="center", fontsize=22)
    fig.text(0.58, 0.915,
             f"骨架峰值温度 {sc_peak:.1f} K | 最大连通液相簇 {args.largest_cluster_cells:,} 个网格单元 | "
             f"{args.largest_cluster_nm3:.2f} nm³ | 压力 {args.pressure_mpa:g} MPa",
             ha="center", va="center", fontsize=17)
    fig.text(0.985, 0.018, "网格单元中心可视化；连续介质计算，不代表原子模型",
             ha="right", va="bottom", fontsize=9, color="#555555")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=125, facecolor="white")
    print(f"wrote {args.output}; t={time_ps:.9f} ps; grid={nx}x{ny}x{nz}; dx={dx_nm:.9f} nm")
    print(f"cells: skeleton={int(skeleton.sum())}, solid Ag={int(solid_ag.sum())}, "
          f"partial Ag={int(partial_ag.sum())}, liquid-dominant Ag={int(liquid_ag.sum())}")


if __name__ == "__main__":
    main()
