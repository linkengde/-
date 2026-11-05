#!/usr/bin/env python3
"""Render P1H fixed-grid thermal/VOF fields as a 3D melt-pool view."""
from __future__ import annotations

import argparse
import gzip
import struct
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HEADER = struct.Struct("<8s6I5d")
EXPECTED_MAGIC = b"P1HSNAP1"
FIELD_COUNT = 20


def read_frame(path: Path):
    raw = gzip.open(path, "rb").read() if path.suffix == ".gz" else path.read_bytes()
    magic, version, nx, ny, nz, nfields, flags, time_ps, dx_m, x0, y0, z0 = HEADER.unpack_from(raw)
    if magic != EXPECTED_MAGIC or version != 1 or nfields != FIELD_COUNT:
        raise ValueError(f"unexpected P1H snapshot header in {path}")
    expected = HEADER.size + 4 * nfields * nx * ny * nz
    if len(raw) != expected:
        raise ValueError(f"incomplete snapshot: expected {expected} bytes, got {len(raw)}")
    fields = np.frombuffer(raw, dtype="<f4", offset=HEADER.size).reshape(nfields, nx, ny, nz)
    if not np.isfinite(fields).all():
        raise ValueError("snapshot contains non-finite values")
    return fields, (nx, ny, nz), flags, time_ps, dx_m, (x0, y0, z0)


def draw_panel(ax, coords, masks, colors, alphas, limits, title, downsample):
    x, y, z = coords
    colors_by_name = colors
    for name in ("solid_ag", "skeleton", "partial_ag", "liquid_ag"):
        mask = masks[name]
        if downsample is not None and name in ("solid_ag", "skeleton"):
            mask = mask & downsample
        if not np.any(mask):
            continue
        ax.scatter(x[mask], y[mask], z[mask], s=4 if name in ("partial_ag", "liquid_ag") else 2.4,
                   c=colors_by_name[name], alpha=alphas[name], depthshade=False,
                   linewidths=0, rasterized=True)
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = limits
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin, zmax)
    ax.set_box_aspect((xmax-xmin, ymax-ymin, zmax-zmin))
    ax.set_xlabel("x (nm)", labelpad=1)
    ax.set_ylabel("y (nm)", labelpad=1)
    ax.set_zlabel("z (nm)", labelpad=1)
    ax.set_title(title, pad=4)
    ax.view_init(elev=25, azim=-58)
    ax.grid(False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--main", required=True, type=Path, help="Compressed or raw P1H thermal/VOF frame")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--pressure-mpa", type=float, default=0.0)
    parser.add_argument("--cluster-volume-nm3", type=float, default=float("nan"))
    parser.add_argument("--cluster-bbox-nm", type=float, nargs=3, metavar=("X", "Y", "Z"), default=(float("nan"),)*3)
    parser.add_argument("--pool-bounds-nm", type=float, nargs=6, metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
                        default=(3.32, 15.09, 3.62, 15.39, 11.77, 15.70))
    args = parser.parse_args()

    fields, (nx, ny, nz), flags, time_ps, dx_m, (x0_m, y0_m, z0_m) = read_frame(args.main)
    vof, ag, skeleton, liquid_fraction = fields[0], fields[2], fields[3], fields[4]
    temperature = fields[5]
    dx_nm = dx_m * 1e9
    x = x0_m*1e9 + np.arange(nx, dtype=np.float32)*dx_nm
    y = y0_m*1e9 + np.arange(ny, dtype=np.float32)*dx_nm
    z = z0_m*1e9 + np.arange(nz, dtype=np.float32)*dx_nm
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")

    sample = vof > 0.5
    sc = sample & (skeleton > 0.5)
    ag_cells = sample & (ag > 0.5)
    solid = ag_cells & (liquid_fraction < 0.01)
    partial = ag_cells & (liquid_fraction >= 0.01) & (liquid_fraction < 0.5)
    liquid = ag_cells & (liquid_fraction >= 0.5)
    masks = {"skeleton": sc, "solid_ag": solid, "partial_ag": partial, "liquid_ag": liquid}
    coords = (xx, yy, zz)
    colors = {"skeleton":"#1c2933", "solid_ag":"#83b6d2", "partial_ag":"#ffe24a", "liquid_ag":"#f27b22"}
    alphas = {"skeleton":0.35, "solid_ag":0.18, "partial_ag":0.88, "liquid_ag":0.56}

    ix, iy, iz = np.indices((nx, ny, nz))
    stride = ((ix % 2) == 0) & ((iy % 2) == 0) & ((iz % 2) == 0)
    full_limits = ((x[0], x[-1]), (y[0], y[-1]), (z[0], min(z[-1], 16.3)))
    bx0,bx1,by0,by1,bz0,bz1=args.pool_bounds_nm
    zoom_limits=((bx0-1.0,bx1+1.0),(by0-1.0,by1+1.0),(bz0-0.9,bz1+0.7))

    fig = plt.figure(figsize=(14, 7.6), facecolor="white")
    ax1=fig.add_subplot(1,2,1,projection="3d")
    ax2=fig.add_subplot(1,2,2,projection="3d")
    draw_panel(ax1,coords,masks,colors,alphas,full_limits,"Ag70 exposed surface region",stride)
    draw_panel(ax2,coords,masks,colors,{"skeleton":0.72,"solid_ag":0.3,"partial_ag":0.95,"liquid_ag":0.66},zoom_limits,"Peak liquid-cluster close-up",None)
    # The solver counts any non-negligible skeleton occupancy when reporting
    # Tmax_SC; visualization uses >0.5 only to assign opaque skeleton markers.
    sc_occupied=skeleton > 1e-8
    max_sc_temp=float(np.max(temperature[sc_occupied])) if np.any(sc_occupied) else float("nan")
    fig.suptitle(f"Ag–Ti₃SiC₂ melt pool at t = {time_ps:.6f} ps | Gaussian pressure = {args.pressure_mpa:g} MPa\n"
                 f"Skeleton peak = {max_sc_temp:.1f} K | largest connected liquid cluster = {args.cluster_volume_nm3:.2f} nm³ | "
                 f"bbox = {args.cluster_bbox_nm[0]:.2f} × {args.cluster_bbox_nm[1]:.2f} × {args.cluster_bbox_nm[2]:.2f} nm",
                 fontsize=14, y=0.98)
    legend=[Line2D([0],[0],marker="o",color="w",label="Ti₃SiC₂ skeleton",markerfacecolor=colors["skeleton"],markersize=8),
            Line2D([0],[0],marker="o",color="w",label="Solid Ag",markerfacecolor=colors["solid_ag"],markersize=8),
            Line2D([0],[0],marker="o",color="w",label="Partially molten Ag",markerfacecolor=colors["partial_ag"],markersize=8),
            Line2D([0],[0],marker="o",color="w",label="Liquid-dominant Ag",markerfacecolor=colors["liquid_ag"],alpha=0.65,markersize=8)]
    fig.legend(handles=legend,loc="lower center",ncol=4,frameon=False,bbox_to_anchor=(0.5,0.015))
    fig.subplots_adjust(top=0.84,bottom=0.11,wspace=0.02)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=220,bbox_inches="tight")
    print(f"wrote {args.output}; shape={nx}x{ny}x{nz}, time={time_ps:.9f} ps, dx={dx_nm:.9f} nm")
    print(f"cells: skeleton={int(sc.sum())}, solid Ag={int(solid.sum())}, partial Ag={int(partial.sum())}, liquid-dominant Ag={int(liquid.sum())}")


if __name__ == "__main__":
    main()
