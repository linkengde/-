#!/usr/bin/env python3
"""Render a real atomistic cutaway with a top surface derived from atom positions."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parents[1]
PARENT = HERE / "audit" / "audit_crater.py"
import importlib.util
spec = importlib.util.spec_from_file_location("audit_crater", PARENT)
ac = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ac)

CENTER = np.array([9.9573898076, 24.2578814220])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=HERE / "structures" / "cooled_after50_0p5ps.data")
    parser.add_argument("--out", type=Path, default=HERE / "images" / "crater_after50_highflux_3d.png")
    args = parser.parse_args()
    atoms = ac.read_data(args.data)
    xyz = atoms["xyz"]
    typ = atoms["typ"]
    a = np.array([atoms["lx"], 0.0])
    b = np.array([atoms["tilt"][0], atoms["ly"]])
    d0 = xyz[:, :2] - CENTER
    candidates = np.stack([d0 + i * a + j * b for i in (-1, 0, 1) for j in (-1, 0, 1)])
    best = np.argmin(np.sum(candidates * candidates, axis=2), axis=0)
    xy = candidates[best, np.arange(len(xyz))]
    radial = np.linalg.norm(xy, axis=1)
    z = xyz[:, 2]

    # Start with the highest atom per fractional XY tile. A sub-angstrom
    # normalized convolution removes single-atom spikes before interpolation.
    nbin = 72
    u = np.clip((atoms["frac"][:, 0] * nbin).astype(int), 0, nbin - 1)
    v = np.clip((atoms["frac"][:, 1] * nbin).astype(int), 0, nbin - 1)
    raw = np.full((nbin, nbin), np.nan)
    for k, (iu, iv) in enumerate(zip(u, v)):
        if 8.0 < z[k] < 52.0 and (not np.isfinite(raw[iv, iu]) or z[k] > raw[iv, iu]):
            raw[iv, iu] = z[k]
    valid = np.isfinite(raw)
    numerator = gaussian_filter(np.where(valid, raw, 0.0), 2.0, mode="wrap")
    denominator = gaussian_filter(valid.astype(float), 2.0, mode="wrap")
    smooth = np.full_like(raw, np.nan)
    good = denominator > 0.12
    smooth[good] = numerator[good] / denominator[good]

    uu = (np.arange(nbin) + 0.5) / nbin
    vv = (np.arange(nbin) + 0.5) / nbin
    vv_grid, uu_grid = np.meshgrid(vv, uu, indexing="ij")
    xlo, ylo = atoms["bounds"]["x"][0], atoms["bounds"]["y"][0]
    Xfrac = xlo + uu_grid * atoms["lx"] + vv_grid * atoms["tilt"][0]
    Yfrac = ylo + vv_grid * atoms["ly"]
    mesh_xy = np.column_stack((Xfrac.ravel(), Yfrac.ravel())) - CENTER
    mesh_images = np.stack([mesh_xy + i * a + j * b for i in (-1, 0, 1) for j in (-1, 0, 1)])
    mesh_best = np.argmin(np.sum(mesh_images * mesh_images, axis=2), axis=0)
    mesh_xy = mesh_images[mesh_best, np.arange(len(mesh_xy))]
    mesh_z = smooth.ravel()
    mesh_keep = np.isfinite(mesh_z) & (np.linalg.norm(mesh_xy, axis=1) < 11.5)
    gx = np.linspace(-14.25, 14.25, 100)
    gy = np.linspace(-14.25, 14.25, 100)
    X, Y = np.meshgrid(gx, gy)
    Z = griddata(mesh_xy[mesh_keep], mesh_z[mesh_keep], (X, Y), method="linear")
    Z[np.hypot(X, Y) > 11.25] = np.nan

    # A visual half-space cut opens the front so the real bicontinuous interior is visible.
    crop = (radial < 14.0) & (z >= 8.0) & (z <= 52.0) & (xy[:, 1] > -2.0)
    tsc = crop & (typ != 4)
    ag = crop & (typ == 4)

    fig = plt.figure(figsize=(13, 9.5), dpi=180)
    ax = fig.add_subplot(111, projection="3d")
    zlo, zhi = np.nanpercentile(Z, [3, 97])
    face = plt.get_cmap("coolwarm")(np.clip((np.nan_to_num(Z, nan=zlo) - zlo) / (zhi - zlo), 0, 1))
    face[..., 3] = np.where(np.isfinite(Z), 0.24, 0.0)
    ax.plot_surface(X, Y, np.ma.masked_invalid(Z), facecolors=face, rstride=1, cstride=1,
                    linewidth=0.0, edgecolor="none", antialiased=True,
                    shade=False, zorder=1)

    # Geometric TSC connectivity links (3.1 A criterion), not a force-field bond list.
    tsc_idx = np.flatnonzero(tsc)
    p = np.column_stack((xy[tsc_idx], z[tsc_idx]))
    tree = cKDTree(p)
    pairs = np.array(list(tree.query_pairs(3.10)), dtype=int)
    if len(pairs):
        segs = np.stack((p[pairs[:, 0]], p[pairs[:, 1]]), axis=1)
        lines = Line3DCollection(segs, colors=(0.64, 0.40, 0.79, 0.22), linewidths=0.40, zorder=3)
        ax.add_collection3d(lines)
    ax.scatter(xy[tsc, 0], xy[tsc, 1], z[tsc], s=8.0, color="#7935a8", alpha=0.88,
               edgecolors="none", depthshade=False, zorder=4)
    ax.scatter(xy[ag, 0], xy[ag, 1], z[ag], s=19.0, color="#f0b51b", alpha=0.95,
               edgecolors="#704b00", linewidths=0.25, depthshade=False, zorder=5)

    ax.set_xlim(-14, 14)
    ax.set_ylim(-14, 14)
    ax.set_zlim(8, 53)
    ax.set_xlabel("x from impact center (Å)", labelpad=5, fontsize=10)
    ax.set_ylabel("y from impact center (Å)", labelpad=5, fontsize=10)
    ax.set_zlabel("z (Å)", labelpad=5, fontsize=10)
    ax.set_title("", pad=0)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", label="Ti₃SiC₂ skeleton atoms", markerfacecolor="#7935a8", markersize=8),
        Line2D([0], [0], marker="o", color="w", label="Ag atoms", markerfacecolor="#f0b51b", markeredgecolor="#704b00", markersize=9),
        Line2D([0], [0], color="#9e6bc0", lw=1.2, label="TSC geometric links, 3.1 Å cutoff"),
        Patch(facecolor="#863850", alpha=0.42, label="Interpolated top-atom surface"),
    ], loc="upper left", bbox_to_anchor=(0.015, 0.90), framealpha=0.92, fontsize=9)
    ax.view_init(elev=29, azim=-58)
    ax.set_box_aspect((28, 28, 45))
    ax.xaxis.pane.set_facecolor((0.96, 0.97, 0.99, 0.8))
    ax.yaxis.pane.set_facecolor((0.96, 0.97, 0.99, 0.8))
    ax.zaxis.pane.set_facecolor((0.96, 0.97, 0.99, 0.8))
    fig.suptitle("Ag–Ti₃SiC₂ crater candidate | 45 prior impacts + 5 additional impacts at 200 fs + 0.5 ps local heat-sink relaxation",
                 y=0.985, fontsize=15)
    fig.text(0.5, 0.945,
             "Gold: Ag atoms  •  Purple: Ti₃SiC₂ atoms and 3.1 Å geometric links  •  Mesh: atom-derived, smoothed top heights",
             ha="center", va="top", fontsize=9.5)
    fig.text(0.5, 0.018,
             "Front half-space (y < −2 Å) hidden for visibility only; geometry is unchanged. No raised rim is resolved yet.",
             ha="center", va="bottom", fontsize=9)
    fig.subplots_adjust(top=0.90, bottom=0.11)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    print(f"saved {args.out} ({args.out.stat().st_size} bytes); displayed TSC={tsc.sum()} Ag={ag.sum()} atoms")


if __name__ == "__main__":
    main()
