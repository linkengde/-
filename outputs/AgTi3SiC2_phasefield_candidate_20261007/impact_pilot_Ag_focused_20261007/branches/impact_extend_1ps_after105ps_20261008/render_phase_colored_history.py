#!/usr/bin/env python3
"""Render the saved 80/230/476 impact snapshots colored by material phase."""
from pathlib import Path
import importlib.util
import numpy as np
from scipy.ndimage import gaussian_filter, gaussian_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
PILOT = ROOT.parents[1]
TO100 = ROOT.parent / "to100ps_200fs_20261008"
parser_path = TO100 / "audit" / "audit_crater.py"
spec = importlib.util.spec_from_file_location("audit_crater", parser_path)
ac = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ac)

CENTER = np.array([9.9573898076, 24.2578814220])
BASE = ac.read_data(PILOT / "structures" / "base_active_local_relaxed.data")
STATES = [
    ("80 次 Ag 轰击", TO100 / "checkpoints" / "after_hit_080.data"),
    ("230 次 Ag 轰击", TO100 / "checkpoints" / "after_hit_230.data"),
    ("476 次 Ag 轰击", ROOT / "checkpoints" / "after_hit_476.data"),
]
TSC_COLOR = "#6841A5"
AG_COLOR = "#E8AD27"


def local_xy(atoms):
    h = atoms["H"]
    origin = np.array([atoms["bounds"]["x"][0], atoms["bounds"]["y"][0]])
    center_frac = np.linalg.solve(h, CENTER - origin)
    df = atoms["frac"] - center_frac
    candidates = np.stack([
        (df + np.array([i, j])) @ h.T
        for i in (-2, -1, 0, 1, 2)
        for j in (-2, -1, 0, 1, 2)
    ])
    ix = np.argmin(np.sum(candidates * candidates, axis=2), axis=0)
    return candidates[ix, np.arange(len(df))]


def radial_surface_profile(atoms, nbin=96):
    u = np.mod(atoms["frac"][:, 0], 1.0)
    v = np.mod(atoms["frac"][:, 1], 1.0)
    z = atoms["xyz"][:, 2]
    grid = np.full((nbin, nbin), np.nan)
    keep = (z > -5.0) & (z < 55.0)
    for i, j, zz in zip((u[keep] * nbin).astype(int), (v[keep] * nbin).astype(int), z[keep]):
        if not np.isfinite(grid[j, i]) or zz > grid[j, i]:
            grid[j, i] = zz
    valid = np.isfinite(grid)
    num = gaussian_filter(np.where(valid, grid, 0.0), 3.0, mode="wrap")
    den = gaussian_filter(valid.astype(float), 3.0, mode="wrap")
    surface = np.full_like(grid, np.nan)
    good = den > 0.08
    surface[good] = num[good] / den[good]

    f = (np.arange(nbin) + 0.5) / nbin
    vv, uu = np.meshgrid(f, f, indexing="ij")
    h = atoms["H"]
    origin = np.array([atoms["bounds"]["x"][0], atoms["bounds"]["y"][0]])
    center_frac = np.linalg.solve(h, CENTER - origin)
    frac = np.stack((uu.ravel(), vv.ravel()), axis=1)
    df = frac - center_frac
    candidates = np.stack([
        (df + np.array([i, j])) @ h.T
        for i in (-1, 0, 1)
        for j in (-1, 0, 1)
    ])
    ix = np.argmin(np.sum(candidates * candidates, axis=2), axis=0)
    closest = candidates[ix, np.arange(len(df))]
    radius = np.linalg.norm(closest, axis=1).reshape(nbin, nbin)
    edges = np.arange(0.0, 20.5001, 0.5)
    centers = (edges[:-1] + edges[1:]) / 2
    med = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (radius >= lo) & (radius < hi) & np.isfinite(surface)
        med.append(float(np.median(surface[mask])) if mask.sum() >= 3 else np.nan)
    med = np.asarray(med)
    g = np.isfinite(med)
    if g.sum() > 2:
        idx = np.arange(len(med))
        med[~g] = np.interp(idx[~g], idx[g], med[g])
        med = gaussian_filter1d(med, 1.0, mode="nearest")
    return centers, med


FONT_PATH = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_manager.fontManager.addfont(FONT_PATH)
FONT_FAMILY = font_manager.FontProperties(fname=FONT_PATH).get_name()
plt.rcParams.update({
    "font.family": FONT_FAMILY,
    "font.size": 9,
    "axes.titleweight": "bold",
    "axes.unicode_minus": False,
})
fig = plt.figure(figsize=(15.5, 9.0), dpi=190)
gs = fig.add_gridspec(2, 3, height_ratios=[1.14, 1.0], left=0.055,
                      right=0.985, top=0.88, bottom=0.16, hspace=0.32, wspace=0.22)
base_r, base_h = radial_surface_profile(BASE)
legend = [
    Line2D([], [], marker="o", linestyle="", markersize=6, color=TSC_COLOR,
           label="Ti3SiC2 骨架（Ti、Si、C）"),
    Line2D([], [], marker="o", linestyle="", markersize=6, color=AG_COLOR,
           label="Ag 银相"),
    Line2D([], [], color="#777777", linestyle="--", label="未轰击表面中位轮廓"),
    Line2D([], [], color="#111111", linestyle="-", label="轰击后表面中位轮廓"),
]

for col, (label, path) in enumerate(STATES):
    state = ac.read_data(path)
    loc = local_xy(state)
    xyz = state["xyz"]
    typ = state["typ"]
    view = ((np.abs(loc[:, 0]) < 20) & (np.abs(loc[:, 1]) < 20) &
            (xyz[:, 2] >= -5) & (xyz[:, 2] < 45))
    tsc = view & np.isin(typ, [1, 2, 3])
    ag = view & (typ == 4)

    ax = fig.add_subplot(gs[0, col], projection="3d")
    # Plot Ag first and the connected ceramic skeleton afterward so its branches
    # remain visible through the dense point cloud.
    ax.scatter(loc[ag, 0], loc[ag, 1], xyz[ag, 2], s=7, c=AG_COLOR,
               alpha=0.58, linewidths=0, depthshade=False)
    ax.scatter(loc[tsc, 0], loc[tsc, 1], xyz[tsc, 2], s=8, c=TSC_COLOR,
               alpha=0.94, linewidths=0, depthshade=False)
    ax.set(xlim=(-20, 20), ylim=(-20, 20), zlim=(-5, 45),
           xlabel="x（Å）", ylabel="y（Å）", zlabel="z（Å）")
    ax.set_title(label, fontsize=12, pad=2)
    ax.view_init(elev=22, azim=-50)
    ax.set_box_aspect((40, 40, 50))
    ax.tick_params(labelsize=7, pad=0)
    ax.xaxis.labelpad = -1
    ax.yaxis.labelpad = -1
    ax.zaxis.labelpad = -1

    ax2 = fig.add_subplot(gs[1, col])
    section = view & (np.abs(loc[:, 1]) < 2.0)
    tsc_sec = section & np.isin(typ, [1, 2, 3])
    ag_sec = section & (typ == 4)
    ax2.scatter(loc[ag_sec, 0], xyz[ag_sec, 2], s=13, c=AG_COLOR,
                alpha=0.68, linewidths=0, zorder=2)
    ax2.scatter(loc[tsc_sec, 0], xyz[tsc_sec, 2], s=14, c=TSC_COLOR,
                alpha=0.94, linewidths=0, zorder=3)
    r, h = radial_surface_profile(state)
    ax2.plot(base_r, base_h, color="#777777", lw=1.25, ls="--", zorder=4)
    ax2.plot(-base_r, base_h, color="#777777", lw=1.25, ls="--", zorder=4)
    ax2.plot(r, h, color="#111111", lw=1.6, zorder=5)
    ax2.plot(-r, h, color="#111111", lw=1.6, zorder=5)
    ax2.set(xlim=(-20, 20), ylim=(-5, 45), xlabel="y≈0 处的 x 坐标（Å）", ylabel="z（Å）")
    ax2.set_title("轰击中心纵向原子截面", fontsize=10, pad=3)
    ax2.grid(alpha=0.18)

fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, 0.935),
           ncol=4, frameon=False, fontsize=9)
fig.suptitle("Ag–Ti3SiC2 轰击形貌：骨架与银相分色", fontsize=15,
             fontweight="bold", y=0.985)
fig.text(0.5, 0.055,
         "紫色：Ti3SiC2 骨架（Ti、Si、C）；金黄色：Ag。颜色表示材料相，不表示温度或位移。",
         ha="center", fontsize=10)
out = ROOT / "images" / "deformation_history_phase_colors_80_230_476.png"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(out)
print(f"bytes={out.stat().st_size}")
