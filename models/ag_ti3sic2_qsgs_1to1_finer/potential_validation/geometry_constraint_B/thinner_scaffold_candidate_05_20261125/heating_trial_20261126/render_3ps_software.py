#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 OVITO 导出的粒子坐标/qbar6 绘制 3.0 ps 三维视图（Matplotlib 软件后端）。"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = np.load(HERE / "frame_3ps_ovito_analysis_arrays.npz")
pos = DATA["positions"].astype(np.float64)
types = DATA["types"]
qbar = DATA["qbar6"]
ag = types == 4
skeleton = ~ag
norm = Normalize(vmin=0.17, vmax=0.62, clip=True)
cmap = plt.get_cmap("coolwarm_r")  # qbar6 低值红色，高值蓝色


def draw_scene(ax, mask, crop=False):
    """绘制一个全模型或顶部热区的三维散点图。"""
    p = pos[mask]
    t = types[mask]
    q = qbar[mask]
    is_ag = t == 4
    ax.scatter(p[~is_ag, 0], p[~is_ag, 1], p[~is_ag, 2],
               c="#103b78" if not crop else "#527eb2",
               alpha=1.0 if not crop else 0.24,
               s=0.23 if not crop else 0.62,
               depthshade=False, linewidths=0, rasterized=True)
    ax.scatter(p[is_ag, 0], p[is_ag, 1], p[is_ag, 2], c=q[is_ag],
               cmap=cmap, norm=norm, s=0.36 if not crop else 5.0,
               depthshade=False, linewidths=0, rasterized=True)
    if crop:
        ax.set_xlim(67, 133)
        ax.set_ylim(67, 133)
        ax.set_zlim(116, 151)
        ax.set_title("Heated core (rxy <= 30 A; top 30 A)", fontsize=16, pad=8)
    else:
        ax.set_xlim(0, 200)
        ax.set_ylim(0, 200)
        ax.set_zlim(0, 150)
        ax.set_title("Full model", fontsize=16, pad=8)
    ax.view_init(elev=23 if not crop else 25, azim=-54)
    ax.set_box_aspect((1, 1, 0.75 if not crop else 0.64))
    ax.set_proj_type("persp", focal_length=1.4)
    ax.set_axis_off()


fig = plt.figure(figsize=(18, 8.6), facecolor="white")
fig.suptitle("Ag / Ti3SiC2 candidate 05 after 3.0 ps TTM heating", fontsize=20, y=0.96)
full_mask = np.ones(len(pos), dtype=bool)
r_xy2 = (pos[:, 0] - 100.0)**2 + (pos[:, 1] - 100.0)**2
hot_mask = (r_xy2 <= 30.0**2) & (pos[:, 2] >= 120.0)
ax1 = fig.add_subplot(1, 2, 1, projection="3d")
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
draw_scene(ax1, full_mask)
draw_scene(ax2, hot_mask, crop=True)
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=[ax1, ax2], fraction=0.026, pad=0.03, shrink=0.76)
cbar.set_label(r"Ag local $\bar{q}_6$ (relative order indicator)", fontsize=12)
cbar.ax.tick_params(labelsize=10)
legend = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#103b78", markersize=9,
           label="Rigid Ti/Si/C support"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#e76f51", markersize=9,
           label="Ag (colored by local qbar6)"),
]
fig.legend(handles=legend, loc="lower center", ncol=2, frameon=False,
           bbox_to_anchor=(0.49, 0.035), fontsize=11)
fig.text(0.5, 0.008,
         "Low qbar6 = more local disorder in this trajectory; the initial PTM classification is mostly 'Other'.",
         ha="center", fontsize=10, color="#434a52")
fig.subplots_adjust(left=0.01, right=0.91, top=0.89, bottom=0.11, wspace=0.02)
out = HERE / "candidate05_heated_3ps_OVITO_data_view.png"
fig.savefig(out, dpi=180, facecolor="white")
print(f"Wrote {out}; hot-region particles={int(hot_mask.sum())}, total={len(pos)}", flush=True)
