#!/usr/bin/env python3
"""Render a direct visual preview of the phase-field morphology."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
mask = np.load(ROOT / "phasefield_mask_tsc.npy").astype(bool)
nx, ny, nz = mask.shape
lx, ly, xy, h = 37.285974861019966, 32.290601434511245, -18.642987430509983, 50.0
purple, gold = "#7851A9", "#E0B84F"

def surface_cells(phase):
    edge = np.zeros_like(phase)
    for axis in (0, 1):
        edge |= phase != np.roll(phase, 1, axis=axis)
        edge |= phase != np.roll(phase, -1, axis=axis)
    edge[:, :, 1:] |= phase[:, :, 1:] != phase[:, :, :-1]
    edge[:, :, :-1] |= phase[:, :, :-1] != phase[:, :, 1:]
    edge[:, :, 0] |= phase[:, :, 0]
    edge[:, :, -1] |= phase[:, :, -1]
    return phase & edge

fig = plt.figure(figsize=(14.5, 5.5), dpi=180, facecolor="white")
ax = fig.add_subplot(1, 3, 1, projection="3d")
for phase, color, label in ((mask, purple, "Ti$_3$SiC$_2$ skeleton"), (~mask, gold, "Ag filling phase")):
    ii, jj, kk = np.where(surface_cells(phase))
    u, v = (ii + 0.5) / nx, (jj + 0.5) / ny
    x, y, z = u * lx + v * xy, v * ly, (kk + 0.5) * h / nz
    ax.scatter(x, y, z, s=2.2, c=color, alpha=0.82, linewidths=0, label=label, depthshade=False)
ax.set_xlim(-19, 38)
ax.set_ylim(0, 33)
ax.set_zlim(0, 70)
ax.set_box_aspect((57, 33, 70))
ax.set_xlabel("x (Å)", labelpad=-2)
ax.set_ylabel("y (Å)", labelpad=-2)
ax.set_zlabel("z (Å)", labelpad=-2)
ax.set_title("3D phase boundaries\n(20 Å top vacuum)", pad=5)
ax.view_init(elev=23, azim=-54)
ax.legend(loc="upper left", fontsize=7, frameon=False, markerscale=2.5)

ax2 = fig.add_subplot(1, 3, 2)
v_idx = ny // 2
ax2.imshow(mask[:, v_idx, :].T, origin="lower", interpolation="nearest", aspect="auto",
           extent=(0, lx, 0, h), cmap=matplotlib.colors.ListedColormap([gold, purple]), vmin=0, vmax=1)
ax2.axhspan(h, h + 20, color="#eeeeee", zorder=0)
ax2.set_xlim(0, lx)
ax2.set_ylim(0, h + 20)
ax2.set_xlabel("x along a (Å)")
ax2.set_ylabel("z (Å)")
ax2.set_title("Vertical slice (one y plane)")
ax2.text(0.04, 0.94, "vacuum", transform=ax2.transAxes, va="top", color="#666666", fontsize=8)

ax3 = fig.add_subplot(1, 3, 3)
top = mask[:, :, -3:].mean(axis=2)
ax3.imshow(top.T, origin="lower", interpolation="nearest", aspect="equal",
           extent=(0, 1, 0, 1), cmap=matplotlib.colors.ListedColormap([gold, purple]), vmin=0, vmax=1)
ax3.set_xlabel("fractional coordinate u")
ax3.set_ylabel("fractional coordinate v")
ax3.set_title("Top layer (upper 4.7 Å)\nTSC purple / Ag gold")

fig.suptitle("Ag–Ti$_3$SiC$_2$ bicontinuous morphology | conserved Cahn–Hilliard phase field",
             fontsize=13, y=0.99)
fig.text(0.5, 0.015, "Morphology-only voxel representation; phase-field parameters are not material-calibrated.",
         ha="center", fontsize=9, color="#555555")
fig.tight_layout(rect=(0, 0.045, 1, 0.94))
fig.savefig(ROOT / "phasefield_model_preview.png", dpi=180, bbox_inches="tight")
