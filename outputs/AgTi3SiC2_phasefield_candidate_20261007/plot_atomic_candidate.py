#!/usr/bin/env python3
"""Plot the locally clearance-screened atomic candidate for morphology review."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "AgTi3SiC2_phasefield_atomic_candidate.data"
lx, ly, xy = 37.285974861019966, 32.290601434511245, -18.642987430509983
purple, gold = "#7851A9", "#E0B84F"


def read_atoms(path):
    atoms = []
    in_atoms = False
    for line in path.read_text().splitlines():
        fields = line.split()
        if line.strip().startswith("Atoms"):
            in_atoms = True
            continue
        if in_atoms and len(fields) >= 5 and fields[0].isdigit() and fields[1].isdigit():
            atoms.append((int(fields[0]), int(fields[1]), *map(float, fields[2:5])))
    arr = np.asarray(atoms, float)
    if len(arr) == 0:
        raise RuntimeError("Atomic candidate contains no atoms")
    return arr


atoms = read_atoms(DATA)
types = atoms[:, 1].astype(int)
xyz = atoms[:, 2:5]
v = (xyz[:, 1] / ly) % 1.0
u = ((xyz[:, 0] - xy * (xyz[:, 1] / ly)) / lx) % 1.0
colors = np.where(types == 4, gold, purple)

fig = plt.figure(figsize=(14.5, 5.6), dpi=180, facecolor="white")
ax = fig.add_subplot(1, 3, 1, projection="3d")
for typ, color, label, size in ((1, purple, "Ti", 4.5), (2, purple, "Si", 4.5),
                                (3, purple, "C", 4.5), (4, gold, "Ag", 5.0)):
    sel = types == typ
    ax.scatter(xyz[sel, 0], xyz[sel, 1], xyz[sel, 2], s=size, c=color,
               alpha=0.75, linewidths=0, label=label, depthshade=False)
ax.set_xlim(-19, 38)
ax.set_ylim(0, 33)
ax.set_zlim(0, 70)
ax.set_box_aspect((57, 33, 70))
ax.set_xlabel("x (Å)", labelpad=-2)
ax.set_ylabel("y (Å)", labelpad=-2)
ax.set_zlabel("z (Å)", labelpad=-2)
ax.set_title("3D atom arrangement\n(20 Å nominal top vacuum)", pad=5)
ax.view_init(elev=22, azim=-54)
ax.legend(loc="upper left", fontsize=7, frameon=False, markerscale=1.5, ncol=2)

side = (np.abs(v - 0.5) < 0.12) & (xyz[:, 2] <= 52.0)
ax2 = fig.add_subplot(1, 3, 2)
ax2.scatter(u[side] * lx, xyz[side, 2], c=colors[side], s=8, alpha=0.85, linewidths=0)
ax2.axhspan(50.4, 70, color="#eeeeee", zorder=0)
ax2.set_xlim(0, lx)
ax2.set_ylim(0, 70)
ax2.set_xlabel("x along a (Å)")
ax2.set_ylabel("z (Å)")
ax2.set_title("Vertical atomic slice\n(central 24% of the y cell)")
ax2.text(0.04, 0.94, "vacuum", transform=ax2.transAxes, va="top", color="#666666", fontsize=8)

top = xyz[:, 2] >= 45.0
ax3 = fig.add_subplot(1, 3, 3)
ax3.scatter(xyz[top, 0], xyz[top, 1], c=colors[top], s=10, alpha=0.9, linewidths=0)
poly = np.array([[0, 0], [lx, 0], [lx + xy, ly], [xy, ly], [0, 0]])
ax3.plot(poly[:, 0], poly[:, 1], color="#333333", lw=0.8)
ax3.set_aspect("equal")
ax3.set_xlim(xy - 1, lx + 1)
ax3.set_ylim(-1, ly + 1)
ax3.set_xlabel("x (Å)")
ax3.set_ylabel("y (Å)")
ax3.set_title("Top atomic layers\n(z ≥ 45 Å)")

fig.suptitle("Ag–Ti$_3$SiC$_2$ phase-field mapped geometry | TSC purple, Ag gold", fontsize=13, y=0.99)
fig.text(0.5, 0.015,
         "Geometry-only local clearance projection; no force-field relaxation or dynamic validation.",
         ha="center", fontsize=9, color="#555555")
fig.tight_layout(rect=(0, 0.045, 1, 0.94))
fig.savefig(ROOT / "phasefield_atomic_candidate_preview.png", dpi=180, bbox_inches="tight")
