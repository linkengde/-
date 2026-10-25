#!/usr/bin/env python3
"""Render a static cutaway PNG from the saved 3-D QSGS morphology."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
DISPLAY_BLOCK_XYZ = (4, 4, 5)


def main() -> None:
    with np.load(ROOT / "qsgs_morphology.npz", allow_pickle=False) as data:
        skeleton = data["skeleton"].astype(bool)
        size_A = data["size_A"].astype(float)
    if any(n % block for n, block in zip(skeleton.shape, DISPLAY_BLOCK_XYZ)):
        raise ValueError("QSGS source dimensions must be divisible by display block size")

    nx, ny, nz = (n // block for n, block in zip(skeleton.shape, DISPLAY_BLOCK_XYZ))
    bx, by, bz = DISPLAY_BLOCK_XYZ
    # Majority projection is only for drawing. The stored QSGS model and atom data
    # remain untouched; the caption discloses this display coarsening.
    display = skeleton.reshape(nx, bx, ny, by, nz, bz).mean(axis=(1, 3, 5)) >= 0.5
    blue = display.copy()
    grey = ~display

    # Remove one front-bottom corner so the image exposes the interior phases.
    cut = tuple(n // 3 for n in display.shape)
    blue[:cut[0], :cut[1], :cut[2]] = False
    grey[:cut[0], :cut[1], :cut[2]] = False

    spacing_A = size_A / np.asarray(display.shape, dtype=float)
    edges = [
        np.arange(n + 1, dtype=float) * spacing_A[axis]
        for axis, n in enumerate(display.shape)
    ]
    x, y, z = np.meshgrid(*edges, indexing="ij")

    fig = plt.figure(figsize=(10, 8), dpi=180)
    ax = fig.add_subplot(111, projection="3d")
    ax.voxels(x, y, z, grey, facecolors="#bcbfc4", edgecolor="none")
    ax.voxels(x, y, z, blue, facecolors="#205dcc", edgecolor="none")
    ax.set_xlim(0.0, size_A[0])
    ax.set_ylim(0.0, size_A[1])
    ax.set_zlim(0.0, size_A[2])
    ax.set_box_aspect(size_A)
    ax.view_init(elev=24, azim=-58)
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")
    ax.set_title(
        "Ag–Ti₃SiC₂ QSGS bicontinuous model\n"
        "Blue: Ti₃SiC₂ skeleton | Grey: Ag | mass ratio 1:1\n"
        "3-D cutaway; display voxels coarsened from the saved QSGS geometry",
        pad=16,
    )
    ax.grid(False)
    fig.tight_layout()
    output = ROOT / "morphology_3d.png"
    fig.savefig(output, dpi=180, bbox_inches="tight", facecolor="white")
    print(output)


if __name__ == "__main__":
    main()
