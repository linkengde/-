#!/usr/bin/env python3
"""Plot matched 5 Å/10 Å near-surface electron-temperature slices."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ag_ttm_matplotlib_config")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--half-width-A", type=float, default=100.0)
    parser.add_argument("--depth-A", type=float, default=3.0)
    args = parser.parse_args()

    inputs = {
        "10 Å": args.directory / "qsgs_source_gridcompare_3x3_10A_dt0p3125fs_fields.npz",
        "5 Å": args.directory / "qsgs_source_gridcompare_3x3_5A_dt0p3125fs_fields.npz",
    }
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.0), constrained_layout=True)
    image = None
    for col, (resolution, path) in enumerate(inputs.items()):
        with np.load(path, allow_pickle=False) as data:
            ag = data["ag_fraction"].astype(float)
            x = data["x_centers_A"].astype(float)
            y = data["y_centers_A"].astype(float)
            z = data["z_centers_A"].astype(float)
            size = data["size_A"].astype(float)
            dx, dy = data["spacing_A"][:2]
            ix = np.flatnonzero(np.abs(x - size[0] / 2.0) <= args.half_width_A)
            iy = np.flatnonzero(np.abs(y - size[1] / 2.0) <= args.half_width_A)
            iz = int(np.argmin(np.abs(size[2] - z - args.depth_A)))
            x_rel = x[ix] - size[0] / 2.0
            y_rel = y[iy] - size[1] / 2.0
            extent = (x_rel[0] - dx / 2.0, x_rel[-1] + dx / 2.0,
                      y_rel[0] - dy / 2.0, y_rel[-1] + dy / 2.0)
            for row, time_ps in enumerate((0.001, 0.005)):
                label = f"t_{time_ps:.9g}ps_electron_temperature_K"
                temperature = data[label].astype(float)[np.ix_(ix, iy, [iz])][:, :, 0].T
                ag_slice = ag[np.ix_(ix, iy, [iz])][:, :, 0].T
                rise = np.ma.masked_where(ag_slice <= 0.0, temperature - 300.0)
                ax = axes[row, col]
                ax.set_facecolor("#c7c9cc")
                image = ax.imshow(
                    rise,
                    origin="lower",
                    extent=extent,
                    interpolation="nearest",
                    cmap="inferno",
                    vmin=0.0,
                    vmax=700.0,
                    alpha=np.clip(ag_slice, 0.0, 1.0),
                )
                ax.add_patch(Circle((0.0, 0.0), 30.0, fill=False,
                                    edgecolor="white", linestyle="--", linewidth=0.8))
                ax.set_aspect("equal")
                ax.set_xlim(-args.half_width_A, args.half_width_A)
                ax.set_ylim(-args.half_width_A, args.half_width_A)
                ax.set_title(f"{resolution}, t = {time_ps:g} ps; depth ≈ {size[2] - z[iz]:.1f} Å")
                ax.set_xlabel("x from source center (Å)")
                ax.set_ylabel("y from source center (Å)")
    fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.88,
                 label="Ag electron temperature rise, Te − 300 K")
    fig.suptitle("Matched QSGS thermal grids, 3×3 lateral repeat, Case A")
    output = args.output or args.directory / "qsgs_source_gridcompare_3x3_5fs.png"
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
