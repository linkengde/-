#!/usr/bin/env python3
"""Plot saved 3-D TTM/enthalpy snapshots as comparable XY and X-depth slices."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ag_ttm_matplotlib_config")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="*_fields.npz written by qsgs_ttm_enthalpy.py")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--xy-depth-A", type=float, default=10.0,
                        help="depth below the heated surface for the XY slice")
    args = parser.parse_args()

    with np.load(args.input, allow_pickle=False) as data:
        spacing = data["spacing_A"].astype(float)
        size = data["size_A"].astype(float)
        times = sorted(
            float(key[2:].split("ps_", 1)[0])
            for key in data.files
            if key.startswith("t_") and key.endswith("ps_electron_temperature_K")
        )
        if not times:
            parser.error("Input file contains no saved temperature snapshots")
        snapshots = {}
        for time in times:
            label = f"t_{time:.9g}ps"
            snapshots[time] = {
                "electron": data[f"{label}_electron_temperature_K"].astype(float),
                "lattice": data[f"{label}_lattice_temperature_K"].astype(float),
                "liquid": data[f"{label}_ag_liquid_fraction"].astype(float),
            }

    dz = spacing[2]
    depth_centers = size[2] - (np.arange(snapshots[times[0]]["electron"].shape[2]) + 0.5) * dz
    z_index = int(np.argmin(np.abs(depth_centers - args.xy_depth_A)))
    y_index = snapshots[times[0]]["electron"].shape[1] // 2
    fields = [
        ("Electron temperature, XY", "electron", "xy", "K"),
        ("Electron temperature, X-depth", "electron", "xz", "K"),
        ("Lattice temperature, X-depth", "lattice", "xz", "K"),
        ("Ag liquid fraction, X-depth", "liquid", "xz", "fraction"),
    ]
    fig, axes = plt.subplots(len(fields), len(times), figsize=(4.1 * len(times), 3.2 * len(fields)),
                             constrained_layout=True, squeeze=False)
    cmap = {"electron": "inferno", "lattice": "inferno", "liquid": "Blues"}
    for row, (title, name, plane, unit) in enumerate(fields):
        values = []
        for time in times:
            field = snapshots[time][name]
            if plane == "xy":
                values.append(field[:, :, z_index].T)
            else:
                # Reverse z so the heated surface (depth=0) appears at the top.
                values.append(field[:, y_index, ::-1].T)
        vmin = 0.0 if name == "liquid" else 300.0
        vmax = 1.0 if name == "liquid" else max(301.0, max(float(np.nanmax(value)) for value in values))
        last_image = None
        for col, (time, value) in enumerate(zip(times, values)):
            ax = axes[row, col]
            if plane == "xy":
                extent = (0.0, size[0], 0.0, size[1])
                xlabel, ylabel = "x (Å)", "y (Å)"
            else:
                extent = (0.0, size[0], size[2], 0.0)
                xlabel, ylabel = "x (Å)", "depth below surface (Å)"
            last_image = ax.imshow(value, origin="upper" if plane == "xz" else "lower",
                                   extent=extent, interpolation="nearest", aspect="auto",
                                   vmin=vmin, vmax=vmax, cmap=cmap[name])
            ax.set_title(f"{title}\n{time:g} ps", fontsize=9)
            ax.set_xlabel(xlabel)
            if col == 0:
                ax.set_ylabel(ylabel)
            if row == 0 and plane == "xy":
                ax.text(0.02, 0.98, f"z depth ≈ {depth_centers[z_index]:.1f} Å",
                        transform=ax.transAxes, va="top", ha="left", color="white", fontsize=8)
        fig.colorbar(last_image, ax=axes[row, :].tolist(), shrink=0.82, label=unit,
                     location="right")

    fig.suptitle("QSGS two-temperature/enthalpy fields: 10 Å, Case A, no bottom sink")
    output = args.output or args.input.with_name(f"{args.input.stem}_slices.png")
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
