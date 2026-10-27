#!/usr/bin/env python3
"""Render XY/XZ temperature, phase mask and liquid fraction for a checkpoint."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("geometry", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with np.load(args.checkpoint, allow_pickle=False) as data:
        temp = np.asarray(data["temperature_K"])
        liquid = np.asarray(data["ag_liquid_fraction"])
        meta = json.loads(str(data["metadata_json"].item()))
    with np.load(args.geometry, allow_pickle=False) as data:
        ag = ~np.asarray(data["skeleton"], dtype=bool)
        size = np.asarray(data["size_A"], dtype=float)
    refine = temp.shape[0] // ag.shape[0]
    if temp.shape != tuple(np.asarray(ag.shape) * refine):
        raise ValueError("Checkpoint shape is not an integer native-grid refinement")
    if refine > 1:
        ag = np.repeat(np.repeat(np.repeat(ag, refine, axis=0), refine, axis=1), refine, axis=2)
    spacing = size / np.asarray(ag.shape)
    hotspot = np.unravel_index(np.argmax(temp), temp.shape)
    z_index = hotspot[2]
    y_index = hotspot[1]
    time_ps = float(meta["time_ps"])
    peak = float(temp.max())
    liquid_vol = float(np.sum(ag * liquid) / ag.sum())
    x_extent = (0.0, float(size[0]))
    y_extent = (0.0, float(size[1]))
    z_extent = (0.0, float(size[2]))

    fig, ax = plt.subplots(2, 2, figsize=(12.4, 8.6), constrained_layout=True)
    im0 = ax[0, 0].imshow(ag[:, :, z_index].T, origin="lower", extent=(*x_extent, *y_extent), aspect="equal", cmap="gray_r", interpolation="nearest")
    ax[0, 0].set_title(f"Ag geometry, XY at z={((z_index + .5) * spacing[2]):.1f} Å")
    ax[0, 0].set_xlabel("x (Å)"); ax[0, 0].set_ylabel("y (Å)")
    im1 = ax[0, 1].imshow(temp[:, :, z_index].T, origin="lower", extent=(*x_extent, *y_extent), aspect="equal", cmap="inferno", interpolation="nearest")
    if np.any(liquid[:, :, z_index] >= 0.5):
        ax[0, 1].contour((np.arange(temp.shape[0]) + .5) * spacing[0], (np.arange(temp.shape[1]) + .5) * spacing[1], (liquid[:, :, z_index] >= .5).T, levels=[.5], colors=["cyan"], linewidths=0.8)
    ax[0, 1].set_title("Temperature, XY through hotspot")
    ax[0, 1].set_xlabel("x (Å)"); ax[0, 1].set_ylabel("y (Å)")
    fig.colorbar(im1, ax=ax[0, 1], label="T (K)")
    im2 = ax[1, 0].imshow(temp[:, y_index, :].T, origin="lower", extent=(*x_extent, *z_extent), aspect="auto", cmap="inferno", interpolation="nearest")
    ax[1, 0].set_title(f"Temperature, XZ at y={((y_index + .5) * spacing[1]):.1f} Å")
    ax[1, 0].set_xlabel("x (Å)"); ax[1, 0].set_ylabel("z (Å)")
    fig.colorbar(im2, ax=ax[1, 0], label="T (K)")
    im3 = ax[1, 1].imshow(liquid[:, y_index, :].T, origin="lower", extent=(*x_extent, *z_extent), aspect="auto", cmap="magma", vmin=0.0, vmax=1.0, interpolation="nearest")
    ax[1, 1].set_title("Ag liquid fraction, same XZ slice")
    ax[1, 1].set_xlabel("x (Å)"); ax[1, 1].set_ylabel("z (Å)")
    fig.colorbar(im3, ax=ax[1, 1], label="liquid fraction")
    fig.suptitle(f"Ag70 Case A · t={time_ps:.4f} ps · peak={peak:.1f} K · Ag liquid vol. fraction={liquid_vol:.5g}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    print(json.dumps({"output": str(args.output), "time_ps": time_ps, "peak_temperature_K": peak, "ag_liquid_volume_fraction": liquid_vol, "grid_shape_xyz": list(temp.shape)}, indent=2))


if __name__ == "__main__":
    main()
