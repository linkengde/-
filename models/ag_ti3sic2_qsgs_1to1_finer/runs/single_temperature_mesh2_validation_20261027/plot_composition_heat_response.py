#!/usr/bin/env python3
"""Plot saved 0.2 ps Ag 41.7 vs 70 vol% states; performs no time stepping."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
CASES = [
    (
        "Ag 41.6654 vol%",
        PROJECT / "qsgs_morphology.npz",
        PROJECT / "runs/single_temperature_restart_track_20261027/checkpoints/caseA_heat_0p200ps.npz",
    ),
    (
        "Ag 69.9996 vol%",
        HERE / "qsgs_morphology_ag70vol.npz",
        HERE / "ag70vol_melt_window/1x_0p2ps/checkpoints/caseA_ag70_0p2ps_0p2ps.npz",
    ),
]


def main() -> None:
    states = []
    for label, geometry, checkpoint in CASES:
        with np.load(geometry, allow_pickle=False) as archive:
            ag = ~np.asarray(archive["skeleton"], dtype=bool)
            size = np.asarray(archive["size_A"], dtype=float)
        with np.load(checkpoint, allow_pickle=False) as archive:
            temperature = np.asarray(archive["temperature_K"], dtype=float)
            liquid = np.asarray(archive["ag_liquid_fraction"], dtype=float)
            metadata = json.loads(str(archive["metadata_json"].item()))
        if temperature.shape != ag.shape or liquid.shape != ag.shape:
            raise ValueError(f"{label}: checkpoint does not match its native geometry")
        hotspot = np.unravel_index(int(np.argmax(temperature)), temperature.shape)
        states.append((label, ag, temperature, liquid, size, metadata, hotspot))

    z_index = 49
    fig, axes = plt.subplots(3, 2, figsize=(13.0, 13.0), constrained_layout=True)
    for col, (label, ag, temperature, liquid, size, metadata, hotspot) in enumerate(states):
        extent = (0.0, float(size[0]), 0.0, float(size[1]))
        t = float(metadata["time_ps"])
        peak = float(temperature.max())
        mean_ag = float(temperature[ag].mean())
        liquid_fraction = float((liquid * ag).sum() / ag.sum())
        axes[0, col].imshow(
            ag[:, :, z_index].T, origin="lower", extent=extent,
            aspect="equal", cmap="gray", vmin=0, vmax=1, interpolation="nearest",
        )
        axes[0, col].set_title(f"{label} · QSGS Ag mask, z={(z_index + .5) * size[2] / ag.shape[2]:.1f} Å")
        temp_slice = np.ma.masked_where(~ag[:, :, z_index].T, temperature[:, :, z_index].T)
        axes[1, col].set_facecolor("0.18")
        im_t = axes[1, col].imshow(
            temp_slice, origin="lower", extent=extent, aspect="equal",
            cmap="inferno", vmin=300, vmax=1300, interpolation="nearest",
        )
        axes[1, col].plot(size[0] / 2.0, size[1] / 2.0, marker="+", color="cyan", ms=8, mew=1.2)
        axes[1, col].set_title(f"T at 0.2 ps · peak {peak:.1f} K · Ag mean {mean_ag:.2f} K")
        liquid_slice = np.ma.masked_where(~ag[:, :, z_index].T, liquid[:, :, z_index].T)
        axes[2, col].set_facecolor("0.18")
        im_l = axes[2, col].imshow(
            liquid_slice, origin="lower", extent=extent, aspect="equal",
            cmap="magma", vmin=0, vmax=1, interpolation="nearest",
        )
        axes[2, col].set_title(f"Ag liquid fraction · global Ag fraction {liquid_fraction:.6g}")
        for row in range(3):
            axes[row, col].set_xlabel("x (Å)")
            axes[row, col].set_ylabel("y (Å)")

    fig.colorbar(im_t, ax=axes[1, :], label="Temperature (K)", shrink=0.82)
    fig.colorbar(im_l, ax=axes[2, :], label="Liquid fraction", shrink=0.82)
    fig.suptitle("Matched 0.2 ps Case A comparison · same source, properties and boundary; different QSGS geometries")
    output = HERE / "ag41_vs_ag70_heat_response_0p2ps.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(json.dumps({"output": str(output), "z_index": z_index, "grid": list(states[0][1].shape)}, indent=2))


if __name__ == "__main__":
    main()
