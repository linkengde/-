#!/usr/bin/env python3
"""Plot CSV outputs from analyze_npt.py (requires optional matplotlib)."""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {name: np.array([float(row[name]) for row in rows]) for name in rows[0]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("series", type=Path)
    parser.add_argument("rdf", type=Path)
    parser.add_argument("--output", type=Path, default=Path("calibration_1600K_npt.png"))
    args = parser.parse_args()
    s = read_csv(args.series)
    r = read_csv(args.rdf)
    t = s["time_ps"]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), constrained_layout=True)
    axes[0, 0].plot(t, s["temperature_K"], lw=1.2)
    axes[0, 0].set(xlabel="Time (ps)", ylabel="Temperature (K)", title="NPT temperature")
    axes[0, 1].plot(t, s["pressure_bar"], lw=1.0)
    axes[0, 1].axhline(0.0, color="black", lw=0.7)
    axes[0, 1].set(xlabel="Time (ps)", ylabel="Pressure (bar)", title="Instantaneous pressure")
    axes[0, 2].plot(t, s["density_g_cm3"], label="Density")
    axes[0, 2].set(xlabel="Time (ps)", ylabel="Density (g cm$^{-3}$)", title="NPT density")
    ax2 = axes[0, 2].twinx()
    ax2.plot(t, s["volume_A3"], color="tab:orange", alpha=0.65, label="Volume")
    ax2.set_ylabel("Volume (Å$^3$)")
    axes[1, 0].plot(t, s["MSD_A2"], lw=1.2)
    axes[1, 0].set(xlabel="Time (ps)", ylabel="MSD (Å$^2$)", title="Affine-corrected MSD")
    axes[1, 1].plot(t, s["potential_energy_eV"], label="Potential")
    axes[1, 1].plot(t, s["total_energy_eV"], label="Total", alpha=0.8)
    axes[1, 1].legend(frameon=False, fontsize=8)
    axes[1, 1].set(xlabel="Time (ps)", ylabel="Energy (eV)", title="Thermostatted energies")
    axes[1, 2].plot(r["r_A"], r["g_r"], lw=1.2)
    axes[1, 2].set(xlabel="r (Å)", ylabel="g(r)", title="Mean RDF, 5–10 ps", xlim=(0.0, 8.0))
    for ax in axes.flat:
        ax.grid(alpha=0.2)
    fig.suptitle("Pure Ag EAM calibration: 1600 K, zero-pressure NPT, 0.5 fs", fontsize=13)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)


if __name__ == "__main__":
    main()
