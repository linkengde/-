#!/usr/bin/env python3
"""Plot the segmented Ag70 Case B heating history from heat_metrics.csv."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "heat_metrics.csv"
FIG_PATH = ROOT / "caseB_heating_curves.png"


def main() -> None:
    with CSV_PATH.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"No heating metrics in {CSV_PATH}")
    time_ps = np.asarray([float(row["time_ps"]) for row in rows])
    ag_peak = np.asarray([float(row["ag_peak_temperature_K"]) for row in rows])
    ag_mean = np.asarray([float(row["ag_mean_temperature_K"]) for row in rows])
    sc_peak = np.asarray([float(row["ti3sic2_peak_temperature_K"]) for row in rows])
    sc_mean = np.asarray([float(row["ti3sic2_mean_temperature_K"]) for row in rows])
    liquid_pct = 100.0 * np.asarray([float(row["ag_liquid_volume_fraction"]) for row in rows])
    scale = 1.0e15
    ag_sensible = scale * np.asarray([float(row["ag_sensible_enthalpy_J"]) for row in rows])
    ag_latent = scale * np.asarray([float(row["ag_latent_enthalpy_J"]) for row in rows])
    sc_energy = scale * np.asarray([float(row["ti3sic2_stored_enthalpy_J"]) for row in rows])
    total = scale * np.asarray([float(row["system_total_enthalpy_J"]) for row in rows])
    injected = scale * np.asarray([float(row["cumulative_external_source_energy_J"]) for row in rows])

    fig, axes = plt.subplots(3, 1, figsize=(10.5, 10.0), sharex=True)
    axes[0].plot(time_ps, ag_peak, color="#b7352d", label="Ag peak")
    axes[0].plot(time_ps, ag_mean, color="#d2746b", linestyle="--", label="Ag mean")
    axes[0].plot(time_ps, sc_peak, color="#245a83", label="Ti₃SiC₂ peak")
    axes[0].plot(time_ps, sc_mean, color="#5d94bb", linestyle="--", label="Ti₃SiC₂ mean")
    axes[0].axhline(1235.08, color="#555555", linestyle=":", label="Ag melting point")
    axes[0].axhline(1473.15, color="#805da3", linestyle="-.", label="Ti₃SiC₂ data ceiling")
    axes[0].set_ylabel("Temperature (K)")
    axes[0].legend(ncol=3, fontsize=8)
    axes[0].grid(alpha=0.25)

    axes[1].plot(time_ps, liquid_pct, color="#e18b00")
    axes[1].set_ylabel("Ag liquid volume (%)")
    axes[1].set_ylim(bottom=0.0)
    axes[1].grid(alpha=0.25)

    axes[2].plot(time_ps, ag_sensible, label="Ag sensible", color="#b7352d")
    axes[2].plot(time_ps, ag_latent, label="Ag latent", color="#e18b00")
    axes[2].plot(time_ps, sc_energy, label="Ti₃SiC₂ stored", color="#245a83")
    axes[2].plot(time_ps, total, label="System total H", color="#333333")
    axes[2].plot(time_ps, injected, label="Input energy", color="#48915b", linestyle="--")
    axes[2].set_ylabel("Energy (10⁻¹⁵ J)")
    axes[2].set_xlabel("Simulation time (ps)")
    axes[2].legend(ncol=3, fontsize=8)
    axes[2].grid(alpha=0.25)

    fig.suptitle("Ag70 1× Case B: segmented heating, ideal contact, insulated outer boundary")
    fig.tight_layout()
    fig.savefig(FIG_PATH, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(FIG_PATH)


if __name__ == "__main__":
    main()
