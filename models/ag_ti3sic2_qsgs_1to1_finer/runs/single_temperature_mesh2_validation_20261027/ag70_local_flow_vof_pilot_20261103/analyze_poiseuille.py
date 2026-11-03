#!/usr/bin/env python3
"""Compare a pressure-driven channel-flow result with planar Poiseuille flow."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
CASE = ROOT / "poiseuille_case"
NX, NY = 120, 40
LENGTH, HEIGHT, NU, DP = 12.0, 1.0, 0.01, 1.0


def read_vector(path: Path) -> np.ndarray:
    text = path.read_text(encoding="ascii")
    match = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not match:
        raise ValueError(f"nonuniform vector field not found in {path}")
    values = np.asarray([
        tuple(map(float, item))
        for item in re.findall(r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)", match.group(2))
    ])
    expected = int(match.group(1))
    if values.shape != (expected, 3) or expected != NX * NY:
        raise ValueError(f"expected {NX * NY} vectors; got {values.shape} in {path}")
    return values.reshape(NY, NX, 3)


def main() -> None:
    times = sorted(
        (p for p in CASE.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?", p.name) and float(p.name) > 0),
        key=lambda p: float(p.name),
    )
    if not times:
        raise RuntimeError("no solver output time directories found")
    final = times[-1]
    velocity = read_vector(final / "U")
    x = (np.arange(NX) + 0.5) * LENGTH / NX
    y = (np.arange(NY) + 0.5) * HEIGHT / NY
    center_band = (x >= 0.4 * LENGTH) & (x <= 0.6 * LENGTH)
    simulated = velocity[:, center_band, 0].mean(axis=1)
    gradient = DP / LENGTH
    analytic = gradient * y * (HEIGHT - y) / (2 * NU)
    u_max_exact = gradient * HEIGHT**2 / (8 * NU)
    flow_exact = gradient * HEIGHT**3 / (12 * NU)
    flow_simulated = float(simulated.sum() * HEIGHT / NY)
    profile_rel_l2 = float(np.sqrt(np.mean((simulated - analytic) ** 2)) / u_max_exact)
    previous = read_vector(times[-2] / "U") if len(times) > 1 else velocity
    previous_profile = previous[:, center_band, 0].mean(axis=1)
    tail_drift = float(np.max(np.abs(simulated - previous_profile)) / u_max_exact)

    log_text = (CASE / "log.simpleFoam").read_text(encoding="utf-8", errors="replace")
    continuity = re.findall(r"cumulative\s*=\s*([-+0-9.eE]+)", log_text)
    metrics = {
        "status": "DIMENSIONLESS_POISEUILLE_BENCHMARK_COMPLETED",
        "solver": "OpenFOAM simpleFoam 2406; laminar, single-phase",
        "grid": [NX, NY, 1],
        "domain_length": LENGTH,
        "channel_height": HEIGHT,
        "kinematic_viscosity_dimensionless": NU,
        "imposed_pressure_drop_dimensionless": DP,
        "sampled_x_band": [float(x[center_band][0]), float(x[center_band][-1])],
        "final_iteration": float(final.name),
        "analytic_centerline_max_velocity": u_max_exact,
        "simulated_centerline_max_velocity": float(simulated.max()),
        "max_velocity_relative_error": float((simulated.max() - u_max_exact) / u_max_exact),
        "velocity_profile_relative_l2_error": profile_rel_l2,
        "analytic_volume_flux_per_unit_span": flow_exact,
        "simulated_volume_flux_per_unit_span": flow_simulated,
        "volume_flux_relative_error": float((flow_simulated - flow_exact) / flow_exact),
        "profile_symmetry_error_normalized": float(np.max(np.abs(simulated - simulated[::-1])) / u_max_exact),
        "max_profile_change_over_last_write_interval_normalized": tail_drift,
        "cumulative_continuity_error_last": float(continuity[-1]) if continuity else None,
        "physical_scope": "dimensionless solver verification only; imposed pressure drop is not an arc load and result is not an Ag flow prediction",
    }
    (ROOT / "poiseuille_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    with (ROOT / "poiseuille_profile.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["y", "u_simulated", "u_analytic"])
        writer.writerows(zip(y, simulated, analytic))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.3, 4.8), constrained_layout=True)
    ax.plot(analytic, y, color="#222222", linewidth=2, label="analytic Poiseuille")
    ax.plot(simulated, y, "o", color="#e87500", markersize=3.5, label="OpenFOAM center-band mean")
    ax.set_xlabel("streamwise velocity (dimensionless)")
    ax.set_ylabel("wall-normal coordinate (dimensionless)")
    ax.set_title("Pressure-driven channel verification (not Ag / arc flow)")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(ROOT / "poiseuille_profile.png", dpi=180, facecolor="white")
    plt.close(fig)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
