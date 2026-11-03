#!/usr/bin/env python3
"""Audit mass, pressure jump, and velocity in the OpenFOAM static-drop test."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
CASE = ROOT / "static_drop_case"
NX = NY = 80
DX = DY = 2.0 / NX
SIGMA = 1.0
RADIUS_INITIAL = 0.4


def scalar_field(path: Path) -> np.ndarray:
    text = path.read_text(encoding="ascii")
    match = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not match:
        raise ValueError(f"nonuniform scalar field not found in {path}")
    n = int(match.group(1))
    values = np.fromstring(match.group(2), sep=" ")
    if values.size != n:
        raise ValueError(f"expected {n} values, found {values.size} in {path}")
    return values.reshape(NY, NX)


def vector_field(path: Path) -> np.ndarray:
    text = path.read_text(encoding="ascii")
    match = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not match:
        raise ValueError(f"nonuniform vector field not found in {path}")
    values = np.asarray([
        tuple(map(float, item))
        for item in re.findall(r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)", match.group(2))
    ])
    n = int(match.group(1))
    if values.shape != (n, 3):
        raise ValueError(f"expected {n} vectors, found {values.shape} in {path}")
    return values.reshape(NY, NX, 3)


def main() -> None:
    global CASE, NX, NY, DX, DY
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=80)
    parser.add_argument("--case-name", default=None)
    parser.add_argument("--output-tag", default=None)
    args = parser.parse_args()
    NX = NY = args.resolution
    DX = DY = 2.0 / NX
    default_case = "static_drop_case" if NX == 80 else f"static_drop_case_{NX}"
    CASE = ROOT / (args.case_name or default_case)
    output_tag = args.output_tag or str(NX)
    initial_alpha = scalar_field(CASE / "0/alpha.liquidAg")
    time_dirs = sorted(
        (p for p in CASE.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?", p.name)),
        key=lambda p: float(p.name),
    )
    output_dirs = [p for p in time_dirs if float(p.name) > 0]
    if not output_dirs:
        raise RuntimeError("no OpenFOAM output time directories found")
    final_dir = output_dirs[-1]
    alpha = scalar_field(final_dir / "alpha.liquidAg")
    pressure = scalar_field(final_dir / "p")
    velocity = vector_field(final_dir / "U")
    speed = np.linalg.norm(velocity, axis=2)

    initial_area = float(initial_alpha.sum() * DX * DY)
    final_area = float(alpha.sum() * DX * DY)
    equivalent_radius = float(np.sqrt(final_area / np.pi))
    liquid_bulk = alpha > 0.99
    gas_bulk = alpha < 0.01
    liquid_pressure = float(pressure[liquid_bulk].mean())
    gas_pressure = float(pressure[gas_bulk].mean())
    measured_jump = liquid_pressure - gas_pressure
    expected_jump = SIGMA / equivalent_radius  # 2-D cylindrical interface.
    continuity_matches = re.findall(r"cumulative\s*=\s*([-+0-9.eE]+)", (CASE / "log.interFoam").read_text(encoding="utf-8", errors="replace"))
    transport_text = (CASE / "constant/transportProperties").read_text(encoding="ascii")
    control_text = (CASE / "system/controlDict").read_text(encoding="ascii")
    nu_match = re.search(r"liquidAg\s*\{[^}]*\bnu\s+([-+0-9.eE]+)", transport_text)
    end_match = re.search(r"\bendTime\s+([-+0-9.eE]+)", control_text)

    history = []
    for directory in output_dirs:
        frame_alpha = scalar_field(directory / "alpha.liquidAg")
        frame_pressure = scalar_field(directory / "p")
        frame_velocity = vector_field(directory / "U")
        frame_liquid = frame_alpha > 0.99
        frame_gas = frame_alpha < 0.01
        frame_speed = np.linalg.norm(frame_velocity, axis=2)
        history.append({
            "time": float(directory.name),
            "liquid_area_per_unit_depth": float(frame_alpha.sum() * DX * DY),
            "relative_area_drift_from_initial": float((frame_alpha.sum() * DX * DY - initial_area) / initial_area),
            "bulk_pressure_jump": float(frame_pressure[frame_liquid].mean() - frame_pressure[frame_gas].mean()),
            "max_speed": float(frame_speed.max()),
        })

    x = -1.0 + (np.arange(NX) + 0.5) * DX
    y = -1.0 + (np.arange(NY) + 0.5) * DY
    np.savez_compressed(
        ROOT / f"static_drop_final_fields_{output_tag}.npz",
        alpha_liquid=alpha,
        pressure=pressure,
        velocity=velocity,
        x=x,
        y=y,
        time=np.asarray(float(final_dir.name)),
    )
    weights = np.clip(alpha, 0, 1)
    total = float(weights.sum())
    xbar = float((weights * x[None, :]).sum() / total)
    ybar = float((weights * y[:, None]).sum() / total)
    xc, yc = x[None, :] - xbar, y[:, None] - ybar
    cov = np.array([
        [(weights * xc * xc).sum(), (weights * xc * yc).sum()],
        [(weights * xc * yc).sum(), (weights * yc * yc).sum()],
    ]) / total
    eig = np.linalg.eigvalsh(cov)
    aspect = float(np.sqrt(eig[-1] / eig[0]))

    metrics = {
        "status": "DIMENSIONLESS_STATIC_DROP_VOF_TEST_COMPLETED",
        "solver": "OpenFOAM interFoam 2406",
        "case_name": CASE.name,
        "mesh_cells": NX * NY,
        "grid_shape_xy": [NX, NY],
        "cell_spacing": [DX, DY],
        "kinematic_viscosity_both_phases_dimensionless": float(nu_match.group(1)) if nu_match else None,
        "requested_end_time": float(end_match.group(1)) if end_match else None,
        "initial_liquid_area_per_unit_depth": initial_area,
        "final_liquid_area_per_unit_depth": final_area,
        "relative_liquid_area_drift": (final_area - initial_area) / initial_area,
        "initial_nominal_radius": RADIUS_INITIAL,
        "final_equivalent_radius": equivalent_radius,
        "expected_2d_laplace_pressure_jump": expected_jump,
        "measured_bulk_pressure_jump": measured_jump,
        "relative_pressure_jump_error": (measured_jump - expected_jump) / expected_jump,
        "liquid_bulk_pressure_mean": liquid_pressure,
        "gas_bulk_pressure_mean": gas_pressure,
        "final_max_speed": float(speed.max()),
        "final_rms_speed": float(np.sqrt(np.mean(speed**2))),
        "final_alpha_min": float(alpha.min()),
        "final_alpha_max": float(alpha.max()),
        "final_drop_centroid": [xbar, ybar],
        "final_drop_covariance_aspect_ratio": aspect,
        "final_time": float(final_dir.name),
        "saved_time_history": history,
        "cumulative_continuity_error_last": float(continuity_matches[-1]) if continuity_matches else None,
        "physical_scope": "dimensionless method verification only; no Ag properties, heating, skeleton, arc forcing, evaporation, or splash prediction",
    }
    (ROOT / f"static_drop_metrics_{output_tag}.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    with (ROOT / f"static_drop_time_history_{output_tag}.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=history[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(history)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    xx, yy = np.meshgrid(x, y)
    levels = np.linspace(0, 1, 11)
    im = axes[0].contourf(xx, yy, np.clip(alpha, 0, 1), levels=levels, cmap="Blues")
    axes[0].contour(xx, yy, alpha, levels=(0.5,), colors=("#d1495b",), linewidths=1.5)
    fig.colorbar(im, ax=axes[0], label="liquid volume fraction")
    axes[0].set_title("VOF phase field at final time")
    axes[0].set_xlabel("x (dimensionless)")
    axes[0].set_ylabel("y (dimensionless)")
    skip = (slice(None, None, 5), slice(None, None, 5))
    q = axes[1].quiver(xx[skip], yy[skip], velocity[:, :, 0][skip], velocity[:, :, 1][skip], speed[skip])
    fig.colorbar(q, ax=axes[1], label="speed (dimensionless)")
    axes[1].contour(xx, yy, alpha, levels=(0.5,), colors=("#d1495b",), linewidths=1.5)
    axes[1].set_title("Residual velocity / interface")
    axes[1].set_xlabel("x (dimensionless)")
    axes[1].set_ylabel("y (dimensionless)")
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlim(-0.65, 0.65)
        ax.set_ylim(-0.65, 0.65)
    fig.suptitle("Static circular-drop VOF benchmark (not an Ag or arc simulation)")
    fig.savefig(ROOT / f"static_drop_benchmark_{output_tag}.png", dpi=180, facecolor="white")
    plt.close(fig)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
