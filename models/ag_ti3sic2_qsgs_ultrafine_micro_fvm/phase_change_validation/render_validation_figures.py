#!/usr/bin/env python3
"""由实际验证CSV/NPZ重绘Stefan、时间步、能量及三维结构结果图。"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/fvm-phase-mplconfig")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np

from material_model import CONFIG


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
OUT = HERE / "results"
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def draw_geometry() -> None:
    with np.load(PROJECT / CONFIG["geometry_source"]) as data:
        sc = np.asarray(data["skeleton"], dtype=bool)
    stride = max(1, sc.shape[0] // 40)
    volume = sc[::stride, ::stride, ::stride]
    fig = plt.figure(figsize=(8.4, 7.0))
    ax = fig.add_subplot(111, projection="3d")
    face = np.empty(volume.shape, dtype=object)
    face[:] = "#12a89d"
    ax.voxels(volume, facecolors=face, edgecolor=(0.05, 0.2, 0.2, 0.10), linewidth=0.12)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=24, azim=36)
    ax.set_xlabel("x (coarse voxels)")
    ax.set_ylabel("y (coarse voxels)")
    ax.set_zlabel("z (coarse voxels)")
    ax.set_title("Ti₃SiC₂ skeleton, QSGS seed 1105 (volume rendering)", pad=22)
    fig.subplots_adjust(left=0.02, right=0.96, bottom=0.02, top=0.90)
    fig.savefig(FIG / "initial_qsgs_three_dimensional.png", dpi=200)
    plt.close(fig)


def draw_stefan() -> None:
    rows = read_csv(OUT / "STEFAN_FRONT_CONVERGENCE.csv")
    grid_rows = [r for r in rows if float(r["dt_s"]) == 2e-6 and int(r["n"]) in (50, 100, 200)]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4), constrained_layout=True)
    styles = {"melt": ("#b21f43", "o", "Melting"), "freeze": ("#1e6b9a", "s", "Freezing")}
    for kind, (color, marker, label) in styles.items():
        subset = sorted([r for r in grid_rows if r["type"] == kind], key=lambda x: int(x["n"]))
        n = np.array([int(r["n"]) for r in subset])
        measured = np.array([float(r["measured_front_m"]) for r in subset]) * 1e6
        exact = float(subset[0]["stefan_front_m"]) * 1e6
        axes[0].plot(n, measured, marker=marker, color=color, lw=1.7, label=f"{label}, FVM")
        axes[0].axhline(exact, color=color, ls="--", lw=1.1, alpha=0.75)
        err = 100.0 * np.array([float(r["front_relative_error"]) for r in subset])
        axes[1].plot(n, err, marker=marker, color=color, lw=1.7, label=label)
    axes[0].set_xlabel("Cells along 0.5 mm domain")
    axes[0].set_ylabel("Front position at 1 ms (µm)")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    axes[1].axhline(2.0, color="#555555", ls=":", label="2% target")
    axes[1].set_xlabel("Cells along domain")
    axes[1].set_ylabel("Absolute front error (%)")
    axes[1].set_yscale("log")
    axes[1].legend()
    axes[1].grid(alpha=0.25, which="both")
    fig.savefig(FIG / "stefan_front_validation.png", dpi=220)
    plt.close(fig)


def draw_time_step() -> None:
    rows = read_csv(OUT / "TIME_STEP_CONVERGENCE.csv")
    rows = sorted(rows, key=lambda r: float(r["dt_factor"]))
    factors = np.array([float(r["dt_factor"]) for r in rows])
    liquid = np.array([float(r["heated_Ag_bulk_liquid_fraction"]) for r in rows]) * 1e4
    tmax = np.array([float(r["heated_max_temperature_K"]) for r in rows])
    fig, ax = plt.subplots(figsize=(7.2, 4.5), constrained_layout=True)
    ax.plot(factors, liquid, marker="o", color="#b21f43", lw=2, label="Ag liquid fraction × 10⁴")
    ax.set_xlabel("Nominal time-step factor (smaller is finer)")
    ax.set_ylabel("Ag bulk liquid fraction × 10⁴", color="#b21f43")
    ax.tick_params(axis="y", labelcolor="#b21f43")
    ax.set_xscale("log", base=2)
    ax.invert_xaxis()
    ax.grid(alpha=0.25)
    ax2 = ax.twinx()
    ax2.plot(factors, tmax, marker="s", color="#e18b22", lw=1.7, label="Maximum temperature")
    ax2.set_ylabel("Maximum temperature (K)", color="#a76412")
    ax2.tick_params(axis="y", labelcolor="#a76412")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], loc="best")
    ax.set_title("50³ seed 1105 QSGS transient time-step check")
    fig.savefig(FIG / "time_step_convergence.png", dpi=220)
    plt.close(fig)


def draw_energy_history() -> None:
    rows = read_csv(OUT / "bicontinuous_qsgs_seed1105_n50_dt025_history.csv")
    time_us = np.array([float(r["physical_time_s"]) for r in rows]) * 1e6
    source = np.array([float(r["cumulative_source_energy_J"]) for r in rows]) * 1e9
    out = np.array([float(r["cumulative_boundary_heat_out_J"]) for r in rows]) * 1e9
    stored = np.array([float(r["stored_enthalpy_change_J"]) for r in rows]) * 1e9
    residual = np.array([float(r["energy_residual_J"]) for r in rows]) * 1e9
    fig, ax = plt.subplots(figsize=(7.6, 4.6), constrained_layout=True)
    ax.plot(time_us, source, lw=1.8, label="Cumulative source input")
    ax.plot(time_us, out, lw=1.8, label="Cumulative heat out")
    ax.plot(time_us, stored, lw=1.8, label="Stored enthalpy change")
    ax.plot(time_us, residual, ls=":", lw=1.2, label="Energy residual")
    ax.set_xlabel("Physical time (µs)")
    ax.set_ylabel("Energy (nJ)")
    ax.grid(alpha=0.25)
    ax.legend()
    ax.set_title("50³ seed 1105 QSGS energy closure")
    fig.savefig(FIG / "energy_balance_time_history.png", dpi=220)
    plt.close(fig)


def draw_3d_hotspot() -> None:
    field = OUT / "bicontinuous_qsgs_seed1105_quarterdt_n50_heated.npz"
    with np.load(field) as data:
        t = data["temperature_K"].astype(float)
        ag = data["ag_fraction"].astype(float)
        liquid = data["ag_liquid_fraction"].astype(float)
    stride = 1
    t, ag, liquid = t[::stride, ::stride, ::stride], ag[::stride, ::stride, ::stride], liquid[::stride, ::stride, ::stride]
    hot = t > 700.0
    phi = (ag > 0.5) & (liquid >= 0.5)
    norm = colors.Normalize(vmin=700.0, vmax=float(t[hot].max()) if np.any(hot) else 1300.0)
    face = cm.inferno(norm(t))
    fig = plt.figure(figsize=(10.0, 4.8))
    ax = fig.add_subplot(121, projection="3d")
    ax.voxels(hot, facecolors=face, edgecolor=(0.2, 0.1, 0.1, 0.18), linewidth=0.12)
    if np.any(phi):
        points = np.where(phi)
        ax.scatter(*points, color="#00c9e0", s=20, label="Ag cells, liquid fraction ≥ 0.5")
        ax.legend(loc="upper left", fontsize=8)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=24, azim=36)
    ax.set_title("Temperature > 700 K")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    bx = fig.add_subplot(122, projection="3d")
    if np.any(phi):
        points = np.where(phi)
        bx.scatter(*points, c=liquid[phi], cmap="magma", vmin=0.0, vmax=1.0, s=35)
    else:
        bx.text(0.5, 0.5, 0.5, "No cells above liquid threshold", ha="center")
    bx.set_box_aspect((1, 1, 1)); bx.view_init(elev=24, azim=36)
    bx.set_title("Connected liquid-bearing Ag cells")
    bx.set_xlabel("x"); bx.set_ylabel("y"); bx.set_zlabel("z")
    fig.suptitle("QSGS seed 1105, 50³, local numerical heating", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "qsgs_50_heated_3d_phase_volume.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    draw_geometry()
    draw_stefan()
    draw_time_step()
    draw_energy_history()
    draw_3d_hotspot()
    print(json.dumps({"figures": sorted(p.name for p in FIG.glob("*.png"))}, ensure_ascii=False))


if __name__ == "__main__":
    main()
