#!/usr/bin/env python3
"""绘制稳态导热求解产生的张量、收敛、种子、残差和场图。"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results" / "cases"
FIGURES = HERE / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)


def read_cases() -> list[dict[str, object]]:
    return [json.loads(path.read_text(encoding="utf-8")) for path in sorted(RESULTS.glob("*.json"))
            if path.name != "visualization_export.json" and "case_id" in json.loads(path.read_text(encoding="utf-8"))]


def tensor(record: dict[str, object]) -> np.ndarray:
    return np.asarray(record["effective_conductivity_tensor_W_m_K"], dtype=float)


def mean_k(record: dict[str, object]) -> float:
    return float(np.mean(np.diag(tensor(record))))


def plot_morphology(cases: list[dict[str, object]]) -> None:
    selected = [r for r in cases if r["geometry_kind"] == "qsgs_morphology"
                and r["size_name"] == "ultrafine_500nm"
                and float(r["interface_resistance_m2_K_W"]) == 0.0
                and int(r["grid_n"]) in (150, 200, 300)]
    fig, ax = plt.subplots(figsize=(7.4, 4.5), constrained_layout=True)
    values_by_n = []
    for n in (150, 200, 300):
        values = [mean_k(r) for r in selected if int(r["grid_n"]) == n]
        values_by_n.append(values)
        ax.scatter([n] * len(values), values, s=36, label=f"{n}³ seeds" if n == 150 else None)
    means = [np.mean(v) for v in values_by_n]
    stds = [np.std(v, ddof=1) if len(v) > 1 else 0.0 for v in values_by_n]
    ax.errorbar([150, 200, 300], means, yerr=stds, color="black", capsize=4,
                marker="D", linewidth=1.2, label="mean ± sample SD")
    ax.set_xlabel("Independently regenerated QSGS grid")
    ax.set_ylabel("Mean diagonal effective conductivity (W m⁻¹ K⁻¹)")
    ax.set_title("QSGS morphology comparison at 300 K, ideal contact")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIGURES / "effective_conductivity_qsgs_morphology.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_seeds(cases: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.5), constrained_layout=True)
    colors = {1105: "#1f77b4", 1136: "#e07b39", 1182: "#388e5b"}
    for seed in (1105, 1136, 1182):
        rows = [r for r in cases if r["geometry_kind"] == "qsgs_morphology"
                and r["size_name"] == "ultrafine_500nm"
                and int(r["seed"]) == seed
                and float(r["interface_resistance_m2_K_W"]) == 0.0
                and int(r["grid_n"]) in (150, 200, 300)]
        rows.sort(key=lambda r: int(r["grid_n"]))
        ax.plot([int(r["grid_n"]) for r in rows], [mean_k(r) for r in rows],
                marker="o", color=colors[seed], label=f"seed {seed}")
    ax.set_xticks([150, 200, 300])
    ax.set_xlabel("Independently regenerated QSGS grid")
    ax.set_ylabel("Mean diagonal effective conductivity (W m⁻¹ K⁻¹)")
    ax.set_title("Three random seeds at ideal contact")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIGURES / "effective_conductivity_seed_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_fixed_convergence(cases: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.5), constrained_layout=True)
    for resistance, label, color in ((0.0, "R = 0", "#1f77b4"), (1e-9, "R = 10⁻⁹ m² K/W", "#c34c4c")):
        base = next(r for r in cases if r["case_id"] == f"qsgs_morphology_n150_seed1105_ultrafine_500nm_R{resistance:.0e}")
        fine = next(r for r in cases if r["case_id"] == f"fixed_nested_300_n300_seed1105_ultrafine_500nm_R{resistance:.0e}")
        a, b = tensor(base), tensor(fine)
        rel = 100.0 * np.abs(np.diag(b) - np.diag(a)) / np.diag(a)
        ax.plot([150, 300], [0.0, float(np.mean(rel))], marker="o", linewidth=2,
                color=color, label=f"{label}; mean diagonal Δ={np.mean(rel):.2f}%")
        for i, direction in enumerate("xyz"):
            ax.annotate(f"{direction.upper()} {rel[i]:.2f}%", (300, rel[i]), xytext=(4, i * 9 - 6),
                        textcoords="offset points", fontsize=8, color=color)
    ax.axhline(3.0, color="black", linestyle="--", linewidth=1, label="3% target")
    ax.set_xticks([150, 300])
    ax.set_xlabel("Grid on fixed 150³ voxel geometry")
    ax.set_ylabel("Relative change in diagonal effective conductivity (%)")
    ax.set_title("150³ base geometry repeated 2× per axis; physical size fixed at 10 μm")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIGURES / "fixed_geometry_convergence.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_interface(cases: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.5), constrained_layout=True)
    for size, label, color in (("ultrafine_500nm", "10 μm RVE", "#1f77b4"),
                               ("micro_10um", "200 μm RVE", "#c34c4c")):
        rows = [r for r in cases if r["geometry_kind"] == "qsgs_morphology"
                and int(r["grid_n"]) == 150 and int(r["seed"]) == 1105 and r["size_name"] == size
                and float(r["interface_resistance_m2_K_W"]) in (0.0, 1e-10, 1e-9, 1e-8)]
        rows.sort(key=lambda r: float(r["interface_resistance_m2_K_W"]))
        ax.semilogx([max(float(r["interface_resistance_m2_K_W"]), 1e-11) for r in rows],
                    [mean_k(r) for r in rows], marker="o", label=label, color=color)
    ax.set_xlabel("Numerical interface resistance R (m² K/W); R=0 shown at 10⁻¹¹")
    ax.set_ylabel("Mean diagonal effective conductivity (W m⁻¹ K⁻¹)")
    ax.set_title("Seed 1105 interface-resistance sensitivity; values are not measured data")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIGURES / "interface_resistance_sensitivity.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_residuals(cases: list[dict[str, object]]) -> None:
    record = next(r for r in cases if r["case_id"] == "performance_n150_seed1105_ultrafine_500nm_R0e+00")
    fig, ax = plt.subplots(figsize=(6.7, 4.3), constrained_layout=True)
    for axis, result in enumerate(record["solver_by_imposed_axis"]):
        ax.semilogy(range(len(result["residual_history"])), result["residual_history"],
                    label=f"{'XYZ'[axis]} gradient")
    ax.axhline(1e-9, color="black", linestyle="--", linewidth=1, label="target")
    ax.set_xlabel("PCG iteration")
    ax.set_ylabel("Relative residual")
    ax.set_title("150³ representative QSGS RVE residual convergence")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIGURES / "solver_residual_convergence.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_fields() -> None:
    for axis in "xyz":
        npz = RESULTS / f"performance_n150_seed1105_ultrafine_500nm_R0e+00_axis{'xyz'.index(axis)}_fields.npz"
        with np.load(npz, allow_pickle=False) as archive:
            temp = np.asarray(archive["temperature_K"])
            q = np.sqrt(sum(np.asarray(archive[f"heat_flux_{c}_W_m2"], dtype=float) ** 2 for c in range(3)))
            material = np.asarray(archive["material_id"], dtype=bool)
        index = temp.shape["xyz".index(axis)] // 2
        ts = np.take(temp, index, axis="xyz".index(axis)).T
        qs = np.take(q, index, axis="xyz".index(axis)).T
        ms = np.take(material, index, axis="xyz".index(axis)).T
        extent = [0, 10, 0, 10]
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
        image = axes[0].imshow(ts, origin="lower", extent=extent, cmap="inferno", interpolation="nearest")
        axes[0].contour(ms.astype(float), levels=[0.5], colors=["#39b5ff"], linewidths=0.3,
                        origin="lower", extent=extent)
        fig.colorbar(image, ax=axes[0], label="Temperature (K)")
        axes[0].set_title(f"{axis.upper()} macroscopic gradient: T field")
        image = axes[1].imshow(qs, origin="lower", extent=extent, cmap="magma", interpolation="nearest")
        axes[1].contour(ms.astype(float), levels=[0.5], colors=["#39b5ff"], linewidths=0.3,
                        origin="lower", extent=extent)
        fig.colorbar(image, ax=axes[1], label="|q| (W/m²)")
        axes[1].set_title("Heat-flow magnitude; blue is Ag/Ti₃SiC₂ interface")
        for ax in axes:
            ax.set_xlabel("μm")
            ax.set_ylabel("μm")
            ax.set_aspect("equal")
        fig.suptitle(f"Computed 150³ RVE, seed 1105; 300–310 K")
        fig.savefig(FIGURES / f"temperature_heatflow_axis_{axis}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    cases = read_cases()
    plot_morphology(cases)
    plot_seeds(cases)
    plot_fixed_convergence(cases)
    plot_interface(cases)
    plot_residuals(cases)
    plot_fields()
    print(json.dumps({"status": "PLOTS_CREATED", "files": sorted(p.name for p in FIGURES.glob("*.png"))}, indent=2))


if __name__ == "__main__":
    main()
