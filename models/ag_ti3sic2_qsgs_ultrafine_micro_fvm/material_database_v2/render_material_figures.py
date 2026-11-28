#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从v2物性表与实际50³结果生成论文级核查图。"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
FIG = HERE / "figures"
RESULTS = HERE / "results" / "same_load_50"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(PROJECT))
from material_database_v2.material_model import PankratzMaterials


def read_csv(name: str) -> list[dict[str, str]]:
    with (HERE / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def plot_silver_conductivity():
    solid = read_csv("Ag_solid_properties.csv")
    liquid = read_csv("Ag_liquid_properties.csv")
    ts = np.array([float(row["temperature_K"]) for row in solid])
    ks = np.array([float(row["thermal_conductivity_W_m_K"]) for row in solid])
    tl = np.array([float(row["temperature_K"]) for row in liquid])
    kl = np.array([float(row["thermal_conductivity_W_m_K"]) for row in liquid])
    fig, ax = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    ax.plot(ts, ks, "o-", label="Solid Ag, Ho recommended values", color="#1668a8")
    ax.plot(tl, kl, "s--", label="Liquid Ag, Ho provisional estimates (≈20% uncertainty)", color="#c74b35")
    ax.axvline(1235.08, color="black", lw=1, ls=":", label="Pankratz Tm = 1235.08 K")
    ax.scatter([1235.08], [358.0], color="#1668a8", s=55, zorder=4)
    ax.scatter([1235.08], [175.0], color="#c74b35", s=55, zorder=4)
    ax.set(xlabel="Temperature (K)", ylabel="Thermal conductivity (W m$^{-1}$ K$^{-1}$)",
           title="Silver thermal conductivity: separate solid and liquid datasets")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "Ag_thermal_conductivity_v2.png", dpi=220)
    plt.close(fig)


def plot_silver_heat_capacity():
    mat = PankratzMaterials(1.0)
    t_s = np.linspace(300.0, mat.Tm_K, 240)
    cp_s = np.asarray(mat.ag_cp_solid_molar_J_mol_K(t_s)) / mat.ag_M_kg_mol
    fig, ax = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    ax.plot(t_s, cp_s, color="#1668a8", lw=2, label="Pankratz solid Cp(T), main set")
    ax.hlines(mat.ag_cp_liquid_J_kg_K, mat.Tm_K, 2000.0, color="#1668a8", lw=2,
              label="Pankratz liquid Cp, 310.30 J kg$^{-1}$ K$^{-1}$")
    ax.hlines(314.6101, 1243.0, 1588.0, color="#d07b17", lw=2, ls="--",
              label="Feber liquid fit, 314.61 (1243–1588 K)")
    ax.hlines(259.8, 1234.0, 2000.0, color="#7a5aa6", lw=2, ls="-.",
              label="Cagran liquid fit, 259.8 (1234–2000 K)")
    ax.hlines(165.6, 1000.0, 1234.0, color="#4d8b57", lw=2, ls=":",
              label="Cagran solid high-T fit, 165.6 (1000–1234 K)")
    ax.axvline(mat.Tm_K, color="black", lw=1, ls=":")
    ax.set(xlim=(300, 2000), ylim=(0, 380), xlabel="Temperature (K)",
           ylabel="Specific heat capacity (J kg$^{-1}$ K$^{-1}$)",
           title="Silver heat capacity sources and stated temperature ranges")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.savefig(FIG / "Ag_heat_capacity_literature_v2.png", dpi=220)
    plt.close(fig)


def plot_silver_enthalpy():
    mat = PankratzMaterials(1.0)
    temperatures = np.linspace(300.0, 1500.0, 1200)
    enthalpy = mat.ag_specific_enthalpy_J_kg(temperatures) / 1000.0
    fig, ax = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    ax.plot(temperatures, enthalpy, color="#a83e31", lw=2)
    ax.axvspan(mat.Ts_K, mat.Tl_K, color="#f0c75e", alpha=0.35,
               label="1 K numerical mushy interval; full latent heat once")
    ax.axvline(mat.Tm_K, color="black", lw=1, ls=":", label="Pankratz Tm")
    ax.set(xlabel="Temperature (K)", ylabel="Specific enthalpy from 300 K (kJ kg$^{-1}$)",
           title="Pankratz Ag enthalpy–temperature relation")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "Ag_enthalpy_temperature_v2.png", dpi=220)
    plt.close(fig)


def plot_tisc_properties():
    rows = read_csv("Ti3SiC2_solid_properties.csv")
    t = np.array([float(row["temperature_K"]) for row in rows])
    cp = np.array([float(row["Cp_mass_J_kg_K"]) for row in rows])
    k = np.array([float(row["thermal_conductivity_W_m_K"]) for row in rows])
    valid = t <= 1473.15
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.6), constrained_layout=True)
    axes[0].plot(t[valid], cp[valid], "o-", color="#1668a8", label="Solver range / fit data")
    axes[0].plot(t[~valid], cp[~valid], "o--", color="#999999", label="Extrapolation only")
    axes[0].axvline(1473.15, color="black", ls=":", lw=1)
    axes[0].set(xlabel="Temperature (K)", ylabel="Cp (J kg$^{-1}$ K$^{-1}$)", title="Ti$_3$SiC$_2$ heat capacity")
    axes[1].plot(t[valid], k[valid], "o-", color="#4d8b57", label="Solver range / fit data")
    axes[1].plot(t[~valid], k[~valid], "o--", color="#999999", label="Extrapolation only")
    axes[1].axvline(1473.15, color="black", ls=":", lw=1)
    axes[1].set(xlabel="Temperature (K)", ylabel="k (W m$^{-1}$ K$^{-1}$)", title="Ti$_3$SiC$_2$ thermal conductivity")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Barsoum et al. (1999); dashed points are extrapolations, not measurements")
    fig.savefig(FIG / "Ti3SiC2_thermal_properties_v2.png", dpi=220)
    plt.close(fig)


def plot_same_load_timeseries():
    summaries = ["old_candidate", "pankratz_v2"]
    labels = {"old_candidate": "Old candidate", "pankratz_v2": "Pankratz v2"}
    colors = {"old_candidate": "#3b6fb6", "pankratz_v2": "#d04a36"}
    fig, axes = plt.subplots(3, 1, figsize=(8.5, 10.4), constrained_layout=True)
    for tag in summaries:
        with (RESULTS / f"{tag}_step_history.csv").open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        time_ns = np.array([float(row["time_s"]) * 1e9 for row in rows])
        phase = np.array([row["phase"] for row in rows])
        max_t = np.array([float(row["max_temperature_K"]) for row in rows])
        bulk = np.array([float(row["ag_bulk_liquid_fraction"]) for row in rows])
        norm = np.array([float(row["normalized_energy_residual"]) for row in rows])
        axes[0].semilogx(time_ns, max_t, color=colors[tag], lw=1.8, label=labels[tag])
        axes[1].semilogx(time_ns, np.where(bulk > 0, bulk, np.nan), color=colors[tag], lw=1.8, label=labels[tag])
        axes[2].semilogx(time_ns, norm, color=colors[tag], lw=1.6, label=labels[tag])
    axes[0].axhline(1235.08, color="black", ls=":", lw=1)
    axes[0].set(ylabel="Maximum temperature (K)", title="Same 50³ QGS numerical load; log physical-time axis")
    axes[1].set(ylabel="Ag bulk liquid fraction", title="Ag-weighted liquid fraction (zeros omitted)")
    axes[2].set(xlabel="Physical time (ns)", ylabel="Normalized cumulative energy residual", title="Cumulative energy balance")
    for ax in axes:
        ax.grid(alpha=0.25, which="both")
        ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "same_load_old_vs_pankratz_50_v2.png", dpi=220)
    plt.close(fig)


def plot_heated_endpoint():
    for tag, title, color in (("old_candidate", "Old candidate", "#e5a443"),
                              ("pankratz_v2", "Pankratz v2", "#4e83c4")):
        data = np.load(RESULTS / f"{tag}_n50_heated_state.npz", allow_pickle=False)
        t = data["temperature_K"]
        liquid = data["liquid_fraction"]
        spacing = data["spacing_m"]
        n = t.shape[0]
        edge_um = n * spacing[0] * 1e6
        extent = (-edge_um / 2, edge_um / 2, -edge_um / 2, edge_um / 2)
        fig = plt.figure(figsize=(12.2, 5.5))
        grid = fig.add_gridspec(1, 2, width_ratios=(1.0, 1.15),
                                left=0.07, right=0.91, bottom=0.12, top=0.84,
                                wspace=0.28)
        ax0 = fig.add_subplot(grid[0, 0])
        ax1 = fig.add_subplot(grid[0, 1], projection="3d")
        z = n // 2
        im0 = ax0.imshow(t[:, :, z].T, origin="lower", extent=extent, cmap="inferno",
                             vmin=300, vmax=float(t.max()), interpolation="nearest")
        ax0.contour(np.linspace(-edge_um/2+spacing[0]*0.5e6, edge_um/2-spacing[0]*0.5e6, n),
                    np.linspace(-edge_um/2+spacing[1]*0.5e6, edge_um/2-spacing[1]*0.5e6, n),
                        liquid[:, :, z].T, levels=[0.5], colors=["#4ee8ff"], linewidths=1.2)
        fig.colorbar(im0, ax=ax0, shrink=0.86, pad=0.035, label="Temperature (K)")
        ax0.set(title="Midplane temperature", xlabel="x (μm)", ylabel="y (μm)")
        ax0.set_aspect("equal")
        ag = data["ag_fraction"]
        coords = np.where((t >= 1225.0) & (liquid > 0.0) & (ag > 1.0e-14))
        if len(coords[0]):
            xyz = [(axis + 0.5 - n / 2) * spacing[i] * 1e6 for i, axis in enumerate(coords)]
            sc = ax1.scatter(*xyz, c=liquid[coords], cmap="cool", vmin=0, vmax=1,
                                 s=18, alpha=0.9, depthshade=False)
            fig.colorbar(sc, ax=ax1, shrink=0.78, pad=0.10, label="Ag liquid fraction")
        ax1.set(xlim=(-edge_um/2, edge_um/2), ylim=(-edge_um/2, edge_um/2),
                    zlim=(-edge_um/2, edge_um/2), xlabel="x (μm)", ylabel="y (μm)", zlabel="z (μm)",
                    title="3-D hot Ag cells (T ≥ 1225 K)")
        ax1.set_box_aspect((1, 1, 1))
        ax1.view_init(elev=23, azim=38)
        ax0.tick_params(labelsize=8)
        ax1.tick_params(labelsize=8)
        fig.suptitle(f"{title} · 50³ · 10 μm RVE · heat end 1 ns · numerical validation pulse",
                     y=0.96, fontsize=12)
        fig.savefig(FIG / f"{tag}_heated_endpoint_3d_slice.png", dpi=220)
        plt.close(fig)


def plot_n100_heated_endpoint():
    data_path = HERE / "results" / "pankratz_v2_n100_dt0.25" / "heated_state.npz"
    if not data_path.exists():
        return
    data = np.load(data_path, allow_pickle=False)
    temperature = data["temperature_K"]
    liquid = data["liquid_fraction"]
    ag = data["ag_fraction"]
    spacing = data["spacing_m"]
    n = temperature.shape[0]
    edge_um = n * spacing[0] * 1e6
    centers = (np.arange(n) + 0.5 - n / 2.0) * spacing[0] * 1e6
    extent = (-edge_um / 2, edge_um / 2, -edge_um / 2, edge_um / 2)
    fig = plt.figure(figsize=(12.2, 5.5))
    grid = fig.add_gridspec(1, 2, width_ratios=(1.0, 1.15),
                            left=0.07, right=0.91, bottom=0.12, top=0.84, wspace=0.28)
    ax0 = fig.add_subplot(grid[0, 0])
    ax1 = fig.add_subplot(grid[0, 1], projection="3d")
    z = n // 2
    im = ax0.imshow(temperature[:, :, z].T, origin="lower", extent=extent,
                    cmap="inferno", vmin=300.0, vmax=float(temperature.max()),
                    interpolation="nearest")
    ax0.contour(centers, centers, liquid[:, :, z].T, levels=[0.5],
                colors=["#4ee8ff"], linewidths=1.0)
    fig.colorbar(im, ax=ax0, shrink=0.86, pad=0.035, label="Temperature (K)")
    ax0.set(title="100³ midplane temperature", xlabel="x (μm)", ylabel="y (μm)")
    ax0.set_aspect("equal")
    hot = (temperature >= 1225.0) & (liquid > 0.0) & (ag > 1e-14)
    coords = np.where(hot)
    if len(coords[0]):
        xyz = [(axis + 0.5 - n / 2) * spacing[i] * 1e6 for i, axis in enumerate(coords)]
        points = ax1.scatter(*xyz, c=liquid[coords], cmap="cool", vmin=0, vmax=1,
                             s=4, alpha=0.9, depthshade=False)
        fig.colorbar(points, ax=ax1, shrink=0.78, pad=0.10, label="Ag liquid fraction")
    ax1.set(xlim=(-edge_um/2, edge_um/2), ylim=(-edge_um/2, edge_um/2),
            zlim=(-edge_um/2, edge_um/2), xlabel="x (μm)", ylabel="y (μm)",
            zlabel="z (μm)", title="3-D hot Ag cells (T ≥ 1225 K)")
    ax1.set_box_aspect((1, 1, 1))
    ax1.view_init(elev=23, azim=38)
    fig.suptitle("Pankratz v2 · 100³ · 10 μm RVE · heat end 1 ns · numerical validation pulse",
                 y=0.96, fontsize=12)
    fig.savefig(FIG / "pankratz_v2_n100_heated_endpoint_3d_slice.png", dpi=220)
    plt.close(fig)


def plot_n100_first_freeze():
    path = (HERE / "results" / "pankratz_v2_n100_freezing_event_cap_10ps_freeze_event_only"
            / "first_frozen_state.npz")
    if not path.exists():
        return
    with np.load(path, allow_pickle=False) as data:
        temperature = data["temperature_K"]
        liquid = data["liquid_fraction"]
        ag = data["ag_fraction"]
        spacing = data["spacing_m"]
        time_ns = float(data["time_s"]) * 1e9
    n = temperature.shape[0]
    edge_um = n * spacing[0] * 1e6
    centers = (np.arange(n) + 0.5 - n / 2.0) * spacing[0] * 1e6
    extent = (-edge_um / 2, edge_um / 2, -edge_um / 2, edge_um / 2)
    fig = plt.figure(figsize=(12.2, 5.5))
    grid = fig.add_gridspec(1, 2, width_ratios=(1.0, 1.15),
                            left=0.07, right=0.91, bottom=0.12, top=0.84, wspace=0.28)
    ax0 = fig.add_subplot(grid[0, 0])
    ax1 = fig.add_subplot(grid[0, 1], projection="3d")
    z = n // 2
    im = ax0.imshow(temperature[:, :, z].T, origin="lower", extent=extent,
                    cmap="inferno", vmin=300.0, vmax=float(temperature.max()),
                    interpolation="nearest")
    fig.colorbar(im, ax=ax0, shrink=0.86, pad=0.035, label="Temperature (K)")
    ax0.set(title="Midplane temperature at complete freeze", xlabel="x (μm)", ylabel="y (μm)")
    ax0.set_aspect("equal")
    hot_solid = (temperature >= 1000.0) & (liquid <= 1e-8) & (ag > 1e-14)
    coords = np.where(hot_solid)
    if len(coords[0]):
        xyz = [(axis + 0.5 - n / 2) * spacing[i] * 1e6 for i, axis in enumerate(coords)]
        points = ax1.scatter(*xyz, c=temperature[coords], cmap="inferno", vmin=1000,
                             vmax=float(temperature.max()), s=4, alpha=0.9, depthshade=False)
        fig.colorbar(points, ax=ax1, shrink=0.78, pad=0.10, label="Temperature (K)")
    ax1.set(xlim=(-edge_um/2, edge_um/2), ylim=(-edge_um/2, edge_um/2),
            zlim=(-edge_um/2, edge_um/2), xlabel="x (μm)", ylabel="y (μm)",
            zlabel="z (μm)", title="3-D hot solid Ag cells (T ≥ 1000 K)")
    ax1.set_box_aspect((1, 1, 1))
    ax1.view_init(elev=23, azim=38)
    fig.suptitle(f"Pankratz v2 · 100³ · fully frozen at {time_ns:.3f} ns · cooling event refinement",
                 y=0.96, fontsize=12)
    fig.savefig(FIG / "pankratz_v2_n100_first_freeze_3d_slice.png", dpi=220)
    plt.close(fig)


def plot_exploratory_melt_endpoint():
    """展示独立数值验证热源下的50³局部Ag熔化区域。"""
    path = HERE / "results" / "scaled_n50_exploratory_target1400K" / "heated_state.npz"
    if not path.exists():
        return
    with np.load(path, allow_pickle=False) as data:
        temperature = data["temperature_K"]
        liquid = data["liquid_fraction"]
        ag = data["ag_fraction"]
        spacing = data["spacing_m"]
    n = temperature.shape[0]
    edge_um = n * spacing[0] * 1e6
    centers = (np.arange(n) + 0.5 - n / 2.0) * spacing[0] * 1e6
    extent = (-edge_um / 2, edge_um / 2, -edge_um / 2, edge_um / 2)
    fig = plt.figure(figsize=(12.2, 5.5))
    grid = fig.add_gridspec(1, 2, width_ratios=(1.0, 1.15),
                            left=0.07, right=0.91, bottom=0.12, top=0.84, wspace=0.28)
    ax0 = fig.add_subplot(grid[0, 0])
    ax1 = fig.add_subplot(grid[0, 1], projection="3d")
    z = n // 2
    im = ax0.imshow(temperature[:, :, z].T, origin="lower", extent=extent,
                    cmap="inferno", vmin=300.0, vmax=float(temperature.max()),
                    interpolation="nearest")
    ax0.contour(centers, centers, liquid[:, :, z].T, levels=[0.5],
                colors=["#4ee8ff"], linewidths=1.2)
    fig.colorbar(im, ax=ax0, shrink=0.86, pad=0.035, label="Temperature (K)")
    ax0.set(title="Midplane temperature / liquid-fraction contour",
            xlabel="x (μm)", ylabel="y (μm)")
    ax0.set_aspect("equal")
    melt = (ag >= 0.5) & (liquid >= 0.5)
    coords = np.where(melt)
    if len(coords[0]):
        xyz = [(axis + 0.5 - n / 2.0) * spacing[i] * 1e6
               for i, axis in enumerate(coords)]
        points = ax1.scatter(*xyz, c=liquid[coords], cmap="cool", vmin=0.5, vmax=1.0,
                             s=30, alpha=0.95, depthshade=False)
        fig.colorbar(points, ax=ax1, shrink=0.78, pad=0.10, label="Ag liquid fraction")
    ax1.set(xlim=(-edge_um/2, edge_um/2), ylim=(-edge_um/2, edge_um/2),
            zlim=(-edge_um/2, edge_um/2), xlabel="x (μm)", ylabel="y (μm)",
            zlabel="z (μm)", title="Connected Ag cells: Ag ≥ 0.5, liquid ≥ 0.5")
    ax1.set_box_aspect((1, 1, 1))
    ax1.view_init(elev=23, azim=38)
    fig.suptitle("Pankratz v2 · 50³ · 10 μm RVE · 1 ns numerical validation pulse",
                 y=0.96, fontsize=12)
    fig.savefig(FIG / "exploratory_50_numerical_melt_region_3d.png", dpi=220)
    plt.close(fig)


def plot_n100_time_refinement():
    colors = {"1": "#365f91", "0.5": "#d18028", "0.25": "#438b67"}
    fig, axes = plt.subplots(3, 1, figsize=(8.8, 10.0), constrained_layout=True)
    for factor, color in colors.items():
        path = HERE / "results" / f"pankratz_v2_n100_dt{factor}" / "accepted_steps.csv"
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        time_ns = np.array([float(row["time_s"]) * 1e9 for row in rows])
        max_t = np.array([float(row["max_temperature_K"]) for row in rows])
        bulk = np.array([float(row["ag_bulk_liquid_fraction"]) for row in rows])
        energy = np.array([float(row["normalized_cumulative_energy_residual"]) for row in rows])
        axes[0].semilogx(time_ns, max_t, color=color, lw=1.7, label=f"dt factor {factor}")
        axes[1].semilogx(time_ns, np.maximum(bulk, 1e-12), color=color, lw=1.7,
                         label=f"dt factor {factor}")
        axes[2].loglog(time_ns, np.maximum(energy, 1e-16), color=color, lw=1.5,
                       label=f"dt factor {factor}")
    for ax in axes:
        ax.axvline(1.0, color="black", ls=":", lw=1)
        ax.grid(alpha=0.25, which="both")
        ax.legend(frameon=False, fontsize=8)
    axes[0].set(ylabel="Maximum temperature (K)", title="100³ QGS Pankratz v2 heat/cool; same numerical pulse")
    axes[1].set(ylabel="Ag bulk liquid fraction", title="Ag liquid fraction; values floored only for log display")
    axes[2].set(xlabel="Physical time (ns)", ylabel="Normalized cumulative energy residual",
                title="Cumulative energy balance")
    fig.savefig(FIG / "n100_time_refinement_energy.png", dpi=220)
    plt.close(fig)


def plot_adjacent_time_convergence():
    path = HERE / "results" / "time_refinement_n100_comparison.json"
    if not path.exists():
        return
    payload = json.loads(path.read_text(encoding="utf-8"))
    changes = payload.get("adjacent_relative_changes", {})
    names = list(changes)
    metrics = ("first_melt_time_s", "heat_end_ag_bulk_liquid_fraction",
               "cooling_time_to_complete_freeze_s")
    labels = ("First melt time", "Heat-end Ag bulk liquid", "Complete-freeze time")
    x = np.arange(len(names))
    width = 0.23
    fig, ax = plt.subplots(figsize=(8.2, 5.2), constrained_layout=True)
    for i, (metric, label, color) in enumerate(zip(metrics, labels,
                                                    ("#4d79a8", "#d17c28", "#53966b"))):
        values = [100.0 * changes[name].get(metric, np.nan) for name in names]
        ax.bar(x + (i - 1) * width, values, width, label=label, color=color)
    ax.axhline(3.0, color="black", ls="--", lw=1, label="3% target")
    ax.set_xticks(x, names)
    ax.set(ylabel="Relative change (%)", title="Adjacent 100³ time-step refinement")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "n100_adjacent_time_convergence.png", dpi=220)
    plt.close(fig)


def plot_local_source_liquid_fraction():
    path = HERE / "results" / "LOCAL_SOURCE_REGION_METRICS.json"
    if not path.exists():
        return
    rows = json.loads(path.read_text(encoding="utf-8"))
    labels = [row["case"].replace("50³_same_load_", "50³ ").replace("100³_Pankratz_v2_", "100³ ")
              for row in rows]
    bulk = [row["Ag_bulk_liquid_fraction"] for row in rows]
    local = [row["source_region_Ag_liquid_fraction"] for row in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), constrained_layout=True)
    axes[0].bar(x, bulk, color="#4d79a8")
    axes[0].set(title="Full-RVE Ag liquid fraction", ylabel="Liquid fraction")
    axes[1].bar(x, local, color="#d17c28")
    axes[1].set(title="Ag liquid fraction in source cells", ylabel="Liquid fraction")
    for ax in axes:
        ax.set_xticks(x, labels, rotation=25, ha="right")
        ax.grid(axis="y", alpha=0.25)
    fig.savefig(FIG / "local_vs_global_ag_liquid_fraction.png", dpi=220)
    plt.close(fig)


def plot_stefan_validation():
    rows = read_csv("results/STEFAN_V2_VALIDATION.csv")
    fig, ax = plt.subplots(figsize=(7.6, 5.0), constrained_layout=True)
    for mode, color, marker in (("melt", "#4d79a8", "o"), ("freeze", "#d17c28", "s")):
        selected = [row for row in rows if row["mode"] == mode]
        dx_um = np.array([float(row["dx_m"]) * 1e6 for row in selected])
        error = 100.0 * np.array([float(row["relative_front_error"]) for row in selected])
        ax.loglog(dx_um, error, marker + "-", color=color, label=f"{mode.capitalize()} front")
    ax.axhline(2.0, color="black", ls="--", lw=1, label="2% target")
    ax.set(xlabel="Grid spacing (μm)", ylabel="Stefan-front error (%)",
           title="Pankratz v2 one-dimensional Stefan checks")
    ax.invert_xaxis()
    ax.grid(alpha=0.25, which="both")
    ax.legend(frameon=False)
    fig.savefig(FIG / "stefan_front_convergence_v2.png", dpi=220)
    plt.close(fig)


def plot_interface_resistance():
    rows = read_csv("results/INTERFACE_RESISTANCE_V2.csv")
    predicted = np.array([float(row["predicted_jump_qR_K"]) for row in rows])
    measured = np.array([float(row["temperature_jump_K"]) for row in rows])
    limit = max(float(np.max(predicted)), 1e-12)
    fig, ax = plt.subplots(figsize=(6.7, 5.0), constrained_layout=True)
    ax.plot([0, limit], [0, limit], color="black", ls="--", label="qR reference")
    ax.scatter(predicted, measured, color="#4d79a8", zorder=3)
    ax.set(xlabel="Predicted interface jump qR (K)", ylabel="Computed face-temperature jump (K)",
           title="Numerical interface-resistance check")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(FIG / "interface_resistance_qR_check.png", dpi=220)
    plt.close(fig)


def plot_mushy_width_sensitivity():
    payload = json.loads((HERE / "results" / "same_load_50" / "mushy_width_sensitivity.json")
                         .read_text(encoding="utf-8"))
    rows = payload["cases"]
    width = np.array([float(row["Ag_mushy_interval_width_K"]) for row in rows])
    peak = np.array([float(row["maximum_temperature_K"]) for row in rows])
    bulk = np.array([float(row["Ag_bulk_liquid_fraction_at_heating_end"]) for row in rows])
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.6), constrained_layout=True)
    axes[0].plot(width, peak, "o-", color="#4d79a8")
    axes[0].set(xlabel="Numerical mushy interval width (K)", ylabel="Peak temperature (K)",
                title="Peak temperature")
    axes[1].plot(width, bulk, "s-", color="#d17c28")
    axes[1].set(xlabel="Numerical mushy interval width (K)", ylabel="Ag bulk liquid fraction",
                title="Ag liquid at heat end")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.savefig(FIG / "mushy_width_sensitivity_v2.png", dpi=220)
    plt.close(fig)


def plot_solver_performance():
    path = HERE / "results" / "time_refinement_n100_comparison.csv"
    if not path.exists():
        return
    rows = read_csv("results/time_refinement_n100_comparison.csv")
    factors = [row["dt_factor"] for row in rows]
    wall = [float(row["run_wall_seconds"]) / 60.0 for row in rows]
    steps = [int(row["accepted_steps"]) for row in rows]
    rss = [float(row["peak_rss_mib"]) for row in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 4.3), constrained_layout=True)
    axes[0].bar(x, wall, color="#4d79a8")
    axes[0].set(ylabel="Wall time (min)", title="100³ run")
    axes[1].bar(x, steps, color="#d17c28")
    axes[1].set(ylabel="Accepted steps", title="Time resolution")
    axes[2].bar(x, rss, color="#53966b")
    axes[2].set(ylabel="Peak RSS (MiB)", title="Memory")
    for ax in axes:
        ax.set_xticks(x, factors)
        ax.set_xlabel("Time-step factor")
        ax.grid(axis="y", alpha=0.25)
    fig.savefig(FIG / "n100_solver_performance.png", dpi=220)
    plt.close(fig)


def plot_scale_similarity():
    micro = HERE / "results" / "scaled_n50_micro200um_similarity" / "accepted_steps.csv"
    fine = RESULTS / "pankratz_v2_step_history.csv"
    if not micro.exists() or not fine.exists():
        return
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7.0), constrained_layout=True)
    for path, scale, label, color in (
            (fine, 1.0, "10 μm ultrafine baseline", "#365f91"),
            (micro, 1.0 / 400.0, "200 μm micro; time rescaled by 1/400", "#d17c28")):
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        # 新旧运行器CSV字段大小写不同，统一按不区分大小写读取。
        columns = {name.lower(): name for name in rows[0]}
        t = np.array([float(row[columns["time_s"]]) for row in rows])
        maximum = np.array([float(row[columns["max_temperature_k"]]) for row in rows])
        liquid = np.array([float(row[columns["ag_bulk_liquid_fraction"]]) for row in rows])
        t_ns = t * scale * 1e9
        axes[0].semilogx(t_ns, maximum, color=color, lw=1.7, label=label)
        axes[1].semilogx(t_ns, np.maximum(liquid, 1e-12), color=color, lw=1.7, label=label)
    for ax in axes:
        ax.axvline(1.0, color="black", ls=":", lw=1)
        ax.grid(alpha=0.25, which="both")
        ax.legend(frameon=False, fontsize=8)
    axes[0].set(ylabel="Maximum temperature (K)", title="Scale-similar temperature histories")
    axes[1].set(xlabel="Rescaled physical time (ns)", ylabel="Ag bulk liquid fraction",
                title="Equal dimensionless pulse, same QGS geometry")
    fig.savefig(FIG / "scale_similarity_10um_vs_200um.png", dpi=220)
    plt.close(fig)


def plot_freezing_event_refinement():
    path = HERE / "results" / "FREEZING_EVENT_REFINEMENT.json"
    if not path.exists():
        return
    data = json.loads(path.read_text(encoding="utf-8"))
    cases = data["cases"]
    order = ("uniform_dt_1", "uniform_dt_0.5", "uniform_dt_0.25",
             "liquid_event_cap_10ps", "liquid_event_cap_5ps")
    labels = ("dt 1", "dt 1/2", "dt 1/4", "liquid cap 10 ps", "liquid cap 5 ps")
    values = [float(cases[key]["cooling_time_to_complete_freeze_s"]) * 1e9 for key in order]
    colors = ("#365f91", "#4d79a8", "#8a9db5", "#d18028", "#53966b")
    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    bars = ax.bar(np.arange(len(order)), values, color=colors)
    ax.set_xticks(np.arange(len(order)), labels, rotation=20, ha="right")
    ax.set(ylabel="Cooling time to complete freeze (ns)",
           title="Uniform steps can cross the last liquid cells")
    ax.grid(axis="y", alpha=0.25)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, value, f"{value:.3f}",
                ha="center", va="bottom", fontsize=8)
    fig.savefig(FIG / "n100_freezing_event_time_refinement.png", dpi=220)
    plt.close(fig)


def main():
    plot_silver_conductivity()
    plot_silver_heat_capacity()
    plot_silver_enthalpy()
    plot_tisc_properties()
    plot_same_load_timeseries()
    plot_heated_endpoint()
    plot_n100_heated_endpoint()
    plot_n100_first_freeze()
    plot_exploratory_melt_endpoint()
    plot_n100_time_refinement()
    plot_adjacent_time_convergence()
    plot_local_source_liquid_fraction()
    plot_stefan_validation()
    plot_interface_resistance()
    plot_mushy_width_sensitivity()
    plot_solver_performance()
    plot_scale_similarity()
    plot_freezing_event_refinement()
    print(f"Wrote figures to {FIG}")


if __name__ == "__main__":
    main()
