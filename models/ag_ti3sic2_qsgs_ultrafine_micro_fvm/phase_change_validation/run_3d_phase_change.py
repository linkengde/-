#!/usr/bin/env python3
"""固定几何源的三维短程Ag加热/熔化/热沉冷却验证；不是电弧工况。"""

from __future__ import annotations

import csv
import json
import math
import os
import resource
import sys
import time
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/fvm-phase-mplconfig")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

from implicit_enthalpy_fvm import BoundarySet, ImplicitEnthalpyFVM
from material_model import CandidateMaterials, CONFIG


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
OUT = HERE / "results"
FIG = HERE / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)


def maxrss_mib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024.0 if sys.platform != "darwin" else value / (1024.0 * 1024.0)


def conservative_axis_weights(n_source: int, n_target: int) -> np.ndarray:
    """返回体积重叠权重；对任意正整数重映射保持体积分数平均值。"""
    weights = np.zeros((n_target, n_source), dtype=np.float64)
    source_width = 1.0 / n_source
    target_width = 1.0 / n_target
    for target in range(n_target):
        lo_t, hi_t = target * target_width, (target + 1) * target_width
        first = max(0, int(math.floor(lo_t / source_width)))
        last = min(n_source - 1, int(math.ceil(hi_t / source_width) - 1))
        for source in range(first, last + 1):
            lo_s, hi_s = source * source_width, (source + 1) * source_width
            overlap = max(0.0, min(hi_t, hi_s) - max(lo_t, lo_s))
            weights[target, source] = overlap / target_width
    weights /= weights.sum(axis=1, keepdims=True)
    return weights


def conservative_remap_3d(field: np.ndarray, n_target: int) -> np.ndarray:
    arr = np.asarray(field, dtype=np.float64)
    weights = conservative_axis_weights(arr.shape[0], n_target)
    for axis in range(3):
        moved = np.moveaxis(arr, axis, 0)
        remapped = np.tensordot(weights, moved, axes=(1, 0))
        arr = np.moveaxis(remapped, 0, axis)
    # 0/1輸入重映射理應落在[0,1]；這裡只消除浮點加總造成的約1e-14越界。
    return np.clip(arr, 0.0, 1.0)


def load_seed_ag() -> np.ndarray:
    source = PROJECT / CONFIG["geometry_source"]
    with np.load(source) as data:
        ag = np.asarray(data["ag"], dtype=np.float64)
    return ag


def make_ag_fraction(kind: str, n: int, seed_ag: np.ndarray) -> np.ndarray:
    if kind == "pure_ag":
        return np.ones((n, n, n), dtype=np.float64)
    if kind == "bilayer":
        ag = np.zeros((n, n, n), dtype=np.float64)
        ag[: n // 2, :, :] = 1.0
        return ag
    if kind == "bicontinuous_qsgs_seed1105":
        return conservative_remap_3d(seed_ag, n)
    raise ValueError(kind)


def source_field(ag: np.ndarray, dx: float, rate_W_m3: float, radius_fraction: float) -> np.ndarray:
    n = ag.shape[0]
    coordinates = (np.arange(n, dtype=np.float64) + 0.5 - 0.5 * n) * dx
    x2 = coordinates[:, None, None] ** 2
    y2 = coordinates[None, :, None] ** 2
    z2 = coordinates[None, None, :] ** 2
    radius = radius_fraction * n * dx
    local = x2 + y2 + z2 <= radius * radius
    return rate_W_m3 * ag * local


def metrics(sim: ImplicitEnthalpyFVM, materials: CandidateMaterials, dx: float,
            with_connectivity: bool = False) -> dict:
    t = sim.temperature_K
    liquid = materials.liquid_fraction(t)
    ag = sim.ag_fraction
    ag_volume = float(np.sum(ag))
    liquid_volume = float(np.sum(ag * liquid))
    if with_connectivity:
        active_liquid = (ag > 0.5) & (liquid >= 0.5)
        labels, count = ndimage.label(active_liquid, structure=ndimage.generate_binary_structure(3, 3))
        if count:
            sizes = np.bincount(labels.ravel())[1:]
            max_region_cells = int(sizes.max())
            max_label = 1 + int(sizes.argmax())
            positions = np.where(labels == max_label)
            max_region_depth_m = float((max(positions[2].max() - positions[2].min() + 1, 1)) * dx)
        else:
            max_region_cells = 0
            max_region_depth_m = 0.0
    else:
        max_region_cells = math.nan
        max_region_depth_m = math.nan
    ag_tmean = float(np.sum(ag * t) / max(ag_volume, 1.0e-30))
    sc = 1.0 - ag
    sc_mean = float(np.sum(sc * t) / max(float(np.sum(sc)), 1.0e-30))
    return {
        "physical_time_s": sim.time_s,
        "maximum_temperature_K": float(t.max()),
        "Ag_mean_temperature_K": ag_tmean,
        "Ti3SiC2_mean_temperature_K": sc_mean,
        "Ag_liquid_volume_m3": liquid_volume * sim.volume_m3,
        "Ag_bulk_liquid_fraction": liquid_volume / max(ag_volume, 1.0e-30),
        "Ag_max_liquid_fraction": float(np.max(liquid)),
        "max_connected_liquid_cells_phi_ge_0p5": max_region_cells,
        "max_connected_liquid_volume_m3_phi_ge_0p5": max_region_cells * sim.volume_m3,
        "max_connected_liquid_depth_m_phi_ge_0p5": max_region_depth_m,
        "cumulative_source_energy_J": sim.source_energy_J,
        "cumulative_boundary_heat_out_J": sim.boundary_energy_out_J,
        "stored_enthalpy_change_J": float(np.sum(sim.enthalpy_density_J_m3) * sim.volume_m3 - sim.initial_stored_enthalpy_J),
        "energy_residual_J": float(np.sum(sim.enthalpy_density_J_m3) * sim.volume_m3 - sim.initial_stored_enthalpy_J
                                    - sim.source_energy_J + sim.boundary_energy_out_J),
    }


def interval(sim: ImplicitEnthalpyFVM, duration_s: float, nominal_dt_s: float,
             source: np.ndarray | float = 0.0, boundary: BoundarySet | None = None,
             phase: str = "heating", deadline_s: float | None = None,
             state_callback=None) -> list[dict]:
    remaining = float(duration_s)
    dt = float(nominal_dt_s)
    records: list[dict] = []
    failures = 0
    min_dt = nominal_dt_s / 256.0
    while remaining > max(1.0e-18, duration_s * 1.0e-14):
        if deadline_s is not None and time.perf_counter() >= deadline_s:
            raise TimeoutError(f"wall budget expired during {phase}; physical_time_s={sim.time_s:.9e}")
        trial_dt = min(dt, remaining)
        try:
            record = sim.advance(trial_dt, source, boundary)
            record["phase"] = phase
            record["adaptive_retries"] = failures
            if state_callback is not None:
                state_callback(record)
            records.append(record)
            remaining -= trial_dt
            dt = min(nominal_dt_s, trial_dt * 1.5)
        except RuntimeError:
            failures += 1
            dt = trial_dt * 0.5
            if dt < min_dt:
                raise
    return records


def save_fields(sim: ImplicitEnthalpyFVM, path: Path) -> None:
    np.savez_compressed(path, temperature_K=sim.temperature_K.astype(np.float32),
                        ag_fraction=sim.ag_fraction.astype(np.float32),
                        ag_liquid_fraction=sim.materials.liquid_fraction(sim.temperature_K).astype(np.float32),
                        spacing_m=sim.spacing_m)


def plot_slice(case: str, n: int, phase: str, sim: ImplicitEnthalpyFVM, materials: CandidateMaterials) -> None:
    center = sim.temperature_K.shape[2] // 2
    t_slice = sim.temperature_K[:, :, center]
    l_slice = (sim.ag_fraction * materials.liquid_fraction(sim.temperature_K))[:, :, center]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    im0 = axes[0].imshow(t_slice.T, origin="lower", cmap="inferno", vmin=300,
                          vmax=max(materials.Tm_K + 20.0, float(t_slice.max())))
    axes[0].set_title(f"T (K), t={sim.time_s:.3e} s")
    fig.colorbar(im0, ax=axes[0], shrink=0.82)
    im1 = axes[1].imshow(l_slice.T, origin="lower", cmap="magma", vmin=0.0, vmax=1.0)
    axes[1].set_title("Ag liquid fraction × Ag voxel fraction")
    fig.colorbar(im1, ax=axes[1], shrink=0.82)
    for ax in axes:
        ax.set_xlabel("x cell")
        ax.set_ylabel("y cell")
    fig.suptitle(f"{case}, {n}³, {phase}; numerical validation load")
    fig.savefig(FIG / f"{case}_n{n}_{phase}_slice.png", dpi=180)
    plt.close(fig)


def plot_3d_volume(case: str, n: int, phase: str, sim: ImplicitEnthalpyFVM,
                   materials: CandidateMaterials) -> None:
    stride = max(1, n // 36)
    sampled_t = sim.temperature_K[::stride, ::stride, ::stride]
    sampled_ag = sim.ag_fraction[::stride, ::stride, ::stride]
    sampled_l = sampled_ag * materials.liquid_fraction(sampled_t)
    grid = np.indices(sampled_t.shape)
    fig = plt.figure(figsize=(8.5, 7.5))
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(grid[0].ravel(), grid[1].ravel(), grid[2].ravel(), c=sampled_t.ravel(),
                    cmap="inferno", s=3.0, alpha=0.16, linewidths=0)
    liquid = sampled_l >= 0.05
    if np.any(liquid):
        ax.scatter(grid[0][liquid], grid[1][liquid], grid[2][liquid],
                   c=sampled_l[liquid], cmap="cool", vmin=0.0, vmax=1.0,
                   s=10.0, alpha=0.9, linewidths=0, label="Ag liquid-bearing cells")
        ax.legend(loc="upper left")
    fig.colorbar(sc, ax=ax, shrink=0.68, pad=0.08, label="Temperature (K)")
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(f"{case} {n}³ {phase}; T volume samples and Ag phase fraction")
    fig.tight_layout()
    fig.savefig(FIG / f"{case}_n{n}_{phase}_3d_volume.png", dpi=180)
    plt.close(fig)


def plot_structure(seed_ag: np.ndarray) -> None:
    n = seed_ag.shape[0]
    step = max(1, n // 36)
    ag = seed_ag[::step, ::step, ::step] >= 0.5
    sc = ~ag
    coords_ag = np.where(ag)
    coords_sc = np.where(sc)
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(*[c[:: max(1, len(c) // 15000)] for c in coords_ag], s=1.2, c="#d5dbe5", alpha=0.16, label="Ag")
    ax.scatter(*[c[:: max(1, len(c) // 15000)] for c in coords_sc], s=1.5, c="#12a89d", alpha=0.72, label="Ti₃SiC₂")
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("x (coarse voxels)")
    ax.set_ylabel("y (coarse voxels)")
    ax.set_zlabel("z (coarse voxels)")
    ax.set_title("QSGS seed 1105; conservative volume-average representation")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "initial_qsgs_three_dimensional.png", dpi=200)
    plt.close(fig)


def run_case(kind: str, n: int, seed_ag: np.ndarray, edge_m: float,
             source_rate: float, heat_duration: float, cool_duration: float,
             temporal_refinement: bool = False, dt_factor: float = 1.0,
             case_label: str | None = None,
             wall_time_budget_s: float | None = None) -> tuple[list[dict], dict]:
    output_label = case_label or (kind if dt_factor == 1.0 else f"{kind}_dt{dt_factor:g}")
    ag = make_ag_fraction(kind, n, seed_ag)
    dx = edge_m / n
    materials = CandidateMaterials()
    sim = ImplicitEnthalpyFVM(ag, (dx, dx, dx), materials=materials,
                              interface_resistance_m2_K_W=0.0, initial_temperature_K=300.0,
                              nonlinear_relative_tolerance=1.0e-8,
                              linear_relative_tolerance=1.0e-8,
                              maximum_nonlinear_iterations=48, maximum_linear_iterations=400)
    q_source = source_field(ag, dx, source_rate, float(CONFIG["numerical_heating_test"]["source_region_radius_fraction_of_edge"]))
    sink = BoundarySet(dirichlet_K={"x_low": 300.0, "x_high": 300.0,
                                    "y_low": 300.0, "y_high": 300.0,
                                    "z_low": 300.0, "z_high": 300.0})
    records: list[dict] = []
    pre = metrics(sim, materials, dx, with_connectivity=True)
    pre.update({"case": kind, "resolution": n, "phase": "initial", "wall_seconds": 0.0,
                "pcg_iterations_total": 0, "nonlinear_iterations": 0, "max_residual": 0.0,
                "maxrss_mib": maxrss_mib(), "ag_fraction_mean": float(ag.mean()), "spacing_m": dx})
    records.append(pre)
    start = time.perf_counter()
    deadline = None if wall_time_budget_s is None else start + wall_time_budget_s
    heat_dt = heat_duration / (20.0 if n >= 100 else 10.0) * dt_factor
    def capture_state(step_record: dict) -> None:
        item = metrics(sim, materials, dx)
        item.update({"case": kind, "resolution": n, "phase": step_record["phase"],
                     "wall_seconds": step_record["step_wall_seconds"],
                     "pcg_iterations_total": step_record["pcg_iterations_total"],
                     "nonlinear_iterations": step_record["nonlinear_iterations"],
                     "max_residual": step_record["nonlinear_relative_residual"],
                     "maxrss_mib": maxrss_mib(), "ag_fraction_mean": float(ag.mean()), "spacing_m": dx,
                     "dt_s": step_record["dt_s"], "adaptive_retries": step_record["adaptive_retries"]})
        records.append(item)

    try:
        heat_records = interval(sim, heat_duration, heat_dt, q_source, None, "heating", deadline,
                                state_callback=capture_state)
    except TimeoutError as exc:
        partial = metrics(sim, materials, dx, with_connectivity=True)
        partial.update({"case": kind, "resolution": n, "phase": "incomplete_heating",
                        "run_status": "TIMEOUT", "status_detail": str(exc),
                        "physical_time_s": sim.time_s, "wall_seconds": time.perf_counter() - start,
                        "maxrss_mib": maxrss_mib(), "ag_fraction_mean": float(ag.mean()), "spacing_m": dx})
        records.append(partial)
        save_fields(sim, OUT / f"{output_label}_n{n}_incomplete_heating.npz")
        scale = max(abs(sim.source_energy_J), abs(partial["stored_enthalpy_change_J"]), 1.0e-30)
        return records, {"case": kind, "resolution": n, "domain_edge_m": edge_m,
                         "spacing_m": dx, "run_status": "TIMEOUT", "timeout_phase": "heating",
                         "physical_time_completed_s": sim.time_s, "wall_seconds": time.perf_counter() - start,
                         "peak_rss_mib": maxrss_mib(), "energy_residual_J": partial["energy_residual_J"],
                         "normalized_energy_residual": abs(partial["energy_residual_J"]) / scale,
                         "source_integral_J": sim.source_energy_J, "melted": partial["Ag_max_liquid_fraction"] > 0.0,
                         "frozen_after_cooling": False, "energy_check_pass": False,
                         "status_detail": str(exc)}
    heated_metrics = metrics(sim, materials, dx, with_connectivity=True)
    plot_slice(output_label, n, "heated", sim, materials)
    if kind == "bicontinuous_qsgs_seed1105":
        plot_3d_volume(output_label, n, "heated", sim, materials)
    save_fields(sim, OUT / f"{output_label}_n{n}_heated.npz")

    cooling_dt = cool_duration / (40.0 if n >= 100 else 20.0) * dt_factor
    try:
        cool_records = interval(sim, cool_duration, cooling_dt, 0.0, sink, "cooling", deadline,
                                state_callback=capture_state)
    except TimeoutError as exc:
        partial = metrics(sim, materials, dx, with_connectivity=True)
        partial.update({"case": kind, "resolution": n, "phase": "incomplete_cooling",
                        "run_status": "TIMEOUT", "status_detail": str(exc),
                        "physical_time_s": sim.time_s, "wall_seconds": time.perf_counter() - start,
                        "maxrss_mib": maxrss_mib(), "ag_fraction_mean": float(ag.mean()), "spacing_m": dx})
        records.append(partial)
        save_fields(sim, OUT / f"{output_label}_n{n}_incomplete_cooling.npz")
        scale = max(abs(sim.source_energy_J), abs(sim.boundary_energy_out_J),
                    abs(partial["stored_enthalpy_change_J"]), 1.0e-30)
        return records, {"case": kind, "resolution": n, "domain_edge_m": edge_m,
                         "spacing_m": dx, "run_status": "TIMEOUT", "timeout_phase": "cooling",
                         "physical_time_completed_s": sim.time_s, "wall_seconds": time.perf_counter() - start,
                         "peak_rss_mib": maxrss_mib(), "energy_residual_J": partial["energy_residual_J"],
                         "normalized_energy_residual": abs(partial["energy_residual_J"]) / scale,
                         "source_integral_J": sim.source_energy_J, "melted": heated_metrics["Ag_max_liquid_fraction"] > 0.0,
                         "frozen_after_cooling": False, "energy_check_pass": False,
                         "status_detail": str(exc)}
    final_metrics = metrics(sim, materials, dx, with_connectivity=True)
    plot_slice(output_label, n, "cooled", sim, materials)
    if kind == "bicontinuous_qsgs_seed1105":
        plot_3d_volume(output_label, n, "cooled", sim, materials)
    save_fields(sim, OUT / f"{output_label}_n{n}_cooled.npz")

    energy_scale = max(abs(sim.source_energy_J), abs(sim.boundary_energy_out_J),
                       abs(final_metrics["stored_enthalpy_change_J"]), 1.0e-30)
    normalized_energy_error = abs(final_metrics["energy_residual_J"]) / energy_scale
    heat_stats = {
        "case": kind,
        "resolution": n,
        "domain_edge_m": edge_m,
        "spacing_m": dx,
        "ag_fraction_mean": float(ag.mean()),
        "ti3sic2_fraction_mean": float((1.0 - ag).mean()),
        "source_rate_W_m3_on_Ag_volume": source_rate,
        "source_integral_J": sim.source_energy_J,
        "heat_duration_s": heat_duration,
        "cool_duration_s": cool_duration,
        "initial_temperature_K": 300.0,
        "hot_sink_boundary_K": 300.0,
        "heated_max_temperature_K": heated_metrics["maximum_temperature_K"],
        "heated_Ag_mean_temperature_K": heated_metrics["Ag_mean_temperature_K"],
        "heated_Ti3SiC2_mean_temperature_K": heated_metrics["Ti3SiC2_mean_temperature_K"],
        "heated_bulk_Ag_liquid_fraction": heated_metrics["Ag_bulk_liquid_fraction"],
        "heated_max_Ag_liquid_fraction": heated_metrics["Ag_max_liquid_fraction"],
        "heated_liquid_volume_m3": heated_metrics["Ag_liquid_volume_m3"],
        "heated_max_connected_liquid_volume_m3_phi_ge_0p5": heated_metrics["max_connected_liquid_volume_m3_phi_ge_0p5"],
        "cooled_max_temperature_K": final_metrics["maximum_temperature_K"],
        "cooled_Ag_mean_temperature_K": final_metrics["Ag_mean_temperature_K"],
        "cooled_Ag_bulk_liquid_fraction": final_metrics["Ag_bulk_liquid_fraction"],
        "cooled_max_Ag_liquid_fraction": final_metrics["Ag_max_liquid_fraction"],
        "sink_energy_out_J": sim.boundary_energy_out_J,
        "energy_residual_J": final_metrics["energy_residual_J"],
        "normalized_energy_residual": normalized_energy_error,
        "melted": heated_metrics["Ag_max_liquid_fraction"] > 0.0,
        "frozen_after_cooling": final_metrics["Ag_bulk_liquid_fraction"] < 1.0e-6,
        "wall_seconds": time.perf_counter() - start,
        "peak_rss_mib": maxrss_mib(),
        "heat_steps": len(heat_records),
        "cooling_steps": len(cool_records),
        "max_pcg_iterations_per_step": int(max((r["pcg_iterations_total"] for r in heat_records + cool_records), default=0)),
        "max_newton_iterations_per_step": int(max((r["nonlinear_iterations"] for r in heat_records + cool_records), default=0)),
        "adaptive_retry_count": int(sum(r["adaptive_retries"] for r in heat_records + cool_records)),
        "time_step_refinement_case": temporal_refinement,
        "nominal_time_step_factor": dt_factor,
    }
    if normalized_energy_error > 1.0e-3:
        heat_stats["energy_check_pass"] = False
    else:
        heat_stats["energy_check_pass"] = True
    return records, heat_stats


def save_checkpoint(cases: list[dict], history: list[dict], status: str) -> None:
    """每个算例完成即落盘，避免后续长算例失败时丢失已完成结果。"""
    (OUT / "3d_phase_change_summary.json").write_text(
        json.dumps({"cases": cases, "status": status,
                    "source_interpretation": "numerical local volumetric heating only; not an arc load",
                    "cooling_boundary": "all six external faces at 300 K",
                    "geometry_map": "conservative volume-overlap remap from same seed 1105"},
                   indent=2, ensure_ascii=False), encoding="utf-8")
    if history:
        fields = list(dict.fromkeys(key for row in history for key in row))
        with (OUT / "3d_phase_change_history.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(history)
    if cases:
        fields = list(dict.fromkeys(key for row in cases for key in row))
        with (OUT / "MELT_FREEZE_METRICS.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(cases)


def main() -> None:
    seed_ag = load_seed_ag()
    plot_structure(seed_ag)
    cfg = CONFIG["numerical_heating_test"]
    edge = float(CONFIG["domain_edge_m"]["ultrafine_500nm"])
    rates = float(cfg["source_rate_W_m3_at_10um"])
    duration = float(cfg["heat_duration_s_at_10um"])
    cool = float(cfg["cool_duration_s_at_10um"])
    cases: list[dict] = []
    history: list[dict] = []
    # 先完成50³三种拓扑及其时间步对照，再升到100³。
    for kind in ("pure_ag", "bilayer", "bicontinuous_qsgs_seed1105"):
        records, summary = run_case(kind, 50, seed_ag, edge, rates, duration, cool)
        cases.append(summary)
        history.extend(records)
        save_checkpoint(cases, history, f"50^3 complete; next={kind}")
        print(json.dumps({"completed": f"{kind} 50^3", **summary}, ensure_ascii=False), flush=True)
    base_50 = next(c for c in cases if c["resolution"] == 50 and c["case"] == "bicontinuous_qsgs_seed1105")
    half_records, half_50 = run_case("bicontinuous_qsgs_seed1105", 50, seed_ag, edge, rates,
                                     duration, cool, temporal_refinement=True, dt_factor=0.5,
                                     case_label="bicontinuous_qsgs_seed1105_halfdt")
    history.extend(half_records)
    cases.append(half_50)
    save_checkpoint(cases, history, "50^3 base and half time step complete")
    quarter_records, refined_50 = run_case("bicontinuous_qsgs_seed1105", 50, seed_ag, edge, rates,
                                           duration, cool, temporal_refinement=True, dt_factor=0.25,
                                           case_label="bicontinuous_qsgs_seed1105_quarterdt")
    history.extend(quarter_records)
    for metric_name in ("heated_maximum_temperature_K", "heated_bulk_Ag_liquid_fraction",
                        "cooled_Ag_bulk_liquid_fraction"):
        baseline_key = {
            "heated_maximum_temperature_K": "heated_max_temperature_K",
            "heated_bulk_Ag_liquid_fraction": "heated_bulk_Ag_liquid_fraction",
            "cooled_Ag_bulk_liquid_fraction": "cooled_Ag_bulk_liquid_fraction",
        }[metric_name]
        half_key = baseline_key
        a, b = half_50[half_key], refined_50[baseline_key]
        refined_50[f"relative_change_{metric_name}"] = abs(b - a) / max(abs(b), 1.0e-12)
    time_changes = [refined_50[k] for k in refined_50 if k.startswith("relative_change_")]
    refined_50["time_step_convergence_pass"] = max(time_changes) <= 0.03
    base_50["time_step_convergence_pass"] = refined_50["time_step_convergence_pass"]
    cases.append(refined_50)
    save_checkpoint(cases, history, "50^3 complete including half time step")
    print(json.dumps({"completed": "bicontinuous_qsgs_seed1105 50^3 dt/2-to-dt/4 convergence",
                      "time_step_convergence_pass": refined_50["time_step_convergence_pass"],
                      "changes": {k: refined_50[k] for k in refined_50 if k.startswith("relative_change_")}},
                     ensure_ascii=False), flush=True)
    timeout_occurred = False
    case_budget = float(cfg.get("wall_time_budget_s_per_case", 3600.0))
    for kind in ("pure_ag", "bilayer", "bicontinuous_qsgs_seed1105"):
        records, summary = run_case(kind, 100, seed_ag, edge, rates, duration, cool,
                                    wall_time_budget_s=case_budget)
        cases.append(summary)
        history.extend(records)
        save_checkpoint(cases, history, f"100^3 checkpoint after {kind}")
        print(json.dumps({"completed": f"{kind} 100^3", **summary}, ensure_ascii=False), flush=True)
        if summary.get("run_status") == "TIMEOUT":
            timeout_occurred = True
            break
    # 150³仅在前六个短算例稳定且能量闭合后运行；作为原分辨率探索性延伸。
    can_start_150 = (not timeout_occurred
                     and all(c.get("energy_check_pass", False) and c.get("melted", False)
                             and c.get("frozen_after_cooling", False) for c in cases)
                     and refined_50["time_step_convergence_pass"])
    if can_start_150:
        for kind in ("bicontinuous_qsgs_seed1105",):
          for label, edge_case, scale in (("ultrafine_10um", edge, 1.0),
                                          ("micro_200um", edge * 20.0, 400.0)):
            records, summary = run_case(kind, 150, seed_ag, edge_case, rates / scale,
                                        duration * scale, cool * scale,
                                        case_label=f"{kind}_{label}", wall_time_budget_s=case_budget)
            summary["scale_model"] = label
            summary["scope_note"] = "150^3 short numerical heating/cooling extension; real material source provenance incomplete"
            cases.append(summary)
            history.extend(records)
            print(json.dumps({"completed": f"{kind} 150^3 {label}", **summary}, ensure_ascii=False), flush=True)
            if summary.get("run_status") == "TIMEOUT":
                break
    save_checkpoint(cases, history, f"complete; 150_started={can_start_150}; peak_rss_mib={maxrss_mib():.1f}")
    print(json.dumps({"cases_completed": len(cases), "150_started": can_start_150,
                      "peak_rss_mib": maxrss_mib()}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
