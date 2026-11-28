#!/usr/bin/env python3
"""第一轮：焓法反算、能量守恒、Stefan前沿和界面热阻验证。"""

from __future__ import annotations

import csv
import json
import math
import resource
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import erf

from implicit_enthalpy_fvm import BoundarySet, ImplicitEnthalpyFVM
from material_model import CandidateMaterials


HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)


class ConstantAgBenchmark:
    """只供解析基准使用的常物性Ag，不代表实际材料参数。"""

    def __init__(self, phase_width: float = 0.1, freezing: bool = False):
        self.Tm_K = 1000.0
        self.width_K = phase_width
        self.Ts_K = self.Tm_K - phase_width if freezing else self.Tm_K
        self.Tl_K = self.Tm_K if freezing else self.Tm_K + phase_width
        self.Lf_J_kg = 2.0e5
        self.rho_ag_solid_kg_m3 = 1000.0
        self.rho_ti_kg_m3 = 4000.0
        self.ti_Tmax_K = 1473.15
        self.reference_temperature_K = 300.0
        self.cp = 1000.0
        self.k = 10.0

    def liquid_fraction(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        value = np.clip((t - self.Ts_K) / (self.Tl_K - self.Ts_K), 0.0, 1.0)
        return float(value) if t.ndim == 0 else value

    def ag_specific_enthalpy_J_kg(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        h = self.cp * (t - self.Tm_K) + self.Lf_J_kg * self.liquid_fraction(t)
        return float(h) if t.ndim == 0 else h

    def ag_apparent_cp_J_kg_K(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        mushy = (t > self.Ts_K) & (t < self.Tl_K)
        cp = self.cp + np.where(mushy, self.Lf_J_kg / (self.Tl_K - self.Ts_K), 0.0)
        return float(cp) if t.ndim == 0 else cp

    def ag_conductivity_W_m_K(self, temperature_K):
        return np.full_like(np.asarray(temperature_K, dtype=float), self.k)

    def enthalpy_density_J_m3(self, temperature_K, ag_fraction):
        t, ag = np.broadcast_arrays(np.asarray(temperature_K, dtype=float), np.asarray(ag_fraction, dtype=float))
        return (ag * self.rho_ag_solid_kg_m3 * self.ag_specific_enthalpy_J_kg(t)
                + (1.0 - ag) * 1.0e6 * (t - self.reference_temperature_K))

    def volumetric_capacity_J_m3_K(self, temperature_K, ag_fraction):
        t, ag = np.broadcast_arrays(np.asarray(temperature_K, dtype=float), np.asarray(ag_fraction, dtype=float))
        return ag * self.rho_ag_solid_kg_m3 * self.ag_apparent_cp_J_kg_K(t) + (1.0 - ag) * 1.0e6

    def conductivity_W_m_K(self, temperature_K, ag_fraction):
        t, ag = np.broadcast_arrays(np.asarray(temperature_K, dtype=float), np.asarray(ag_fraction, dtype=float))
        return ag * self.k + (1.0 - ag) * 2.0


def stefan_lambda(ste: float) -> float:
    def residual(lam):
        return math.sqrt(math.pi) * lam * math.exp(lam * lam) * math.erf(lam) - ste
    return brentq(residual, 1.0e-12, 5.0)


def maxrss_mib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024.0 if sys.platform != "darwin" else value / (1024.0 * 1024.0)


def result_row(name: str, passed: bool, metric: float, criterion: str, details: str = ""):
    return {"case": name, "pass": bool(passed), "metric": float(metric), "criterion": criterion, "details": details}


def run_material_roundtrip(rows: list[dict]) -> dict:
    materials = CandidateMaterials()
    temperatures = np.linspace(300.0, materials.Tm_K + 30.0, 600)
    fractions = np.ones_like(temperatures)
    h = materials.enthalpy_density_J_m3(temperatures, fractions)
    solver = ImplicitEnthalpyFVM(np.ones(temperatures.size), (1.0e-6,))
    recovered = solver.temperature_from_enthalpy_density(h, temperatures - 20.0)
    h_back = materials.enthalpy_density_J_m3(recovered, fractions)
    max_relative_h_error = float(np.max(np.abs(h_back - h) / np.maximum(np.abs(h), 1.0e6)))
    max_temperature_error = float(np.max(np.abs(recovered - temperatures)))
    max_fraction_error = float(max(np.max(np.maximum(-materials.liquid_fraction(recovered), 0.0)),
                                   np.max(np.maximum(materials.liquid_fraction(recovered) - 1.0, 0.0))))
    rows.append(result_row("A_enthalpy_roundtrip", max_relative_h_error <= 1.0e-8 and max_fraction_error == 0.0,
                           max_relative_h_error, "relative enthalpy error <= 1e-8; liquid fraction in [0,1]",
                           f"maximum_temperature_error_K={max_temperature_error:.6e}; width_K={materials.width_K:g}"))
    return {"max_relative_enthalpy_error": max_relative_h_error,
            "max_temperature_error_K": max_temperature_error,
            "liquid_fraction_range": [float(materials.liquid_fraction(recovered).min()),
                                      float(materials.liquid_fraction(recovered).max())]}


def run_single_cell(rows: list[dict]) -> dict:
    mats = CandidateMaterials()
    dt = 1.0e-3
    target_temperature = mats.Tm_K + 5.0
    ag = np.ones((1,), dtype=float)
    sim = ImplicitEnthalpyFVM(ag, (1.0e-6,), materials=mats, initial_temperature_K=300.0)
    target_h = mats.enthalpy_density_J_m3(np.array([target_temperature]), ag)[0]
    q = target_h / dt
    heated = sim.advance(dt, q)
    latent_expected = mats.rho_ag_solid_kg_m3 * mats.Lf_J_kg * sim.volume_m3
    latent_measured = (mats.rho_ag_solid_kg_m3 * mats.Lf_J_kg
                       * sim.volume_m3 * mats.liquid_fraction(sim.temperature_K)[0])
    latent_error = abs(latent_measured - latent_expected) / latent_expected
    deposited = heated["cumulative_source_energy_J"]
    storage_error = abs(heated["storage_change_step_J"] - deposited) / max(abs(deposited), 1.0e-30)
    rows.append(result_row("B_single_cell_adiabatic_melt", latent_error <= 1.0e-3 and storage_error <= 1.0e-3,
                           max(latent_error, storage_error), "latent and energy errors <= 0.1%",
                           f"target_T_K={target_temperature:g}; final_T_K={sim.temperature_K[0]:.9f}; "
                           f"liquid_fraction={mats.liquid_fraction(sim.temperature_K)[0]:.9f}"))

    # 固定温度热沉双端冷却：散热必须从有限体积边界流出，而非只在域内重分布。
    sink_case = ImplicitEnthalpyFVM(ag, (1.0e-3,), materials=mats, initial_temperature_K=target_temperature)
    sink = BoundarySet(dirichlet_K={"x_low": 300.0, "x_high": 300.0})
    freeze_steps = []
    for _ in range(100):
        freeze_steps.append(sink_case.advance(1.0e-3, 0.0, sink))
    froze = float(mats.liquid_fraction(sink_case.temperature_K)[0]) < 1.0e-8
    sink_loss = sink_case.boundary_energy_out_J
    sink_balance = abs(sink_case.source_energy_J - sink_loss
                       - (float(np.sum(sink_case.enthalpy_density_J_m3) * sink_case.volume_m3)
                          - float(np.sum(mats.enthalpy_density_J_m3(np.array([target_temperature]), ag))
                                   * sink_case.volume_m3)))
    sink_norm = sink_balance / max(abs(sink_loss), 1.0e-30)
    rows.append(result_row("C_single_cell_cooling_freeze", froze and sink_loss > 0.0 and sink_norm <= 1.0e-3,
                           sink_norm, "Ag freezes; positive sink outflow; normalized balance <= 0.1%",
                           f"freeze_time_s={next((x['time_s'] for x in freeze_steps if x['max_ag_liquid_fraction'] < 1e-8), float('nan')):.6e}; "
                           f"final_liquid_fraction={mats.liquid_fraction(sink_case.temperature_K)[0]:.6e}; "
                           f"sink_out_J={sink_loss:.6e}"))
    return {"melt": heated, "cooling": freeze_steps[-1], "latent_expected_J": latent_expected,
            "latent_measured_J": latent_measured, "cooling_sink_out_J": sink_loss, "cooling_balance_rel": sink_norm}


def run_stefan_case(freezing: bool, n: int, dt_s: float, duration_s: float) -> dict:
    mats = ConstantAgBenchmark(phase_width=0.1, freezing=freezing)
    dx = 5.0e-4 / n
    start_temperature = mats.Tm_K if freezing else mats.Tm_K
    boundary_temperature = mats.Tm_K - 100.0 if freezing else mats.Tm_K + 100.0
    sim = ImplicitEnthalpyFVM(np.ones(n), (dx,), materials=mats, initial_temperature_K=start_temperature,
                              nonlinear_relative_tolerance=2.0e-9, linear_relative_tolerance=1.0e-10,
                              maximum_nonlinear_iterations=30, maximum_linear_iterations=800)
    boundary = BoundarySet(dirichlet_K={"x_low": boundary_temperature})
    steps = int(round(duration_s / dt_s))
    last = None
    start = time.perf_counter()
    for _ in range(steps):
        last = sim.advance(dt_s, 0.0, boundary)
    elapsed = time.perf_counter() - start
    phi = mats.liquid_fraction(sim.temperature_K)
    measured_length = float(np.sum((1.0 - phi if freezing else phi)) * dx)
    ste = mats.cp * 100.0 / mats.Lf_J_kg
    lam = stefan_lambda(ste)
    alpha = mats.k / (mats.rho_ag_solid_kg_m3 * mats.cp)
    exact_length = 2.0 * lam * math.sqrt(alpha * duration_s)
    relative_error = abs(measured_length - exact_length) / exact_length
    residual = abs(float(last["energy_residual_J"]))
    energy_scale = max(abs(float(last["cumulative_boundary_out_J"])),
                       abs(float(last["total_stored_enthalpy_J"])), 1.0e-30)
    return {"type": "freeze" if freezing else "melt", "n": n, "dt_s": dt_s, "duration_s": duration_s,
            "steps": steps, "measured_front_m": measured_length, "stefan_front_m": exact_length,
            "front_relative_error": relative_error, "lambda": lam, "stefan_number": ste,
            "normalized_energy_residual": residual / energy_scale, "wall_seconds": elapsed,
            "maxrss_mib": maxrss_mib(), "max_liquid_fraction": float(phi.max()),
            "min_liquid_fraction": float(phi.min()), "final_max_T_K": float(sim.temperature_K.max()),
            "final_min_T_K": float(sim.temperature_K.min()), "pcg_iterations_last": last["pcg_iterations_total"],
            "nonlinear_iterations_last": last["nonlinear_iterations"]}


def run_stefan(rows: list[dict]) -> dict:
    duration = 1.0e-3
    results = []
    # 空间收敛：同一物理时刻和相同时间步，逐级加密。
    for freezing in (False, True):
        family = [run_stefan_case(freezing, n, 2.0e-6, duration) for n in (50, 100, 200)]
        # 时间步收敛：固定N=100并把dt连续减半。
        time_family = [run_stefan_case(freezing, 100, dt, duration) for dt in (4.0e-6, 2.0e-6, 1.0e-6)]
        results.extend(family)
        results.extend(time_family)
        label = "freezing" if freezing else "melting"
        best = family[-1]
        space_change = abs(family[-1]["measured_front_m"] - family[-2]["measured_front_m"]) / max(
            family[-1]["measured_front_m"], 1.0e-30)
        time_change = abs(time_family[-1]["measured_front_m"] - time_family[-2]["measured_front_m"]) / max(
            time_family[-1]["measured_front_m"], 1.0e-30)
        rows.append(result_row(f"D_E_stefan_{label}", best["front_relative_error"] <= 0.02
                               and space_change <= 0.03 and time_change <= 0.03
                               and best["normalized_energy_residual"] <= 1.0e-3,
                               best["front_relative_error"],
                               "front error <= 2%; last space/time changes <= 3%; energy <= 0.1%",
                               f"space_change_100_to_200={space_change:.5%}; time_change_dt2e-6_to1e-6={time_change:.5%}"))
    return {"cases": results}


def run_interface_and_uniform(rows: list[dict]) -> dict:
    # G：双层导热（恒温边界、串联界面热阻）解析解。
    n = 100
    dx = 1.0e-4
    ag = np.zeros(n)
    ag[: n // 2] = 1.0
    mats = ConstantAgBenchmark(phase_width=1.0)
    mats.k = 10.0
    # scaffold 的常导热率由 ConstantAgBenchmark 的混合函数给出为2 W/(m K)。
    r_int = 5.0e-4
    sim = ImplicitEnthalpyFVM(ag, (dx,), materials=mats, interface_resistance_m2_K_W=r_int,
                              initial_temperature_K=300.0, linear_relative_tolerance=1.0e-12,
                              maximum_nonlinear_iterations=40, maximum_linear_iterations=1000)
    boundary = BoundarySet(dirichlet_K={"x_low": 400.0, "x_high": 300.0})
    last = None
    for _ in range(20):
        last = sim.advance(5.0, 0.0, boundary)
    # 总厚度按半域分成两层；边界半控制体和内部中心距离合起来正好是层厚。
    length_ag = (n // 2) * dx
    length_sc = length_ag
    q_exact = 100.0 / (length_ag / 10.0 + r_int + length_sc / 2.0)
    k_field = mats.conductivity_W_m_K(sim.temperature_K, ag)
    high_face = sim._boundary_face_conductance_density(k_field, boundary)["x_high"]
    q_numeric = float(np.sum(high_face * (sim.temperature_K - 300.0)) * sim.volume_m3)
    flux_error = abs(q_numeric - q_exact) / q_exact
    # 界面两侧线性拟合到相界面中心，比较解析跳温。
    jump_exact = q_exact * r_int
    left_face = sim.temperature_K[n // 2 - 1] - q_numeric * (0.5 * dx) / 10.0
    right_face = sim.temperature_K[n // 2] + q_numeric * (0.5 * dx) / 2.0
    jump_numeric = left_face - right_face
    jump_error = abs(jump_numeric - jump_exact) / max(abs(jump_exact), 1.0e-12)
    rows.append(result_row("F_interface_resistance_steady", flux_error <= 0.01 and jump_error <= 0.02,
                           max(flux_error, jump_error), "flux error <= 1%; interface jump error <= 2%",
                           f"q_exact_W_m2={q_exact:.8g}; q_numeric_W_m2={q_numeric:.8g}; "
                           f"dT_interface_exact_K={jump_exact:.8g}; dT_interface_numeric_K={jump_numeric:.8g}"))

    # H1：300 K绝热均匀场不变。
    uniform = ImplicitEnthalpyFVM(np.full((6, 5, 4), 0.7), (1.0e-6,) * 3, initial_temperature_K=300.0)
    for _ in range(5):
        uniform.advance(1.0e-5)
    uniform_drift = float(np.max(np.abs(uniform.temperature_K - 300.0)))
    rows.append(result_row("G_uniform_300K_hold", uniform_drift <= 1.0e-10, uniform_drift,
                           "maximum drift <= 1e-10 K"))

    # H2：加热后停源，再由冷却边界带走能量；用整体能量平衡检验边界散热项。
    pulse = ImplicitEnthalpyFVM(np.ones(20), (5.0e-5,), initial_temperature_K=300.0)
    localized = np.zeros(20)
    localized[8:12] = 3.0e15
    pulse.advance(2.0e-7, localized)
    source_after_pulse = pulse.source_energy_J
    sink = BoundarySet(dirichlet_K={"x_low": 300.0, "x_high": 300.0})
    cooling_history = []
    for _ in range(100):
        cooling_history.append(pulse.advance(2.0e-7, 0.0, sink))
    final_storage = float(np.sum(pulse.enthalpy_density_J_m3) * pulse.volume_m3)
    total_balance = abs(final_storage - pulse.source_energy_J + pulse.boundary_energy_out_J)
    normalized = total_balance / max(source_after_pulse, 1.0e-30)
    sink_positive = pulse.boundary_energy_out_J > 0.0
    rows.append(result_row("H_source_off_boundary_energy_balance", sink_positive and normalized <= 1.0e-3,
                           normalized, "positive outflow; source-off energy residual <= 0.1%",
                           f"source_J={source_after_pulse:.8e}; sink_out_J={pulse.boundary_energy_out_J:.8e}; "
                           f"final_max_T_K={pulse.temperature_K.max():.8f}"))
    return {"interface_flux_error": flux_error, "interface_temperature_jump_error": jump_error,
            "uniform_temperature_drift_K": uniform_drift, "source_off_energy_residual_rel": normalized,
            "source_J": source_after_pulse, "sink_out_J": pulse.boundary_energy_out_J,
            "cooling_steps": len(cooling_history), "cooling_duration_s": pulse.time_s - 2.0e-7}


def main() -> None:
    rows: list[dict] = []
    started = time.perf_counter()
    results = {
        "material_roundtrip": run_material_roundtrip(rows),
        "single_control_volume": run_single_cell(rows),
        "stefan": run_stefan(rows),
        "interface_uniform_and_source_off": run_interface_and_uniform(rows),
    }
    results["checks"] = rows
    results["all_small_checks_pass"] = all(r["pass"] for r in rows)
    results["wall_seconds"] = time.perf_counter() - started
    results["maxrss_mib"] = maxrss_mib()
    (OUT / "small_benchmark_results.json").write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    with (OUT / "SMALL_BENCHMARKS.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    with (OUT / "STEFAN_FRONT_CONVERGENCE.csv").open("w", newline="", encoding="utf-8") as stream:
        records = results["stefan"]["cases"]
        writer = csv.DictWriter(stream, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps({"all_small_checks_pass": results["all_small_checks_pass"],
                      "wall_seconds": results["wall_seconds"], "maxrss_mib": results["maxrss_mib"],
                      "checks": rows}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
