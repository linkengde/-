#!/usr/bin/env python3
"""运行纯相、串并联、界面热阻及离散守恒解析基准。"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import solve_banded


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
sys.path.insert(0, str(HERE))
from steady_periodic_fvm import PeriodicConductivityFVM  # noqa: E402

CONFIG = json.loads((HERE / "VALIDATION_CONFIG.json").read_text(encoding="utf-8"))

K_AG = float(CONFIG["materials_at_temperature_K"]["Ag"]["thermal_conductivity_W_m_K"])
K_TI3SIC2 = float(CONFIG["materials_at_temperature_K"]["Ti3SiC2"]["thermal_conductivity_W_m_K"])
N = int(CONFIG["validation_runs"]["analytic_periodic_grid_n"])
LENGTH_M = float(CONFIG["rves"]["ultrafine_500nm"]["domain_edge_m"])
DX = LENGTH_M / N
T_REF_K = float(CONFIG["boundary_and_solver"]["reference_temperature_K"])
DELTA_T_K = float(CONFIG["boundary_and_solver"]["imposed_macro_temperature_difference_K"])
FINITE_INTERFACE_R_M2_K_W = 1.0e-9  # 数值验证参数，不是实测Ag/Ti3SiC2界面数据


def relative_error(value: float, exact: float) -> float:
    return abs(value - exact) / abs(exact)


def tensor_for(mask: np.ndarray, resistance: float = 0.0) -> np.ndarray:
    solver = PeriodicConductivityFVM(
        mask,
        (DX, DX, DX),
        K_AG,
        K_TI3SIC2,
        resistance,
        fft_workers=2,
    )
    result = solver.solve_effective_tensor(relative_tolerance=1.0e-11)
    return np.asarray(result["effective_conductivity_tensor_W_m_K"])


def solve_1d_dirichlet(
    conductivity: np.ndarray,
    length_m: float,
    temperature_hot_K: float,
    temperature_cold_K: float,
    interface_resistance_m2_K_W: float = 0.0,
) -> dict[str, object]:
    """一维单元中心FVM，边界面Dirichlet，界面面热阻串联。"""
    k = np.asarray(conductivity, dtype=float)
    n = k.size
    dx = length_m / n
    left_g = np.empty(n, dtype=float)
    right_g = np.empty(n, dtype=float)
    left_g[0] = 2.0 * k[0] / dx
    right_g[-1] = 2.0 * k[-1] / dx
    for i in range(1, n):
        resistance = 0.5 * dx / k[i - 1] + 0.5 * dx / k[i]
        if not np.isclose(k[i - 1], k[i], rtol=0.0, atol=0.0):
            resistance += interface_resistance_m2_K_W
        face_g = 1.0 / resistance
        right_g[i - 1] = face_g
        left_g[i] = face_g
    ab = np.zeros((3, n), dtype=float)
    ab[1] = left_g + right_g
    ab[0, 1:] = -left_g[1:]
    ab[2, :-1] = -right_g[:-1]
    rhs = np.zeros(n, dtype=float)
    rhs[0] += left_g[0] * temperature_hot_K
    rhs[-1] += right_g[-1] * temperature_cold_K
    temperature = solve_banded((1, 1), ab, rhs, overwrite_ab=True, overwrite_b=True)
    q_left = left_g[0] * (temperature_hot_K - temperature[0])
    q_right = right_g[-1] * (temperature[-1] - temperature_cold_K)
    q_internal = right_g[:-1] * (temperature[:-1] - temperature[1:])
    all_flux = np.concatenate(([q_left], q_internal, [q_right]))
    flux_scale = max(float(np.mean(np.abs(all_flux))), 1.0)
    max_flux_spread = float((all_flux.max() - all_flux.min()) / flux_scale)
    total_resistance = (
        0.5 * dx / k[0]
        + sum(0.5 * dx / k[i - 1] + 0.5 * dx / k[i]
              + (interface_resistance_m2_K_W if not np.isclose(k[i - 1], k[i], rtol=0.0, atol=0.0) else 0.0)
              for i in range(1, n))
        + 0.5 * dx / k[-1]
    )
    exact_flux = (temperature_hot_K - temperature_cold_K) / total_resistance
    interface_jump = []
    expected_jump = []
    for face in np.flatnonzero(k[:-1] != k[1:]):
        face_flux = q_internal[face]
        hot_side_face_temperature = temperature[face] - face_flux * (0.5 * dx) / k[face]
        cold_side_face_temperature = temperature[face + 1] + face_flux * (0.5 * dx) / k[face + 1]
        interface_jump.append(float(hot_side_face_temperature - cold_side_face_temperature))
        expected_jump.append(float(face_flux * interface_resistance_m2_K_W))
    jump_errors = [abs(actual - expected) for actual, expected in zip(interface_jump, expected_jump)]
    return {
        "mean_heat_flux_W_m2": float(np.mean(all_flux)),
        "exact_heat_flux_W_m2": float(exact_flux),
        "relative_flux_error": abs(float(np.mean(all_flux)) - exact_flux) / abs(exact_flux),
        "normalized_flux_nonuniformity": max_flux_spread,
        "temperature_min_K": float(temperature.min()),
        "temperature_max_K": float(temperature.max()),
        "interface_temperature_jump_K": interface_jump,
        "expected_interface_temperature_jump_qR_K": expected_jump,
        "maximum_interface_temperature_jump_error_K": max(jump_errors, default=0.0),
        "interface_resistance_m2_K_W": interface_resistance_m2_K_W,
    }


def main() -> None:
    started = time.perf_counter()
    rows: list[dict[str, object]] = []

    one_dimensional_cases = []
    for name, conductivity, resistance in (
        ("pure_Ag_1D_Dirichlet", np.full(N, K_AG), 0.0),
        ("pure_Ti3SiC2_1D_Dirichlet", np.full(N, K_TI3SIC2), 0.0),
        ("Ag_Ti3SiC2_series_1D_Dirichlet_R0", np.r_[np.full(N // 2, K_AG), np.full(N // 2, K_TI3SIC2)], 0.0),
        ("Ag_Ti3SiC2_series_1D_Dirichlet_finite_R", np.r_[np.full(N // 2, K_AG), np.full(N // 2, K_TI3SIC2)], FINITE_INTERFACE_R_M2_K_W),
    ):
        result_1d = solve_1d_dirichlet(conductivity, LENGTH_M, T_REF_K + DELTA_T_K, T_REF_K, resistance)
        result_1d["case"] = name
        result_1d["status"] = "PASS" if result_1d["relative_flux_error"] <= 0.01 and result_1d["normalized_flux_nonuniformity"] <= 1e-10 else "FAIL"
        one_dimensional_cases.append(result_1d)
    rows.extend(one_dimensional_cases)

    for name, skeleton, expected in (
        ("pure_Ag_1D", np.zeros((N,) * 3, dtype=bool), K_AG),
        ("pure_Ti3SiC2_1D", np.ones((N,) * 3, dtype=bool), K_TI3SIC2),
    ):
        tensor = tensor_for(skeleton)
        errors = [relative_error(float(tensor[i, i]), expected) for i in range(3)]
        rows.append({
            "case": name,
            "expected_W_m_K": expected,
            "tensor_W_m_K": tensor.tolist(),
            "maximum_relative_diagonal_error": max(errors),
            "status": "PASS" if max(errors) <= 0.01 else "FAIL",
        })

    layers_x = np.broadcast_to(
        (np.arange(N)[:, None, None] >= N // 2), (N, N, N)
    ).copy()
    tensor_series = tensor_for(layers_x)
    expected_series = 1.0 / (0.5 / K_AG + 0.5 / K_TI3SIC2)
    expected_parallel = 0.5 * (K_AG + K_TI3SIC2)
    rows.append({
        "case": "Ag_Ti3SiC2_series_and_parallel_laminate_R0",
        "expected_series_x_W_m_K": expected_series,
        "expected_parallel_yz_W_m_K": expected_parallel,
        "tensor_W_m_K": tensor_series.tolist(),
        "series_relative_error": relative_error(float(tensor_series[0, 0]), expected_series),
        "parallel_relative_error": max(
            relative_error(float(tensor_series[1, 1]), expected_parallel),
            relative_error(float(tensor_series[2, 2]), expected_parallel),
        ),
        "maximum_off_diagonal_W_m_K": float(
            np.max(np.abs(tensor_series - np.diag(np.diag(tensor_series))))
        ),
        "status": "PASS" if max(
            relative_error(float(tensor_series[0, 0]), expected_series),
            relative_error(float(tensor_series[1, 1]), expected_parallel),
            relative_error(float(tensor_series[2, 2]), expected_parallel),
        ) <= 0.01 else "FAIL",
    })

    tensor_y_layers = tensor_for(np.swapaxes(layers_x, 0, 1))
    rows.append({
        "case": "Ag_Ti3SiC2_parallel_x_and_series_y_laminate_R0",
        "expected_series_y_W_m_K": expected_series,
        "expected_parallel_xz_W_m_K": expected_parallel,
        "tensor_W_m_K": tensor_y_layers.tolist(),
        "maximum_diagonal_relative_error": max(
            relative_error(float(tensor_y_layers[1, 1]), expected_series),
            relative_error(float(tensor_y_layers[0, 0]), expected_parallel),
            relative_error(float(tensor_y_layers[2, 2]), expected_parallel),
        ),
        "status": "PASS" if max(
            relative_error(float(tensor_y_layers[1, 1]), expected_series),
            relative_error(float(tensor_y_layers[0, 0]), expected_parallel),
            relative_error(float(tensor_y_layers[2, 2]), expected_parallel),
        ) <= 0.01 else "FAIL",
    })

    tensor_finite = tensor_for(layers_x, FINITE_INTERFACE_R_M2_K_W)
    expected_finite_series = 1.0 / (
        0.5 / K_AG + 0.5 / K_TI3SIC2 + 2.0 * FINITE_INTERFACE_R_M2_K_W / LENGTH_M
    )
    rows.append({
        "case": "Ag_Ti3SiC2_series_laminate_finite_R",
        "interface_resistance_m2_K_W": FINITE_INTERFACE_R_M2_K_W,
        "expected_series_x_W_m_K": expected_finite_series,
        "expected_parallel_yz_W_m_K": expected_parallel,
        "tensor_W_m_K": tensor_finite.tolist(),
        "maximum_relative_error": max(
            relative_error(float(tensor_finite[0, 0]), expected_finite_series),
            relative_error(float(tensor_finite[1, 1]), expected_parallel),
            relative_error(float(tensor_finite[2, 2]), expected_parallel),
        ),
        "status": "PASS" if max(
            relative_error(float(tensor_finite[0, 0]), expected_finite_series),
            relative_error(float(tensor_finite[1, 1]), expected_parallel),
            relative_error(float(tensor_finite[2, 2]), expected_parallel),
        ) <= 0.01 else "FAIL",
    })

    conservation_solver = PeriodicConductivityFVM(
        layers_x, (DX, DX, DX), K_AG, K_TI3SIC2, fft_workers=2
    )
    uniform_300 = np.full((N,) * 3, T_REF_K)
    uniform_residual = float(np.max(np.abs(conservation_solver.apply_operator(uniform_300))))
    random_field = np.random.default_rng(20261128).normal(size=(N,) * 3)
    divergence = conservation_solver.apply_operator(random_field)
    global_balance = float(np.sum(divergence))
    normalized_balance = abs(global_balance) / max(float(np.sum(np.abs(divergence))), 1.0)
    rows.extend([
        {
            "case": "no_source_uniform_300K_hold",
            "maximum_operator_residual_W": uniform_residual,
            "status": "PASS" if uniform_residual == 0.0 else "FAIL",
        },
        {
            "case": "periodic_adiabatic_global_energy_balance",
            "sum_divergence_W": global_balance,
            "normalized_global_balance": normalized_balance,
            "status": "PASS" if normalized_balance <= 1.0e-12 else "FAIL",
        },
    ])

    result = {
        "status": "PASS" if all(row["status"] == "PASS" for row in rows) else "FAIL",
        "solver": "periodic cell-centred conservative FV with FFT-preconditioned CG",
        "grid_shape": [N, N, N],
        "domain_length_m": [LENGTH_M] * 3,
        "temperature_reference_K": T_REF_K,
        "temperature_delta_K_for_linear_cases": DELTA_T_K,
        "thermal_conductivity_inputs_W_m_K": {"Ag": K_AG, "Ti3SiC2": K_TI3SIC2},
        "finite_interface_resistance_test_parameter_m2_K_W": FINITE_INTERFACE_R_M2_K_W,
        "finite_R_parameter_is_measured_data": False,
        "acceptance_relative_error": 0.01,
        "cases": rows,
        "wall_seconds": time.perf_counter() - started,
    }
    output = HERE / "benchmarks" / "analytic_benchmark_results.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), **result}, indent=2))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
