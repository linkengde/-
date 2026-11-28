#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用Pankratz v2物性验证Stefan前沿、界面热阻及能量守恒。"""
from __future__ import annotations

import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import erf

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
PHASE = PROJECT / "phase_change_validation"
SOLVER = PROJECT / "solver_optimization"
for _path in (str(SOLVER), str(PHASE), str(PROJECT)):
    if _path in sys.path:
        sys.path.remove(_path)
    sys.path.insert(0, _path)

from implicit_enthalpy_fvm import BoundarySet, ImplicitEnthalpyFVM
from material_database_v2.material_model import PankratzMaterials

OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)


def stefan_lambda(ste: float) -> float:
    """Solve the one-phase semi-infinite Stefan similarity relation."""
    return float(brentq(lambda lam: math.sqrt(math.pi) * lam * math.exp(lam * lam) * erf(lam) - ste,
                        1e-12, 3.0))


def run_stefan_case(mode: str, n: int, dt_s: float, duration_s: float = 1.0e-5) -> dict:
    materials = PankratzMaterials(1.0)
    dx = 128.0e-6 / n
    tm = materials.Tm_K
    initial = materials.Ts_K - 1.0e-6 if mode == "melt" else materials.Tl_K + 1.0e-6
    boundary_temperature = tm + 50.0 if mode == "melt" else 300.0
    cp = materials.ag_cp_liquid_J_kg_K if mode == "melt" else materials.ag_cp_solid_mid
    k = materials.ag_k_liquid_W_m_K if mode == "melt" else materials.ag_k[-1]
    delta_t = boundary_temperature - tm if mode == "melt" else tm - boundary_temperature
    ste = cp * delta_t / materials.Lf_J_kg
    lam = stefan_lambda(ste)
    alpha = k / (materials.rho_ag_solid_kg_m3 * cp)
    sim = ImplicitEnthalpyFVM(
        np.ones(n, dtype=np.float64), (dx,), materials=materials,
        initial_temperature_K=initial, nonlinear_relative_tolerance=1e-8,
        linear_relative_tolerance=1e-9, maximum_nonlinear_iterations=40,
        maximum_linear_iterations=400)
    boundary = BoundarySet(dirichlet_K={"x_low": boundary_temperature})
    started = time.perf_counter()
    retries = 0
    accepted = 0
    nominal_dt = dt_s
    trial_dt = nominal_dt
    while sim.time_s < duration_s - max(1e-18, duration_s * 1e-14):
        trial_dt = min(trial_dt, duration_s - sim.time_s)
        try:
            sim.advance(trial_dt, 0.0, boundary)
            accepted += 1
            trial_dt = min(nominal_dt, trial_dt * 1.4)
        except RuntimeError:
            retries += 1
            trial_dt *= 0.5
            if trial_dt < nominal_dt / 4096.0:
                raise
    liquid = materials.liquid_fraction(sim.temperature_K)
    if mode == "melt":
        active = np.flatnonzero(liquid >= 0.5)
        if active.size and active[-1] + 1 < n:
            j = int(active[-1])
            left_x, right_x = (j + 0.5) * dx, (j + 1.5) * dx
            f0, f1 = float(liquid[j]), float(liquid[j + 1])
            measured = left_x + (f0 - 0.5) / max(f0 - f1, 1e-30) * (right_x - left_x)
        else:
            measured = float((active[-1] + 0.5) * dx) if active.size else 0.0
    else:
        active = np.flatnonzero(liquid >= 0.5)
        if active.size and active[0] > 0:
            j = int(active[0])
            left_x, right_x = (j - 0.5) * dx, (j + 0.5) * dx
            f0, f1 = float(liquid[j - 1]), float(liquid[j])
            measured = left_x + (0.5 - f0) / max(f1 - f0, 1e-30) * (right_x - left_x)
        else:
            measured = float(active[0] * dx) if active.size else n * dx
    analytic = 2.0 * lam * math.sqrt(alpha * duration_s)
    energy_change = float(np.sum(sim.enthalpy_density_J_m3) * sim.volume_m3
                          - sim.initial_stored_enthalpy_J)
    residual = energy_change - sim.source_energy_J + sim.boundary_energy_out_J
    scale = max(abs(energy_change), abs(sim.boundary_energy_out_J), 1e-30)
    return {
        "mode": mode, "n": n, "dx_m": dx, "dt_nominal_s": nominal_dt,
        "duration_s": duration_s, "accepted_steps": accepted, "retries": retries,
        "Stefan_number": ste, "lambda": lam, "analytic_front_m": analytic,
        "measured_front_m": measured,
        "relative_front_error": abs(measured - analytic) / max(analytic, 1e-30),
        "liquid_fraction_min": float(liquid.min()), "liquid_fraction_max": float(liquid.max()),
        "max_temperature_K": float(sim.temperature_K.max()),
        "normalized_energy_residual": abs(residual) / scale,
        "wall_seconds": time.perf_counter() - started,
    }


def run_interface_case(resistance_m2_K_W: float) -> dict:
    """Solve a 1-D Ag/Ti3SiC2 bilayer and recover the interface temperature jump."""
    n = 100
    dx = 1.0e-6
    ag = np.zeros(n, dtype=np.float64)
    ag[: n // 2] = 1.0
    materials = PankratzMaterials(1.0)
    sim = ImplicitEnthalpyFVM(
        ag, (dx,), materials=materials, interface_resistance_m2_K_W=resistance_m2_K_W,
        initial_temperature_K=300.0, nonlinear_relative_tolerance=1e-9,
        linear_relative_tolerance=1e-10, maximum_nonlinear_iterations=40,
        maximum_linear_iterations=400)
    boundary = BoundarySet(dirichlet_K={"x_low": 400.0, "x_high": 300.0})
    started = time.perf_counter()
    records = []
    for _ in range(12):
        records.append(sim.advance(5.0e-3, 0.0, boundary))
    conductivity = materials.conductivity_W_m_K(sim.temperature_K, ag)
    faces, _ = sim._face_conductances_density(conductivity, boundary)
    interface_flux_W_m2 = float(faces[0][n // 2 - 1]
                                 * (sim.temperature_K[n // 2 - 1] - sim.temperature_K[n // 2]) * dx)
    t_ag_face = float(sim.temperature_K[n // 2 - 1]
                       - interface_flux_W_m2 * (0.5 * dx) / conductivity[n // 2 - 1])
    t_ti_face = float(sim.temperature_K[n // 2]
                       + interface_flux_W_m2 * (0.5 * dx) / conductivity[n // 2])
    jump = t_ag_face - t_ti_face
    expected_jump = interface_flux_W_m2 * resistance_m2_K_W
    total_storage = float(np.sum(sim.enthalpy_density_J_m3) * sim.volume_m3
                          - sim.initial_stored_enthalpy_J)
    residual = total_storage - sim.source_energy_J + sim.boundary_energy_out_J
    scale = max(abs(total_storage), abs(sim.boundary_energy_out_J), 1e-30)
    return {
        "interface_resistance_m2_K_W": resistance_m2_K_W,
        "interface_flux_W_m2": interface_flux_W_m2,
        "Ag_side_face_temperature_K": t_ag_face,
        "Ti3SiC2_side_face_temperature_K": t_ti_face,
        "temperature_jump_K": jump, "predicted_jump_qR_K": expected_jump,
        "jump_relative_error": abs(jump - expected_jump) / max(abs(expected_jump), 1e-30)
            if resistance_m2_K_W else abs(jump),
        "last_step_max_temperature_change_K": float(np.max(np.abs(
            records[-1]["temperature_change_K"]))) if "temperature_change_K" in records[-1] else None,
        "normalized_energy_residual": abs(residual) / scale,
        "wall_seconds": time.perf_counter() - started,
    }


def write_csv(path: Path, rows: list[dict]):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    started = time.perf_counter()
    stefan_rows = []
    for mode in ("melt", "freeze"):
        for n in (64, 128, 256):
            stefan_rows.append(run_stefan_case(mode, n, 2.0e-7))
        for dt in (4.0e-7, 2.0e-7, 1.0e-7):
            stefan_rows.append(run_stefan_case(mode, 128, dt))
    interface_rows = [run_interface_case(value) for value in (0.0, 1e-10, 1e-9, 1e-8)]
    width_rows = []
    for width in (0.5, 1.0, 2.0):
        materials = PankratzMaterials(width)
        sensible = 0.5 * (materials.ag_cp_solid_mid + materials.ag_cp_liquid_J_kg_K) * width
        jump = materials.ag_specific_enthalpy_J_kg(materials.Tl_K) - materials.ag_specific_enthalpy_J_kg(materials.Ts_K)
        width_rows.append({"mushy_width_K": width, "latent_heat_J_kg": materials.Lf_J_kg,
                           "enthalpy_jump_J_kg": jump, "sensible_in_interval_J_kg": sensible,
                           "latent_integral_error_J_kg": jump - sensible - materials.Lf_J_kg})
    payload = {
        "status": "COMPLETED",
        "stefan_cases": stefan_rows,
        "interface_resistance_cases": interface_rows,
        "mushy_width_cases": width_rows,
        "wall_seconds": time.perf_counter() - started,
        "interpretation": "One-dimensional analytical/numerical verification with v2 materials; not a microstructure or arc prediction.",
        "acceptance_targets": {"normalized_energy_residual_max": 1e-3,
                               "liquid_fraction_range": [0.0, 1.0],
                               "interface_jump_matches_qR": "within numerical roundoff"},
    }
    (OUT / "material_validation_suite.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    write_csv(OUT / "STEFAN_V2_VALIDATION.csv", stefan_rows)
    write_csv(OUT / "INTERFACE_RESISTANCE_V2.csv", interface_rows)
    write_csv(OUT / "MUSHY_WIDTH_V2.csv", width_rows)
    print(json.dumps({"wall_seconds": payload["wall_seconds"],
                      "max_energy_residual": max(r["normalized_energy_residual"] for r in stefan_rows + interface_rows),
                      "stefan": stefan_rows, "interface": interface_rows}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
