#!/usr/bin/env python3
"""1-D one-phase Stefan melting benchmark for the Ag enthalpy curve."""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from thermophysical_model import AgEnthalpy


@dataclass(frozen=True)
class StefanResult:
    dx_m: float
    cell_count: int
    end_time_s: float
    time_steps: int
    boundary_temperature_K: float
    analytical_lambda: float
    analytical_interface_m: float
    numerical_liquid_equivalent_length_m: float
    relative_interface_error: float
    boundary_energy_input_J_m2: float
    enthalpy_change_J_m2: float
    energy_balance_residual_J_m2: float


def analytical_lambda(stefan_number: float) -> float:
    """Solve Ste = sqrt(pi)*lambda*exp(lambda^2)*erf(lambda)."""
    if not np.isfinite(stefan_number) or stefan_number <= 0.0:
        raise ValueError("Stefan number must be finite and positive")

    def residual(value: float) -> float:
        return math.sqrt(math.pi) * value * math.exp(value * value) * math.erf(value) - stefan_number

    low, high = 0.0, 2.0
    while residual(high) < 0.0:
        high *= 1.5
        if high > 8.0:
            raise RuntimeError("Could not bracket the Stefan similarity parameter")
    for _ in range(80):
        middle = 0.5 * (low + high)
        if residual(middle) < 0.0:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def run_one_phase_melting(
    dx_m: float = 40.0e-6,
    end_time_s: float = 0.01,
    boundary_superheat_K: float = 100.0,
    cell_count: int = 250,
    enthalpy: AgEnthalpy | None = None,
) -> StefanResult:
    """Melt semi-infinite Ag initially at Tm using a fixed hot boundary.

    The far end is adiabatic and must remain far from the melting front. The
    reported front is the sum of per-cell liquid fractions times dx. Density is
    the first-pass fixed solid-reference density on both sides of the front.
    """
    if enthalpy is None:
        enthalpy = AgEnthalpy()
    if not np.isfinite(dx_m) or dx_m <= 0.0 or cell_count < 2:
        raise ValueError("dx_m and cell_count must be positive")
    if not np.isfinite(end_time_s) or end_time_s <= 0.0:
        raise ValueError("end_time_s must be finite and positive")
    if not np.isfinite(boundary_superheat_K) or boundary_superheat_K <= 0.0:
        raise ValueError("boundary_superheat_K must be finite and positive")

    rho = enthalpy.reference_density_kg_m3
    tm = enthalpy.melting_temperature_K
    tb = tm + boundary_superheat_K
    h_initial = float(enthalpy.specific_enthalpy_J_kg(tm, 0.0))
    energy = np.full(cell_count, rho * h_initial, dtype=float)  # J/m^3 relative to 300 K

    cp_min = min(
        float(enthalpy.cp_solid_J_kg_K(tm)),
        float(enthalpy.cp_liquid_J_kg_K(tm)),
    )
    k_max = float(max(enthalpy.solid_k_W_m_K.max(), enthalpy.liquid_k_W_m_K.max()))
    dt_stable = 0.25 * rho * cp_min * dx_m**2 / k_max
    steps = int(math.ceil(end_time_s / dt_stable))
    dt = end_time_s / steps
    boundary_input = 0.0
    time = 0.0

    while time < end_time_s - 1.0e-15:
        step = min(dt, end_time_s - time)
        state = enthalpy.state_from_specific_enthalpy(energy / rho)
        temperature = np.asarray(state.temperature_K, dtype=float)
        liquid = np.asarray(state.liquid_fraction, dtype=float)
        conductivity = np.asarray(enthalpy.conductivity_W_m_K(temperature, liquid), dtype=float)
        face_k = 2.0 * conductivity[:-1] * conductivity[1:] / (
            conductivity[:-1] + conductivity[1:]
        )
        flow_right = np.empty(cell_count + 1, dtype=float)  # W/m^2, positive toward +x
        flow_right[0] = 2.0 * conductivity[0] * (tb - temperature[0]) / dx_m
        flow_right[1:-1] = face_k * (temperature[:-1] - temperature[1:]) / dx_m
        flow_right[-1] = 0.0
        energy += step * (flow_right[:-1] - flow_right[1:]) / dx_m
        boundary_input += step * flow_right[0]
        time += step

    final_state = enthalpy.state_from_specific_enthalpy(energy / rho)
    liquid_length = float(np.sum(np.asarray(final_state.liquid_fraction)) * dx_m)
    cp_liquid = float(enthalpy.cp_liquid_J_kg_K(tm))
    stefan_number = cp_liquid * boundary_superheat_K / enthalpy.latent_heat_J_kg
    lam = analytical_lambda(stefan_number)
    alpha_liquid = float(enthalpy.liquid_k_W_m_K[0]) / (rho * cp_liquid)
    analytical_front = 2.0 * lam * math.sqrt(alpha_liquid * end_time_s)
    delta_energy = float((energy - rho * h_initial).sum() * dx_m)
    residual = delta_energy - boundary_input
    return StefanResult(
        dx_m=dx_m,
        cell_count=cell_count,
        end_time_s=end_time_s,
        time_steps=steps,
        boundary_temperature_K=tb,
        analytical_lambda=lam,
        analytical_interface_m=analytical_front,
        numerical_liquid_equivalent_length_m=liquid_length,
        relative_interface_error=(liquid_length - analytical_front) / analytical_front,
        boundary_energy_input_J_m2=boundary_input,
        enthalpy_change_J_m2=delta_energy,
        energy_balance_residual_J_m2=residual,
    )


def main() -> None:
    enthalpy = AgEnthalpy()
    results = [
        run_one_phase_melting(40.0e-6, 0.01, 100.0, 250, enthalpy),
        run_one_phase_melting(20.0e-6, 0.01, 100.0, 500, enthalpy),
    ]
    payload = {
        "status": "AG_ENTHALPY_ONE_PHASE_STEFAN_BENCHMARK",
        "model": "1D finite-volume enthalpy conduction; fixed-temperature boundary at Tm+100 K; initial solid at Tm; insulated far boundary",
        "fixed_density_kg_m3": enthalpy.reference_density_kg_m3,
        "results": [asdict(result) for result in results],
        "limitations": [
            "The benchmark validates latent-heat bookkeeping and front propagation, not the QSGS composite or its scaffold boundary.",
            "It uses the first-pass candidate Cp/k interpolation and fixed solid-reference density.",
            "The one-phase similarity solution assumes the solid ahead of the front remains at Tm.",
        ],
    }
    report = Path(__file__).resolve().parent / "stefan_benchmark_report.json"
    report.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
