#!/usr/bin/env python3
"""对两个RVE做10步局部表面高斯热流冒烟测试；不是正式电弧工况。"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parents[1]
PHASE_FIELD = REPOSITORY / "models" / "ag_ti3sic2_qsgs_1to1_finer" / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_FIELD))
from qsgs_enthalpy_solver import (  # noqa: E402
    EV_PER_PS_ANGSTROM2_TO_W_M2,
    QSGSEnthalpyGrid,
)

CASES = {
    "ultrafine_500nm": ROOT / "ultrafine_500nm" / "fvm_inputs",
    "micro_10um": ROOT / "micro_10um" / "fvm_inputs",
}
SEED = 1105
BOUNDARY = "B"
PEAK_FLUX_W_M2 = 1.0e12
STEP_S = 0.5e-12
STEP_COUNT = 10


def main() -> None:
    summaries = []
    for model_name, input_dir in CASES.items():
        started = time.perf_counter()
        projected = input_dir / f"qsgs_seed_{SEED}_fvm.npz"
        source_geometry = input_dir / f"qsgs_seed_{SEED}_native.npz"
        with np.load(projected, allow_pickle=False) as f:
            size_A = np.asarray(f["size_A"], dtype=float)
        gaussian_radius_A = float(size_A[0] * 0.10)
        absorption_depth_A = float(size_A[2] * 0.01)
        solver = QSGSEnthalpyGrid(
            projected,
            BOUNDARY,
            peak_flux_eV_ps_A2=PEAK_FLUX_W_M2 / EV_PER_PS_ANGSTROM2_TO_W_M2,
            gaussian_radius_A=gaussian_radius_A,
            absorption_depth_A=absorption_depth_A,
            source_geometry_npz=source_geometry,
        )
        stable_dt_s = solver.stable_timestep_s()
        if STEP_S >= stable_dt_s:
            raise RuntimeError(f"{model_name}: test step is not below the explicit stability limit")
        if solver.actual_source_power_W <= 0.0:
            raise RuntimeError(f"{model_name}: positive boundary flux produced no heat-source power")
        rows = [solver.advance(STEP_S, source_on=True, enforce_stability=True) for _ in range(STEP_COUNT)]
        temperature = solver.state.temperature_K
        liquid = solver.state.ag_liquid_fraction
        enthalpy_J = float(solver.enthalpy_density_J_m3.sum() * solver.cell_volume_m3)
        residual_J = enthalpy_J - solver.injected_energy_J
        if not np.isfinite(temperature).all() or np.any(temperature < 300.0):
            raise RuntimeError(f"{model_name}: non-finite or below-initial temperatures after heating")
        if not np.isfinite(liquid).all() or np.any((liquid < 0.0) | (liquid > 1.0)):
            raise RuntimeError(f"{model_name}: liquid fraction escaped [0, 1]")
        if float(temperature.max()) <= 300.0:
            raise RuntimeError(f"{model_name}: heating failed to raise temperature")
        if abs(residual_J) > max(1.0e-25, abs(solver.injected_energy_J) * 1.0e-8):
            raise RuntimeError(f"{model_name}: global energy balance failed")

        step_csv = ROOT / "validation" / f"short_heating_{model_name}.csv"
        with step_csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        row = {
            "model": model_name,
            "seed": SEED,
            "boundary_case": BOUNDARY,
            "grid_shape_xyz": list(solver.temperature_shape),
            "size_A": solver.size_A.tolist(),
            "size_m": (solver.size_A * 1.0e-10).tolist(),
            "initial_temperature_K": 300.0,
            "peak_surface_gaussian_flux_W_m2": PEAK_FLUX_W_M2,
            "source_footprint_definition": "Gaussian radial radius = 0.10 of RVE x length; exponential depth = 0.01 of RVE z length",
            "gaussian_radius_A": gaussian_radius_A,
            "absorption_depth_A": absorption_depth_A,
            "source_geometry_integration": solver.source_integration_method,
            "source_power_W": solver.actual_source_power_W,
            "stable_timestep_limit_s": stable_dt_s,
            "test_timestep_s": STEP_S,
            "test_step_count": STEP_COUNT,
            "test_duration_s": solver.time_s,
            "injected_energy_J": solver.injected_energy_J,
            "final_enthalpy_J": enthalpy_J,
            "energy_residual_J": residual_J,
            "temperature_min_K": float(temperature.min()),
            "temperature_mean_K": float(temperature.mean()),
            "temperature_max_K": float(temperature.max()),
            "max_ag_liquid_fraction": float(liquid.max()),
            "liquid_fraction_range_valid": bool(np.all((liquid >= 0.0) & (liquid <= 1.0))),
            "elapsed_seconds": time.perf_counter() - started,
            "process_peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
            "step_history_csv": step_csv.name,
            "interpretation": "synthetic local-source code smoke only; flux and footprint are numerical test inputs, not measured arc conditions",
        }
        summaries.append(row)
        print(json.dumps(row), flush=True)
    payload = {
        "status": "SHORT_LOCAL_HEATING_SMOKE_PASSED",
        "heat_sources_applied": True,
        "flow_or_evaporation_run": False,
        "models": summaries,
        "limitations": [
            "Only 10 stable explicit steps were advanced from a uniform 300 K state.",
            "The Gaussian footprint is scaled with each RVE and must not be interpreted as the same physical arc source.",
            "No mesh convergence, measured heat boundary, finite interface resistance, flow, or full melting/cooling cycle was tested.",
        ],
    }
    (ROOT / "validation" / "short_heating_smoke.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
