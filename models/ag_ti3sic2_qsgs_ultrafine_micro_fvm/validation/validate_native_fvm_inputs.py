#!/usr/bin/env python3
"""用 single-temperature-study 分支原生网格入口检查无热源输入。"""

from __future__ import annotations

import json
from pathlib import Path
import resource
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOLVER_SOURCE = (
    ROOT / "geometry_generation" / "source_extracts" / "single-temperature-study" / "solver"
)
sys.path.insert(0, str(SOLVER_SOURCE))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid  # noqa: E402


MODEL_DIRS = (
    ROOT / "ultrafine_500nm" / "fvm_inputs",
    ROOT / "micro_10um" / "fvm_inputs",
)
OUTPUT = ROOT / "validation" / "native_fvm_input_validation.json"


def main() -> None:
    rows = []
    for model_dir in MODEL_DIRS:
        for geometry_path in sorted(model_dir.glob("*_native.npz")):
            with np.load(geometry_path, allow_pickle=False) as archive:
                skeleton = np.asarray(archive["skeleton"], dtype=bool)
                size_A = np.asarray(archive["size_A"], dtype=np.float64)
                size_m = np.asarray(archive["size_m"], dtype=np.float64)
                spacing_m = np.asarray(archive["spacing_m"], dtype=np.float64)
                shape_xyz = tuple(int(value) for value in archive["shape_xyz"])
                seed = int(archive["seed"])
            if skeleton.shape != shape_xyz or shape_xyz != (150, 150, 150):
                raise ValueError(f"{geometry_path}: native skeleton shape mismatch")
            if not np.allclose(size_A * 1.0e-10, size_m, rtol=1.0e-14, atol=0.0):
                raise ValueError(f"{geometry_path}: native size_A conversion failed")
            if not np.allclose(size_m / np.asarray(shape_xyz), spacing_m, rtol=1.0e-14, atol=0.0):
                raise ValueError(f"{geometry_path}: native cell spacing mismatch")
            for boundary in ("A", "B"):
                started = time.perf_counter()
                model = QSGSEnthalpyGrid(
                    None,
                    boundary,
                    peak_flux_eV_ps_A2=0.0,
                    native_geometry_npz=geometry_path,
                )
                if model.temperature_shape != shape_xyz:
                    raise ValueError(f"{geometry_path}: native solver shape mismatch")
                if not np.array_equal(model.sc_fraction, skeleton.astype(float)):
                    raise ValueError(f"{geometry_path}: native solver changed skeleton phase")
                if not np.allclose(model.size_A, size_A, rtol=0.0, atol=1.0e-12):
                    raise ValueError(f"{geometry_path}: native solver changed physical size")
                if not np.allclose(model.spacing_A * 1.0e-10, spacing_m, rtol=1.0e-14, atol=0.0):
                    raise ValueError(f"{geometry_path}: native solver SI spacing mismatch")
                if model.actual_source_power_W != 0.0 or np.any(model.source_W_m3 != 0.0):
                    raise ValueError(f"{geometry_path}: zero-flux native initialization added heat")
                result = model.advance(1.0e-18, source_on=False, enforce_stability=True)
                if not np.allclose(model.state.temperature_K, 300.0, rtol=0.0, atol=1.0e-10):
                    raise ValueError(f"{geometry_path}: uniform zero-source state changed")
                if np.any(model.state.ag_liquid_fraction != 0.0):
                    raise ValueError(f"{geometry_path}: liquid appeared without heating")
                if result["injected_energy_J"] != 0.0 or result["energy_balance_error_J"] != 0.0:
                    raise ValueError(f"{geometry_path}: zero-source energy ledger is not closed")
                row = {
                    "native_geometry": str(geometry_path.relative_to(ROOT)),
                    "seed": seed,
                    "boundary_case": boundary,
                    "shape_xyz": shape_xyz,
                    "ag_fraction": float(model.ag_fraction.mean()),
                    "skeleton_fraction": float(model.sc_fraction.mean()),
                    "size_m": size_m.tolist(),
                    "spacing_m": spacing_m.tolist(),
                    "size_A": size_A.tolist(),
                    "spacing_A": model.spacing_A.tolist(),
                    "source_power_W": model.actual_source_power_W,
                    "final_temperature_K": float(model.state.temperature_K.max()),
                    "liquid_fraction_max": float(model.state.ag_liquid_fraction.max()),
                    "energy_residual_J": result["energy_balance_error_J"],
                    "stable_timestep_s": model.stable_timestep_s(),
                    "elapsed_seconds": time.perf_counter() - started,
                    "process_peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
                }
                rows.append(row)
                print(json.dumps(row), flush=True)

    payload = {
        "status": "NATIVE_GEOMETRY_ZERO_SOURCE_CHECKS_PASSED",
        "solver_branch": "codex/single-temperature-study",
        "solver_source_commit": "703b63ccf840dc0f8d0e85349221c8457163417f",
        "solver_source_path": str(SOLVER_SOURCE.relative_to(ROOT)),
        "cases": rows,
        "heat_sources_applied": False,
        "flow_or_phase_change_run": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
