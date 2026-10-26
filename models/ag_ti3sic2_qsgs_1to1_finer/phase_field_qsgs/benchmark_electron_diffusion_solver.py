#!/usr/bin/env python3
"""Compare strict and optimized CG tolerances on one QSGS diffusion step."""

from __future__ import annotations

import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy
from scipy.sparse.linalg import cg as scipy_cg

HERE = Path(__file__).resolve().parent
MODEL_DIR = HERE.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import qsgs_ttm_enthalpy as ttm_module  # noqa: E402
from qsgs_ttm_enthalpy import QSGSDoubleTemperatureGrid  # noqa: E402


GRID_FILE = HERE / "qsgs_phase_fractions_10A_nested.npz"
CONTACT_FILE = HERE / "qsgs_ag_face_contacts_10A_nested.npz"
GEOMETRY_FILE = MODEL_DIR / "qsgs_morphology.npz"
TILE_XY = 3
DT_PS = 0.0003125


def build_model() -> QSGSDoubleTemperatureGrid:
    model = QSGSDoubleTemperatureGrid(
        GRID_FILE,
        "A",
        peak_flux_eV_ps_A2=0.0,
        face_contacts_npz=CONTACT_FILE,
        tile_xy=TILE_XY,
        source_geometry_npz=GEOMETRY_FILE,
    )
    x = (np.arange(model.electron.temperature_K.shape[0]) + 0.5) * model.lattice.spacing_A[0]
    y = (np.arange(model.electron.temperature_K.shape[1]) + 0.5) * model.lattice.spacing_A[1]
    z = (np.arange(model.electron.temperature_K.shape[2]) + 0.5) * model.lattice.spacing_A[2]
    xx, yy = np.meshgrid(x, y, indexing="ij")
    r2 = (xx - model.lattice.size_A[0] / 2.0) ** 2 + (
        yy - model.lattice.size_A[1] / 2.0
    ) ** 2
    field = 300.0 + 800.0 * np.exp(-r2[:, :, None] / 30.0**2) * np.exp(
        -z[None, None, :] / 20.0
    )
    active = model.lattice.ag_fraction > 0.0
    model.electron.temperature_K[active] = field[active]
    return model


def run_case(rtol: float) -> tuple[dict, np.ndarray]:
    model = build_model()
    iterations: list[int] = []
    original_cg = ttm_module.cg
    target_rtol = float(rtol)

    def counted_cg(A, b, M=None, rtol=1.0e-5, atol=0.0, maxiter=None, **kwargs):
        count = [0]
        solution, info = scipy_cg(
            A,
            b,
            M=M,
            rtol=target_rtol,
            atol=atol,
            maxiter=maxiter,
            callback=lambda _x: count.__setitem__(0, count[0] + 1),
            **kwargs,
        )
        iterations.append(count[0])
        return solution, info

    ttm_module.cg = counted_cg
    start = time.perf_counter()
    try:
        net_energy_change_J = model._implicit_electron_diffusion(DT_PS)
    finally:
        ttm_module.cg = original_cg
    elapsed = time.perf_counter() - start
    return (
        {
            "rtol": rtol,
            "elapsed_s": elapsed,
            "cg_iterations_per_newton_solve": iterations,
            "cg_iterations_total": int(sum(iterations)),
            "net_diffusion_energy_change_J": net_energy_change_J,
            "active_ag_cells": int(np.count_nonzero(model.lattice.ag_fraction > 0.0)),
        },
        model.electron.temperature_K.copy(),
    )


def main() -> None:
    strict, strict_te = run_case(1.0e-10)
    optimized, optimized_te = run_case(ttm_module.ELECTRON_DIFFUSION_CG_RTOL)
    report = {
        "status": "QSGS_ELECTRON_DIFFUSION_CG_TOLERANCE_COMPARISON",
        "grid_file": GRID_FILE.name,
        "face_contact_file": CONTACT_FILE.name,
        "geometry_file": GEOMETRY_FILE.name,
        "boundary_case": "A",
        "xy_tile": TILE_XY,
        "dt_ps": DT_PS,
        "grid_shape_xyz": list(strict_te.shape),
        "active_ag_cells": strict["active_ag_cells"],
        "temperature_initialization": "300 K plus 800 K radial/depth exponential profile; diffusion only",
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "strict_reference": strict,
        "optimized": optimized,
        "comparison": {
            "max_abs_electron_temperature_difference_K": float(
                np.max(np.abs(strict_te - optimized_te))
            ),
            "cg_iteration_reduction_fraction": 1.0
            - optimized["cg_iterations_total"] / strict["cg_iterations_total"],
            "net_diffusion_energy_difference_J": optimized[
                "net_diffusion_energy_change_J"
            ]
            - strict["net_diffusion_energy_change_J"],
        },
        "note": "Timing is a single local run and is environment-sensitive; compare iteration counts and fields as the reproducible numerical measures.",
    }
    output = HERE / "nested_grid_diagnostics" / "electron_diffusion_cg_tolerance_comparison.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
