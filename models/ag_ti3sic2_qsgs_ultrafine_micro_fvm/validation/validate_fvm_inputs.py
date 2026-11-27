#!/usr/bin/env python3
"""检查 NPZ/VTI 输入字段，并做零热源均匀温度初始化一致性步。"""

from __future__ import annotations

import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parents[1]
PHASE_FIELD = REPOSITORY / "models" / "ag_ti3sic2_qsgs_1to1_finer" / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_FIELD))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid  # noqa: E402


MODEL_DIRS = (
    ROOT / "ultrafine_500nm" / "fvm_inputs",
    ROOT / "micro_10um" / "fvm_inputs",
)
OUTPUT = ROOT / "validation" / "fvm_input_validation.json"


def vti_field(path: Path, name: str) -> np.ndarray:
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    image = reader.GetOutput()
    array = image.GetCellData().GetArray(name)
    if array is None:
        raise ValueError(f"{path}: missing VTI field {name}")
    return vtk_to_numpy(array)


def main() -> None:
    results = []
    for model_dir in MODEL_DIRS:
        for npz_path in sorted(model_dir.glob("*_fvm.npz")):
            started = time.perf_counter()
            with np.load(npz_path, allow_pickle=False) as archive:
                shape = tuple(int(value) for value in archive["shape_xyz"])
                ag_fraction = np.asarray(archive["ag_fraction"], dtype=np.float64)
                skeleton_fraction = np.asarray(archive["skeleton_fraction"], dtype=np.float64)
                ag_mask = np.asarray(archive["ag_mask"], dtype=bool)
                skeleton_mask = np.asarray(archive["skeleton_mask"], dtype=bool)
                size_A = np.asarray(archive["size_A"], dtype=np.float64)
                spacing_A = np.asarray(archive["spacing_A"], dtype=np.float64)
                size_m = np.asarray(archive["size_m"], dtype=np.float64)
                spacing_m = np.asarray(archive["spacing_m"], dtype=np.float64)
                seed = int(archive["seed"])
            if shape != (150, 150, 150) or ag_fraction.shape != shape:
                raise ValueError(f"{npz_path}: unexpected shape")
            if not np.array_equal(ag_mask, ~skeleton_mask):
                raise ValueError(f"{npz_path}: binary phases do not complement each other")
            if not np.allclose(ag_fraction + skeleton_fraction, 1.0, atol=0.0, rtol=0.0):
                raise ValueError(f"{npz_path}: phase fractions do not sum to one")
            if not np.array_equal(ag_fraction, ag_mask.astype(float)):
                raise ValueError(f"{npz_path}: Ag fractions differ from binary mask")
            if not np.array_equal(skeleton_fraction, skeleton_mask.astype(float)):
                raise ValueError(f"{npz_path}: skeleton fractions differ from binary mask")
            if not np.allclose(size_A * 1.0e-10, size_m, rtol=1.0e-14, atol=0.0):
                raise ValueError(f"{npz_path}: size_A to SI conversion failed")
            if not np.allclose(spacing_A * 1.0e-10, spacing_m, rtol=1.0e-14, atol=0.0):
                raise ValueError(f"{npz_path}: spacing_A to SI conversion failed")
            if not np.allclose(size_m / np.asarray(shape), spacing_m, rtol=1.0e-14, atol=0.0):
                raise ValueError(f"{npz_path}: physical size and voxel pitch disagree")

            vti_path = npz_path.with_name(npz_path.name.replace("_fvm.npz", ".vti"))
            vtk_ag = vti_field(vti_path, "ag_fraction")
            vtk_sc = vti_field(vti_path, "skeleton_fraction")
            if vtk_ag.size != int(np.prod(shape)) or vtk_sc.size != int(np.prod(shape)):
                raise ValueError(f"{vti_path}: VTI cell count mismatch")
            if not np.array_equal(vtk_ag.astype(bool), ag_fraction.T.ravel(order="C")):
                raise ValueError(f"{vti_path}: cell order or Ag values mismatch")
            if not np.array_equal(vtk_sc.astype(bool), skeleton_fraction.T.ravel(order="C")):
                raise ValueError(f"{vti_path}: cell order or skeleton values mismatch")

            for boundary in ("A", "B"):
                solver = QSGSEnthalpyGrid(npz_path, boundary, peak_flux_eV_ps_A2=0.0)
                if solver.temperature_shape != shape:
                    raise ValueError(f"{npz_path}: solver read a different shape")
                if not np.allclose(solver.size_A, size_A, rtol=0.0, atol=1.0e-12):
                    raise ValueError(f"{npz_path}: solver changed the physical dimensions")
                if not np.allclose(solver.spacing_A, spacing_A, rtol=0.0, atol=1.0e-12):
                    raise ValueError(f"{npz_path}: solver changed the grid spacing")
                if solver.actual_source_power_W != 0.0 or np.any(solver.source_W_m3 != 0.0):
                    raise ValueError(f"{npz_path}: zero-source initialization created heat input")
                dt_s = 1.0e-18
                row = solver.advance(dt_s, source_on=False, enforce_stability=True)
                if np.any(np.abs(solver.state.temperature_K - 300.0) > 1.0e-10):
                    raise ValueError(f"{npz_path}: uniform 300 K state changed")
                if np.any(solver.state.ag_liquid_fraction != 0.0):
                    raise ValueError(f"{npz_path}: zero-source step created liquid Ag")
                if abs(row["energy_balance_error_J"]) > 1.0e-28:
                    raise ValueError(f"{npz_path}: zero-source step did not conserve energy")
                results.append(
                    {
                        "file": str(npz_path.relative_to(ROOT)),
                        "seed": seed,
                        "shape_xyz": shape,
                        "boundary_case": boundary,
                        "ag_fraction": float(ag_fraction.mean()),
                        "skeleton_fraction": float(skeleton_fraction.mean()),
                        "size_m": size_m.tolist(),
                        "spacing_m": spacing_m.tolist(),
                        "size_A": size_A.tolist(),
                        "spacing_A": spacing_A.tolist(),
                        "fvm_zero_source_power_W": solver.actual_source_power_W,
                        "uniform_state_final_temperature_K": float(solver.state.temperature_K.max()),
                        "liquid_fraction_max": float(solver.state.ag_liquid_fraction.max()),
                        "energy_residual_J": row["energy_balance_error_J"],
                        "stable_timestep_s": solver.stable_timestep_s(),
                    }
                )
            results[-1]["archive_check_seconds_including_two_boundary_cases"] = time.perf_counter() - started
            results[-1]["process_peak_rss_mib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
            print(json.dumps(results[-1]), flush=True)

    payload = {
        "status": "ZERO_SOURCE_FVM_INPUT_AND_INITIALIZATION_CHECKS_PASSED",
        "solver_source": str(PHASE_FIELD.relative_to(REPOSITORY)),
        "cases": results,
        "heat_sources_applied": False,
        "flow_or_phase_change_run": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
