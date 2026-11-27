#!/usr/bin/env python3
"""把二值 QSGS 几何写成现有单温焓法 NPZ 和 ParaView VTI。"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "geometry_generation"
MODELS = {
    "ultrafine_500nm": {"size_m": 10.0e-6, "directory": ROOT / "ultrafine_500nm" / "fvm_inputs"},
    "micro_10um": {"size_m": 200.0e-6, "directory": ROOT / "micro_10um" / "fvm_inputs"},
}
SHAPE = (150, 150, 150)
SEEDS = (1105, 1136, 1182)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def make_vti(path: Path, skeleton: np.ndarray, size_m: float, seed: int) -> None:
    spacing_m = size_m / SHAPE[0]
    image = vtk.vtkImageData()
    image.SetDimensions(*(n + 1 for n in SHAPE))
    image.SetOrigin(0.0, 0.0, 0.0)
    image.SetSpacing(spacing_m, spacing_m, spacing_m)
    vtk_order = np.transpose(skeleton, (2, 1, 0)).ravel(order="C")
    material_id = numpy_to_vtk(vtk_order.astype(np.uint8), deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
    material_id.SetName("material_id")
    image.GetCellData().AddArray(material_id)
    ag_fraction = numpy_to_vtk((1 - vtk_order).astype(np.float32), deep=True, array_type=vtk.VTK_FLOAT)
    ag_fraction.SetName("ag_fraction")
    image.GetCellData().AddArray(ag_fraction)
    skeleton_fraction = numpy_to_vtk(vtk_order.astype(np.float32), deep=True, array_type=vtk.VTK_FLOAT)
    skeleton_fraction.SetName("skeleton_fraction")
    image.GetCellData().AddArray(skeleton_fraction)

    field = vtk.vtkIntArray()
    field.SetName("seed")
    field.InsertNextValue(seed)
    image.GetFieldData().AddArray(field)
    shape_field = vtk.vtkIntArray()
    shape_field.SetName("shape_xyz")
    for value in SHAPE:
        shape_field.InsertNextValue(value)
    image.GetFieldData().AddArray(shape_field)
    size_field = vtk.vtkDoubleArray()
    size_field.SetName("size_m")
    spacing_field = vtk.vtkDoubleArray()
    spacing_field.SetName("spacing_m")
    for _ in range(3):
        size_field.InsertNextValue(size_m)
        spacing_field.InsertNextValue(spacing_m)
    image.GetFieldData().AddArray(size_field)
    image.GetFieldData().AddArray(spacing_field)

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(image)
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    if writer.Write() != 1:
        raise RuntimeError(f"VTK XML writer failed for {path}")


def main() -> None:
    source_files: dict[int, tuple[Path, np.ndarray, np.ndarray]] = {}
    for seed in SEEDS:
        path = GEOMETRY / f"qsgs_seed_{seed}.npz"
        with np.load(path, allow_pickle=False) as archive:
            skeleton = np.asarray(archive["skeleton"], dtype=bool)
            ag = np.asarray(archive["ag"], dtype=bool)
            if skeleton.shape != SHAPE or not np.array_equal(ag, ~skeleton):
                raise ValueError(f"seed {seed}: unexpected geometry shape or phase complement")
            source_files[seed] = path, skeleton, ag

    for model_name, model in MODELS.items():
        size_m = float(model["size_m"])
        spacing_m = size_m / SHAPE[0]
        size_A = np.full(3, size_m * 1.0e10, dtype=np.float64)
        spacing_A = np.full(3, spacing_m * 1.0e10, dtype=np.float64)
        output_dir = Path(model["directory"])
        output_dir.mkdir(parents=True, exist_ok=True)
        for seed, (source_path, skeleton, ag) in source_files.items():
            started = time.perf_counter()
            npz_path = output_dir / f"qsgs_seed_{seed}_fvm.npz"
            np.savez_compressed(
                npz_path,
                ag_fraction=ag.astype(np.float64),
                skeleton_fraction=skeleton.astype(np.float64),
                ag_mask=ag,
                skeleton_mask=skeleton,
                shape_xyz=np.asarray(SHAPE, dtype=np.int32),
                size_A=size_A,
                spacing_A=spacing_A,
                size_m=np.full(3, size_m, dtype=np.float64),
                spacing_m=np.full(3, spacing_m, dtype=np.float64),
                seed=np.asarray(seed, dtype=np.int32),
                qsgs_p0=np.asarray(0.05, dtype=np.float64),
                qsgs_growth=np.full(3, 0.11, dtype=np.float64),
                correlation_sigma_voxels=np.asarray(4.0, dtype=np.float64),
                source_geometry_sha256=np.asarray(sha256(source_path)),
            )
            native_path = output_dir / f"qsgs_seed_{seed}_native.npz"
            np.savez_compressed(
                native_path,
                skeleton=skeleton,
                size_A=size_A,
                shape_xyz=np.asarray(SHAPE, dtype=np.int32),
                size_m=np.full(3, size_m, dtype=np.float64),
                spacing_m=np.full(3, spacing_m, dtype=np.float64),
                seed=np.asarray(seed, dtype=np.int32),
                qsgs_p0=np.asarray(0.05, dtype=np.float64),
                qsgs_growth=np.full(3, 0.11, dtype=np.float64),
                correlation_sigma_voxels=np.asarray(4.0, dtype=np.float64),
                source_geometry_sha256=np.asarray(sha256(source_path)),
            )
            vti_path = output_dir / f"qsgs_seed_{seed}.vti"
            make_vti(vti_path, skeleton, size_m, seed)
            meta = {
                "model": model_name,
                "seed": seed,
                "shape_xyz": SHAPE,
                "size_m": [size_m] * 3,
                "spacing_m": [spacing_m] * 3,
                "size_A": size_A.tolist(),
                "spacing_A": spacing_A.tolist(),
                "ag_fraction": float(ag.mean()),
                "skeleton_fraction": float(skeleton.mean()),
                "boundary_cases_supported_by_current_solver": {
                    "A": "adiabatic Ti3SiC2 scaffold limit; scaffold heat capacity and conductivity omitted",
                    "B": "ideal thermal contact; shared temperature and volume-weighted properties, zero added interface resistance",
                },
                "no_heat_source_applied": True,
                "source_geometry_file": str(source_path.relative_to(ROOT)),
                "source_geometry_sha256": sha256(source_path),
                "fvm_npz": npz_path.name,
                "fvm_npz_sha256": sha256(npz_path),
                "native_solver_npz": native_path.name,
                "native_solver_npz_sha256": sha256(native_path),
                "vti": vti_path.name,
                "vti_sha256": sha256(vti_path),
                "vtk_version": vtk.vtkVersion.GetVTKVersion(),
                "export_wall_seconds": time.perf_counter() - started,
                "fvm_npz_bytes": npz_path.stat().st_size,
                "native_solver_npz_bytes": native_path.stat().st_size,
                "vti_bytes": vti_path.stat().st_size,
            }
            (output_dir / f"qsgs_seed_{seed}_metadata.json").write_text(
                json.dumps(meta, indent=2) + "\n", encoding="utf-8"
            )
            print(json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()
