#!/usr/bin/env python3
"""导出真实150³周期导热场为VTI并绘制中心截面及高热流位置。"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
CASE_ID = "performance_n150_seed1105_ultrafine_500nm_R0e+00"
CASE_DIR = HERE / "results" / "cases"
FIGURES = HERE / "figures"
N = 150
SIZE_M = 10.0e-6
SPACING_M = SIZE_M / N


def vtk_array(array: np.ndarray, name: str, array_type: int) -> vtk.vtkDataArray:
    order = np.transpose(array, (2, 1, 0)).ravel(order="C")
    result = numpy_to_vtk(order, deep=True, array_type=array_type)
    result.SetName(name)
    return result


def export_vti(axis: int, data: dict[str, np.ndarray]) -> Path:
    image = vtk.vtkImageData()
    image.SetDimensions(N + 1, N + 1, N + 1)
    image.SetOrigin(0.0, 0.0, 0.0)
    image.SetSpacing(SPACING_M, SPACING_M, SPACING_M)
    image.GetCellData().AddArray(vtk_array(data["material_id"].astype(np.uint8), "material_id", vtk.VTK_UNSIGNED_CHAR))
    image.GetCellData().AddArray(vtk_array(data["temperature_K"].astype(np.float32), "temperature_K", vtk.VTK_FLOAT))
    for component in range(3):
        image.GetCellData().AddArray(vtk_array(data[f"heat_flux_{component}_W_m2"].astype(np.float32), f"heat_flux_{component}_W_m2", vtk.VTK_FLOAT))
    magnitude = np.sqrt(sum(data[f"heat_flux_{component}_W_m2"].astype(np.float64) ** 2 for component in range(3)))
    image.GetCellData().AddArray(vtk_array(magnitude.astype(np.float32), "heat_flux_magnitude_W_m2", vtk.VTK_FLOAT))
    path = HERE / "results" / f"seed1105_150_{'xyz'[axis]}_gradient_fields.vti"
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(image)
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    if writer.Write() != 1:
        raise RuntimeError(f"VTK writer failed for {path}")
    return path


def plot_axis(axis: int, data: dict[str, np.ndarray]) -> None:
    mask = data["material_id"].astype(bool)
    temp = data["temperature_K"]
    q = np.sqrt(sum(data[f"heat_flux_{c}_W_m2"].astype(np.float64) ** 2 for c in range(3)))
    cut = N // 2
    temperature_slice = np.take(temp, cut, axis=axis).T
    flux_slice = np.take(q, cut, axis=axis).T
    phase_slice = np.take(mask, cut, axis=axis).T
    axis_names = "xyz"
    plane_names = ("yz", "xz", "xy")
    extent = [0.0, SIZE_M * 1e6, 0.0, SIZE_M * 1e6]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    image = axes[0].imshow(temperature_slice, origin="lower", extent=extent, cmap="inferno", interpolation="nearest")
    axes[0].contour(phase_slice.astype(float), levels=[0.5], colors=["#36b8ff"], linewidths=0.25,
                    origin="lower", extent=extent)
    fig.colorbar(image, ax=axes[0], label="Temperature (K)")
    axes[0].set_title(f"{axis_names[axis].upper()} gradient · {plane_names[axis].upper()} center plane")
    flux_image = axes[1].imshow(flux_slice, origin="lower", extent=extent, cmap="magma", interpolation="nearest")
    axes[1].contour(phase_slice.astype(float), levels=[0.5], colors=["#43baff"], linewidths=0.25,
                    origin="lower", extent=extent)
    fig.colorbar(flux_image, ax=axes[1], label="|q| (W/m²)")
    axes[1].set_title("Heat-flux concentration; blue line is phase interface")
    for ax in axes:
        ax.set_xlabel("μm")
        ax.set_ylabel("μm")
        ax.set_aspect("equal")
    fig.suptitle("Computed RVE field · seed 1105 · 150³ · 300–310 K affine periodic gradient")
    path = FIGURES / f"temperature_heatflow_axis_{axis_names[axis]}.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_flux_points(data: dict[str, np.ndarray]) -> Path:
    q = np.stack([data[f"heat_flux_{c}_W_m2"] for c in range(3)], axis=-1).astype(np.float64)
    magnitude = np.linalg.norm(q, axis=-1)
    count = min(100, magnitude.size)
    flat = np.argpartition(magnitude.ravel(), -count)[-count:]
    flat = flat[np.argsort(magnitude.ravel()[flat])[::-1]]
    coords = np.column_stack(np.unravel_index(flat, magnitude.shape))
    path = HERE / "results" / "LOCAL_HEAT_FLOW_POINTS.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["rank", "x_um", "y_um", "z_um", "temperature_K", "material", "qx_W_m2", "qy_W_m2", "qz_W_m2", "q_magnitude_W_m2"])
        for rank, (i, j, k) in enumerate(coords, start=1):
            phase = "Ti3SiC2" if data["material_id"][i, j, k] else "Ag"
            writer.writerow([
                rank,
                (i + 0.5) * SPACING_M * 1e6,
                (j + 0.5) * SPACING_M * 1e6,
                (k + 0.5) * SPACING_M * 1e6,
                float(data["temperature_K"][i, j, k]),
                phase,
                *[float(q[i, j, k, c]) for c in range(3)],
                float(magnitude[i, j, k]),
            ])
    return path


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    paths = []
    for axis in range(3):
        snapshot = CASE_DIR / f"{CASE_ID}_axis{axis}_fields.npz"
        with np.load(snapshot, allow_pickle=False) as archive:
            data = {key: np.asarray(archive[key]) for key in archive.files}
        paths.append(str(export_vti(axis, data).relative_to(HERE)))
        plot_axis(axis, data)
        if axis == 0:
            flux_points = str(write_flux_points(data).relative_to(HERE))
    summary = {
        "status": "COMPUTED_FIELDS_EXPORTED",
        "source_case": CASE_ID,
        "grid": [N, N, N],
        "domain_size_m": [SIZE_M] * 3,
        "delta_temperature_K": 10.0,
        "vti_files": paths,
        "plots": [f"figures/temperature_heatflow_axis_{axis}.png" for axis in "xyz"],
        "high_flux_points_csv": flux_points,
        "vtk_version": vtk.vtkVersion.GetVTKVersion(),
        "note": "VTI fields are computed cell-centred temperature and heat-flux values; flux-point CSV ranks the 100 largest local |q| cells for the X-gradient case.",
    }
    (HERE / "results" / "visualization_export.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
