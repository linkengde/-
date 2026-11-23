#!/usr/bin/env python3
"""Build OVITO-ready display meshes from the model's QSGS skeleton mask.

This script creates visualization geometry only. It does not edit or relax
the atomistic model.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import ndimage
from vtkmodules.vtkCommonDataModel import vtkImageData
from vtkmodules.vtkFiltersCore import (
    vtkPolyDataConnectivityFilter,
    vtkFlyingEdges3D,
    vtkTriangleFilter,
    vtkWindowedSincPolyDataFilter,
)
from vtkmodules.vtkImagingGeneral import vtkImageGaussianSmooth
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy

HERE = Path(__file__).resolve().parent
MODEL_DIR = HERE.parents[1]
MASK_FILE = MODEL_DIR / "qsgs_morphology.npz"
ATOM_FILE = MODEL_DIR / "ag_ti3sic2_qsgs.extxyz"
WHOLE_MESH = HERE / "Ti3SiC2_whole_surface_mesh.npz"
LOCAL_MESH = HERE / "Ti3SiC2_local_branch_surface_mesh.npz"
LOCAL_ATOMS = HERE / "Ti3SiC2_local_branch_atoms.extxyz"
REPORT = HERE / "mesh_build_report.json"

VOXEL_A = 2.5
LOCAL_START = (27, 56, 41)  # (x, y, z) voxel indices
LOCAL_SIZE = (10, 10, 10)   # 25 x 25 x 25 A


def build_mesh(
    volume: np.ndarray,
    first_voxel_index: tuple[int, int, int],
    path: Path,
    *,
    iso_value: float = 0.55,
    keep_largest_surface: bool = False,
) -> dict:
    """Create a smoothed isosurface, preserving the original voxel geometry."""
    padded = np.pad(volume.astype(np.float32), 1, mode="constant")
    # Image samples sit at voxel centers. The padded origin is one sample before
    # the first real voxel center, so the isosurface stays in model coordinates.
    origin = tuple((np.asarray(first_voxel_index, dtype=float) * VOXEL_A - 1.25).tolist())

    image = vtkImageData()
    image.SetDimensions(*padded.shape)
    image.SetOrigin(*origin)
    image.SetSpacing(VOXEL_A, VOXEL_A, VOXEL_A)
    scalars = numpy_to_vtk(padded.ravel(order="F"), deep=True)
    image.GetPointData().SetScalars(scalars)

    gaussian = vtkImageGaussianSmooth()
    gaussian.SetInputData(image)
    gaussian.SetDimensionality(3)
    gaussian.SetStandardDeviation(0.55, 0.55, 0.55)
    gaussian.SetRadiusFactors(2.0, 2.0, 2.0)

    contour = vtkFlyingEdges3D()
    contour.SetInputConnection(gaussian.GetOutputPort())
    # A raised display threshold contracts the phase surface and makes the
    # branching visible at this scale. This is a render-only change.
    contour.SetValue(0, iso_value)
    contour.ComputeNormalsOff()
    contour.ComputeGradientsOff()

    input_port = contour.GetOutputPort()
    if keep_largest_surface:
        connected = vtkPolyDataConnectivityFilter()
        connected.SetInputConnection(input_port)
        connected.SetExtractionModeToLargestRegion()
        input_port = connected.GetOutputPort()

    smooth = vtkWindowedSincPolyDataFilter()
    smooth.SetInputConnection(input_port)
    smooth.SetNumberOfIterations(8)
    smooth.SetPassBand(0.12)
    smooth.BoundarySmoothingOn()
    smooth.FeatureEdgeSmoothingOff()
    smooth.NonManifoldSmoothingOn()
    smooth.NormalizeCoordinatesOn()

    triangles = vtkTriangleFilter()
    triangle_input = smooth.GetOutputPort()
    if keep_largest_surface:
        final_connected = vtkPolyDataConnectivityFilter()
        final_connected.SetInputConnection(triangle_input)
        final_connected.SetExtractionModeToLargestRegion()
        triangle_input = final_connected.GetOutputPort()
    triangles.SetInputConnection(triangle_input)
    triangles.Update()

    poly = triangles.GetOutput()
    vertices = vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float32, copy=True)
    faces = vtk_to_numpy(poly.GetPolys().GetData()).reshape(-1, 4)[:, 1:].astype(np.int32, copy=True)
    np.savez_compressed(path, vertices=vertices, faces=faces)
    return {
        "file": path.name,
        "vertex_count": int(len(vertices)),
        "triangle_count": int(len(faces)),
        "bounds_A": [vertices.min(axis=0).tolist(), vertices.max(axis=0).tolist()],
    }


def export_local_atoms(component: np.ndarray, path: Path) -> dict[str, int]:
    """Export actual source atoms that lie within the selected mask component."""
    x0, y0, z0 = (np.asarray(LOCAL_START, dtype=float) * VOXEL_A).tolist()
    nx, ny, nz = LOCAL_SIZE
    x1, y1, z1 = x0 + nx * VOXEL_A, y0 + ny * VOXEL_A, z0 + nz * VOXEL_A
    counts: Counter[str] = Counter()
    lines: list[tuple[str, float, float, float]] = []

    with ATOM_FILE.open() as source:
        int(next(source))  # atom count
        next(source)       # extxyz comment
        for row in source:
            fields = row.split()
            if fields[4] != "1":  # phase 1 is Ti3SiC2; phase 2 is Ag
                continue
            xyz = np.asarray([float(v) for v in fields[1:4]], dtype=float)
            if not (x0 <= xyz[0] < x1 and y0 <= xyz[1] < y1 and z0 <= xyz[2] < z1):
                continue
            voxel = np.floor(xyz / VOXEL_A).astype(int) - np.asarray(LOCAL_START)
            if not component[tuple(voxel)]:
                continue
            symbol = fields[0]
            lines.append((symbol, *(xyz - np.asarray([x0, y0, z0])).tolist()))
            counts[symbol] += 1

    extent = np.asarray(LOCAL_SIZE, dtype=float) * VOXEL_A
    comment = (
        f'Lattice="{extent[0]:g} 0 0 0 {extent[1]:g} 0 0 0 {extent[2]:g}" '
        'Properties=species:S:1:pos:R:3 pbc="F F F" '
        'source="unchanged atom positions from Ag-Ti3SiC2 QSGS model" '
        f'crop_origin_A="{x0:g} {y0:g} {z0:g}" extent_A="{extent[0]:g} {extent[1]:g} {extent[2]:g}" '
        'mask_component="largest face-connected Ti3SiC2 component in ROI"'
    )
    with path.open("w") as target:
        target.write(f"{len(lines)}\n{comment}\n")
        for symbol, x, y, z in lines:
            target.write(f"{symbol} {x:.8f} {y:.8f} {z:.8f}\n")
    return dict(counts)


def main() -> None:
    archive = np.load(MASK_FILE)
    mask = np.asarray(archive["skeleton"], dtype=bool)
    size_a = np.asarray(archive["size_A"], dtype=float).tolist()
    whole = build_mesh(mask, (0, 0, 0), WHOLE_MESH)

    x, y, z = LOCAL_START
    sx, sy, sz = LOCAL_SIZE
    crop = mask[x : x + sx, y : y + sy, z : z + sz]
    labels, component_count = ndimage.label(
        crop, structure=ndimage.generate_binary_structure(3, 1)
    )
    component_sizes = np.bincount(labels.ravel())
    largest_id = int(np.argmax(component_sizes[1:]) + 1)
    largest = labels == largest_id
    local = build_mesh(largest, LOCAL_START, LOCAL_MESH, iso_value=0.60, keep_largest_surface=True)
    atom_counts = export_local_atoms(largest, LOCAL_ATOMS)

    report = {
        "status": "OVITO_WHOLE_SKELETON_AND_LOCAL_BRANCH_MESHES_BUILT",
        "source_model": str(MODEL_DIR / "ag_ti3sic2_qsgs.extxyz"),
        "source_mask": str(MASK_FILE),
        "source_coordinates_modified": False,
        "material_extent_A": size_a,
        "voxel_spacing_A": VOXEL_A,
        "whole_mesh": whole,
        "local_roi_A": {
            "x": [LOCAL_START[0] * VOXEL_A, (LOCAL_START[0] + LOCAL_SIZE[0]) * VOXEL_A],
            "y": [LOCAL_START[1] * VOXEL_A, (LOCAL_START[1] + LOCAL_SIZE[1]) * VOXEL_A],
            "z": [LOCAL_START[2] * VOXEL_A, (LOCAL_START[2] + LOCAL_SIZE[2]) * VOXEL_A],
        },
        "local_roi_voxel_start_xyz": list(LOCAL_START),
        "local_roi_voxel_extent_xyz": list(LOCAL_SIZE),
        "local_roi_occupied_voxels": int(crop.sum()),
        "local_roi_face_connected_components": int(component_count),
        "local_roi_largest_component_voxels": int(largest.sum()),
        "local_atom_counts_in_largest_component": atom_counts,
        "local_mesh": local,
        "display_processing": {
            "method": "Gaussian smoothed QSGS phase-mask isosurface, followed by windowed-sinc mesh smoothing",
            "gaussian_sigma_voxels": 0.55,
            "whole_isosurface_level": 0.55,
            "local_isosurface_level": 0.60,
            "windowed_sinc_iterations": 8,
            "windowed_sinc_passband": 0.12,
            "local_crop_filter": "largest face-connected mask component, then largest resulting surface component; display-only removal of small disconnected fragments",
            "note": "Only the display mesh is smoothed and cropped. The atomistic model and its coordinates remain unchanged.",
        },
        "outputs": [WHOLE_MESH.name, LOCAL_MESH.name, LOCAL_ATOMS.name],
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
