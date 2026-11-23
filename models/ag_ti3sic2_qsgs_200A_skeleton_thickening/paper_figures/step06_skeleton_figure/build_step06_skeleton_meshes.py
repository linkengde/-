#!/usr/bin/env python3
"""从 step06 的 Ti3SiC2 相场提取整体骨架和局部枝状网格用于论文配图。"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read
from scipy import ndimage as ndi
from skimage.morphology import skeletonize

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1]
REPO = HERE.parents[3]
STEP06 = MODEL / "output/step06_coarsened_networks"
FIELD_FILE = STEP06 / "coarsened_step06_phase_fields.npz"
ATOM_FILE = STEP06 / "ag_ti3sic2_coarsened_networks_step06.extxyz"
OUT = HERE

# 与参考图一致的 25 × 25 × 25 Å 局部取景，避开盒壁。
# 该局部窗口含一个 6 邻接连通的分叉骨架，保留该连通体的全部枝条，尺寸为 25 Å。
# 该 ROI 的局部中心线有 4 个末端和 1 个分叉簇，能直接显示为连通树枝。
LOCAL_START = (20, 52, 6)
LOCAL_SIZE = (10, 10, 10)
VOXEL_A = 2.5
REPORT = OUT / "mesh_build_report.json"


def load_surface_builder():
    """重用已经核验过的 VTK marching-cubes 和表面平滑程序。"""
    path = HERE / "step06_surface_mesh_utils.py"
    spec = importlib.util.spec_from_file_location("step06_mesh_builder", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"缺少 VTK 网格工具：{path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_local_atom_crop(atoms, selected: np.ndarray, path: Path) -> dict:
    """导出当前 step06 原子模型中对应局部 mask 的实际 Ti3SiC2 原子。"""
    phase = np.asarray(atoms.arrays["phase"], dtype=np.int32)
    symbols = np.asarray(atoms.get_chemical_symbols())
    pos = np.asarray(atoms.positions)
    origin_vox = np.asarray(LOCAL_START)
    size_vox = np.asarray(LOCAL_SIZE)
    inside_roi = np.all(pos >= origin_vox * VOXEL_A, axis=1) & np.all(
        pos < (origin_vox + size_vox) * VOXEL_A, axis=1
    )
    ids = np.flatnonzero((phase == 1) & inside_roi)
    vox = np.floor(pos[ids] / VOXEL_A).astype(int) - origin_vox
    inside_mask = selected[vox[:, 0], vox[:, 1], vox[:, 2]]
    ids = ids[inside_mask]
    local_pos = pos[ids] - origin_vox * VOXEL_A
    extent = size_vox * VOXEL_A
    comment = (
        f'Lattice="{extent[0]:g} 0 0 0 {extent[1]:g} 0 0 0 {extent[2]:g}" '
        'Properties=species:S:1:pos:R:3 pbc="F F F" '
        'source="step06 coarsened atom model; coordinates unchanged" '
        f'crop_origin_A="{(origin_vox * VOXEL_A).tolist()}" '
        'phase="Ti3SiC2"'
    )
    with path.open("w") as f:
        f.write(f"{len(ids)}\n{comment}\n")
        for symbol, xyz in zip(symbols[ids], local_pos):
            f.write(f"{symbol} {xyz[0]:.8f} {xyz[1]:.8f} {xyz[2]:.8f}\n")
    return {s: int(np.count_nonzero(symbols[ids] == s)) for s in ("Ti", "Si", "C")}


def count_surface_components(faces: np.ndarray, vertex_count: int) -> int:
    """统计三角网格的面连通分量数，确认局部表面没有断开。"""
    parent = np.arange(vertex_count, dtype=np.int64)

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = int(parent[index])
        return index

    for face in faces:
        root0 = find(int(face[0]))
        for vertex in face[1:]:
            root = find(int(vertex))
            if root != root0:
                parent[root] = root0
    return int(len({find(i) for i in range(vertex_count)}))


def main() -> None:
    if not FIELD_FILE.is_file() or not ATOM_FILE.is_file():
        raise FileNotFoundError("step06 原子文件或平滑相场缺失")
    mesh_builder = load_surface_builder()
    field_data = np.load(FIELD_FILE)
    mask = np.asarray(field_data["step06_equal_volume_Ti3SiC2_mask"], dtype=bool)
    if mask.shape != (80, 80, 80):
        raise ValueError(f"step06 Ti3SiC2 相场尺寸异常：{mask.shape}")
    sx, sy, sz = LOCAL_START
    nx, ny, nz = LOCAL_SIZE
    crop = mask[sx:sx + nx, sy:sy + ny, sz:sz + nz]
    labels, n_components = ndi.label(crop, structure=ndi.generate_binary_structure(3, 1))
    counts = np.bincount(labels.ravel())
    if n_components != 1:
        raise RuntimeError(f"局部 ROI 必须是一个 6 邻接连通骨架，实际连通分量数：{n_components}")
    selected = labels == 1

    whole_mesh_path = OUT / "Ti3SiC2_step06_whole_surface_mesh.npz"
    local_mesh_path = OUT / "Ti3SiC2_step06_local_branch_surface_mesh.npz"
    whole_info = mesh_builder.build_mesh(mask, (0, 0, 0), whole_mesh_path, iso_value=0.55)
    local_info = mesh_builder.build_mesh(
        selected, LOCAL_START, local_mesh_path, iso_value=0.60, keep_largest_surface=True
    )
    local_mesh = np.load(local_mesh_path)
    surface_components = count_surface_components(local_mesh["faces"], len(local_mesh["vertices"]))
    if surface_components != 1:
        raise RuntimeError(f"局部三角网格没有保持连续，表面连通分量数：{surface_components}")
    centerline = skeletonize(selected)
    neighbor_count = ndi.convolve(centerline.astype(np.uint8), np.ones((3, 3, 3), dtype=np.uint8), mode="constant") - centerline
    branch_voxels = centerline & (neighbor_count > 3)
    branch_clusters, branch_cluster_count = ndi.label(
        branch_voxels, structure=np.ones((3, 3, 3), dtype=bool)
    )
    endpoint_count = int(np.count_nonzero(centerline & (neighbor_count == 1)))
    atoms = read(ATOM_FILE)
    local_atom_counts = write_local_atom_crop(
        atoms, selected, OUT / "Ti3SiC2_step06_local_branch_atoms.extxyz"
    )

    report = {
        "status": "STEP06_WHOLE_AND_LOCAL_Ti3SiC2_SKELETON_MESHES_BUILT",
        "source_atom_model": str(ATOM_FILE.relative_to(REPO)),
        "source_phase_field": str(FIELD_FILE.relative_to(REPO)),
        "source_mask_key": "step06_equal_volume_Ti3SiC2_mask",
        "material_region_A": [200.0, 200.0, 200.0],
        "mesh_voxel_spacing_A": VOXEL_A,
        "local_roi_voxel_start_xyz": list(LOCAL_START),
        "local_roi_voxel_extent_xyz": list(LOCAL_SIZE),
        "local_roi_A": {
            "x": [LOCAL_START[0] * VOXEL_A, (LOCAL_START[0] + LOCAL_SIZE[0]) * VOXEL_A],
            "y": [LOCAL_START[1] * VOXEL_A, (LOCAL_START[1] + LOCAL_SIZE[1]) * VOXEL_A],
            "z": [LOCAL_START[2] * VOXEL_A, (LOCAL_START[2] + LOCAL_SIZE[2]) * VOXEL_A],
        },
        "local_roi_occupied_voxels": int(crop.sum()),
        "local_roi_6_neighbor_connected_components": int(n_components),
        "retained_connected_branch_voxels": int(selected.sum()),
        "local_surface_connected_components": surface_components,
        "skeleton_topology": {
            "centerline_voxels": int(centerline.sum()),
            "terminal_branches": endpoint_count,
            "branch_junction_clusters": int(branch_cluster_count),
            "branch_junction_voxels": int(branch_voxels.sum()),
        },
        "local_atom_counts_by_element": local_atom_counts,
        "whole_mesh": whole_info,
        "local_mesh": local_info,
        "display_processing": {
            "method": "step06 equal-volume smoothed Ti3SiC2 occupancy mask rendered as a Gaussian-smoothed isosurface",
            "isosurface_level_whole": 0.55,
            "isosurface_level_local": 0.60,
            "mesh_smoothing": "windowed-sinc display smoothing only",
            "coordinate_system": "Å",
            "atom_model_modified": False,
        },
        "outputs": [whole_mesh_path.name, local_mesh_path.name,
                    "Ti3SiC2_step06_local_branch_atoms.extxyz"],
    }
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
