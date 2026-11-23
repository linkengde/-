#!/usr/bin/env python3
"""从已生成的 Ag/Ti3SiC2 原子模型重建 Ti 相并提取整体和局部骨架。"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch
import numpy as np
from mpl_toolkits.mplot3d import proj3d
from scipy import ndimage as ndi
from skimage.morphology import skeletonize
from vtk.util.numpy_support import vtk_to_numpy
from ase.io import read, write

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
GEOMETRY_TOOLS = REPO / "models/ag_ti3sic2_qsgs_200A_spacing10x"
sys.path.insert(0, str(GEOMETRY_TOOLS))
import spacing10x_geometry_utils as geom

OUT = HERE / "output"
FIG = HERE / "visualization"
ATOMS_FILE = OUT / "ag_ti3sic2_spacing10x.extxyz"
MODEL_REPORT = OUT / "atomic_model_report.json"
SIZE_A = np.array([200.0, 200.0, 200.0])
BOX_A = np.array([200.0, 200.0, 300.0])
GRID = np.array([80, 80, 80])
BOX_GRID = (80, 80, 120)
VOXEL_A = 2.5
TUBE_RADIUS_A = 7.5
DISPLAY_SIGMA = 0.7
LOCAL_SIZE = np.array([52, 52, 52])
FONT_CJK = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
if FONT_CJK.exists():
    font_manager.fontManager.addfont(str(FONT_CJK))
    plt.rcParams["font.family"] = [font_manager.FontProperties(fname=str(FONT_CJK)).get_name()]
plt.rcParams["axes.unicode_minus"] = False


def largest_periodic_component(mask: np.ndarray) -> tuple[np.ndarray, int, float]:
    """保留 26 邻接下最大的 x/y 周期、z 开放连通相。"""
    structure = ndi.generate_binary_structure(3, 3)
    labels, count = ndi.label(mask, structure=structure)
    parent = np.arange(count + 1, dtype=np.int64)

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = int(parent[i])
        return i

    for axis in (0, 1):
        low = np.take(labels, 0, axis=axis)
        high = np.take(labels, -1, axis=axis)
        for a, b in zip(low.ravel(), high.ravel()):
            if a and b:
                ra, rb = root(int(a)), root(int(b))
                if ra != rb:
                    parent[rb] = ra
    mapping = np.fromiter((root(i) for i in range(count + 1)), dtype=np.int64)
    labels = mapping[labels]
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    sizes[0] = 0
    selected = int(np.argmax(sizes))
    result = labels == selected
    return result, int(np.count_nonzero(sizes)), float(result.sum() / mask.sum()) if mask.any() else 0.0


def atomic_voxel_mask(ti_positions: np.ndarray, target_fraction: float) -> tuple[np.ndarray, dict]:
    """以最终 Ti 原子位置重建相场；体积分数仅用于体素采样的阈值校准。"""
    # 将最终 Ti/Si/C 原子计数沉积到 2.5 Å 网格，再轻微平滑晶格尺度的点阵起伏。
    # 这样恢复连续相界面，不会把单个晶格空位误认成枝杈或孔洞。
    voxel = np.floor(ti_positions / VOXEL_A).astype(np.int32)
    voxel[:, :2] %= GRID[:2]
    voxel[:, 2] = np.clip(voxel[:, 2], 0, GRID[2] - 1)
    atom_density = np.zeros(tuple(GRID), dtype=np.float32)
    np.add.at(atom_density, tuple(voxel.T), 1.0)
    atom_density = ndi.gaussian_filter(
        atom_density, sigma=0.8, mode=("wrap", "wrap", "nearest")
    )
    n_select = int(round(float(target_fraction) * atom_density.size))
    if not 0 < n_select < atom_density.size:
        raise ValueError(f"原子体积分数非法：{target_fraction}")
    chosen = np.argpartition(atom_density.ravel(), atom_density.size - n_select)[-n_select:]
    raw = np.zeros(atom_density.size, dtype=bool)
    raw[chosen] = True
    raw = raw.reshape(tuple(GRID))
    mask, components, largest_share = largest_periodic_component(raw)
    threshold = float(np.min(atom_density.ravel()[chosen]))
    return mask, {
        "method": "deposit final Ti/Si/C atoms onto 2.5 A voxels, Gaussian smooth sigma=0.8 voxel, threshold to final Ti phase volume fraction",
        "target_volume_fraction_from_final_atoms": float(target_fraction),
        "selected_voxels_before_connectivity_cleanup": int(n_select),
        "smoothed_atom_density_threshold_atoms_per_voxel": threshold,
        "raw_component_count": int(components),
        "largest_component_share_before_cleanup": largest_share,
        "extracted_voxels": int(mask.sum()),
        "extracted_volume_fraction": float(mask.mean()),
    }


def periodic_tube(centerline: np.ndarray) -> np.ndarray:
    """将中心线扩成便于观察的显示管，并正确处理 x/y 周期边界。"""
    tiled = np.tile(centerline, (3, 3, 1))
    distance = ndi.distance_transform_edt(~tiled, sampling=(VOXEL_A,) * 3)
    central = distance[GRID[0]:2*GRID[0], GRID[1]:2*GRID[1], :]
    return central <= TUBE_RADIUS_A


def mesh_arrays(poly) -> tuple[np.ndarray, np.ndarray]:
    vertices = vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float32, copy=True)
    faces = vtk_to_numpy(poly.GetPolys().GetData()).reshape(-1, 4)[:, 1:].astype(np.int32, copy=True)
    return vertices, faces


def make_mesh(volume: np.ndarray, periodic_xy: bool = True):
    mode = ("wrap", "wrap", "nearest") if periodic_xy else "constant"
    field = ndi.gaussian_filter(volume.astype(np.float32), sigma=DISPLAY_SIGMA, mode=mode)
    vertices, faces = geom.mesh_from_field(field, np.zeros(3))
    return geom.mesh_polydata(vertices, faces)


def pick_display_shift(gap_points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """平移周期显示原点，让测量间隙和局部截取框处在模型中央。"""
    midpoint = gap_points.mean(axis=0)
    shift_xy = np.mod(SIZE_A[:2] / 2.0 - midpoint[:2], SIZE_A[:2])
    shift_voxels = np.rint(shift_xy / VOXEL_A).astype(int) % GRID[:2]
    shown = gap_points.copy()
    shown[:, :2] = (shown[:, :2] + shift_voxels * VOXEL_A) % SIZE_A[:2]
    shown_midpoint = (midpoint[:2] + shift_voxels * VOXEL_A) % SIZE_A[:2]
    for row in shown:
        row[:2] += SIZE_A[:2] * np.rint((shown_midpoint - row[:2]) / SIZE_A[:2])
    return shift_voxels, shown


def draw_roi(ax, start: np.ndarray, size: np.ndarray) -> None:
    lo = start.astype(float)
    hi = lo + size
    corners = np.array([
        [lo[0],lo[1],lo[2]],[hi[0],lo[1],lo[2]],[hi[0],hi[1],lo[2]],[lo[0],hi[1],lo[2]],
        [lo[0],lo[1],hi[2]],[hi[0],lo[1],hi[2]],[hi[0],hi[1],hi[2]],[lo[0],hi[1],hi[2]],
    ])
    edges = ((0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7))
    for i, j in edges:
        ax.plot(*corners[[i,j]].T, color="#2370c8", linestyle=(0,(3,3)), linewidth=1.0)


def render_figure(whole_poly, local_poly, gap_points_local: np.ndarray,
                  gap_nm: float, roi_start: np.ndarray) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    whole_v, whole_f = mesh_arrays(whole_poly)
    local_v, local_f = mesh_arrays(local_poly)
    fig = plt.figure(figsize=(15.5, 9.0), facecolor="white")
    fig.text(0.055, 0.944, r"Ti$_3$SiC$_2$ 原子模型骨架与局部分叉",
             fontsize=23, fontweight="bold", color="#1e2d39")
    fig.text(0.056, 0.908, "从完整 Ag/Ti3SiC2 原子模型的 Ti/Si/C 原子位置提取；蓝色虚框标出局部截取范围",
             fontsize=11.5, color="#657684")

    left = fig.add_subplot(1, 2, 1, projection="3d")
    left.add_collection3d(geom.shaded_collection(whole_v, whole_f, "#8299aa"))
    geom.draw_box(left, tuple(BOX_A), color="#1c2730", linewidth=1.0)
    plane = np.array([[0,0,200],[200,0,200],[200,200,200],[0,200,200],[0,0,200]])
    left.plot(*plane.T, color="#536574", linestyle="--", linewidth=0.8)
    draw_roi(left, roi_start, LOCAL_SIZE * VOXEL_A)
    geom.set_3d_axes(left, ((0,200),(0,200),(0,300)), camera=(22,-54))
    left.set_title("整体骨架", loc="left", pad=2, fontsize=15, fontweight="bold", color="#273846")
    left.text2D(0.03, 0.91, "材料区 200 × 200 × 200 Å · 20 nm\n模拟盒 200 × 200 × 300 Å · 30 nm",
                transform=left.transAxes, fontsize=9.5, color="#526574")

    right = fig.add_subplot(1, 2, 2, projection="3d")
    right.add_collection3d(geom.shaded_collection(local_v, local_f, "#8299aa"))
    local_extent = LOCAL_SIZE * VOXEL_A
    geom.draw_box(right, tuple(local_extent), color="#2370c8", linewidth=1.2)
    p0, p1 = gap_points_local
    geom.set_3d_axes(right, ((0,local_extent[0]),(0,local_extent[1]),(0,local_extent[2])), camera=(42,63))
    right.set_title("局部分叉与枝间距", loc="left", pad=2,
                    fontsize=15, fontweight="bold", color="#273846")
    right.text2D(0.03, 0.91, "截取范围 13 × 13 × 13 nm\n橙色标尺：两枝表面间距",
                 transform=right.transAxes, fontsize=9.5, color="#526574")

    fig.canvas.draw()
    projected = [proj3d.proj_transform(*point, right.get_proj())[:2] for point in (p0, p1)]
    figure_points = [fig.transFigure.inverted().transform(right.transData.transform(p)) for p in projected]
    fig.add_artist(FancyArrowPatch(figure_points[0], figure_points[1], transform=fig.transFigure,
                                   arrowstyle="<->", mutation_scale=15, linewidth=2.4,
                                   color="#ef7c32", zorder=1000))
    midpoint = np.mean(figure_points, axis=0)
    fig.text(float(midpoint[0]), float(midpoint[1] + 0.022), f"{gap_nm:.2f} nm",
             color="#d96525", fontsize=11, fontweight="bold", ha="center", va="center",
             zorder=1001, bbox={"facecolor":"white", "edgecolor":"none", "alpha":0.96, "pad":2.2})
    fig.subplots_adjust(left=0.02, right=0.985, top=0.87, bottom=0.045, wspace=0.06)
    prefix = FIG / "Ti3SiC2_spacing10x_atomic_whole_and_local"
    fig.savefig(prefix.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(prefix.with_suffix(".pdf"), dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    report = json.loads(MODEL_REPORT.read_text())
    atoms = read(ATOMS_FILE)
    phases = atoms.arrays["phase"]
    ti_positions = atoms.positions[phases == 1].copy()
    target_fraction = report["source_qsgs"]["mask_volume_fraction_after_two_phase_connectivity_cleanup"]

    # 第一步：由最终原子坐标重建相场；第二步：在此相场中提取中轴骨架。
    atomic_mask, reconstruction = atomic_voxel_mask(ti_positions, target_fraction)
    centerline = skeletonize(atomic_mask, method="lee")
    tube = periodic_tube(centerline)
    skeleton_metrics = {
        "centerline_voxels": int(centerline.sum()),
        "branchpoint_voxels_degree_ge_4": 0,
        "branchpoint_clusters": 0,
        "endpoint_voxels": 0,
        "display_tube_radius_A": TUBE_RADIUS_A,
        "display_tube_diameter_nm": 2.0 * TUBE_RADIUS_A / 10.0,
        "display_tube_fraction": float(tube.mean()),
    }
    neighbors = np.ones((3,3,3), dtype=np.uint8)
    neighbors[1,1,1] = 0
    degree = ndi.convolve(centerline.astype(np.uint8), neighbors, mode="constant")
    branch = centerline & (degree >= 4)
    _, branch_clusters = ndi.label(branch, structure=ndi.generate_binary_structure(3,3))
    skeleton_metrics["branchpoint_voxels_degree_ge_4"] = int(branch.sum())
    skeleton_metrics["branchpoint_clusters"] = int(branch_clusters)
    skeleton_metrics["endpoint_voxels"] = int((centerline & (degree == 1)).sum())

    # 枝间距从原子重建的 Ti3SiC2 相表面测量，而非从人为膨胀的中轴显示管测量。
    gap = geom.find_display_gap(atomic_mask)
    points = np.asarray(gap["surface_points_xyz_A"], dtype=float)
    shift_xy, displayed_points = pick_display_shift(points)
    shift_voxels = np.array([shift_xy[0], shift_xy[1], 0], dtype=int)
    phase_display = np.roll(np.roll(atomic_mask, shift_xy[0], axis=0), shift_xy[1], axis=1)
    center_display = np.roll(np.roll(centerline, shift_xy[0], axis=0), shift_xy[1], axis=1)
    center = displayed_points.mean(axis=0)
    start_voxel = np.rint(center / VOXEL_A - LOCAL_SIZE / 2).astype(int)
    start_voxel = np.clip(start_voxel, [0,0,0], GRID - LOCAL_SIZE)
    roi_start = start_voxel * VOXEL_A
    sx, sy, sz = start_voxel
    nx, ny, nz = LOCAL_SIZE
    local_mask = phase_display[sx:sx+nx, sy:sy+ny, sz:sz+nz]
    if local_mask.shape != tuple(LOCAL_SIZE):
        raise RuntimeError("局部骨架区域超出完整相场")
    local_points = displayed_points - roi_start
    local_poly = make_mesh(local_mask, periodic_xy=False)
    whole_poly = make_mesh(phase_display, periodic_xy=True)

    atomfile = OUT / "ag_ti3sic2_spacing10x.extxyz"
    local_atoms = read(atomfile)
    positions = local_atoms.positions.copy()
    positions[:, :2] = (positions[:, :2] + shift_xy * VOXEL_A) % SIZE_A[:2]
    in_roi = ((positions[:, 0] >= roi_start[0]) & (positions[:, 0] < roi_start[0] + LOCAL_SIZE[0]*VOXEL_A) &
              (positions[:, 1] >= roi_start[1]) & (positions[:, 1] < roi_start[1] + LOCAL_SIZE[1]*VOXEL_A) &
              (positions[:, 2] >= roi_start[2]) & (positions[:, 2] < roi_start[2] + LOCAL_SIZE[2]*VOXEL_A))
    local_atoms = local_atoms[in_roi]
    local_positions = positions[in_roi] - roi_start
    local_atoms.set_positions(local_positions)
    local_atoms.set_cell(LOCAL_SIZE * VOXEL_A)
    local_atoms.set_pbc([False, False, False])
    write(OUT / "local_branch_atoms.extxyz", local_atoms, format="extxyz")

    geom.write_stl(OUT / "Ti3SiC2_skeleton_from_atoms.stl", whole_poly)
    geom.write_stl(OUT / "Ti3SiC2_local_branch_from_atoms.stl", local_poly)
    np.savez_compressed(
        OUT / "atomic_skeleton_extraction.npz",
        Ti3SiC2_voxel_mask=atomic_mask.astype(np.uint8),
        skeleton_centerline=center_display.astype(np.uint8),
        skeleton_display_tube=np.roll(np.roll(tube, shift_xy[0], axis=0), shift_xy[1], axis=1).astype(np.uint8),
        voxel_spacing_A=np.asarray([VOXEL_A]*3),
        display_shift_voxels=shift_voxels,
        local_roi_start_A=roi_start,
        local_roi_size_A=LOCAL_SIZE*VOXEL_A,
    )
    render_figure(whole_poly, local_poly, local_points, float(gap["gap_nm"]), roi_start)

    gap["surface_points_display_A"] = displayed_points.tolist()
    gap["ratio_to_source_gap"] = float(gap["gap_nm"] / 0.93719)
    final = {
        "status": "SKELETON_EXTRACTED_FROM_FINAL_ATOMISTIC_MODEL",
        "input_atomic_model": str(ATOMS_FILE.name),
        "input_atom_count": int(len(atoms)),
        "Ti3SiC2_atom_count": int(len(ti_positions)),
        "atomic_to_voxel_reconstruction": reconstruction,
        "skeleton_metrics": skeleton_metrics,
        "measured_local_branch_surface_gap": gap,
        "local_roi": {
            "start_voxel_xyz": start_voxel.tolist(),
            "size_voxels_xyz": LOCAL_SIZE.tolist(),
            "size_A": (LOCAL_SIZE*VOXEL_A).tolist(),
            "local_atom_count_including_Ag": int(len(local_atoms)),
        },
        "outputs": {
            "whole_skeleton_stl": "Ti3SiC2_skeleton_from_atoms.stl",
            "local_skeleton_stl": "Ti3SiC2_local_branch_from_atoms.stl",
            "local_atomic_crop": "local_branch_atoms.extxyz",
            "extraction_npz": "atomic_skeleton_extraction.npz",
            "paper_figure_png": "../visualization/Ti3SiC2_spacing10x_atomic_whole_and_local.png",
        },
        "method_note": "骨架相表面与中轴均由最终 extxyz 中 phase=1 的 Ti/Si/C 原子坐标重建；枝间距在平滑后的实际 Ti3SiC2 相表面间测量。中轴显示管半径仅用于中轴拓扑检查，不用于枝间距标注。",
    }
    (OUT / "skeleton_extraction_report.json").write_text(json.dumps(final, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({
        "status": final["status"],
        "atom_count": len(atoms),
        "Ti3SiC2_voxel_fraction": reconstruction["extracted_volume_fraction"],
        "branch_gap_nm": gap["gap_nm"],
        "ratio_to_source": gap["ratio_to_source_gap"],
        "centerline_voxels": skeleton_metrics["centerline_voxels"],
        "local_atoms": len(local_atoms),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
