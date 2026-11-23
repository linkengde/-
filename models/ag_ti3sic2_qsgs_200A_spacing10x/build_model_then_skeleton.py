#!/usr/bin/env python3
"""先构建固定尺寸的双相粗化样本，再提取骨架和局部分叉。"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch
import numpy as np
from mpl_toolkits.mplot3d import proj3d
from scipy import ndimage as ndi
from skimage.morphology import skeletonize
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import spacing10x_geometry_utils as geom

FONT_CJK = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
if FONT_CJK.exists():
    font_manager.fontManager.addfont(str(FONT_CJK))
    plt.rcParams["font.family"] = [font_manager.FontProperties(fname=str(FONT_CJK)).get_name()]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SOURCE = REPO / "models/ag_ti3sic2_qsgs_200A_ttm_finer/qsgs_morphology.npz"
OUT = HERE / "output"
FIG = HERE / "visualization"

VOXEL_A = 2.5
MATERIAL_GRID = (80, 80, 80)
BOX_GRID = (80, 80, 120)
MATERIAL_A = (200.0, 200.0, 200.0)
BOX_A = (200.0, 200.0, 300.0)

# 从原始 QSGS 完整场中选取 20 Å 立方采样窗；将其映射到 200 Å 材料区。
# 新模型始终先完整生成，再从 Ti3SiC2 相抽取中心线骨架。
SOURCE_CROP_START = (24, 53, 41)
SOURCE_CROP_VOXELS = 8
SPATIAL_SCALE = 10
MODEL_THRESHOLD = 0.5

# 骨架只用于可视化，中心线膨胀为半径 7.5 Å 的显示管径。
SKELETON_RADIUS_A = 7.5
DISPLAY_SIGMA = 0.7

# 为了将目标分叉放在视图中央，仅滚动周期显示坐标，不改变模型尺寸。
DISPLAY_SHIFT_X_A = 50.0
DISPLAY_SHIFT_Y_A = 100.0
LOCAL_START = (16, 0, 24)
LOCAL_SIZE = (52, 52, 52)  # 130 Å = 13 nm


def largest_component(mask: np.ndarray) -> tuple[np.ndarray, int, float]:
    """保留 26 邻接下最大连通区，并返回连通区数及占比。"""
    labels, count = ndi.label(mask, structure=ndi.generate_binary_structure(3, 3))
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    selected = int(np.argmax(sizes))
    result = labels == selected
    share = float(result.sum() / mask.sum()) if mask.any() else 0.0
    return result, int(count), share


def build_full_model(source: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    """从完整 QSGS 场缩放出目标全盒的 Ti3SiC2/Ag 双相体素模型。"""
    sx, sy, sz = SOURCE_CROP_START
    n = SOURCE_CROP_VOXELS
    window = source[sx:sx+n, sy:sy+n, sz:sz+n].astype(np.float32)
    if window.shape != (n, n, n):
        raise ValueError(f"QSGS 源窗口尺寸异常：{window.shape}")
    scaled = ndi.zoom(window, SPATIAL_SCALE, order=1, mode="nearest", prefilter=False)
    if scaled.shape != MATERIAL_GRID:
        raise ValueError(f"缩放后的材料区尺寸异常：{scaled.shape}")
    raw_ti = scaled >= MODEL_THRESHOLD
    ti, ti_count, ti_share = largest_component(raw_ti)
    raw_ag = ~raw_ti
    _, ag_count, ag_share = largest_component(raw_ag)
    # 两相分区互补。仅移除缩放插值产生的非主 Ti 相小孤岛。
    ag = ~ti
    return ti, ag, {
        "source_crop_ti_fraction": float(window.mean()),
        "scaled_ti_fraction_before_cleanup": float((scaled >= MODEL_THRESHOLD).mean()),
        "ti_components_before_cleanup": ti_count,
        "ti_largest_component_share": ti_share,
        "ag_components_before_cleanup": ag_count,
        "ag_largest_component_share_before_partition": ag_share,
    }


def extract_centerline_and_tubes(ti_phase: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    """从已建好的 Ti3SiC2 相提取中轴，再生成便于看清的枝状表面。"""
    centerline = skeletonize(ti_phase, method="lee")
    distances = ndi.distance_transform_edt(~centerline, sampling=(VOXEL_A,) * 3)
    tube = distances <= SKELETON_RADIUS_A
    tube, components, largest_share = largest_component(tube)
    neighbors = np.ones((3, 3, 3), dtype=np.uint8)
    neighbors[1, 1, 1] = 0
    degree = ndi.convolve(centerline.astype(np.uint8), neighbors, mode="constant")
    branch_voxels = centerline & (degree >= 4)
    branch_labels, branch_count = ndi.label(
        branch_voxels, structure=ndi.generate_binary_structure(3, 3)
    )
    endpoint_voxels = centerline & (degree == 1)
    metrics = {
        "centerline_voxels": int(centerline.sum()),
        "centerline_branchpoint_voxels_26_neighbor_degree_ge_4": int(branch_voxels.sum()),
        "centerline_branchpoint_clusters": int(branch_count),
        "centerline_endpoint_voxels": int(endpoint_voxels.sum()),
        "display_tube_radius_A": SKELETON_RADIUS_A,
        "display_tube_diameter_nm": 2.0 * SKELETON_RADIUS_A / 10.0,
        "display_tube_volume_fraction": float(tube.mean()),
        "display_tube_components": components,
        "display_tube_largest_component_share": largest_share,
    }
    return centerline, tube, metrics


def make_mesh(volume: np.ndarray, sigma: float = DISPLAY_SIGMA,
              periodic_xy: bool = True) -> vtk.vtkPolyData:
    """从体素相场生成平滑显示用表面。"""
    mode = ("wrap", "wrap", "nearest") if periodic_xy else "constant"
    field = ndi.gaussian_filter(volume.astype(np.float32), sigma=sigma, mode=mode)
    vertices, faces = geom.mesh_from_field(field, np.zeros(3))
    return geom.mesh_polydata(vertices, faces)


def poly_arrays(poly: vtk.vtkPolyData) -> tuple[np.ndarray, np.ndarray]:
    vertices = vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float32, copy=True)
    faces = vtk_to_numpy(poly.GetPolys().GetData()).reshape(-1, 4)[:, 1:].astype(np.int32, copy=True)
    return vertices, faces


def transform_points_for_display(points: np.ndarray) -> np.ndarray:
    result = np.asarray(points, dtype=float).copy()
    result[:, 0] = (result[:, 0] + DISPLAY_SHIFT_X_A) % MATERIAL_A[0]
    result[:, 1] = (result[:, 1] + DISPLAY_SHIFT_Y_A) % MATERIAL_A[1]
    return result


def draw_roi(ax, start_a: np.ndarray, size_a: np.ndarray) -> None:
    lo = start_a
    hi = start_a + size_a
    corners = np.array([
        [lo[0],lo[1],lo[2]],[hi[0],lo[1],lo[2]],[hi[0],hi[1],lo[2]],[lo[0],hi[1],lo[2]],
        [lo[0],lo[1],hi[2]],[hi[0],lo[1],hi[2]],[hi[0],hi[1],hi[2]],[lo[0],hi[1],hi[2]],
    ])
    edges = ((0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7))
    for i, j in edges:
        ax.plot(*corners[[i,j]].T, color="#2370c8", linestyle=(0,(3,3)), linewidth=1.0)


def render_figure(whole_poly: vtk.vtkPolyData, local_poly: vtk.vtkPolyData,
                  gap: dict, roi_start: np.ndarray, roi_size: np.ndarray) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    whole_v, whole_f = poly_arrays(whole_poly)
    local_v, local_f = poly_arrays(local_poly)
    fig = plt.figure(figsize=(15.5, 9.0), facecolor="white")
    fig.text(0.055, 0.944, r"Ti$_3$SiC$_2$ 粗化骨架与局部枝状结构",
             fontsize=23, fontweight="bold", color="#1e2d39")
    fig.text(0.056, 0.908, "完整双相模型内提取 Ti3SiC2 骨架；蓝色虚线框为局部截取范围",
             fontsize=11.5, color="#657684")

    left = fig.add_subplot(1, 2, 1, projection="3d")
    left.add_collection3d(geom.shaded_collection(whole_v, whole_f, "#8299aa"))
    geom.draw_box(left, BOX_A, color="#1c2730", linewidth=1.0)
    # 材料顶部与 z 向真空区的分界。
    plane = np.array([[0,0,200],[200,0,200],[200,200,200],[0,200,200],[0,0,200]])
    left.plot(*plane.T, color="#536574", linestyle="--", linewidth=0.8)
    draw_roi(left, roi_start, roi_size)
    geom.set_3d_axes(left, ((0,200),(0,200),(0,300)), camera=(22,-54))
    left.set_title("整体骨架", loc="left", pad=2, fontsize=15, fontweight="bold", color="#273846")
    left.text2D(0.03, 0.91, "材料区 200 × 200 × 200 Å · 20 nm\n模拟盒 200 × 200 × 300 Å · 30 nm",
                transform=left.transAxes, fontsize=9.5, color="#526574")

    right = fig.add_subplot(1, 2, 2, projection="3d")
    right.add_collection3d(geom.shaded_collection(local_v, local_f, "#8299aa"))
    local_extent = roi_size
    geom.draw_box(right, tuple(local_extent), color="#2370c8", linewidth=1.2)
    points = np.asarray(gap["surface_points_display_A"], dtype=float).T
    points -= roi_start[:, None]
    p0, p1 = points[:,0], points[:,1]
    delta = p1 - p0
    distance = float(np.linalg.norm(delta))
    midpoint = (p0 + p1) / 2.0
    geom.set_3d_axes(right, ((0,local_extent[0]),(0,local_extent[1]),(0,local_extent[2])), camera=(42,63))
    # 将三维测量端点投影到图面后覆盖绘制，确保标尺和数值不会被骨架遮挡。
    fig.canvas.draw()
    projected = [proj3d.proj_transform(*point, right.get_proj())[:2] for point in (p0, p1)]
    figure_points = [fig.transFigure.inverted().transform(right.transData.transform(point))
                     for point in projected]
    marker = FancyArrowPatch(figure_points[0], figure_points[1], transform=fig.transFigure,
                             arrowstyle="<->", mutation_scale=15, linewidth=2.4,
                             color="#ef7c32", zorder=1000)
    fig.add_artist(marker)
    figure_midpoint = np.mean(figure_points, axis=0)
    fig.text(float(figure_midpoint[0]), float(figure_midpoint[1] + 0.022),
             f"{distance/10.0:.2f} nm", color="#d96525", fontsize=11,
             fontweight="bold", ha="center", va="center", zorder=1001,
             bbox={"facecolor":"white", "edgecolor":"none", "alpha":0.96, "pad":2.2})
    right.set_title("局部分叉与枝间距", loc="left", pad=2,
                    fontsize=15, fontweight="bold", color="#273846")
    right.text2D(0.03, 0.91, "截取范围 13 × 13 × 13 nm\n橙色标尺：两枝表面间距",
                 transform=right.transAxes, fontsize=9.5, color="#526574")

    fig.subplots_adjust(left=0.02, right=0.985, top=0.87, bottom=0.045, wspace=0.06)
    prefix = FIG / "Ti3SiC2_spacing10x_whole_and_local"
    fig.savefig(prefix.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(prefix.with_suffix(".pdf"), dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    archive = np.load(SOURCE)
    source_mask = np.asarray(archive["skeleton"], dtype=bool)
    if source_mask.shape != (80,80,80):
        raise ValueError(f"源 QSGS 网格异常：{source_mask.shape}")

    # 1) 先生成完整 20 nm 双相材料区；再补上 10 nm 真空层成为 30 nm 模拟盒。
    ti_phase, ag_phase, phase_metrics = build_full_model(source_mask)
    # 2) 从已建好的 Ti3SiC2 相提取中心线，再构成可视化骨架表面。
    centerline, skeleton_tube, skeleton_metrics = extract_centerline_and_tubes(ti_phase)
    # 3) 在完整模型坐标中滚动周期显示，之后从同一整副骨架截取局部。
    shift_x = int(round(DISPLAY_SHIFT_X_A / VOXEL_A))
    shift_y = int(round(DISPLAY_SHIFT_Y_A / VOXEL_A))
    ti_display = np.roll(np.roll(ti_phase, shift_x, axis=0), shift_y, axis=1)
    ag_display = np.roll(np.roll(ag_phase, shift_x, axis=0), shift_y, axis=1)
    center_display = np.roll(np.roll(centerline, shift_x, axis=0), shift_y, axis=1)
    skeleton_display = np.roll(np.roll(skeleton_tube, shift_x, axis=0), shift_y, axis=1)

    gap = geom.find_display_gap(skeleton_tube)
    gap["surface_points_display_A"] = transform_points_for_display(
        np.asarray(gap["surface_points_xyz_A"], dtype=float)
    ).tolist()
    gap["ratio_to_original_spacing"] = float(gap["gap_nm"] / (geom.SOURCE_GAP_A / 10.0))

    # 13 nm 截取框围住两处相对枝面，整块 ROI 直接来自已生成的全模型骨架。
    roi_start = np.asarray(LOCAL_START, dtype=float) * VOXEL_A
    roi_size = np.asarray(LOCAL_SIZE, dtype=float) * VOXEL_A
    sx, sy, sz = LOCAL_START
    nx, ny, nz = LOCAL_SIZE
    local_mask = skeleton_display[sx:sx+nx, sy:sy+ny, sz:sz+nz]
    if local_mask.shape != tuple(LOCAL_SIZE):
        raise ValueError("局部骨架截取超出完整材料区")

    # 完整骨架单独构网；局部掩膜来自同一完整骨架的体素裁切，再独立构成显示表面。
    whole_poly = make_mesh(skeleton_display)
    local_poly = make_mesh(local_mask, periodic_xy=False)

    # 输出包含完整 Ti3SiC2/Ag 相场、抽取中心线及展示用骨架管。
    labels = np.zeros(BOX_GRID, dtype=np.uint8)
    labels[:,:,:80] = np.where(ti_display, 1, 2)
    ti_box = np.zeros(BOX_GRID, dtype=np.uint8); ti_box[:,:,:80] = ti_display.astype(np.uint8)
    ag_box = np.zeros(BOX_GRID, dtype=np.uint8); ag_box[:,:,:80] = ag_display.astype(np.uint8)
    skel_box = np.zeros(BOX_GRID, dtype=np.uint8); skel_box[:,:,:80] = skeleton_display.astype(np.uint8)
    np.savez_compressed(
        OUT / "spacing10x_phase_field.npz",
        phase_label=labels,
        model_ti3sic2=ti_box,
        model_ag=ag_box,
        skeleton_centerline=np.pad(center_display, ((0,0),(0,0),(0,40))),
        ti3sic2_skeleton=skel_box,
        voxel_spacing_A=np.asarray([VOXEL_A]*3),
        material_dimensions_A=np.asarray(MATERIAL_A),
        simulation_box_A=np.asarray(BOX_A),
        labels=np.asarray(["0=vacuum", "1=Ti3SiC2", "2=Ag"]),
    )
    geom.write_vtk_field(OUT / "spacing10x_model_ti3sic2.vtk", ti_box, "Ti3SiC2")
    geom.write_vtk_field(OUT / "spacing10x_model_ag.vtk", ag_box, "Ag")
    geom.write_vtk_field(OUT / "spacing10x_composite.vtk", labels, "Phase")
    geom.write_vtk_field(OUT / "spacing10x_extracted_skeleton.vtk", skel_box, "Ti3SiC2_skeleton")
    geom.write_stl(OUT / "spacing10x_ti3sic2_skeleton.stl", whole_poly)
    geom.write_stl(OUT / "spacing10x_local_skeleton.stl", local_poly)

    # 将两相分别导出；NPZ/VTK 保留完整相标签，便于后续分析。
    ti_poly = make_mesh(ti_display)
    ag_poly = make_mesh(ag_display)
    geom.write_stl(OUT / "spacing10x_model_ti3sic2.stl", ti_poly)
    geom.write_stl(OUT / "spacing10x_model_ag.stl", ag_poly)

    render_figure(whole_poly, local_poly, gap, roi_start, roi_size)

    ti_components, ti_sizes = geom.periodic_labels(ti_display)
    ag_components, ag_sizes = geom.periodic_labels(ag_display)
    tube_components, tube_sizes = geom.periodic_labels(skeleton_display)
    report = {
        "status": "FULL_TWO_PHASE_GEOMETRY_BUILT_THEN_SKELETON_EXTRACTED",
        "source": {
            "file": str(SOURCE.relative_to(REPO)),
            "sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            "source_grid": [80,80,80],
            "source_material_A": [200.0,200.0,200.0],
        },
        "construction": {
            "order": ["build_scaled_two_phase_model", "extract_Ti3SiC2_centerline", "crop_local_skeleton_from_full_model"],
            "source_crop_start_voxel_xyz": list(SOURCE_CROP_START),
            "source_crop_extent_voxels_xyz": [SOURCE_CROP_VOXELS]*3,
            "source_crop_extent_A": [20.0]*3,
            "spatial_scale_factor": SPATIAL_SCALE,
            "target_material_region_A": list(MATERIAL_A),
            "target_simulation_box_A": list(BOX_A),
            "grid_xyz": list(BOX_GRID),
            "voxel_spacing_A": VOXEL_A,
            "nominal_pbc": [True,True,False],
            "display_coordinate_shift_A": [DISPLAY_SHIFT_X_A,DISPLAY_SHIFT_Y_A,0.0],
            "local_roi_start_voxel_xyz": list(LOCAL_START),
            "local_roi_extent_voxels_xyz": list(LOCAL_SIZE),
            "local_roi_extent_nm_xyz": (roi_size/10.0).tolist(),
        },
        "full_two_phase_model": {
            **phase_metrics,
            "Ti3SiC2_volume_fraction": float(ti_display.mean()),
            "Ag_volume_fraction": float(ag_display.mean()),
            "Ti3SiC2_periodic_xy_components": int(np.count_nonzero(ti_sizes)),
            "Ag_periodic_xy_components": int(np.count_nonzero(ag_sizes)),
            "Ti3SiC2_largest_periodic_xy_share": float(ti_sizes.max()/ti_display.sum()),
            "Ag_largest_periodic_xy_share": float(ag_sizes.max()/ag_display.sum()),
        },
        "extracted_skeleton": {
            **skeleton_metrics,
            "periodic_xy_components": int(np.count_nonzero(tube_sizes)),
            "largest_periodic_xy_share": float(tube_sizes.max()/skeleton_display.sum()),
            "measured_local_branch_gap": gap,
        },
        "method_note": "这是几何形貌模型。枝间距在从完整粗化双相模型提取的骨架显示表面上测量；骨架管径是可视化厚度，不等同于物理原子结构。未进行原子重建或 MD 弛豫。",
        "outputs": {
            "phase_field_npz": "output/spacing10x_phase_field.npz",
            "full_composite_vtk": "output/spacing10x_composite.vtk",
            "extracted_skeleton_stl": "output/spacing10x_ti3sic2_skeleton.stl",
            "local_skeleton_stl": "output/spacing10x_local_skeleton.stl",
            "extracted_skeleton_vtk": "output/spacing10x_extracted_skeleton.vtk",
            "paper_figure_png": "visualization/Ti3SiC2_spacing10x_whole_and_local.png",
            "paper_figure_pdf": "visualization/Ti3SiC2_spacing10x_whole_and_local.pdf",
        },
    }
    (OUT / "topology_metrics.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "material_A": MATERIAL_A,
        "box_A": BOX_A,
        "Ti_fraction": report["full_two_phase_model"]["Ti3SiC2_volume_fraction"],
        "skeleton_fraction": skeleton_display.mean(),
        "branch_gap_nm": gap["gap_nm"],
        "gap_ratio": gap["ratio_to_original_spacing"],
        "centerline_branch_clusters": skeleton_metrics["centerline_branchpoint_clusters"],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
