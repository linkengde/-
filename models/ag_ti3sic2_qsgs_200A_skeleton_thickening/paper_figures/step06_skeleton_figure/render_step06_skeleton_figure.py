#!/usr/bin/env python3
"""Render a paper-style whole/local view of the accepted step06 skeleton."""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from ovito import scene


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
HELPER_PATH = HERE / "skeleton_render_helper.py"
WHOLE_MESH = HERE / "Ti3SiC2_step06_whole_surface_mesh.npz"
LOCAL_MESH = HERE / "Ti3SiC2_step06_local_branch_surface_mesh.npz"
WHOLE_PNG = HERE / "Ti3SiC2_step06_whole_skeleton_OVITO.png"
LOCAL_PNG = HERE / "Ti3SiC2_step06_local_branch_spacing.png"
COMPOSITE_PNG = HERE / "Ti3SiC2_step06_whole_skeleton_with_local_branch.png"
COMPOSITE_TIF = HERE / "Ti3SiC2_step06_whole_skeleton_with_local_branch.tif"
COMPOSITE_PDF = HERE / "Ti3SiC2_step06_whole_skeleton_with_local_branch.pdf"
REPORT = HERE / "step06_figure_report.json"

# ROI 来自 step06 局部网格裁剪：起点 (40, 25, 63)，每个方向 10 个 2.5 Å 体素。
ROI = {"x": (100.0, 125.0), "y": (62.5, 87.5), "z": (157.5, 182.5)}
WHOLE_SIZE = (1600, 1600)
LOCAL_SIZE = (1100, 1100)
CANVAS_SIZE = (2920, 1900)
MAIN_ORIGIN = (20, 175)
DETAIL_ORIGIN = (1780, 455)
DETAIL_SIZE = (1000, 1000)


def load_helper():
    """导入已有 OVITO 渲染与论文版式辅助函数。"""
    spec = importlib.util.spec_from_file_location("skeleton_render_helper", HELPER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法载入渲染辅助脚本：{HELPER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.ROI = ROI
    module.WHOLE_SIZE = WHOLE_SIZE
    module.LOCAL_SIZE = LOCAL_SIZE
    return module


def mesh_component_vertex_groups(faces: np.ndarray, vertex_count: int) -> list[np.ndarray]:
    """按三角面连通性划分网格顶点，供两条相邻骨架枝条测距。"""
    parent = np.arange(vertex_count, dtype=np.int64)

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = int(parent[index])
        return index

    def union(first: int, second: int) -> None:
        root_a, root_b = find(first), find(second)
        if root_a != root_b:
            parent[root_b] = root_a

    for face in faces:
        union(int(face[0]), int(face[1]))
        union(int(face[0]), int(face[2]))
    roots = np.fromiter((find(i) for i in range(vertex_count)), dtype=np.int64)
    unique, counts = np.unique(roots, return_counts=True)
    order = np.argsort(counts)[::-1]
    return [np.flatnonzero(roots == unique[index]) for index in order]


def nearest_surface_vertices(vertices: np.ndarray, faces: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, list[int]]:
    """在局部两个主表面网格之间寻找最短顶点间距。"""
    groups = mesh_component_vertex_groups(faces, len(vertices))
    if len(groups) < 2:
        raise RuntimeError("局部裁剪没有包含两条独立骨架枝条，无法标注枝间距。")
    first_ids, second_ids = groups[0], groups[1]
    first, second = vertices[first_ids], vertices[second_ids]
    # 当前局部裁剪仅有数百个顶点，直接计算两表面顶点距离矩阵即可。
    delta = first[:, None, :] - second[None, :, :]
    distance = np.linalg.norm(delta, axis=2)
    row, column = np.unravel_index(int(np.argmin(distance)), distance.shape)
    point_a, point_b = first[row], second[column]
    return point_a, point_b, float(distance[row, column]), [len(first_ids), len(second_ids)]


def make_paper_canvas(
    helper,
    whole_image: Image.Image,
    local_image: Image.Image,
    whole_viewport,
    dimension_label: str,
) -> Image.Image:
    """排版整体模型、局部放大、选区虚线及实际枝间距标注。"""
    font_path = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
    title_font = ImageFont.truetype(font_path, 58, index=0)
    panel_font = ImageFont.truetype(font_path, 42, index=0)
    subtitle_font = ImageFont.truetype(font_path, 29, index=0)

    canvas = Image.new("RGB", CANVAS_SIZE, "white")
    whole_with_roi = whole_image.copy()
    helper.draw_roi_box(whole_with_roi, whole_viewport, WHOLE_SIZE)
    canvas.paste(whole_with_roi, MAIN_ORIGIN)
    canvas.paste(local_image.resize(DETAIL_SIZE, Image.Resampling.LANCZOS), DETAIL_ORIGIN)
    draw = ImageDraw.Draw(canvas)

    blue = helper.BLUE
    draw.rectangle(
        (
            DETAIL_ORIGIN[0] - 4,
            DETAIL_ORIGIN[1] - 4,
            DETAIL_ORIGIN[0] + DETAIL_SIZE[0] + 3,
            DETAIL_ORIGIN[1] + DETAIL_SIZE[1] + 3,
        ),
        outline=blue,
        width=5,
    )

    # 两条引线从整体图中 ROI 的同侧角点延伸到局部放大框。
    point_a = helper.projected_point((ROI["x"][1], ROI["y"][1], ROI["z"][1]), whole_viewport, WHOLE_SIZE)
    point_b = helper.projected_point((ROI["x"][1], ROI["y"][0], ROI["z"][0]), whole_viewport, WHOLE_SIZE)
    point_a = (point_a[0] + MAIN_ORIGIN[0], point_a[1] + MAIN_ORIGIN[1])
    point_b = (point_b[0] + MAIN_ORIGIN[0], point_b[1] + MAIN_ORIGIN[1])
    helper.dashed_line(draw, point_a, (DETAIL_ORIGIN[0] - 6, DETAIL_ORIGIN[1] + 120), fill=blue, width=2, dash=14, gap=10, halo=False)
    helper.dashed_line(draw, point_b, (DETAIL_ORIGIN[0] - 6, DETAIL_ORIGIN[1] + 880), fill=blue, width=2, dash=14, gap=10, halo=False)

    draw.text((70, 32), "Ti3SiC2 连通骨架与局部枝状结构", font=title_font, fill=(28, 42, 55), stroke_width=1)
    draw.text((72, 108), "蓝色虚线框标示从整体模型中截取的局部范围", font=subtitle_font, fill=(75, 91, 106))
    draw.text((95, 185), "整体骨架", font=panel_font, fill=(28, 42, 55), stroke_width=1)
    draw.text((95, 238), "200 × 200 × 200 Å  ·  20 nm", font=subtitle_font, fill=(75, 91, 106))
    draw.text((1785, 325), "局部枝状骨架", font=panel_font, fill=(28, 42, 55), stroke_width=1)
    draw.text((1788, 380), "25 × 25 × 25 Å  ·  2.5 nm", font=subtitle_font, fill=(75, 91, 106))
    return canvas


def main() -> None:
    helper = load_helper()
    renderer = helper.make_renderer()
    whole_data = np.load(WHOLE_MESH)
    local_data = np.load(LOCAL_MESH)
    whole_vertices = whole_data["vertices"].astype(float)
    local_vertices = local_data["vertices"].astype(float)
    local_faces = local_data["faces"].astype(int)

    point_a, point_b, vertex_gap_A, component_sizes = nearest_surface_vertices(local_vertices, local_faces)

    # 先渲染整体模型，并在最终视口上投影 ROI。
    whole_pipeline, _ = helper.load_mesh(WHOLE_MESH, "step06 whole Ti3SiC2 skeleton")
    whole_pipeline.add_to_scene()
    whole_viewport = helper.make_viewport((100.0, 100.0, 100.0), (-1.0, -1.08, -0.74), WHOLE_SIZE)
    scene.viewports.active_vp = whole_viewport
    whole_image = helper.render_view(whole_viewport, WHOLE_PNG, WHOLE_SIZE, renderer)
    whole_pipeline.remove_from_scene()

    # 相机沿枝间隙法向观察，使最短枝间距方向落在图像平面中且两端可见。
    local_center = tuple(sum(ROI[axis]) / 2 for axis in ("x", "y", "z"))
    gap_direction = (point_b - point_a) / np.linalg.norm(point_b - point_a)
    horizontal = np.cross(gap_direction, np.asarray([0.0, 0.0, 1.0]))
    if np.linalg.norm(horizontal) < 1.0e-8:
        horizontal = np.cross(gap_direction, np.asarray([0.0, 1.0, 0.0]))
    horizontal /= np.linalg.norm(horizontal)
    vertical = np.cross(gap_direction, horizontal)
    vertical /= np.linalg.norm(vertical)
    # 沿箭头轴选择可见面，使两端均能直视，且实际间距箭头水平显示。
    local_direction = -horizontal
    local_up = -vertical
    local_pipeline, _ = helper.load_mesh(LOCAL_MESH, "step06 local Ti3SiC2 branches")
    local_pipeline.add_to_scene()
    local_viewport = helper.make_viewport(local_center, local_direction, LOCAL_SIZE)
    local_viewport.camera_up = local_up
    scene.viewports.active_vp = local_viewport
    local_image = helper.render_view(local_viewport, LOCAL_PNG, LOCAL_SIZE, renderer)

    start_pixel = helper.projected_point(tuple(point_a), local_viewport, LOCAL_SIZE)
    end_pixel = helper.projected_point(tuple(point_b), local_viewport, LOCAL_SIZE)
    # 将箭头尖端映射回最终显示的三角网格表面，并按真实 3D 坐标计算标注值。
    measured_a = helper.surface_point_on_camera_ray(local_vertices, local_faces, start_pixel, local_viewport, LOCAL_SIZE)
    measured_b = helper.surface_point_on_camera_ray(local_vertices, local_faces, end_pixel, local_viewport, LOCAL_SIZE)
    gap_A = float(np.linalg.norm(measured_b - measured_a))
    gap_nm = gap_A / 10.0
    label = f"{gap_A:.2f} Å"
    helper.draw_branch_gap_annotation(local_image, start_pixel, end_pixel, label, ImageFont.truetype("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", 36, index=0))
    local_image.save(LOCAL_PNG, "PNG", dpi=(300, 300), optimize=True)
    local_pipeline.remove_from_scene()

    canvas = make_paper_canvas(helper, whole_image, local_image, whole_viewport, label)
    canvas.save(COMPOSITE_PNG, format="PNG", dpi=(300, 300), optimize=True)
    canvas.save(COMPOSITE_TIF, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
    canvas.save(COMPOSITE_PDF, format="PDF", resolution=300.0)

    report = {
        "status": "STEP06_PAPER_STYLE_SKELETON_FIGURE_RENDERED",
        "source_model": "step06 coarsened Ag-Ti3SiC2 atomic model; Ti3SiC2 display mesh extracted from its step06 equal-volume phase field",
        "source_phase_field": "output/step06_coarsened_networks/coarsened_step06_phase_fields.npz",
        "source_mask_key": "step06_equal_volume_Ti3SiC2_mask",
        "material_region_A": [200.0, 200.0, 200.0],
        "local_roi_A": {axis: list(limits) for axis, limits in ROI.items()},
        "local_roi_extent_A": [25.0, 25.0, 25.0],
        "whole_image_px": list(WHOLE_SIZE),
        "local_image_px": list(LOCAL_SIZE),
        "composite_image_px": list(CANVAS_SIZE),
        "local_camera_direction": list(local_direction),
        "local_camera_up": list(local_up),
        "local_surface_component_vertex_counts": component_sizes,
        "branch_spacing": {
            "quantity": "3D Euclidean distance between the two nearest visible branch surface points",
            "distance_A": round(gap_A, 6),
            "distance_nm": round(gap_nm, 6),
            "label": label,
            "vertex_seed_distance_A": round(vertex_gap_A, 6),
            "endpoint_coordinates_A": [measured_a.tolist(), measured_b.tolist()],
            "endpoint_pixels_local": [list(map(float, start_pixel)), list(map(float, end_pixel))],
            "method": "nearest vertices on the two largest face-connected surface components; final endpoints are obtained by intersecting the rendered camera rays with the actual triangle mesh, then measuring their 3D Euclidean distance",
        },
        "display_note": "The surface mesh is a smoothed visualization of the step06 Ti3SiC2 phase field; no atomic coordinates or model dimensions were changed.",
        "outputs": [WHOLE_PNG.name, LOCAL_PNG.name, COMPOSITE_PNG.name, COMPOSITE_TIF.name, COMPOSITE_PDF.name, REPORT.name],
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
