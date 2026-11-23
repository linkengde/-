#!/usr/bin/env python3
"""从归档的 Round7 原子轨迹制作图2、图4并核验数据。

本程序只做已有轨迹后处理，不运行 MD。图2读取五帧全体系 dump；图4将更正径向区
CSV、真实 1.60 ps 坐标、候选原子 dump 与最大连通分量 dump 按原子 ID 对齐，并重算配位数和连通分量。
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/round7_mplconfig")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from ovito import scene
from ovito.io import import_file
from ovito.modifiers import DeleteSelectedModifier, ExpressionSelectionModifier
from ovito.vis import ParticlesVis, StandardRenderer, Viewport


FRAME_TIMES_PS = [0.60, 0.85, 1.10, 1.35, 1.60]
TYPE_NAMES = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}
TYPE_COLORS = {
    1: (0.34, 0.50, 0.66),  # Ti：蓝灰
    2: (0.91, 0.61, 0.20),  # Si：琥珀色
    3: (0.16, 0.20, 0.25),  # C：深灰
    4: (0.49, 0.76, 0.84),  # Ag：浅青
}
ORANGE = "#e87532"
ORANGE_LIGHT = "#f4b183"
INK = "#1b2935"
MUTED = "#5f7180"
FONT_CJK = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"


def load_font(size: int) -> ImageFont.FreeTypeFont:
    """加载有中日韩字形的字体，确保图中文字可读。"""
    if Path(FONT_CJK).exists():
        return ImageFont.truetype(FONT_CJK, size, index=0)
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", size)


def save_pil_figure(image: Image.Image, stem: Path) -> None:
    """保存论文常用的 PNG、TIFF 和 PDF 三种格式。"""
    image.save(stem.with_suffix(".png"), format="PNG", dpi=(300, 300), optimize=True)
    image.save(stem.with_suffix(".tif"), format="TIFF", dpi=(300, 300), compression="tiff_lzw")
    image.save(stem.with_suffix(".pdf"), format="PDF", resolution=300.0)


def draw_scale_bar(draw: ImageDraw.ImageDraw, xy: tuple[int, int], length_px: int, label: str) -> None:
    """在白底小框内绘制经正交视口比例换算的长度标尺。"""
    x, y = xy
    draw.rounded_rectangle((x - 12, y - 38, x + length_px + 82, y + 24), radius=8,
                           fill="white", outline="#9aa9b4", width=2)
    draw.line((x, y, x + length_px, y), fill=INK, width=6)
    draw.line((x, y - 9, x, y + 9), fill=INK, width=3)
    draw.line((x + length_px, y - 9, x + length_px, y + 9), fill=INK, width=3)
    draw.text((x + length_px + 12, y - 19), label, font=load_font(20), fill=INK)


def draw_projected_roi(
    image: Image.Image,
    roi: dict[str, list[float]],
    center: np.ndarray,
    direction: np.ndarray,
    viewport_half_height_A: float,
) -> None:
    """把真实 Å 坐标定义的局部裁切盒投影成全体系面板上的蓝色虚线框。"""
    width, height = image.size
    up = np.array((0.0, 0.0, 1.0))
    up -= np.dot(up, direction) * direction
    up /= np.linalg.norm(up)
    right = np.cross(direction, up)
    right /= np.linalg.norm(right)
    px_per_A = height / (2.0 * viewport_half_height_A)
    xs, ys, zs = roi["x_A"], roi["y_A"], roi["z_A"]
    points = {}
    for ix, x in enumerate(xs):
        for iy, y in enumerate(ys):
            for iz, z in enumerate(zs):
                offset = np.array((x, y, z)) - center
                screen_x = width / 2.0 + float(np.dot(offset, right)) * px_per_A
                screen_y = height / 2.0 - float(np.dot(offset, up)) * px_per_A
                points[(ix, iy, iz)] = (round(screen_x), round(screen_y))

    draw = ImageDraw.Draw(image)
    edges = []
    for key, point in points.items():
        for axis in range(3):
            if key[axis] == 0:
                other = list(key)
                other[axis] = 1
                edges.append((point, points[tuple(other)]))
    # 白色底线提升对比度，蓝色虚线标示真正参与裁切的坐标范围。
    for start, end in edges:
        draw.line((*start, *end), fill="white", width=7)
    for start, end in edges:
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = math.hypot(dx, dy)
        if length == 0:
            continue
        ux, uy = dx / length, dy / length
        position = 0.0
        while position < length:
            stop = min(position + 10.0, length)
            draw.line((start[0] + ux * position, start[1] + uy * position,
                       start[0] + ux * stop, start[1] + uy * stop),
                      fill="#1769aa", width=4)
            position += 17.0


def render_figure2(trajectory: Path, out_dir: Path, run_summary: dict) -> dict:
    """渲染真实五帧、独立热源邻域放大图和物理长度标尺。"""
    pipeline = import_file(str(trajectory))
    if pipeline.source.num_frames != len(FRAME_TIMES_PS):
        raise ValueError(f"期望 5 帧，实际读取到 {pipeline.source.num_frames} 帧")

    frame_data = [pipeline.compute(i) for i in range(pipeline.source.num_frames)]
    atom_counts = [int(data.particles.count) for data in frame_data]
    if len(set(atom_counts)) != 1:
        raise ValueError(f"五帧原子数不一致：{atom_counts}")

    first_ids = np.asarray(frame_data[0].particles["Particle Identifier"], dtype=np.int64)
    first_order = np.argsort(first_ids)
    sorted_ids = first_ids[first_order]
    first_types = np.asarray(frame_data[0].particles["Particle Type"], dtype=np.int32)[first_order]
    first_positions = np.asarray(frame_data[0].particles["Position"], dtype=float)[first_order]
    rms_displacement = []
    timestep_records = []
    expected_steps = [6000, 8500, 11000, 13500, 16000]
    for frame_index, data in enumerate(frame_data):
        frame_ids = np.asarray(data.particles["Particle Identifier"], dtype=np.int64)
        frame_order = np.argsort(frame_ids)
        if not np.array_equal(frame_ids[frame_order], sorted_ids):
            raise ValueError(f"第 {frame_index} 帧原子 ID 集合不一致")
        types = np.asarray(data.particles["Particle Type"], dtype=np.int32)[frame_order]
        if not np.array_equal(types, first_types):
            raise ValueError(f"第 {frame_index} 帧元素类型与原子 ID 映射不一致")
        step = int(data.attributes["Timestep"])
        if step != expected_steps[frame_index]:
            raise ValueError(f"第 {frame_index} 帧 timestep 不匹配：{step}")
        timestep_records.append(step)
        positions = np.asarray(data.particles["Position"], dtype=float)[frame_order]
        displacement = positions - first_positions
        displacement[:, :2] -= np.rint(displacement[:, :2] / 200.0) * 200.0
        rms_displacement.append(float(np.sqrt(np.mean(np.einsum("ij,ij->i", displacement, displacement)))))
        cell = np.asarray(data.cell.matrix, dtype=float)
        lengths = np.linalg.norm(cell[:, :3], axis=1)
        if not np.allclose(lengths, (200.0, 200.0, 250.0), atol=1e-6):
            raise ValueError(f"第 {frame_index} 帧盒子尺寸异常：{lengths.tolist()}")
        if tuple(data.cell.pbc) != (True, True, False):
            raise ValueError(f"第 {frame_index} 帧边界条件异常：{data.cell.pbc}")

    type_counts = {int(t): int(np.count_nonzero(first_types == t)) for t in np.unique(first_types)}
    if set(type_counts) != {1, 2, 3, 4}:
        raise ValueError(f"原子类型映射与 Round7 预期不符：{type_counts}")

    # 四种元素各用独立颜色，避免把 Ti₃SiC₂ 三种元素合并成同一色。
    source = pipeline.source.data
    for particle_type in source.particles_.particle_types_.types_:
        particle_type.name = TYPE_NAMES[int(particle_type.id)]
        particle_type.color = TYPE_COLORS[int(particle_type.id)]
    source.particles_.vis.radius = 0.8
    source.particles_.vis.shape = ParticlesVis.Shape.Sphere
    source.cell_.vis.enabled = True
    source.cell_.vis.render_cell = True
    source.cell_.vis.line_width = 0.8
    source.cell_.vis.rendering_color = (0.14, 0.18, 0.22)

    pipeline.add_to_scene()
    viewport = Viewport(type=Viewport.Type.Ortho)
    scene.viewports.active_vp = viewport
    full_center = np.array((100.0, 100.0, 125.0))
    direction = np.array((-350.0, 640.0, -325.0), dtype=float)
    direction /= np.linalg.norm(direction)
    viewport.camera_pos = full_center - direction * 1000.0
    viewport.camera_dir = direction
    viewport.camera_up = (0.0, 0.0, 1.0)
    viewport.fov = 180.0  # 正交投影下是可见画面半高，Å

    renderer = StandardRenderer()
    renderer.antialiasing_level = 2
    renderer.ambient_occlusion = True
    renderer.ambient_occlusion_intensity = 0.50
    renderer.ambient_occlusion_radius = 4.0
    renderer.shadows = True
    renderer.shadow_intensity = 0.18

    panel_size = (1240, 900)
    panels: list[Image.Image] = []
    rendered_frame_info = []
    temporary_panel_dir = out_dir / ".round7_render_panels"
    temporary_panel_dir.mkdir(exist_ok=True)
    for frame_index, time_ps in enumerate(FRAME_TIMES_PS):
        scene.anim.current_frame = frame_index
        rendered = viewport.render_image(size=panel_size, background=(1, 1, 1),
                                         renderer=renderer, stop_on_error=True)
        panel_path = temporary_panel_dir / f"full_{frame_index:02d}.png"
        rendered.save(str(panel_path), "PNG")
        with Image.open(panel_path) as panel_image:
            panels.append(panel_image.convert("RGB").copy())
        panel_path.unlink(missing_ok=True)
        rendered_frame_info.append({
            "frame": frame_index,
            "time_ps": time_ps,
            "timestep": timestep_records[frame_index],
            "atoms": atom_counts[frame_index],
            "rms_net_displacement_from_0p60ps_A": round(rms_displacement[frame_index], 6),
        })

    # 局部框选只用于显示，不改变全体系面板的盒子尺寸或原子坐标。
    roi = {"x_A": [70.0, 130.0], "y_A": [70.0, 130.0], "z_A": [120.0, 175.0]}
    outside_roi = (
        "Position.X < 70 || Position.X > 130 || Position.Y < 70 || Position.Y > 130 || "
        "Position.Z < 120 || Position.Z > 175"
    )
    pipeline.modifiers.append(ExpressionSelectionModifier(expression=outside_roi))
    pipeline.modifiers.append(DeleteSelectedModifier())
    source.cell_.vis.enabled = False
    scene.anim.current_frame = 4
    local_data = pipeline.compute(4)
    local_atom_count = int(local_data.particles.count)
    if local_atom_count == 0:
        raise ValueError("热源邻域局部裁剪结果为空")
    local_center = np.array((100.0, 100.0, 147.5))
    viewport.camera_pos = local_center - direction * 400.0
    viewport.camera_dir = direction
    viewport.camera_up = (0.0, 0.0, 1.0)
    viewport.fov = 52.0
    local_rendered = viewport.render_image(size=panel_size, background=(1, 1, 1),
                                           renderer=renderer, stop_on_error=True)
    local_path = temporary_panel_dir / "local_1p60ps.png"
    local_rendered.save(str(local_path), "PNG")
    with Image.open(local_path) as panel_image:
        panels.append(panel_image.convert("RGB").copy())
    local_path.unlink(missing_ok=True)
    pipeline.remove_from_scene()
    temporary_panel_dir.rmdir()

    # 五个完整模型面板加一个同帧局部放大面板。
    canvas = Image.new("RGB", (3980, 2320), "white")
    draw = ImageDraw.Draw(canvas)
    title_font, subtitle_font = load_font(59), load_font(27)
    panel_font, legend_font, foot_font = load_font(29), load_font(22), load_font(20)
    draw.text((78, 26), "Ag/Ti3SiC2 五帧原子结构演化与热源邻域", font=title_font, fill=INK, stroke_width=1)
    draw.text((82, 99), "Round7 已存真实快照｜相同正交视角、盒体比例和原子半径；未插值补帧", font=subtitle_font, fill=MUTED)

    # 元素图例及原子计数。
    legend_x = 86
    for particle_type in (1, 2, 3, 4):
        color = tuple(round(c * 255) for c in TYPE_COLORS[particle_type])
        draw.ellipse((legend_x, 157, legend_x + 22, 179), fill=color, outline="#45515a")
        draw.text((legend_x + 31, 153), f"{TYPE_NAMES[particle_type]}  {type_counts[particle_type]:,}",
                  font=legend_font, fill=INK)
        legend_x += 255
    draw.text((1120, 153), "完整盒：200 × 200 × 250 Å；x/y 周期，z 自由表面", font=legend_font, fill=INK)
    draw.line((2740, 164, 2810, 164), fill="#1769aa", width=4)
    draw.text((2823, 153), "蓝虚线：局部裁切范围", font=legend_font, fill=INK)

    col_x = [70, 1370, 2670]
    row_label_y = [207, 1190]
    image_y = [250, 1233]
    panel_text = [f"({chr(97 + i)})  {t:.2f} ps" for i, t in enumerate(FRAME_TIMES_PS)]
    panel_text[0] += "  (RMS净位移 0.00 Å)"
    for i in range(1, len(FRAME_TIMES_PS)):
        panel_text[i] += f"  (RMS净位移 {rms_displacement[i]:.2f} Å)"
    panel_text.append(f"(f)  热源邻域局部放大｜1.60 ps｜{local_atom_count:,} atoms")
    for i, image in enumerate(panels):
        row, col = divmod(i, 3)
        x, y = col_x[col], image_y[row]
        draw.text((x + 10, row_label_y[row]), panel_text[i], font=panel_font, fill=INK)
        if i < 5:
            draw_projected_roi(image, roi, full_center, direction, 180.0)
        canvas.paste(image, (x, y))
        if i < 5:
            bar_px = round(50.0 * panel_size[1] / (2.0 * 180.0))
            draw_scale_bar(draw, (x + 45, y + panel_size[1] - 53), bar_px, "50 Å")
        else:
            bar_px = round(20.0 * panel_size[1] / (2.0 * 52.0))
            draw_scale_bar(draw, (x + 45, y + panel_size[1] - 53), bar_px, "20 Å")

    source_info = run_summary["source"]
    grid_shape = run_summary["electron_grid_shape"]
    grid_spacing = run_summary["electron_grid_spacing_A"]
    footer1 = (
        f"热源：Ag-only 参考 TTM；中心 ({source_info['center_xy_A'][0]:.0f},"
        f"{source_info['center_xy_A'][1]:.0f}) Å，R={source_info['R_A']:.0f} Å，"
        f"作用深度={source_info['skin_A']:.0f} Å，源面 z={source_info['surface_z_A']:.0f} Å；"
        f"网格 {grid_shape[0]}×{grid_shape[1]}×{grid_shape[2]}，"
        f"间距 {grid_spacing[0]:.0f}×{grid_spacing[1]:.0f}×{grid_spacing[2]:.0f} Å。"
    )
    draw.text((80, 2180), footer1, font=foot_font, fill=MUTED)
    draw.text((80, 2213), "开发版 MTP 与 Python 参考 TTM 结果，仅作探索性结构可视化；候选区不能据此认定为已验证熔池。",
              font=foot_font, fill="#8b4a2d")

    stem = out_dir / "Figure2_round7_fiveframe_structure_evolution"
    save_pil_figure(canvas, stem)
    return {
        "files": [stem.with_suffix(ext).name for ext in (".png", ".tif", ".pdf")],
        "resolution_px": list(canvas.size),
        "panel_resolution_px": list(panel_size),
        "frame_records": rendered_frame_info,
        "type_counts_each_frame": {str(k): v for k, v in type_counts.items()},
        "box_A": [200.0, 200.0, 250.0],
        "boundary": ["periodic", "periodic", "free"],
        "local_roi_A": roi,
        "local_roi_atoms_at_1p60ps": local_atom_count,
        "orthographic_scale_bar_lengths_A": {"full_model": 50.0, "local_view": 20.0},
        "scale_conversion": "bar length = physical length × panel pixel height / (2 × orthographic viewport fov)",
        "element_colors": {TYPE_NAMES[k]: [round(c, 3) for c in TYPE_COLORS[k]] for k in TYPE_COLORS},
    }


def extract_csv_and_positions(csv_path: Path, trajectory: Path) -> tuple[dict, np.ndarray, np.ndarray]:
    """将更正径向筛选 CSV 与 1.60 ps 原子坐标按 atom ID 对齐。"""
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    ids = np.asarray([int(row["atom_id"]) for row in rows], dtype=np.int64)
    q6 = np.asarray([float(row["mean_q6_1ps"]) for row in rows], dtype=float)
    nonfcc = np.asarray([float(row["ptm_nonfcc_fraction_1ps"]) for row in rows], dtype=float)
    candidate = np.asarray([int(row["joint_candidate"]) == 1 for row in rows], dtype=bool)

    pipeline = import_file(str(trajectory))
    data = pipeline.compute(len(FRAME_TIMES_PS) - 1)
    all_ids = np.asarray(data.particles["Particle Identifier"], dtype=np.int64)
    all_types = np.asarray(data.particles["Particle Type"], dtype=np.int32)
    all_positions = np.asarray(data.particles["Position"], dtype=float)
    order = np.argsort(all_ids)
    sorted_ids = all_ids[order]
    indices = np.searchsorted(sorted_ids, ids)
    if np.any(indices >= len(sorted_ids)) or not np.array_equal(sorted_ids[indices], ids):
        raise ValueError("筛选 CSV 中存在无法在 1.60 ps 快照匹配的 atom ID")
    selected_indices = order[indices]
    if not np.all(all_types[selected_indices] == 4):
        raise ValueError("径向核心筛选 CSV 中匹配到非 Ag 原子")
    core_positions = all_positions[selected_indices]
    ag_positions = all_positions[all_types == 4]
    return (
        {"rows": rows, "ids": ids, "candidate": candidate, "q6": q6, "nonfcc": nonfcc},
        core_positions,
        ag_positions,
    )


def calculate_coordination(
    core_positions: np.ndarray,
    ag_positions: np.ndarray,
    candidate: np.ndarray,
) -> tuple[np.ndarray, float]:
    """按 Ag-only 4 Å 邻居数计算配位，x/y 周期、z 开放。"""
    cutoff = 4.0
    cell_width = cutoff
    nxy = int(round(200.0 / cell_width))
    bins: dict[tuple[int, int, int], list[int]] = {}
    for atom_index, position in enumerate(ag_positions):
        key = (
            int(np.floor((position[0] % 200.0) / cell_width)) % nxy,
            int(np.floor((position[1] % 200.0) / cell_width)) % nxy,
            int(np.floor(position[2] / cell_width)),
        )
        bins.setdefault(key, []).append(atom_index)

    coordination = np.zeros(len(core_positions), dtype=np.int32)
    for center_index, center in enumerate(core_positions):
        cell_x = int(np.floor((center[0] % 200.0) / cell_width)) % nxy
        cell_y = int(np.floor((center[1] % 200.0) / cell_width)) % nxy
        cell_z = int(np.floor(center[2] / cell_width))
        count_with_self = 0
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    bucket = bins.get(((cell_x + dx) % nxy, (cell_y + dy) % nxy, cell_z + dz))
                    if not bucket:
                        continue
                    delta = ag_positions[bucket] - center
                    delta[:, 0] -= np.rint(delta[:, 0] / 200.0) * 200.0
                    delta[:, 1] -= np.rint(delta[:, 1] / 200.0) * 200.0
                    count_with_self += int(np.count_nonzero(np.einsum("ij,ij->i", delta, delta) <= cutoff**2))
        coordination[center_index] = count_with_self - 1
    candidate_mean = float(np.mean(coordination[candidate]))
    return coordination, candidate_mean


def read_lammps_atom_table(path: Path) -> tuple[list[str], np.ndarray]:
    """读取小型候选原子 dump 的 ATOMS 表，仅取已保存的数据。"""
    lines = path.read_text(encoding="utf-8").splitlines()
    header = next((i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS")), None)
    if header is None:
        raise ValueError(f"找不到 LAMMPS ATOMS 表：{path}")
    columns = lines[header].split()[2:]
    table = np.loadtxt(lines[header + 1 :], dtype=float, ndmin=2)
    if table.shape[1] != len(columns):
        raise ValueError(f"{path.name} 的列数与 ATOMS 表头不匹配")
    return columns, table


def recompute_connected_components(
    positions: np.ndarray,
    cutoff_A: float,
    box_xy_A: float = 200.0,
) -> list[list[int]]:
    """按 cutoff 重算候选原子连通分量；只对 x/y 使用最小镜像。"""
    nxy = int(math.floor(box_xy_A / cutoff_A))
    width_xy = box_xy_A / nxy
    bins: dict[tuple[int, int, int], list[int]] = {}
    for i, p in enumerate(positions):
        key = (int(p[0] // width_xy) % nxy,
               int(p[1] // width_xy) % nxy,
               int(math.floor(p[2] / cutoff_A)))
        bins.setdefault(key, []).append(i)

    parent = list(range(len(positions)))
    sizes = [1] * len(positions)

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri == rj:
            return
        if sizes[ri] < sizes[rj]:
            ri, rj = rj, ri
        parent[rj] = ri
        sizes[ri] += sizes[rj]

    for i, p in enumerate(positions):
        bx = int(p[0] // width_xy) % nxy
        by = int(p[1] // width_xy) % nxy
        bz = int(math.floor(p[2] / cutoff_A))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    candidates = bins.get(((bx + dx) % nxy, (by + dy) % nxy, bz + dz), ())
                    for j in candidates:
                        if j <= i:
                            continue
                        delta = positions[j] - p
                        delta[:2] -= np.rint(delta[:2] / box_xy_A) * box_xy_A
                        if float(np.dot(delta, delta)) <= cutoff_A**2 + 1e-12:
                            union(i, j)

    components: dict[int, list[int]] = {}
    for i in range(len(positions)):
        components.setdefault(find(i), []).append(i)
    return sorted(components.values(), key=len, reverse=True)


def style_3d_axis(ax, limits: tuple[tuple[float, float], ...]) -> None:
    """统一四个 3D 面板的空间范围、投影和坐标轴。"""
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_zlim(*limits[2])
    ax.set_box_aspect((limits[0][1] - limits[0][0], limits[1][1] - limits[1][0], limits[2][1] - limits[2][0]))
    ax.set_proj_type("ortho")
    ax.view_init(elev=24, azim=-48)
    ax.set_xlabel("x (Å)", labelpad=1)
    ax.set_ylabel("y (Å)", labelpad=1)
    ax.set_zlabel("z (Å)", labelpad=1)
    ax.tick_params(axis="both", which="major", labelsize=8, pad=0)
    ax.tick_params(axis="z", which="major", labelsize=8, pad=0)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor("#dce3e8")
    ax.grid(True, color="#dce3e8", linewidth=0.45, alpha=0.6)


def render_figure4(
    csv_path: Path,
    diagnostic_path: Path,
    candidates_dump: Path,
    largest_dump: Path,
    screen_path: Path,
    trajectory: Path,
    out_dir: Path,
) -> dict:
    """制作 Q6、PTM 非 FCC 持续比例、重算配位数和最大候选连通分量图。"""
    aligned, core_positions, ag_positions = extract_csv_and_positions(csv_path, trajectory)
    core_ids = aligned["ids"]
    candidate = aligned["candidate"]
    q6 = aligned["q6"]
    nonfcc = aligned["nonfcc"]

    with screen_path.open(encoding="utf-8") as handle:
        screen = json.load(handle)
    with diagnostic_path.open(encoding="utf-8") as handle:
        diagnostic = json.load(handle)
    source_expected = int(screen["source_core_ag_atoms"])
    candidate_expected = int(screen["joint_candidate_atoms"])
    cluster_summary = screen["cluster_summary"]
    ranked = screen["clusters_ranked"]
    q6_cutoff = float(screen["rule"]["q6"].split("<")[-1].split(";")[0].strip())
    nonfcc_cutoff = float(screen["rule"]["ptm"].split(">=")[-1].split(";")[0].strip())
    displacement_cutoff = float(screen["rule"]["msd_i"].split(">=")[-1].split()[0])
    connectivity_cutoff = float(cluster_summary["connectivity_cutoff_A"])

    if len(core_positions) != source_expected or len(np.unique(core_ids)) != source_expected:
        raise ValueError(f"source-core CSV 计数/ID 唯一性不符：{len(core_ids)} vs {source_expected}")
    if int(candidate.sum()) != candidate_expected:
        raise ValueError(f"CSV 候选标志计数不符：{candidate.sum()} vs {candidate_expected}")
    if not np.isfinite(q6).all() or not np.isfinite(nonfcc).all():
        raise ValueError("Q6 或 PTM 非 FCC 比例含非有限值")
    if np.any((q6 < 0) | (q6 > 1)) or np.any((nonfcc < 0) | (nonfcc > 1)):
        raise ValueError("Q6 / PTM 比例超出 [0,1] 范围")

    candidate_columns, candidate_table = read_lammps_atom_table(candidates_dump)
    largest_columns, largest_table = read_lammps_atom_table(largest_dump)
    candidate_dump_ids = candidate_table[:, candidate_columns.index("id")].astype(np.int64)
    largest_dump_ids = largest_table[:, largest_columns.index("id")].astype(np.int64)
    if len(np.unique(candidate_dump_ids)) != len(candidate_dump_ids):
        raise ValueError("联合候选 dump 中存在重复原子 ID")
    if set(candidate_dump_ids.tolist()) != set(core_ids[candidate].tolist()):
        raise ValueError("联合候选 dump 与 CSV 的候选原子 ID 集合不一致")
    if len(np.unique(largest_dump_ids)) != len(largest_dump_ids):
        raise ValueError("最大连通团簇 dump 中存在重复原子 ID")
    if set(largest_dump_ids.tolist()) != set(map(int, ranked[0]["atom_ids"])):
        raise ValueError("最大连通团簇 dump 与筛选 JSON 第一名 ID 不一致")
    if not set(largest_dump_ids.tolist()).issubset(set(candidate_dump_ids.tolist())):
        raise ValueError("最大连通团簇包含不在联合候选集中的原子")

    # 检查候选 dump 的坐标与真实 1.60 ps 全体系帧逐 ID 对齐。
    last = import_file(str(trajectory)).compute(len(FRAME_TIMES_PS) - 1)
    last_ids = np.asarray(last.particles["Particle Identifier"], dtype=np.int64)
    last_positions = np.asarray(last.particles["Position"], dtype=float)
    last_order = np.argsort(last_ids)
    sorted_last_ids = last_ids[last_order]
    sorted_last_positions = last_positions[last_order]
    candidate_indices = np.searchsorted(sorted_last_ids, candidate_dump_ids)
    dump_position_columns = [candidate_columns.index(c) for c in ("x", "y", "z")]
    candidate_dump_positions = candidate_table[:, dump_position_columns]
    delta = candidate_dump_positions - sorted_last_positions[candidate_indices]
    delta[:, :2] -= np.rint(delta[:, :2] / 200.0) * 200.0
    candidate_position_max_abs_error = float(np.max(np.abs(delta)))
    if candidate_position_max_abs_error > 1e-5:
        raise ValueError(f"候选 dump 坐标与真实轨迹不一致：max Δ={candidate_position_max_abs_error:.3g} Å")

    coordination, candidate_coordination_mean = calculate_coordination(core_positions, ag_positions, candidate)
    stored_mean = float(diagnostic["frames"][-1]["ag_coordination"]["joint_candidates"]["mean"])
    if abs(candidate_coordination_mean - stored_mean) > 0.05:
        raise ValueError(
            f"配位数复算与诊断均值偏差过大：复算={candidate_coordination_mean:.3f}, 报告={stored_mean:.3f}"
        )

    candidate_positions = core_positions[candidate]
    candidate_ids_ordered = core_ids[candidate]
    largest_mask = np.isin(candidate_ids_ordered, largest_dump_ids)
    if int(largest_mask.sum()) != len(largest_dump_ids):
        raise ValueError("最大连通分量 ID 无法完整映射回 CSV / 轨迹坐标")
    components = recompute_connected_components(candidate_positions, connectivity_cutoff)
    computed_largest_ids = set(candidate_ids_ordered[components[0]].tolist())
    if len(components) != int(cluster_summary["component_count"]):
        raise ValueError(f"连通分量数复算不符：{len(components)} vs {cluster_summary['component_count']}")
    if len(components[0]) != len(largest_dump_ids) or computed_largest_ids != set(largest_dump_ids.tolist()):
        raise ValueError("最大连通分量从轨迹独立复算后与归档成员不一致")

    plt.rcParams.update({
        "font.family": "Noto Sans CJK JP",
        "font.sans-serif": ["Noto Sans CJK JP", "Noto Sans CJK SC", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "font.size": 10,
    })
    fig = plt.figure(figsize=(14.6, 10.4), facecolor="white")
    fig.text(0.055, 0.963, "Ag 局部结构失序与迁移候选区", fontsize=21, fontweight="bold", color=INK)
    fig.text(0.057, 0.925,
             "更正径向源区 1,787 个 Ag｜Q6 与 PTM 为 0.60–1.60 ps 五帧统计｜坐标取 1.60 ps",
             fontsize=10.5, color=MUTED)
    grid = fig.add_gridspec(2, 2, left=0.055, right=0.95, bottom=0.18, top=0.88,
                            wspace=0.20, hspace=0.24)
    axes = [fig.add_subplot(grid[r, c], projection="3d") for r in range(2) for c in range(2)]

    # 所有子图使用同一真实空间范围，范围由筛选数据坐标加 5 Å 留白得出。
    lower = np.floor((core_positions.min(axis=0) - 5.0) / 10.0) * 10.0
    upper = np.ceil((core_positions.max(axis=0) + 5.0) / 10.0) * 10.0
    limits = tuple((float(lower[i]), float(upper[i])) for i in range(3))
    labels = ["(a) 五帧平均 Ag Q6", "(b) 五帧 PTM 非 FCC 比例",
              "(c) 1.60 ps Ag 配位数（4 Å）", "(d) 联合候选及最大连通分量"]
    for ax, label in zip(axes, labels):
        style_3d_axis(ax, limits)
        ax.set_title(label, loc="left", pad=4, fontsize=11, color=INK, fontweight="bold")

    x, y, z = core_positions.T
    scatter_a = axes[0].scatter(x, y, z, c=q6, cmap="viridis", vmin=0.0, vmax=1.0,
                                s=9, depthshade=False, linewidths=0)
    cb_a = fig.colorbar(scatter_a, ax=axes[0], shrink=0.62, pad=0.08, fraction=0.035)
    cb_a.set_label("Q6（五帧均值，无量纲）", fontsize=8)
    cb_a.ax.tick_params(labelsize=7)

    scatter_b = axes[1].scatter(x, y, z, c=nonfcc, cmap="magma", vmin=0.0, vmax=1.0,
                                s=9, depthshade=False, linewidths=0)
    cb_b = fig.colorbar(scatter_b, ax=axes[1], shrink=0.62, pad=0.08, fraction=0.035)
    cb_b.set_label("PTM 非 FCC 帧比例（0–1）", fontsize=8)
    cb_b.ax.tick_params(labelsize=7)

    scatter_c = axes[2].scatter(x, y, z, c=coordination, cmap="plasma",
                                vmin=0, vmax=max(14, int(coordination.max())),
                                s=9, depthshade=False, linewidths=0)
    cb_c = fig.colorbar(scatter_c, ax=axes[2], shrink=0.62, pad=0.08, fraction=0.035)
    cb_c.set_label("Ag 邻居数（4 Å 内）", fontsize=8)
    cb_c.ax.tick_params(labelsize=7)

    others = ~candidate
    nonlargest = candidate & ~np.isin(core_ids, largest_dump_ids)
    largest_core = candidate & np.isin(core_ids, largest_dump_ids)
    axes[3].scatter(x[others], y[others], z[others], color="#b8c5ce", s=8, alpha=0.42,
                    depthshade=False, linewidths=0, label=f"未通过至少一项：{int(others.sum()):,}")
    axes[3].scatter(x[nonlargest], y[nonlargest], z[nonlargest], color=ORANGE_LIGHT, s=12,
                    alpha=0.80, depthshade=False, linewidths=0,
                    label=f"其他候选分量：{int(nonlargest.sum()):,}")
    axes[3].scatter(x[largest_core], y[largest_core], z[largest_core], color=ORANGE, s=15,
                    alpha=0.98, depthshade=False, linewidths=0,
                    label=f"最大连通分量：{int(largest_core.sum()):,}")
    axes[3].legend(loc="upper left", fontsize=7.5, frameon=True, framealpha=0.94)

    fig.text(0.058, 0.115,
             f"筛选：Q6 < {q6_cutoff:.4f}；PTM 非 FCC 比例 ≥ {nonfcc_cutoff:.2f}；"
             f"1.00 ps 净位移平方 ≥ {displacement_cutoff:.4f} Å²。",
             fontsize=9.2, color=INK)
    fig.text(0.058, 0.074,
             f"连通定义：1.60 ps Ag–Ag 距离 ≤ {connectivity_cutoff:.1f} Å，x/y 周期、z 开放；"
             f"复算 {len(components)} 个分量，最大 {len(components[0])} 个原子。",
             fontsize=9.0, color=INK)
    fig.text(0.058, 0.035,
             f"探索性候选，不等同已确认液态/熔池；阈值来自 Round5，未针对 Round7 校准。"
             f" 候选平均 Ag 配位数 {candidate_coordination_mean:.2f}。",
             fontsize=8.8, color="#8b4a2d")

    stem = out_dir / "Figure4_round7_ag_disorder_candidate_region"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tif"), dpi=300, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    candidate_coordination = coordination[candidate]
    q6_anomalies = []
    for row_index in np.flatnonzero(q6 <= 1e-8):
        q6_anomalies.append({
            "atom_id": int(core_ids[row_index]),
            "q6": float(q6[row_index]),
            "ptm_nonfcc_fraction": float(nonfcc[row_index]),
            "ag_coordination_at_4A": int(coordination[row_index]),
            "joint_candidate": bool(candidate[row_index]),
        })
    return {
        "files": [stem.with_suffix(ext).name for ext in (".png", ".tif", ".pdf")],
        "source_core_atoms": int(len(core_positions)),
        "source_core_unique_ids": int(len(np.unique(core_ids))),
        "joint_candidate_atoms": int(candidate.sum()),
        "noncandidate_atoms": int((~candidate).sum()),
        "candidate_dump_ids_match_csv": True,
        "largest_component_dump_matches_json": True,
        "candidate_dump_position_max_abs_error_A": candidate_position_max_abs_error,
        "recomputed_component_count": len(components),
        "recomputed_largest_component_atoms": len(components[0]),
        "recomputed_largest_component_ids_match_dump": True,
        "candidate_atoms_outside_largest_component": int(nonlargest.sum()),
        "q6_range": [float(q6.min()), float(q6.max())],
        "q6_zero_or_near_zero_count_le_1e-8": int(np.count_nonzero(q6 <= 1e-8)),
        "q6_zero_or_near_zero_atoms": q6_anomalies,
        "ptm_nonfcc_fraction_range": [float(nonfcc.min()), float(nonfcc.max())],
        "coordination_cutoff_A": 4.0,
        "candidate_coordination_mean_recomputed": round(candidate_coordination_mean, 4),
        "candidate_coordination_mean_from_diagnostic": round(stored_mean, 4),
        "candidate_coordination_min_median_max": [int(candidate_coordination.min()),
                                                    float(np.median(candidate_coordination)),
                                                    int(candidate_coordination.max())],
        "candidate_coordination_counts_below_8_10_12": {
            str(cutoff): int(np.count_nonzero(candidate_coordination < cutoff)) for cutoff in (8, 10, 12)
        },
        "spatial_bounds_A": [[v[0], v[1]] for v in limits],
        "metric_window_ps": [0.6, 1.6],
        "positions_frame_ps": 1.6,
        "limitations": [
            "Q6/PTM/displacement thresholds are Round5-derived and not calibrated for Round7.",
            "The 874 joint candidates are not confirmed liquid atoms or a validated melt pool.",
            "PTM non-FCC includes Other; low coordination can confound classification.",
            "The stored displacement is a single-origin net displacement squared, not a time-averaged MSD.",
            "One joint candidate (atom ID 303221) has Q6=0 and Ag-only CN=0 at 4 Å, consistent with an isolated low-coordination atom; it is not evidence of a liquid cluster.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", type=Path, required=True, help="五帧全体系 dump.gz")
    parser.add_argument("--source-core-csv", type=Path, required=True, help="更正径向区逐 Ag CSV")
    parser.add_argument("--ptm-diagnostic", type=Path, required=True, help="Ag 配位诊断 JSON")
    parser.add_argument("--candidates-dump", type=Path, required=True, help="874 个联合候选原子 dump")
    parser.add_argument("--largest-dump", type=Path, required=True, help="最大候选连通分量 dump")
    parser.add_argument("--screen-json", type=Path, required=True, help="联合筛选及连通分量汇总 JSON")
    parser.add_argument("--run-summary", type=Path, required=True, help="Round7 1.60 ps 运行汇总 JSON")
    parser.add_argument("--trajectory-manifest", type=Path, required=True, help="五帧源数据 SHA256 清单")
    parser.add_argument("--analysis-manifest", type=Path, required=True, help="更正径向分析 SHA256 清单")
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # 对原始归档文件做 SHA256 校验，避免误用同名但不同内容的数据。
    def verify_manifest(manifest_path: Path, files: list[Path]) -> dict:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        entries = {entry["file"]: entry for entry in manifest}
        checked = {}
        for path in files:
            if path.name not in entries:
                raise ValueError(f"SHA 清单中没有文件记录：{path.name}")
            entry = entries[path.name]
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if path.stat().st_size != int(entry["bytes"]) or digest != entry["sha256"]:
                raise ValueError(f"源文件 SHA/字节数与归档清单不符：{path.name}")
            checked[path.name] = {"bytes": path.stat().st_size, "sha256": digest, "manifest_match": True}
        return checked

    trajectory_sources = verify_manifest(args.trajectory_manifest, [args.trajectory, args.run_summary])
    analysis_sources = verify_manifest(args.analysis_manifest, [
        args.source_core_csv, args.ptm_diagnostic, args.candidates_dump, args.largest_dump, args.screen_json
    ])
    run_summary = json.loads(args.run_summary.read_text(encoding="utf-8"))
    figure2 = render_figure2(args.trajectory, args.out_dir, run_summary)
    figure4 = render_figure4(
        args.source_core_csv, args.ptm_diagnostic, args.candidates_dump, args.largest_dump,
        args.screen_json, args.trajectory, args.out_dir
    )

    checkpoints = run_summary.get("checkpoints", [])
    force_points = [(float(x["max_force_eV_A"]), float(x["time_ps"]))
                    for x in checkpoints if x.get("max_force_eV_A") is not None]
    residual_points = [(float(x["combined_source_corrected_residual_eV"]), float(x["time_ps"]))
                       for x in checkpoints if x.get("combined_source_corrected_residual_eV") is not None]
    max_force, max_force_time = max(force_points, key=lambda item: item[0]) if force_points else (None, None)
    max_residual_abs, max_residual_time = max(residual_points, key=lambda item: abs(item[0])) if residual_points else (None, None)
    quality_report = {
        "status": "ROUND7_STAGE1_FIGURES_2_AND_4_COMPLETE",
        "branch_review": {
            "main": {"commit": "026a1411ebaab92dabfb9ac6df51060907db3032",
                     "finding": "包含上一版图件；指定完整五帧原始轨迹不在此分支。"},
            "codex/round7-main-pool-figure-20261121": {
                "commit": "f64a0b48953a04162eba1afd703d5be3ed96bddb",
                "finding": "包含 0.380 ps、458,248 原子候选图与说明；README 明确未包含原始 dump/运行日志。"
            },
            "codex/round7-v12m-cold-diagnostics-20261119": {
                "commit": "df890b0bf62f6b4e4c60b9f6a14b69276873dd1d",
                "finding": "包含 v12m 等冷态/势函数诊断归档，不含本次所需五帧 TTM-MD 全体系轨迹。"
            },
            "codex/upload-round7-fiveframe-trajectory-1p6ps-20261109": {
                "commit": "878328378869ba4e5860a4081420fcdc59966c52",
                "finding": "包含本次使用的 0.60–1.60 ps 五帧全体系 dump、运行汇总及原始五帧分析。"
            },
            "figure4_corrected_analysis": {
                "branch": "codex/upload-round7-fiveframe-trajectory-1p6ps-20261109",
                "commit": "878328378869ba4e5860a4081420fcdc59966c52",
                "path": "models/ag_ti3sic2_qsgs_200x200x150_md_finer_skeleton/runs/round7_corrected_radial_core_analysis_1p6ps_20261109/",
                "finding": "采用更正径向掩膜版本，避免旧版 6,448 原子 x-z 带状掩膜遗漏 y 条件的问题。"
            },
        },
        "verified_source_hashes": {**trajectory_sources, **analysis_sources},
        "run_diagnostics": {
            "atoms": int(run_summary["atoms"]),
            "duration_ps": float(run_summary["duration_ps"]),
            "dt_fs": float(run_summary["dt_fs"]),
            "potential_classification": run_summary.get("potential_classification"),
            "max_checkpoint_force_eV_A": max_force,
            "max_checkpoint_force_time_ps": max_force_time,
            "max_abs_combined_source_corrected_residual_eV": abs(max_residual_abs) if max_residual_abs is not None else None,
            "signed_combined_source_corrected_residual_at_max_abs_eV": max_residual_abs,
            "time_of_max_abs_residual_ps": max_residual_time,
            "final_source_corrected_combined_energy_residual_eV": float(
                run_summary["source_corrected_combined_energy_residual_eV"]
            ),
            "final_cumulative_ttm_grid_residual_eV": float(run_summary["cumulative_ttm_grid_residual_eV"]),
            "limitations": run_summary.get("limitations", []),
        },
        "figure2": figure2,
        "figure4": figure4,
        "scientific_scope": (
            "只对归档轨迹进行后处理，不运行新 MD。Round7 开发版 MTP 与 Ag-only Python 参考 TTM 仍属探索性；"
            "874 个筛选原子不等同于已确认液态/熔池，图件不能单独证明熔焊或界面传热机制。"
        ),
    }
    report_path = args.out_dir / "round7_stage1_quality_report.json"
    report_path.write_text(json.dumps(quality_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(quality_report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
