#!/usr/bin/env python3
"""Render Round7 paper figures 2 and 4 from the archived trajectory.

No molecular-dynamics steps are run here. Figure 2 uses the five stored
full-system snapshots. Figure 4 combines the corrected radial-core CSV with
the 1.60 ps coordinates and recomputes Ag-only coordination at 4 Å.
"""
from __future__ import annotations

import argparse
import csv
import json
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
from ovito.vis import ParticlesVis, StandardRenderer, Viewport


FRAME_TIMES_PS = [0.60, 0.85, 1.10, 1.35, 1.60]
TYPE_NAMES = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}
FRAMEWORK_RGB = (0.075, 0.105, 0.145)
AG_RGB = (0.60, 0.79, 0.89)
ORANGE = "#e87532"
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


def render_figure2(trajectory: Path, out_dir: Path) -> dict:
    """用同一 OVITO 相机和显示样式渲染五个全体系快照。"""
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
    first_positions = np.asarray(frame_data[0].particles["Position"], dtype=float)[first_order]
    rms_displacement = []
    for data in frame_data:
        frame_ids = np.asarray(data.particles["Particle Identifier"], dtype=np.int64)
        frame_order = np.argsort(frame_ids)
        if not np.array_equal(frame_ids[frame_order], sorted_ids):
            raise ValueError("五帧中的原子 ID 集合不一致")
        positions = np.asarray(data.particles["Position"], dtype=float)[frame_order]
        displacement = positions - first_positions
        displacement[:, :2] -= np.rint(displacement[:, :2] / 200.0) * 200.0
        rms_displacement.append(float(np.sqrt(np.mean(np.einsum("ij,ij->i", displacement, displacement)))))

    first_types = np.asarray(frame_data[0].particles["Particle Type"], dtype=int)
    type_counts = {int(t): int(np.count_nonzero(first_types == t)) for t in np.unique(first_types)}
    if set(type_counts) != {1, 2, 3, 4}:
        raise ValueError(f"原子类型与 Round7 预期不符：{type_counts}")
    framework_count = sum(type_counts[t] for t in (1, 2, 3))
    ag_count = type_counts[4]

    # 给 OVITO 源数据中的类型赋予清晰的元素名和论文配色。
    source = pipeline.source.data
    for particle_type in source.particles_.particle_types_.types_:
        particle_type.name = TYPE_NAMES[int(particle_type.id)]
        particle_type.color = FRAMEWORK_RGB if particle_type.id in (1, 2, 3) else AG_RGB
    source.particles_.vis.radius = 0.8
    source.particles_.vis.shape = ParticlesVis.Shape.Sphere
    source.cell_.vis.enabled = True
    source.cell_.vis.render_cell = True
    source.cell_.vis.line_width = 0.8
    source.cell_.vis.rendering_color = (0.14, 0.18, 0.22)

    pipeline.add_to_scene()
    viewport = Viewport(type=Viewport.Type.Perspective)
    scene.viewports.active_vp = viewport
    center = np.array((100.0, 100.0, 125.0))
    direction = np.array((-350.0, 640.0, -325.0))
    direction /= np.linalg.norm(direction)
    viewport.camera_pos = center - direction * 820.0
    viewport.camera_dir = direction
    viewport.camera_up = (0.0, 0.0, 1.0)
    # 窄视场使 200 × 200 × 250 Å 盒体在每个面板内清晰可辨。
    viewport.fov = 12.0

    renderer = StandardRenderer()
    renderer.antialiasing_level = 2
    renderer.ambient_occlusion = True
    renderer.ambient_occlusion_intensity = 0.55
    renderer.ambient_occlusion_radius = 5.0
    renderer.shadows = True
    renderer.shadow_intensity = 0.22

    panel_size = (1240, 930)
    panels: list[Image.Image] = []
    rendered_frame_info = []
    temporary_panel_dir = out_dir / ".round7_render_panels"
    temporary_panel_dir.mkdir(exist_ok=True)
    for frame_index, time_ps in enumerate(FRAME_TIMES_PS):
        scene.anim.current_frame = frame_index
        rendered = viewport.render_image(
            size=panel_size,
            background=(1.0, 1.0, 1.0),
            renderer=renderer,
            stop_on_error=True,
        )
        panel_path = temporary_panel_dir / f"panel_{frame_index:02d}.png"
        rendered.save(str(panel_path), "PNG")
        with Image.open(panel_path) as panel_image:
            panels.append(panel_image.convert("RGB").copy())
        panel_path.unlink(missing_ok=True)
        rendered_frame_info.append({
            "frame": frame_index,
            "time_ps": time_ps,
            "timestep": int(frame_data[frame_index].attributes["Timestep"]),
            "atoms": atom_counts[frame_index],
            "rms_displacement_from_0p60ps_A": round(rms_displacement[frame_index], 6),
        })
    temporary_panel_dir.rmdir()

    # 采用 3 × 2 面板；最后一格放模型参数、颜色图例和数据限制。
    canvas = Image.new("RGB", (3980, 2390), "white")
    draw = ImageDraw.Draw(canvas)
    title_font = load_font(66)
    subtitle_font = load_font(31)
    panel_font = load_font(34)
    body_font = load_font(27)
    small_font = load_font(23)
    draw.text((90, 35), "局域热输入下 Ag/Ti3SiC2 原子结构演化", font=title_font, fill=INK, stroke_width=1)
    draw.text(
        (96, 122),
        "Round7 已存全体系轨迹｜统一观察方向、相颜色与原子显示尺度",
        font=subtitle_font,
        fill=MUTED,
    )

    column_x = [80, 1370, 2660]
    row_y = [220, 1280]
    panel_labels = [
        f"({chr(ord('a') + i)})  {FRAME_TIMES_PS[i]:.2f} ps   RMS 位移 {rms_displacement[i]:.2f} Å"
        for i in range(len(FRAME_TIMES_PS))
    ]
    for i, image in enumerate(panels):
        row, col = divmod(i, 3)
        x, y = column_x[col], row_y[row]
        draw.text((x + 12, y), panel_labels[i], font=panel_font, fill=INK)
        canvas.paste(image, (x, y + 58))

    card_x, card_y = column_x[2], row_y[1] + 58
    card_w, card_h = panel_size
    draw.rounded_rectangle(
        (card_x, card_y, card_x + card_w, card_y + card_h),
        radius=24,
        fill="#f3f6f8",
        outline="#d5dfe6",
        width=3,
    )
    draw.text((card_x + 48, card_y + 42), "模型与热源设置", font=panel_font, fill=INK, stroke_width=1)
    draw.rectangle((card_x + 52, card_y + 120, card_x + 88, card_y + 156), fill="#141d29")
    draw.text((card_x + 108, card_y + 116), f"Ti3SiC2 骨架：{framework_count:,} atoms", font=body_font, fill=INK)
    draw.rectangle((card_x + 52, card_y + 180, card_x + 88, card_y + 216), fill="#99c9e3")
    draw.text((card_x + 108, card_y + 176), f"Ag 银相：{ag_count:,} atoms", font=body_font, fill=INK)
    draw.text((card_x + 52, card_y + 264), f"全体系：{atom_counts[0]:,} atoms", font=body_font, fill=INK)
    draw.text((card_x + 52, card_y + 326), "模拟盒：200 × 200 × 250 Å", font=body_font, fill=INK)
    draw.text((card_x + 52, card_y + 374), "边界：x/y 周期；z 自由表面", font=body_font, fill=INK)
    draw.text((card_x + 52, card_y + 458), "Ag-only TTM 热源（探索性参考耦合）", font=body_font, fill=INK)
    draw.text((card_x + 52, card_y + 512), "I0 = 5.0 eV/(ps·Å^2)；R = 30 Å", font=small_font, fill=INK)
    draw.text((card_x + 52, card_y + 554), "表面作用深度 20 Å；源面 z = 150 Å", font=small_font, fill=INK)
    draw.text((card_x + 52, card_y + 596), "中心 (x,y) = (100,100) Å", font=small_font, fill=INK)
    draw.text((card_x + 52, card_y + 672), "TTM 网格：20 × 20 × 15；间距 10 Å", font=small_font, fill=INK)
    draw.text((card_x + 52, card_y + 736), "时间帧间隔 0.25 ps；无独立 1.50 ps 帧", font=small_font, fill=INK)
    draw.text(
        (card_x + 52, card_y + 812),
        "开发版 MTP 与 Python 参考 TTM；图示为探索性结果，",
        font=small_font,
        fill="#8b4a2d",
    )
    draw.text(
        (card_x + 52, card_y + 850),
        "不能单独作为熔化或熔焊机制已验证的证据。",
        font=small_font,
        fill="#8b4a2d",
    )
    draw.text(
        (92, 2312),
        "RMS 位移按 x/y 周期最小镜像、z 开放计算；结构来自已保存的五帧，未补帧或重跑 MD。",
        font=small_font,
        fill=MUTED,
    )

    stem = out_dir / "Figure2_round7_fiveframe_structure_evolution"
    save_pil_figure(canvas, stem)
    return {
        "files": [stem.with_suffix(ext).name for ext in (".png", ".tif", ".pdf")],
        "resolution_px": list(canvas.size),
        "panel_resolution_px": list(panel_size),
        "frame_records": rendered_frame_info,
        "type_counts_first_frame": {str(k): v for k, v in type_counts.items()},
        "phase_counts_first_frame": {
            "Ti3SiC2_framework_types_1_to_3": framework_count,
            "Ag_type_4": ag_count,
        },
        "box_A": [200.0, 200.0, 250.0],
        "boundary": ["periodic", "periodic", "free"],
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


def render_figure4(csv_path: Path, diagnostic_path: Path, trajectory: Path, out_dir: Path) -> dict:
    """制作 Q6、PTM、配位数和联合候选区的四面板空间图。"""
    aligned, core_positions, ag_positions = extract_csv_and_positions(csv_path, trajectory)
    candidate = aligned["candidate"]
    q6 = aligned["q6"]
    nonfcc = aligned["nonfcc"]
    coordination, candidate_coordination_mean = calculate_coordination(core_positions, ag_positions, candidate)
    if len(core_positions) != 1787 or int(candidate.sum()) != 874:
        raise ValueError(f"更正径向候选数据计数不符：core={len(core_positions)}, candidate={candidate.sum()}")

    with diagnostic_path.open(encoding="utf-8") as handle:
        diagnostic = json.load(handle)
    stored_mean = float(diagnostic["frames"][-1]["ag_coordination"]["joint_candidates"]["mean"])
    if abs(candidate_coordination_mean - stored_mean) > 0.05:
        raise ValueError(
            f"配位数复算与 Round7 诊断均值偏差过大：复算={candidate_coordination_mean:.3f}, "
            f"报告={stored_mean:.3f}"
        )

    plt.rcParams.update({
        "font.family": "Noto Sans CJK JP",
        "font.sans-serif": ["Noto Sans CJK JP", "Noto Sans CJK SC", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "font.size": 10,
    })
    fig = plt.figure(figsize=(14.4, 10.2), facecolor="white")
    fig.text(0.055, 0.963, "Ag 局部结构失序与迁移候选区域", fontsize=21, fontweight="bold", color=INK)
    fig.text(
        0.057,
        0.925,
        "Round7 更正径向源区：1,787 个 Ag｜坐标取 1.60 ps；Q6/PTM 为 0.60–1.60 ps 五帧统计",
        fontsize=10.5,
        color=MUTED,
    )
    grid = fig.add_gridspec(2, 2, left=0.055, right=0.95, bottom=0.17, top=0.88, wspace=0.20, hspace=0.24)
    axes = [
        fig.add_subplot(grid[0, 0], projection="3d"),
        fig.add_subplot(grid[0, 1], projection="3d"),
        fig.add_subplot(grid[1, 0], projection="3d"),
        fig.add_subplot(grid[1, 1], projection="3d"),
    ]
    limits = ((50.0, 150.0), (50.0, 150.0), (100.0, 165.0))
    labels = [
        "(a) 五帧平均 Ag Q6",
        "(b) 五帧 PTM 非 FCC 比例",
        "(c) 1.60 ps Ag 配位数（4 Å）",
        "(d) 联合筛选候选原子",
    ]
    for ax, label in zip(axes, labels):
        style_3d_axis(ax, limits)
        ax.set_title(label, loc="left", pad=4, fontsize=11, color=INK, fontweight="bold")

    x, y, z = core_positions.T
    scatter_a = axes[0].scatter(x, y, z, c=q6, cmap="viridis", vmin=0.0, vmax=1.0,
                                s=9, depthshade=False, linewidths=0)
    cb_a = fig.colorbar(scatter_a, ax=axes[0], shrink=0.62, pad=0.08, fraction=0.035)
    cb_a.set_label("mean Q6", fontsize=8)
    cb_a.ax.tick_params(labelsize=7)

    scatter_b = axes[1].scatter(x, y, z, c=nonfcc, cmap="magma", vmin=0.0, vmax=1.0,
                                s=9, depthshade=False, linewidths=0)
    cb_b = fig.colorbar(scatter_b, ax=axes[1], shrink=0.62, pad=0.08, fraction=0.035)
    cb_b.set_label("非 FCC 帧比例", fontsize=8)
    cb_b.ax.tick_params(labelsize=7)

    scatter_c = axes[2].scatter(x, y, z, c=coordination, cmap="plasma",
                                vmin=0, vmax=max(14, int(coordination.max())),
                                s=9, depthshade=False, linewidths=0)
    cb_c = fig.colorbar(scatter_c, ax=axes[2], shrink=0.62, pad=0.08, fraction=0.035)
    cb_c.set_label("Ag 邻居数", fontsize=8)
    cb_c.ax.tick_params(labelsize=7)

    axes[3].scatter(x[~candidate], y[~candidate], z[~candidate], color="#a9bbc8",
                    s=8, alpha=0.48, depthshade=False, linewidths=0, label=f"未通过至少一项：{int((~candidate).sum()):,}")
    axes[3].scatter(x[candidate], y[candidate], z[candidate], color=ORANGE,
                    s=11, alpha=0.95, depthshade=False, linewidths=0, label=f"联合候选：{int(candidate.sum()):,}")
    axes[3].legend(loc="upper left", fontsize=8, frameon=True, framealpha=0.93)

    fig.text(
        0.058,
        0.108,
        "筛选条件：mean Q6 < 0.4081；PTM 非 FCC 比例 ≥ 0.80；1 ps 漂移修正净位移平方 ≥ 0.5172 Å²。",
        fontsize=9.5,
        color=INK,
    )
    fig.text(
        0.058,
        0.068,
        f"探索性候选，不等同于已确认液态/熔池；阈值来自 Round5，未用 Round7 专门校准。"
        f"  874 个候选分成 41 个连通分量，最大 693 个；候选区平均配位数 {candidate_coordination_mean:.2f}。",
        fontsize=9.1,
        color="#8b4a2d",
    )

    stem = out_dir / "Figure4_round7_ag_disorder_candidate_region"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(
        stem.with_suffix(".tif"),
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    return {
        "files": [stem.with_suffix(ext).name for ext in (".png", ".tif", ".pdf")],
        "subset_atoms": int(len(core_positions)),
        "joint_candidate_atoms": int(candidate.sum()),
        "other_subset_atoms": int((~candidate).sum()),
        "components": 41,
        "largest_component_atoms": 693,
        "q6_cutoff": 0.408106928,
        "ptm_nonfcc_fraction_cutoff": 0.8,
        "displacement_threshold_A2": 0.517167930,
        "coordination_cutoff_A": 4.0,
        "connectivity_cutoff_A": 3.5,
        "candidate_coordination_mean_recomputed": round(candidate_coordination_mean, 4),
        "candidate_coordination_mean_from_diagnostic": round(stored_mean, 4),
        "positions_frame_ps": 1.6,
        "metric_window_ps": [0.6, 1.6],
        "limitations": [
            "Round5-derived screen thresholds are not calibrated for Round7.",
            "Joint candidates are not confirmed liquid atoms or a validated melt pool.",
            "PTM non-FCC includes Other; low coordination can confound classification.",
            "The stored displacement statistic is a single-origin drift-corrected net displacement squared, not a time-averaged MSD.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", type=Path, required=True, help="Round7 five-frame full-system .dump.gz")
    parser.add_argument("--source-core-csv", type=Path, required=True, help="corrected radial source-core CSV")
    parser.add_argument("--ptm-diagnostic", type=Path, required=True, help="framewise PTM/coordination diagnostic JSON")
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--data-branch-commit", default="878328378869ba4e5860a4081420fcdc59966c52")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    figure2 = render_figure2(args.trajectory, args.out_dir)
    figure4 = render_figure4(args.source_core_csv, args.ptm_diagnostic, args.trajectory, args.out_dir)
    import hashlib
    report = {
        "status": "ROUND7_FIGURES_2_AND_4_RENDERED",
        "data_branch": "codex/upload-round7-fiveframe-trajectory-1p6ps-20261109",
        "data_branch_commit": args.data_branch_commit,
        "trajectory": str(args.trajectory),
        "trajectory_frames_ps": FRAME_TIMES_PS,
        "trajectory_sha256": hashlib.sha256(args.trajectory.read_bytes()).hexdigest(),
        "figure2": figure2,
        "figure4": figure4,
        "scientific_scope": (
            "These are post-processing figures from the archived Round7 development MTP and Ag-only "
            "Python reference TTM run. No new MD was run. They do not establish calibrated melting, "
            "welding, or interfacial heat-transfer behavior."
        ),
    }
    report_path = args.out_dir / "round7_figure_generation_report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
