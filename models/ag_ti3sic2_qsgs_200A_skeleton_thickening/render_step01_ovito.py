#!/usr/bin/env python3
"""用 OVITO StandardRenderer 并排渲染原始与第一档加粗原子模型。"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from PySide6.QtCore import Qt
from ovito.io import import_file
from ovito.modifiers import ColorCodingModifier
from ovito.vis import ParticlesVis, StandardRenderer, Viewport

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASE_MODEL = REPO / "models/ag_ti3sic2_qsgs_200A_ttm_finer/ag_ti3sic2_qsgs.extxyz"
NEW_MODEL = HERE / "output/ag_ti3sic2_skeleton_thickened_step01.extxyz"
FIG = HERE / "visualization"
PNG = FIG / "Ag_Ti3SiC2_skeleton_thickened_step01_comparison.png"
PDF = FIG / "Ag_Ti3SiC2_skeleton_thickened_step01_comparison.pdf"
TIFF = FIG / "Ag_Ti3SiC2_skeleton_thickened_step01_comparison.tif"
REPORT = FIG / "render_report.json"
TITLE = "Ag–Ti3SiC2 原子模型：骨架加粗初版"
SUBTITLE = "同一纳米模型尺寸与模拟盒；比较 Ti3SiC2 骨架第一档加粗前后"
LEFT_LABEL = "原始模型"
RIGHT_LABEL = "骨架外扩 1 Å/侧"

COLORS = {
    1: np.asarray((0.055, 0.078, 0.102)),  # Ti3SiC2 骨架：深灰蓝
    2: np.asarray((0.663, 0.831, 0.910)),  # Ag：浅蓝银色
}
CANVAS = (3500, 1900)
PANEL = (1660, 1540)


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc" if bold else "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def render_model(model: Path, destination: Path) -> tuple[int, dict[int, int]]:
    pipeline = import_file(str(model))
    data = pipeline.compute()
    phase = np.asarray(data.particles["phase"], dtype=np.int32)
    counts = {int(p): int(np.count_nonzero(phase == p)) for p in np.unique(phase)}
    if set(counts) != {1, 2}:
        raise ValueError(f"模型相标签异常：{counts}")

    color = ColorCodingModifier(operate_on="particles", property="phase")
    color.start_value = 1.0
    color.end_value = 2.0
    color.gradient = ColorCodingModifier.Gradient(np.asarray((COLORS[1], COLORS[2])))
    pipeline.modifiers.append(color)
    particles = data.particles.vis
    particles.radius = 0.95
    particles.shape = ParticlesVis.Shape.Sphere
    cell_vis = data.cell.vis
    cell_vis.enabled = True
    cell_vis.render_cell = True
    cell_vis.line_width = 1.25
    cell_vis.rendering_color = (0.035, 0.045, 0.055)

    pipeline.add_to_scene()
    viewport = Viewport(type=Viewport.Type.Perspective)
    viewport.camera_pos = (690.0, -760.0, 755.0)
    viewport.camera_dir = (-590.0, 860.0, -605.0)
    viewport.camera_up = (0.0, 0.0, 1.0)
    viewport.fov = 38.0
    viewport.zoom_all(size=PANEL)
    renderer = StandardRenderer()
    renderer.antialiasing_level = 4
    renderer.ambient_occlusion = True
    renderer.ambient_occlusion_intensity = 0.62
    renderer.ambient_occlusion_radius = 5.0
    renderer.shadows = True
    renderer.shadow_intensity = 0.22
    renderer.order_independent_transparency = True
    image = viewport.render_image(size=PANEL, background=(1.0, 1.0, 1.0), alpha=False,
                                  renderer=renderer, stop_on_error=True)
    image.save(str(destination))
    pipeline.remove_from_scene()
    return int(len(phase)), counts


def main() -> None:
    if not BASE_MODEL.exists() or not NEW_MODEL.exists():
        raise FileNotFoundError("基准或加粗原子模型文件不存在")
    FIG.mkdir(parents=True, exist_ok=True)
    left_path = FIG / "baseline_render.png"
    right_path = FIG / "step01_render.png"
    base_n, base_counts = render_model(BASE_MODEL, left_path)
    new_n, new_counts = render_model(NEW_MODEL, right_path)

    canvas = Image.new("RGB", CANVAS, "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((90, 55), TITLE, font=font(52, True), fill="#1e2d39")
    draw.text((94, 126), SUBTITLE, font=font(28), fill="#607382")
    draw.text((105, 215), LEFT_LABEL, font=font(32, True), fill="#273846")
    draw.text((1720, 215), RIGHT_LABEL, font=font(32, True), fill="#273846")
    draw.text((105, 270), f"Ti3SiC2 {base_counts[1]:,}  ·  Ag {base_counts[2]:,}  ·  总计 {base_n:,}", font=font(23), fill="#526574")
    draw.text((1720, 270), f"Ti3SiC2 {new_counts[1]:,}  ·  Ag {new_counts[2]:,}  ·  总计 {new_n:,}", font=font(23), fill="#526574")
    with Image.open(left_path) as left, Image.open(right_path) as right:
        canvas.paste(left.convert("RGB"), (70, 330))
        canvas.paste(right.convert("RGB"), (1770, 330))
    draw.rectangle((70, 330, 70 + PANEL[0], 330 + PANEL[1]), outline="#b8c2c9", width=2)
    draw.rectangle((1770, 330, 1770 + PANEL[0], 330 + PANEL[1]), outline="#b8c2c9", width=2)
    # 相图例：标签与实际 OVITO 相颜色一致。
    legend_y = 1815
    draw.ellipse((105, legend_y, 135, legend_y + 30), fill="#0e141a")
    draw.text((148, legend_y - 2), "Ti3SiC2 骨架", font=font(22), fill="#354654")
    draw.ellipse((365, legend_y, 395, legend_y + 30), fill="#a9d4e8")
    draw.text((408, legend_y - 2), "Ag 银相", font=font(22), fill="#354654")
    draw.text((950, legend_y - 1), "材料区 200 × 200 × 200 Å（20 nm）  ·  模拟盒 200 × 200 × 300 Å（30 nm）",
              font=font(22), fill="#526574")
    canvas.save(PNG, dpi=(300, 300))
    canvas.save(TIFF, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
    canvas.save(PDF, format="PDF", resolution=300.0)
    REPORT.write_text(json.dumps({
        "status": "OVITO_STANDARD_RENDERER_COMPARISON_COMPLETE",
        "renderer": "OVITO StandardRenderer",
        "resolution_px": list(CANVAS),
        "baseline_atoms": base_n,
        "baseline_phase_counts": {"Ti3SiC2": base_counts[1], "Ag": base_counts[2]},
        "step01_atoms": new_n,
        "step01_phase_counts": {"Ti3SiC2": new_counts[1], "Ag": new_counts[2]},
        "outputs": [PNG.name, PDF.name, TIFF.name],
    }, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"png": str(PNG), "baseline_atoms": base_n, "step01_atoms": new_n,
                      "phase_counts": {"baseline": base_counts, "step01": new_counts}},
                     ensure_ascii=False, indent=2))
    left_path.unlink(missing_ok=True)
    right_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
