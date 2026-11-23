#!/usr/bin/env python3
"""使用 OVITO 为十倍枝间距的 Ag/Ti3SiC2 原子模型生成三维图。"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PIL import Image
from PySide6.QtCore import Qt
from ovito.io import import_file
from ovito.modifiers import ColorCodingModifier
from ovito.vis import ParticlesVis, StandardRenderer, TextLabelOverlay, Viewport

HERE = Path(__file__).resolve().parent
MODEL = HERE / "output/ag_ti3sic2_spacing10x.extxyz"
FIG = HERE / "visualization"
PNG = FIG / "Ag_Ti3SiC2_spacing10x_atomic_OVITO.png"
PDF = FIG / "Ag_Ti3SiC2_spacing10x_atomic_OVITO.pdf"
TIFF = FIG / "Ag_Ti3SiC2_spacing10x_atomic_OVITO.tif"
REPORT = FIG / "Ag_Ti3SiC2_spacing10x_atomic_OVITO.json"

COLORS = {
    1: np.asarray((0.055, 0.078, 0.102)),  # Ti3SiC2 骨架：深灰
    2: np.asarray((0.663, 0.831, 0.910)),  # Ag：浅蓝银色
}


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    pipeline = import_file(str(MODEL))
    data = pipeline.source.data
    phase = np.asarray(data.particles["phase"], dtype=np.int32)
    counts = {int(p): int(np.count_nonzero(phase == p)) for p in np.unique(phase)}
    if set(counts) != {1, 2}:
        raise ValueError(f"相标签异常：{counts}")

    color = ColorCodingModifier(operate_on="particles", property="phase")
    color.start_value = 1.0
    color.end_value = 2.0
    color.gradient = ColorCodingModifier.Gradient(np.asarray((COLORS[1], COLORS[2])))
    color.title = "Phase"
    pipeline.modifiers.append(color)

    particle_vis = data.particles.vis
    particle_vis.radius = 0.95
    particle_vis.shape = ParticlesVis.Shape.Sphere
    cell_vis = data.cell.vis
    cell_vis.enabled = True
    cell_vis.render_cell = True
    cell_vis.line_width = 1.0
    cell_vis.rendering_color = (0.035, 0.045, 0.055)

    label = TextLabelOverlay()
    label.text = (
        "<div style='background:#ffffff; border:1px solid #73818b; padding:8px 11px;'>"
        "<b>Ag–Ti₃SiC₂ 原子模型（枝间距约放大 10 倍，未加热）</b><br/>"
        f"<span style='color:#0E141A;'>■</span> Ti₃SiC₂：{counts[1]:,}<br/>"
        f"<span style='color:#A9D4E8;'>■</span> Ag：{counts[2]:,}<br/>"
        f"总原子数：{len(phase):,}<br/>"
        "材料区：200 × 200 × 200 Å（20 nm）<br/>"
        "模拟盒：200 × 200 × 300 Å（20 × 20 × 30 nm）"
        "</div>"
    )
    label.font_size = 0.024
    label.font_family = "Noto Sans CJK SC"
    label.text_color = (0.08, 0.10, 0.12)
    label.outline_enabled = False
    label.alignment = Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop
    label.offset_x = 0.025
    label.offset_y = -0.025

    pipeline.add_to_scene()
    viewport = Viewport(type=Viewport.Type.Perspective)
    viewport.overlays.append(label)
    viewport.camera_pos = (690.0, -760.0, 755.0)
    viewport.camera_dir = (-590.0, 860.0, -605.0)
    viewport.camera_up = (0.0, 0.0, 1.0)
    viewport.fov = 38.0
    viewport.zoom_all(size=(2200, 1650))

    renderer = StandardRenderer()
    renderer.antialiasing_level = 4
    renderer.ambient_occlusion = True
    renderer.ambient_occlusion_intensity = 0.6
    renderer.ambient_occlusion_radius = 5.0
    renderer.shadows = True
    renderer.shadow_intensity = 0.22
    renderer.order_independent_transparency = True

    image = viewport.render_image(
        size=(2200, 1650), filename=str(PNG), background=(1.0, 1.0, 1.0),
        alpha=False, renderer=renderer, stop_on_error=True,
    )
    image.save(str(PNG), "PNG")
    with Image.open(PNG) as rendered:
        rendered.save(TIFF, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
        rendered.save(PDF, format="PDF", resolution=300.0)

    report = {
        "status": "OVITO_ATOMISTIC_MODEL_RENDERED",
        "model": str(MODEL.relative_to(HERE.parent)),
        "atoms": int(len(phase)),
        "phase_counts": {"Ti3SiC2": counts[1], "Ag": counts[2]},
        "material_region_A": [200, 200, 200],
        "simulation_box_A": [200, 200, 300],
        "pbc": [True, True, False],
        "renderer": "OVITO StandardRenderer",
        "resolution_px": [2200, 1650],
        "outputs": [PNG.name, PDF.name, TIFF.name],
    }
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
