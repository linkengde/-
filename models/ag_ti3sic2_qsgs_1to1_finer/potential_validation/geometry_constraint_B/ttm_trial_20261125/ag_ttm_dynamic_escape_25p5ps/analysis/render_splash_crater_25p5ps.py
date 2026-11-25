#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读读取保存的 25.5 ps dump，用 OVITO Tachyon 导出喷溅全貌和热区局部 3D 图。"""
from pathlib import Path
import sys
import numpy as np
from PySide6.QtCore import QPoint, QRect, Qt
from PySide6.QtGui import QColor, QFont, QImage, QPainter, QPen
from ovito.io import import_file
from ovito.modifiers import DeleteSelectedModifier, PythonModifier
from ovito.vis import TachyonRenderer, Viewport

HERE = Path(__file__).resolve().parent
TRIAL = HERE.parent
FRAME = TRIAL / "checkpoint_25p5ps/frame_25.500ps_step51000.dump.gz"
ANALYSIS_DIR = HERE / "q6_ptm_reference"
sys.path.insert(0, str(ANALYSIS_DIR))
from analyze_q6_ptm_msd import read_ag  # noqa: E402

OUTPUT_OVERVIEW = HERE / "splash_crater_25p5ps_overview.png"
OUTPUT_CLOSEUP = HERE / "splash_crater_25p5ps_closeup.png"
OUTPUT_PANEL = HERE / "splash_crater_25p5ps_3D.png"


def make_atom_style(frame_path):
    """按 Ag-only Qbar6 与 PTM 结果为粒子着色；红色仅标低序结构筛选区。"""
    ids, pos, _q6, qbar6, ptm, _names = read_ag(frame_path)
    qbar_by_id = {int(i): float(q) for i, q in zip(ids, qbar6)}
    ptm_by_id = {int(i): int(t) for i, t in zip(ids, ptm)}

    def apply(_frame, data):
        particle_ids = np.asarray(data.particles["Particle Identifier"].array, dtype=np.int64)
        types = np.asarray(data.particles["Particle Type"].array, dtype=np.int32)
        xyz = np.asarray(data.particles["Position"].array, dtype=np.float64)
        is_ag = types == 4
        qbar = np.full(len(types), np.nan, dtype=np.float32)
        ptm_type = np.full(len(types), -1, dtype=np.int32)
        ag_ids = particle_ids[is_ag]
        qbar[is_ag] = np.fromiter((qbar_by_id[int(i)] for i in ag_ids), dtype=np.float32, count=len(ag_ids))
        ptm_type[is_ag] = np.fromiter((ptm_by_id[int(i)] for i in ag_ids), dtype=np.int32, count=len(ag_ids))
        r2 = (xyz[:, 0] - 100.0) ** 2 + (xyz[:, 1] - 100.0) ** 2
        hot_core = is_ag & (r2 <= 60.0**2) & (xyz[:, 2] >= 100.0) & (xyz[:, 2] <= 150.0)
        disordered_screen = hot_core & (qbar < 0.45) & (ptm_type == 0)
        above_top = is_ag & (xyz[:, 2] > 150.0)

        colors = np.empty((len(types), 3), dtype=np.float32)
        colors[~is_ag] = (0.035, 0.045, 0.060)  # 固定 Ti/Si/C 骨架
        colors[is_ag] = (0.68, 0.81, 0.91)       # Ag
        colors[disordered_screen] = (0.94, 0.08, 0.055)  # Qbar6/PTM 低序筛选区
        colors[above_top] = (1.0, 0.63, 0.02)    # 高于初始材料顶部的 Ag
        radii = np.where(is_ag, 1.05, 0.92).astype(np.float32)
        data.particles_.create_property("Color", data=colors)
        data.particles_.create_property("Radius", data=radii)
        data.particles_.create_property("Selection", data=np.zeros(len(types), dtype=np.int32))
        if _frame < 0:  # 保留处理函数对帧参数的显式引用，避免静态分析误报。
            return

    return apply


def render(filename, crop=False, size=(1440, 1040)):
    """OVITO Tachyon 3D 渲染；局部图通过 Selection 删除裁剪区域外粒子。"""
    pipeline = import_file(str(FRAME))
    pipeline.source.data.cell.vis.enabled = False
    pipeline.modifiers.append(PythonModifier(function=make_atom_style(FRAME)))
    if crop:
        def select_outside(_frame, data):
            xyz = np.asarray(data.particles["Position"].array, dtype=np.float64)
            r2 = (xyz[:, 0] - 100.0) ** 2 + (xyz[:, 1] - 100.0) ** 2
            outside = (r2 > 55.0**2) | (xyz[:, 2] < 95.0) | (xyz[:, 2] > 330.0)
            data.particles_.create_property("Selection", data=outside.astype(np.int32))
        pipeline.modifiers.append(PythonModifier(function=select_outside))
        pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.add_to_scene()
    viewport = Viewport(type=Viewport.Type.Perspective)
    viewport.camera_dir = (-1.0, -1.0, -0.62) if not crop else (-0.18, -0.18, -1.0)
    viewport.zoom_all(size=size)
    viewport.render_image(filename=str(filename), size=size, background=(1.0, 1.0, 1.0), renderer=TachyonRenderer())
    pipeline.remove_from_scene()


def compose_panel():
    """组合全貌及热源局部，并在图注中说明颜色筛选的含义和限制。"""
    panel_w, panel_h, margin = 1440, 1040, 48
    header_h, footer_h = 150, 105
    canvas = QImage(panel_w * 2 + margin * 3, header_h + panel_h + footer_h, QImage.Format.Format_RGB32)
    canvas.fill(QColor("white"))
    painter = QPainter(canvas)
    painter.setPen(QColor(20, 28, 36))
    painter.setFont(QFont("DejaVu Sans", 34, QFont.Weight.Bold))
    painter.drawText(QPoint(margin, 58), "Ag/Ti3SiC2 | dynamic escape branch | 25.5 ps | T_Ag=4478 K")
    painter.setPen(QColor(65, 73, 82))
    painter.setFont(QFont("DejaVu Sans", 19))
    painter.drawText(QPoint(margin, 106), "OVITO Tachyon | T_Ag is mean of TTM-coupled Ag | z box extended to 400 A")
    y = header_h
    for j, file in enumerate((OUTPUT_OVERVIEW, OUTPUT_CLOSEUP)):
        x = margin + j * (panel_w + margin)
        painter.drawImage(QPoint(x, y), QImage(str(file)))
        painter.setPen(QPen(QColor(190, 197, 204), 2))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawRect(QRect(x, y, panel_w, panel_h))
        painter.setPen(QColor(25, 33, 40))
        painter.setFont(QFont("DejaVu Sans", 23, QFont.Weight.Bold))
        painter.drawText(QPoint(x + 18, y + 40), "Full model" if j == 0 else "Heated zone and ejecta close-up")
    painter.setPen(QColor(31, 38, 45))
    painter.setFont(QFont("DejaVu Sans", 17))
    legend_y = header_h + panel_h + 32
    entries = [((9, 12, 15), "Rigid Ti/Si/C scaffold"), ((173, 207, 232), "Ag"),
               ((240, 20, 14), "Qbar6<0.45 + PTM Other screen"), ((255, 160, 5), "Ag above initial z=150 A")]
    x_positions = [margin, margin + 410, margin + 700, margin + 1210]
    for x, (rgb, label) in zip(x_positions, entries):
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(*rgb))
        painter.drawEllipse(QRect(x, legend_y - 17, 20, 20))
        painter.setPen(QColor(31, 38, 45))
        painter.setFont(QFont("DejaVu Sans", 15))
        painter.drawText(QPoint(x + 27, legend_y), label)
    painter.setPen(QColor(90, 95, 100))
    painter.setFont(QFont("DejaVu Sans", 13))
    painter.drawText(QPoint(margin, legend_y + 36),
                     "Red marks a local disorder screen, not a standalone proof of liquid phase; the image reads saved dump data only.")
    painter.end()
    canvas.save(str(OUTPUT_PANEL), "PNG")


def main():
    if not FRAME.exists():
        raise FileNotFoundError(FRAME)
    print("OVITO: rendering full 3D structure", flush=True)
    render(OUTPUT_OVERVIEW)
    print("OVITO: rendering central pool close-up", flush=True)
    render(OUTPUT_CLOSEUP, crop=True)
    compose_panel()
    print(f"Wrote {OUTPUT_PANEL}", flush=True)


if __name__ == "__main__":
    main()
