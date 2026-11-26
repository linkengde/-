#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 OVITO Tachyon 将 20.6 ps 原子帧渲染成全模型与热区双视图。"""
from pathlib import Path
import sys
import numpy as np

from PySide6.QtCore import QPoint, QRect, Qt
from PySide6.QtGui import QColor, QFont, QImage, QPainter, QPen
from ovito.io import import_file
from ovito.modifiers import DeleteSelectedModifier, PythonModifier
from ovito.vis import TachyonRenderer, Viewport

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FRAME = ROOT / "cooling_20p55_to_20p60ps_z500_v6" / "frame_20.600ps_step41200.dump.gz"
sys.path.insert(0, str(HERE))
from analyze_q6_ptm_msd import read_ag  # noqa: E402

OUT_FULL = HERE / "candidate05_20p6ps_OVITO_Tachyon_full.png"
OUT_CLOSE = HERE / "candidate05_20p6ps_OVITO_Tachyon_hotspot.png"
OUT_PANEL = HERE / "candidate05_20p6ps_OVITO_Tachyon.png"

# 0: scaffold (black), 1: normal Ag (light blue), 2: local disorder screen (red), 3: Ag above initial top (orange)
COLORS = {
    0: np.array([0.025, 0.035, 0.045], dtype=np.float32),
    1: np.array([0.63, 0.82, 0.93], dtype=np.float32),
    2: np.array([0.92, 0.055, 0.045], dtype=np.float32),
    3: np.array([1.00, 0.58, 0.015], dtype=np.float32),
}


def make_style_modifier(qbar_by_id, ptm_by_id, crop):
    """应用示意图的四类配色，并为右侧视图建立中心热区选择。"""
    def apply(_frame, data):
        ids = np.asarray(data.particles["Particle Identifier"].array, dtype=np.int64)
        types = np.asarray(data.particles["Particle Type"].array, dtype=np.int32)
        pos = np.asarray(data.particles["Position"].array, dtype=np.float64)
        is_ag = types == 4
        qbar = np.full(len(ids), np.nan, dtype=np.float32)
        ptm = np.full(len(ids), -1, dtype=np.int32)
        if np.any(is_ag):
            qbar[is_ag] = np.fromiter((qbar_by_id[int(i)] for i in ids[is_ag]),
                                      dtype=np.float32, count=int(is_ag.sum()))
            ptm[is_ag] = np.fromiter((ptm_by_id[int(i)] for i in ids[is_ag]),
                                     dtype=np.int32, count=int(is_ag.sum()))
        colors = np.empty((len(ids), 3), dtype=np.float32)
        colors[:] = COLORS[0]  # Ti/Si/C 刚性骨架
        colors[is_ag] = COLORS[1]
        disorder = is_ag & (qbar < 0.45) & (ptm == 0)
        colors[disorder] = COLORS[2]
        elevated = is_ag & (pos[:, 2] > 150.0)
        colors[elevated] = COLORS[3]  # 优先标记初始材料顶面 z=150 Å 以上的 Ag
        # 轻微缩小球半径，避免高密度视图中过度遮挡；仍由 OVITO/Tachyon 绘制三维球体。
        radii = np.where(is_ag, 1.00, 0.82).astype(np.float32)
        data.particles_.create_property("Color", data=colors)
        data.particles_.create_property("Radius", data=radii)
        if crop:
            r2 = (pos[:, 0] - 100.0)**2 + (pos[:, 1] - 100.0)**2
            # 热区上下文：保留中心半径 72 Å、z=112–315 Å 范围中的原子。
            outside = (r2 > 72.0**2) | (pos[:, 2] < 112.0) | (pos[:, 2] > 315.0)
            data.particles_.create_property("Selection", data=outside.astype(np.int32))
    return apply


def render_panel(output, crop=False, size=(1250, 820)):
    """使用 Tachyon 对单个 OVITO Pipeline 视口进行离屏光线追踪。"""
    pipeline = import_file(str(FRAME))
    pipeline.source.data.cell.vis.enabled = False
    ids, _pos, _q6, qbar6, ptm, _names = read_ag(FRAME)
    qbar_by_id = {int(i): float(q) for i, q in zip(ids, qbar6)}
    ptm_by_id = {int(i): int(v) for i, v in zip(ids, ptm)}
    pipeline.modifiers.append(PythonModifier(function=make_style_modifier(qbar_by_id, ptm_by_id, crop)))
    if crop:
        pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.add_to_scene()
    viewport = Viewport(type=Viewport.Type.Perspective)
    viewport.camera_dir = (-1.0, -1.0, -0.52) if not crop else (-0.24, -1.0, -0.42)
    viewport.camera_up = (0.0, 0.0, 1.0)
    viewport.zoom_all(size=size)
    viewport.render_image(filename=str(output), size=size, background=(1.0, 1.0, 1.0),
                          renderer=TachyonRenderer())
    pipeline.remove_from_scene()
    print(f"RENDERED {output}", flush=True)


def compose_panel():
    """按参考图排成双面板，并加入帧时间、温度与颜色图例。"""
    panel_w, panel_h, margin, header, footer = 1250, 820, 38, 110, 100
    canvas = QImage(panel_w * 2 + margin * 3, header + panel_h + footer,
                    QImage.Format.Format_RGB32)
    canvas.fill(QColor("white"))
    painter = QPainter(canvas)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setPen(QColor(28, 36, 45))
    painter.setFont(QFont("DejaVu Sans", 29, QFont.Weight.Bold))
    painter.drawText(QPoint(margin, 43), "Ag/Ti3SiC2 | 20 ps heating + ideal-sink cooling | 20.6 ps | T_Ag ≈ 4172 K")
    painter.setPen(QColor(75, 83, 92))
    painter.setFont(QFont("DejaVu Sans", 16))
    painter.drawText(QPoint(margin, 78), "OVITO Tachyon | T_Ag is mean TTM lattice temperature | z box extended to 500 Å to retain ejected Ag")
    y = header
    titles = ((OUT_FULL, "Full model"), (OUT_CLOSE, "Heated zone and elevated Ag close-up"))
    for j, (img_path, title) in enumerate(titles):
        x = margin + j * (panel_w + margin)
        painter.drawImage(QPoint(x, y), QImage(str(img_path)))
        painter.setPen(QPen(QColor(190, 198, 206), 2))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawRect(QRect(x, y, panel_w, panel_h))
        painter.setPen(QColor(28, 36, 45))
        painter.setFont(QFont("DejaVu Sans", 20, QFont.Weight.Bold))
        painter.drawText(QPoint(x + 16, y + 34), title)

    # 底部图例，颜色与渲染属性使用完全相同的 RGB 值。
    legend_y = header + panel_h + 30
    legend = [
        (COLORS[0], "Rigid Ti/Si/C scaffold"),
        (COLORS[1], "Ag"),
        (COLORS[2], "qbar6 < 0.45 + PTM Other screen"),
        (COLORS[3], "Ag above initial z = 150 Å"),
    ]
    x = margin
    painter.setFont(QFont("DejaVu Sans", 14))
    for color, label in legend:
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor.fromRgbF(float(color[0]), float(color[1]), float(color[2])))
        painter.drawEllipse(QPoint(x + 9, legend_y - 5), 8, 8)
        painter.setPen(QColor(45, 53, 62))
        painter.drawText(QPoint(x + 26, legend_y), label)
        x += 26 + painter.fontMetrics().horizontalAdvance(label) + 30
    painter.setPen(QColor(78, 85, 93))
    painter.setFont(QFont("DejaVu Sans", 12))
    painter.drawText(QPoint(margin, header + panel_h + 58),
        "Red is a local disorder screen, not standalone proof of liquid phase; orange marks Ag above the initial material top.")
    painter.end()
    if not canvas.save(str(OUT_PANEL), "PNG"):
        raise RuntimeError(f"Could not save {OUT_PANEL}")
    print(f"COMPOSED {OUT_PANEL}", flush=True)


if __name__ == "__main__":
    if not FRAME.exists():
        raise FileNotFoundError(FRAME)
    render_panel(OUT_FULL, crop=False)
    render_panel(OUT_CLOSE, crop=True)
    compose_panel()
