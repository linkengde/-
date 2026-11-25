#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 OVITO Tachyon 渲染细化骨架原子模型的全貌和剖切图。"""

from pathlib import Path
import json

import numpy as np
from PySide6.QtCore import QPoint, QRect, Qt
from PySide6.QtGui import QColor, QFont, QImage, QPainter, QPen
from ovito.io import import_file
from ovito.modifiers import DeleteSelectedModifier, PythonModifier
from ovito.vis import TachyonRenderer, Viewport


HERE = Path(__file__).resolve().parent
STRUCTURE = HERE / "ag_ti3sic2_qsgs.extxyz"


def render(cut=False, filename="preview.png", size=(1500, 1050)):
    """渲染完整模型；剖切视图保留 y>=100 Å 的一半以显示内部形貌。"""
    pipeline = import_file(str(STRUCTURE))
    pipeline.source.data.cell.vis.enabled = False

    def set_particle_style(_frame, data):
        type_ids = np.asarray(data.particles["Particle Type"].array, dtype=np.int32)
        pos = np.asarray(data.particles["Position"].array, dtype=np.float64)
        ag = type_ids == 4
        color = np.empty((len(type_ids), 3), dtype=np.float32)
        color[ag] = (0.69, 0.76, 0.82)
        color[~ag] = (0.10, 0.32, 0.66)
        radius = np.where(ag, 1.03, 0.90).astype(np.float32)
        if cut:
            radius *= 1.05
            remove = pos[:, 1] < 100.0
            data.particles_.create_property("Selection", data=remove.astype(np.int32))
        data.particles_.create_property("Color", data=color)
        data.particles_.create_property("Radius", data=radius)

    pipeline.modifiers.append(PythonModifier(function=set_particle_style))
    if cut:
        pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.add_to_scene()

    viewport = Viewport(type=Viewport.Type.Perspective)
    viewport.camera_dir = (-1.0, -1.0, -0.62)
    viewport.zoom_all(size=size)
    viewport.render_image(
        filename=str(filename),
        size=size,
        background=(1.0, 1.0, 1.0),
        renderer=TachyonRenderer(),
    )
    pipeline.remove_from_scene()


def main():
    if not STRUCTURE.exists():
        raise FileNotFoundError(STRUCTURE)
    full = HERE / "thinner_scaffold_OVITO_full.png"
    section = HERE / "thinner_scaffold_OVITO_center_cut.png"
    render(False, full)
    render(True, section)

    panel_w, panel_h, margin = 1500, 1050, 56
    header_h, footer_h = 140, 100
    canvas = QImage(panel_w * 2 + margin * 3, header_h + panel_h + footer_h,
                    QImage.Format.Format_RGB32)
    canvas.fill(QColor("white"))
    painter = QPainter(canvas)
    painter.setPen(QColor(25, 33, 43))
    painter.setFont(QFont("DejaVu Sans", 37, QFont.Weight.Bold))
    geometry = json.loads((HERE / "geometry_report.json").read_text())
    atoms_k = geometry["atoms"] / 1000.0
    painter.drawText(QPoint(margin, 62), f"QSGS thinner scaffold candidate 05 | {atoms_k:.0f}k atoms | seed 20261125")
    painter.setPen(QColor(66, 76, 88))
    painter.setFont(QFont("DejaVu Sans", 23))
    painter.drawText(QPoint(margin, 108),
                     "Same 200 x 200 x 150 A size and 1:1 phase mass | QSGS core probability 0.064")

    left, right = QImage(str(full)), QImage(str(section))
    x1, x2, y = margin, margin * 2 + panel_w, header_h
    painter.drawImage(QPoint(x1, y), left)
    painter.drawImage(QPoint(x2, y), right)
    painter.setPen(QPen(QColor(190, 196, 202), 2))
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawRect(QRect(x1, y, panel_w, panel_h))
    painter.drawRect(QRect(x2, y, panel_w, panel_h))
    painter.setPen(QColor(30, 37, 44))
    painter.setFont(QFont("DejaVu Sans", 27, QFont.Weight.Bold))
    painter.drawText(QPoint(x1 + 15, y + 42), "A  Full atomic model")
    painter.drawText(QPoint(x2 + 15, y + 42), "B  Center cut (y < 100 A removed)")

    legend_y = y + panel_h + 34
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(QColor(26, 82, 168))
    painter.drawEllipse(QRect(margin, legend_y, 26, 26))
    painter.setPen(QColor(35, 43, 51))
    painter.setFont(QFont("DejaVu Sans", 22))
    painter.drawText(QPoint(margin + 38, legend_y + 22), "Ti3SiC2 skeleton (Ti/Si/C)")
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(QColor(176, 194, 210))
    painter.drawEllipse(QRect(margin + 620, legend_y, 26, 26))
    painter.setPen(QColor(35, 43, 51))
    painter.drawText(QPoint(margin + 658, legend_y + 22), "Ag")
    painter.setPen(QColor(80, 87, 94))
    painter.setFont(QFont("DejaVu Sans", 18))
    painter.drawText(QPoint(margin + 810, legend_y + 22),
                     "New stochastic geometry only; not relaxed or connected to prior heating checkpoints.")
    painter.end()
    out = HERE / "thinner_scaffold_OVITO.png"
    canvas.save(str(out), "PNG")
    print(f"Wrote {out}", flush=True)


if __name__ == "__main__":
    main()
