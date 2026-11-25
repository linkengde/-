#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 OVITO 将 5.5 ps 热区液态样 Ag 熔池与固定骨架渲染成双视图图片。"""

from pathlib import Path
import sys

import numpy as np
from PySide6.QtCore import QPoint, QRect, Qt
from PySide6.QtGui import QColor, QFont, QImage, QPainter, QPen
from ovito.io import import_file
from ovito.modifiers import DeleteSelectedModifier, PythonModifier
from ovito.vis import TachyonRenderer, Viewport


HERE = Path(__file__).resolve().parent
FRAME = HERE / "frame_05.500ps_step11000.dump.gz"
OUTPUT = HERE / "Ag_melt_pool_5p5ps_OVITO.png"
ANALYSIS_DIR = HERE
sys.path.insert(0, str(ANALYSIS_DIR))
from analyze_q6_ptm_msd import read_ag  # noqa: E402


def load_ag_qbar(frame_path):
    """使用与相态分析相同的 Ag-only 邻居图得到逐原子 qbar6。"""
    ids, _pos, _q6, qbar6, _ptm, _names = read_ag(frame_path)
    return {int(atom_id): float(qbar) for atom_id, qbar in zip(ids, qbar6)}


def render_view(closeup, filename, size=(1320, 980)):
    """为全貌或热区局部建立单独的 OVITO scene 并用 Tachyon 渲染。"""
    pipeline = import_file(str(FRAME))
    pipeline.source.data.cell.vis.enabled = False

    qbar_by_id = load_ag_qbar(FRAME)
    r_cut_hot = 30.0
    z_hot_min = 120.0
    z_top = 150.0

    def assign_figure_properties(frame, data):
        ids = np.asarray(data.particles["Particle Identifier"].array, dtype=np.int64)
        types = np.asarray(data.particles["Particle Type"].array, dtype=np.int32)
        pos = np.asarray(data.particles["Position"].array, dtype=np.float64)
        qbar = np.full(len(ids), np.nan, dtype=np.float32)

        ag_mask = types == 4
        qbar[ag_mask] = np.fromiter(
            (qbar_by_id[int(atom_id)] for atom_id in ids[ag_mask]),
            dtype=np.float32,
            count=int(ag_mask.sum()),
        )

        # 只把预定义热源核心内、qbar6<0.45 的 Ag 染红；其余 Ag 为银蓝色。
        r2 = (pos[:, 0] - 100.0) ** 2 + (pos[:, 1] - 100.0) ** 2
        hot_core = (
            ag_mask
            & (r2 <= r_cut_hot**2)
            & (pos[:, 2] >= z_hot_min)
            & (pos[:, 2] <= z_top)
        )
        liquid_like = hot_core & (qbar < 0.45)

        # Ti/Si/C 合并显示为深灰色支撑骨架，Ag 使用冷色，液态样热区用红色。
        colors = np.empty((len(ids), 3), dtype=np.float32)
        colors[~ag_mask] = (0.16, 0.20, 0.24)
        colors[ag_mask] = (0.66, 0.79, 0.88)
        colors[liquid_like] = (0.92, 0.12, 0.08)

        radii = np.where(ag_mask, 1.05, 0.95).astype(np.float32)
        if closeup:
            radii *= 1.08

        data.particles_.create_property("Color", data=colors)
        data.particles_.create_property("Radius", data=radii)

        if closeup:
            # 局部图显示熔池及相邻骨架/Ag，范围为热源中心附近。
            outside = (
                (r2 > 46.0**2)
                | (pos[:, 2] < 88.0)
                | (pos[:, 2] > z_top)
            )
            data.particles_.create_property(
                "Selection", data=outside.astype(np.int32)
            )

    pipeline.modifiers.append(PythonModifier(function=assign_figure_properties))
    if closeup:
        pipeline.modifiers.append(DeleteSelectedModifier())

    pipeline.add_to_scene()
    vp = Viewport(type=Viewport.Type.Perspective)
    vp.camera_dir = (-1.0, -1.0, -0.58)
    vp.zoom_all(size=size)
    vp.render_image(
        filename=str(filename),
        size=size,
        background=(1.0, 1.0, 1.0),
        renderer=TachyonRenderer(),
    )
    pipeline.remove_from_scene()


def main():
    if not FRAME.exists():
        raise FileNotFoundError(f"找不到 5.5 ps 轨迹帧：{FRAME}")

    overview = HERE / "_overview_tmp.png"
    detail = HERE / "_detail_tmp.png"
    print("正在用 OVITO 计算 Ag qbar6 并渲染全模型视图……", flush=True)
    render_view(False, overview)
    print("正在渲染热源区局部视图……", flush=True)
    render_view(True, detail)

    # 使用 Qt 绘图添加标题、相态图例和判读口径，得到双视图图版。
    panel_w, panel_h = 1320, 980
    margin = 64
    header_h = 145
    footer_h = 150
    canvas = QImage(
        panel_w * 2 + margin * 3,
        header_h + panel_h + footer_h,
        QImage.Format.Format_RGB32,
    )
    canvas.fill(QColor("white"))
    painter = QPainter(canvas)
    painter.setPen(QColor(24, 31, 38))
    painter.setFont(QFont("DejaVu Sans", 38, QFont.Weight.Bold))
    painter.drawText(
        QPoint(margin, 72),
        "Ag local melt-like region in the 347k Ag/Ti3SiC2 model",
    )
    painter.setPen(QColor(64, 73, 82))
    painter.setFont(QFont("DejaVu Sans", 25))
    painter.drawText(
        QPoint(margin, 118),
        "OVITO snapshot at 5.5 ps | red = hot-core Ag with qbar6 < 0.45",
    )

    left = QImage(str(overview))
    right = QImage(str(detail))
    x_left = margin
    x_right = margin * 2 + panel_w
    y_top = header_h
    painter.drawImage(QPoint(x_left, y_top), left)
    painter.drawImage(QPoint(x_right, y_top), right)
    painter.setPen(QPen(QColor(190, 196, 202), 2))
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawRect(QRect(x_left, y_top, panel_w, panel_h))
    painter.drawRect(QRect(x_right, y_top, panel_w, panel_h))
    painter.setPen(QColor(30, 37, 44))
    painter.setFont(QFont("DejaVu Sans", 27, QFont.Weight.Bold))
    painter.drawText(QPoint(x_left + 14, y_top + 42), "A  Full composite")
    painter.drawText(QPoint(x_right + 14, y_top + 42), "B  Hot-zone close-up")

    legend_y = y_top + panel_h + 40
    legend = [
        ((41, 51, 61), "Ti3SiC2 rigid scaffold"),
        ((168, 201, 224), "Ag outside hot-core liquid-like mask"),
        ((235, 31, 20), "Ag in hot core: qbar6 < 0.45"),
    ]
    legend_x = [margin, margin + 850, margin + 1700]
    for x, (color, label) in zip(legend_x, legend):
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(*color))
        painter.drawEllipse(QRect(x, legend_y + 2, 28, 28))
        painter.setPen(QColor(38, 44, 50))
        painter.setFont(QFont("DejaVu Sans", 22))
        painter.drawText(QPoint(x + 40, legend_y + 24), label)

    painter.setPen(QColor(86, 91, 97))
    painter.setFont(QFont("DejaVu Sans", 19))
    painter.drawText(
        QPoint(margin, legend_y + 66),
        "Exploratory structural classification; hot-core mask: r_xy ≤ 30 Å, z ≥ 120 Å. Not the entire Ag phase.",
    )
    painter.end()
    canvas.save(str(OUTPUT), "PNG")
    overview.unlink(missing_ok=True)
    detail.unlink(missing_ok=True)
    print(f"Wrote {OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
