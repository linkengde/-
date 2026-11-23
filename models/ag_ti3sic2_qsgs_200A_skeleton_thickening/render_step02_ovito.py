#!/usr/bin/env python3
"""用 OVITO 并排渲染第一档和第二档骨架加粗原子模型。"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
STEP1 = HERE / "output/ag_ti3sic2_skeleton_thickened_step01.extxyz"
STEP2 = HERE / "output/step02/ag_ti3sic2_skeleton_thickened_step02.extxyz"
FIG = HERE / "visualization"

renderer_path = HERE / "render_step01_ovito.py"
spec = importlib.util.spec_from_file_location("ovito_step_renderer", renderer_path)
if spec is None or spec.loader is None:
    raise RuntimeError(f"无法加载 OVITO 渲染工具：{renderer_path}")
renderer = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = renderer
spec.loader.exec_module(renderer)

renderer.BASE_MODEL = STEP1
renderer.NEW_MODEL = STEP2
renderer.PNG = FIG / "Ag_Ti3SiC2_skeleton_thickened_step01_vs_step02.png"
renderer.PDF = FIG / "Ag_Ti3SiC2_skeleton_thickened_step01_vs_step02.pdf"
renderer.TIFF = FIG / "Ag_Ti3SiC2_skeleton_thickened_step01_vs_step02.tif"
renderer.REPORT = FIG / "step02_render_report.json"
renderer.TITLE = "Ag–Ti3SiC2 原子模型：骨架继续加粗"
renderer.SUBTITLE = "比较累计外扩 1 Å/侧与 2 Å/侧；材料区和模拟盒保持不变"
renderer.LEFT_LABEL = "第一档：外扩 1 Å/侧"
renderer.RIGHT_LABEL = "第二档：外扩 2 Å/侧"


if __name__ == "__main__":
    renderer.main()
