#!/usr/bin/env python3
"""OVITO 并排渲染原始模型与保持 1:1 质量比的加粗版本。"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASELINE = REPO / "models/ag_ti3sic2_qsgs_200A_ttm_finer/ag_ti3sic2_qsgs.extxyz"
STEP4 = HERE / "output/step04/ag_ti3sic2_skeleton_thickened_step04.extxyz"
FIG = HERE / "visualization"

renderer_path = HERE / "render_step01_ovito.py"
spec = importlib.util.spec_from_file_location("ovito_step_renderer", renderer_path)
if spec is None or spec.loader is None:
    raise RuntimeError(f"无法加载 OVITO 渲染工具：{renderer_path}")
renderer = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = renderer
spec.loader.exec_module(renderer)

renderer.BASE_MODEL = BASELINE
renderer.NEW_MODEL = STEP4
renderer.PNG = FIG / "Ag_Ti3SiC2_skeleton_thickened_equal_mass_1to1.png"
renderer.PDF = FIG / "Ag_Ti3SiC2_skeleton_thickened_equal_mass_1to1.pdf"
renderer.TIFF = FIG / "Ag_Ti3SiC2_skeleton_thickened_equal_mass_1to1.tif"
renderer.REPORT = FIG / "equal_mass_1to1_render_report.json"
renderer.TITLE = "Ag–Ti3SiC2 原子模型：保持 1:1 质量比"
renderer.SUBTITLE = "原模型与可行的小幅骨架加粗；材料区和模拟盒不变"
renderer.LEFT_LABEL = "原始纳米模型：Ag/Ti3SiC2 = 1:1"
renderer.RIGHT_LABEL = "加粗 0.4656 Å/侧：Ag/Ti3SiC2 = 1:1"


if __name__ == "__main__":
    renderer.main()
