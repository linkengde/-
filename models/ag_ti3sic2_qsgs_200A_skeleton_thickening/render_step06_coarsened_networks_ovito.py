#!/usr/bin/env python3
"""用 OVITO StandardRenderer 对比 step05 与两相网络粗化后的 step06。"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
STEP5 = HERE / "output/step05_both_phases/ag_ti3sic2_both_phases_thickened_step05.extxyz"
STEP6 = HERE / "output/step06_coarsened_networks/ag_ti3sic2_coarsened_networks_step06.extxyz"
FIG = HERE / "visualization"

renderer_path = HERE / "render_step01_ovito.py"
spec = importlib.util.spec_from_file_location("ovito_step_renderer", renderer_path)
if spec is None or spec.loader is None:
    raise RuntimeError(f"无法加载 OVITO 渲染工具：{renderer_path}")
renderer = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = renderer
spec.loader.exec_module(renderer)

renderer.BASE_MODEL = STEP5
renderer.NEW_MODEL = STEP6
renderer.PNG = FIG / "Ag_Ti3SiC2_network_coarsened_step06.png"
renderer.PDF = FIG / "Ag_Ti3SiC2_network_coarsened_step06.pdf"
renderer.TIFF = FIG / "Ag_Ti3SiC2_network_coarsened_step06.tif"
renderer.REPORT = FIG / "coarsened_step06_render_report.json"
renderer.TITLE = "Ag–Ti3SiC2 原子模型：两相网络继续加粗"
renderer.SUBTITLE = "在原子总数、1:1 质量比和盒子尺寸不变时，粗化细枝并重分配相内原子"
renderer.LEFT_LABEL = "step05：同步扩展后的枝状网络"
renderer.RIGHT_LABEL = "step06：细枝平滑、主枝加粗；相质量不变"


if __name__ == "__main__":
    renderer.main()
    # 底层共用的渲染器历史字段名为 step01；本渲染结果按当前步骤重命名。
    report = json.loads(renderer.REPORT.read_text())
    report["step06_atoms"] = report.pop("step01_atoms")
    report["step06_phase_counts"] = report.pop("step01_phase_counts")
    renderer.REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
