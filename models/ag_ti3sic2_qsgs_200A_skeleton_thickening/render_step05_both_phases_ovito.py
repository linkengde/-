#!/usr/bin/env python3
"""用 OVITO StandardRenderer 对比原模型与两相同步加粗模型。"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASELINE = REPO / "models/ag_ti3sic2_qsgs_200A_ttm_finer/ag_ti3sic2_qsgs.extxyz"
STEP5 = HERE / "output/step05_both_phases/ag_ti3sic2_both_phases_thickened_step05.extxyz"
FIG = HERE / "visualization"

renderer_path = HERE / "render_step01_ovito.py"
spec = importlib.util.spec_from_file_location("ovito_step_renderer", renderer_path)
if spec is None or spec.loader is None:
    raise RuntimeError(f"无法加载 OVITO 渲染工具：{renderer_path}")
renderer = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = renderer
spec.loader.exec_module(renderer)

renderer.BASE_MODEL = BASELINE
renderer.NEW_MODEL = STEP5
renderer.PNG = FIG / "Ag_Ti3SiC2_both_phases_thickened_step05.png"
renderer.PDF = FIG / "Ag_Ti3SiC2_both_phases_thickened_step05.pdf"
renderer.TIFF = FIG / "Ag_Ti3SiC2_both_phases_thickened_step05.tif"
renderer.REPORT = FIG / "both_phases_step05_render_report.json"
renderer.TITLE = "Ag–Ti3SiC2 原子模型：骨架与银网络同步加粗"
renderer.SUBTITLE = "Ag/Ti3SiC2 质量比保持 1:1；材料区和模拟盒尺寸不变"
renderer.LEFT_LABEL = "原始纳米原子模型：质量比 1:1"
renderer.RIGHT_LABEL = "同步加粗：Ti3SiC2 +0.80 Å/侧；Ag 网络扩展"


if __name__ == "__main__":
    renderer.main()
    # 底层共用的渲染器历史字段名为 step01；本渲染结果按当前步骤重命名。
    report = json.loads(renderer.REPORT.read_text())
    report["step05_atoms"] = report.pop("step01_atoms")
    report["step05_phase_counts"] = report.pop("step01_phase_counts")
    renderer.REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
