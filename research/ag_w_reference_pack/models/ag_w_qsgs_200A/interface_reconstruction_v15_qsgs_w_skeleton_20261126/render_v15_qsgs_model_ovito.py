#!/usr/bin/env python3
"""Render the V15 QSGS W-backbone / Ag-complement model in OVITO."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from ovito.io import import_file
from ovito.modifiers import DeleteSelectedModifier, ExpressionSelectionModifier, PythonScriptModifier
from ovito.vis import TachyonRenderer, Viewport


HERE = Path(__file__).resolve().parent
DATA = HERE / "atomic_candidate" / "ag_w_qsgs.data.gz"
FULL = HERE / "V15_full_box_black_W_pale_blue_Ag.png"
CUT = HERE / "V15_cutaway_black_W_pale_blue_Ag.png"
PREVIEW = HERE / "V15_QSGS_model_preview_CN.png"


def color_particles(frame: int, data) -> None:
    types = np.asarray(data.particles["Particle Type"], dtype=np.int32)
    colors = np.empty((len(types), 3), dtype=np.float32)
    # Type 1 is W, type 2 is Ag (see specorder in the geometry builder).
    colors[types == 1] = (0.025, 0.030, 0.035)  # near-black W scaffold
    colors[types == 2] = (0.58, 0.79, 0.89)  # pale-blue Ag matrix
    # Half of the first-neighbor spacing for bcc W and fcc Ag, respectively.
    radii = np.where(types == 1, 1.368, 1.446).astype(np.float32)
    data.particles_.create_property("Color", data=colors)
    data.particles_.create_property("Radius", data=radii)


def render(cutaway: bool, filename: Path) -> None:
    pipeline = import_file(str(DATA), atom_style="atomic")
    # Keep the full 200 x 200 x 300 A box visible; only the material occupies z=0..200 A.
    pipeline.source.data.cell.vis.enabled = True
    pipeline.modifiers.append(PythonScriptModifier(function=color_particles))
    pipeline.modifiers.append(ExpressionSelectionModifier(expression="Position.Z < 0 || Position.Z > 200"))
    pipeline.modifiers.append(DeleteSelectedModifier())
    if cutaway:
        pipeline.modifiers.append(ExpressionSelectionModifier(expression="Position.X < 100"))
        pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.add_to_scene()

    viewport = Viewport(type=Viewport.Type.Perspective)
    viewport.camera_dir = (1.0, -1.0, -0.72)
    viewport.camera_up = (0.0, 0.0, 1.0)
    viewport.zoom_all(size=(1500, 1400))
    renderer = TachyonRenderer()
    renderer.antialiasing = True
    renderer.antialiasing_samples = 8
    renderer.ambient_occlusion = True
    viewport.render_image(filename=str(filename), size=(1500, 1400), background=(0.98, 0.985, 0.99), renderer=renderer)
    pipeline.remove_from_scene()


def font(path: str, size: int) -> ImageFont.ImageFont:
    try:
        return ImageFont.truetype(path, size)
    except OSError:
        return ImageFont.load_default()


def main() -> None:
    render(False, FULL)
    render(True, CUT)
    left = Image.open(FULL).convert("RGB")
    right = Image.open(CUT).convert("RGB")
    header, gap = 148, 34
    canvas = Image.new("RGB", (left.width + right.width + gap, left.height + header), (250, 251, 252))
    draw = ImageDraw.Draw(canvas)
    cjk = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
    title_font = font(cjk, 34)
    body_font = font(cjk, 22)
    caption_font = font(cjk, 20)
    draw.text((30, 14), "Ag–W 四参数 QSGS：黑色 W 骨架 + 浅蓝色 Ag 互补相", fill=(18, 28, 36), font=title_font)
    draw.text((32, 65), "seed=20261141 | 材料 200×200×200 Å | 计算盒 200×200×300 Å | W≈35.41 vol%", fill=(54, 68, 80), font=body_font)
    draw.text((32, 105), "初始几何预览；Ag/W 原子界面重叠尚未消除，不能作为 MD 初态", fill=(132, 50, 35), font=body_font)
    draw.text((24, header + 12), "完整计算盒（顶部为真空）", fill=(30, 44, 55), font=caption_font)
    draw.text((left.width + gap + 24, header + 12), "半剖视图（移除 x<100 Å）", fill=(30, 44, 55), font=caption_font)
    canvas.paste(left, (0, header))
    canvas.paste(right, (left.width + gap, header))
    canvas.save(PREVIEW, dpi=(300, 300))
    print(PREVIEW)


if __name__ == "__main__":
    main()
