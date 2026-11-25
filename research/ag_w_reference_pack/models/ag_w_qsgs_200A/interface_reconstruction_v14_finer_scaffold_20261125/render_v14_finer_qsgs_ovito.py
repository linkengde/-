#!/usr/bin/env python3
"""OVITO rendering for the V14 finer QSGS W-backbone candidate."""
from __future__ import annotations
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from ovito.io import import_file
from ovito.modifiers import DeleteSelectedModifier, ExpressionSelectionModifier, PythonScriptModifier
from ovito.vis import TachyonRenderer, Viewport

HERE = Path(__file__).resolve().parent
DATA = HERE / "qsgs_v14_finer_candidate/ag_w_qsgs.data.gz"
FULL = HERE / "v14_finer_qsgs_full_model_ovito.png"
CUT = HERE / "v14_finer_qsgs_cutaway_ovito.png"
COMPOSITE = HERE / "v14_finer_qsgs_morphology_preview.png"

def appearance(frame: int, data) -> None:
    typ = np.asarray(data.particles["Particle Type"], dtype=np.int32)
    color = np.empty((len(typ), 3), dtype=np.float32)
    color[typ == 1] = (0.10, 0.25, 0.36)
    color[typ == 2] = (0.59, 0.76, 0.85)
    data.particles_.create_property("Color", data=color)
    data.particles_.create_property("Radius", data=np.where(typ == 1, 1.10, 1.16).astype(np.float32))

def render(cutaway: bool, output: Path) -> None:
    pipeline = import_file(str(DATA), atom_style="atomic")
    pipeline.source.data.cell.vis.enabled = False
    pipeline.modifiers.append(PythonScriptModifier(function=appearance))
    pipeline.modifiers.append(ExpressionSelectionModifier(expression="Position.Z < 0 || Position.Z > 200"))
    pipeline.modifiers.append(DeleteSelectedModifier())
    if cutaway:
        pipeline.modifiers.append(ExpressionSelectionModifier(expression="Position.X < 100"))
        pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.add_to_scene()
    vp = Viewport(type=Viewport.Type.Perspective)
    vp.camera_dir = (1.0, -1.0, -0.6)
    vp.camera_up = (0.0, 0.0, 1.0)
    vp.zoom_all(size=(1600, 1350))
    renderer = TachyonRenderer()
    renderer.antialiasing = True
    renderer.antialiasing_samples = 8
    renderer.ambient_occlusion = True
    vp.render_image(filename=str(output), size=(1600, 1350), background=(0.97, 0.98, 0.99), renderer=renderer)
    pipeline.remove_from_scene()

def main() -> None:
    render(False, FULL)
    render(True, CUT)
    a, b = Image.open(FULL).convert("RGB"), Image.open(CUT).convert("RGB")
    header, gutter = 100, 36
    out = Image.new("RGB", (a.width + b.width + gutter, a.height + header), (247, 249, 251))
    d = ImageDraw.Draw(out)
    try:
        title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 28)
        small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 18)
    except OSError:
        title = small = ImageFont.load_default()
    d.text((28, 12), "Ag-W V14 | finer QSGS W-backbone preview", fill=(20, 35, 49), font=title)
    d.text((28, 56), "Dark blue: W scaffold   Pale blue: Ag matrix   Geometry only; Ag-W overlaps remain", fill=(58, 76, 92), font=small)
    d.text((24, header + 14), "Full 200 × 200 × 200 Å material", fill=(30, 45, 58), font=small)
    d.text((a.width + gutter + 24, header + 14), "Cutaway: x < 100 Å hidden", fill=(30, 45, 58), font=small)
    out.paste(a, (0, header))
    out.paste(b, (a.width + gutter, header))
    out.save(COMPOSITE)
    print(COMPOSITE)

if __name__ == "__main__":
    main()
