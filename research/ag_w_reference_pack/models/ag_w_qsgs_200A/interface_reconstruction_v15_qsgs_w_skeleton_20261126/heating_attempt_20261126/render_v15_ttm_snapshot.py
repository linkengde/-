#!/usr/bin/env python3
"""Render a full 3D V15 Ag-W TTM snapshot with an explicitly tentative disorder overlay."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier, SteinhardtParametersModifier
from ovito.pipeline import Pipeline, StaticSource
from ovito.vis import ParticlesVis, TachyonRenderer, Viewport


def fields(path: Path):
    pipeline = import_file(str(path))
    pipeline.modifiers.append(PolyhedralTemplateMatchingModifier())
    q6 = SteinhardtParametersModifier()
    q6.max_l = 6
    q6.neighbor_mode = SteinhardtParametersModifier.NeighborMode.Voronoi
    pipeline.modifiers.append(q6)
    data = pipeline.compute()
    p = data.particles
    return pipeline, data, {
        "id": np.asarray(p["Particle Identifier"], dtype=np.int64),
        "type": np.asarray(p["Particle Type"], dtype=np.int32),
        "position": np.asarray(p["Position"], dtype=float),
        "ptm": np.asarray(p["Structure Type"], dtype=np.int32),
        "q6": np.asarray(p["Steinhardt"], dtype=float)[:, 4],
    }


def font(size: int, bold: bool = False):
    paths = [
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc" if bold else "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for item in paths:
        if Path(item).exists():
            return ImageFont.truetype(item, size)
    return ImageFont.load_default()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--start", type=Path, required=True)
    ap.add_argument("--snapshot", type=Path, required=True)
    ap.add_argument("--time-ps", type=float, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    final_pipe, data, end = fields(args.snapshot)
    start_pipe, _, start = fields(args.start)
    i0, i1 = np.argsort(start["id"]), np.argsort(end["id"])
    for f in (start, end):
        ix = np.argsort(f["id"])
        for key in f:
            f[key] = f[key][ix]
    if not np.array_equal(start["id"], end["id"]):
        raise RuntimeError("Particle IDs differ between the initial and selected snapshots")
    delta = end["position"] - start["position"]
    delta[:, :2] -= 200.0 * np.rint(delta[:, :2] / 200.0)
    displacement = np.linalg.norm(delta, axis=1)
    # The input Ag structure is already disordered, so orange marks a combined
    # structure/displacement candidate and must not be interpreted as liquid.
    selected_ag = (start["type"] == 2) & (start["position"][:, 2] >= 40.0) & (start["position"][:, 2] <= 200.0)
    ptm_other = end["ptm"] == int(PolyhedralTemplateMatchingModifier.Type.OTHER)
    candidate = selected_ag & ptm_other & (end["q6"] < 0.4) & (displacement > 2.0)
    colors = np.empty((len(end["id"]), 3), dtype=float)
    colors[end["type"] == 1] = (0.035, 0.045, 0.055)  # near-black W skeleton
    colors[end["type"] == 2] = (0.66, 0.82, 0.89)      # pale-blue Ag
    colors[candidate] = (0.96, 0.39, 0.08)             # tentative disorder candidates
    data.particles_.create_property("Color", data=colors)
    data.particles_.vis.radius = 0.88
    data.particles_.vis.shape = ParticlesVis.Shape.Sphere
    data.cell.vis.enabled = True
    data.cell.vis.render_cell = True
    data.cell.vis.line_width = 1.2
    data.cell.vis.rendering_color = (0.18, 0.20, 0.22)

    scene = Pipeline(source=StaticSource(data=data))
    scene.add_to_scene()
    size = (1700, 1250)
    vp = Viewport(type=Viewport.Type.Perspective)
    vp.camera_pos = (229.1, -122.8, 260.0)
    vp.camera_dir = (-0.4609, 0.7957, -0.3930)
    vp.camera_up = (0.0, 0.0, 1.0)
    vp.fov = 36.0
    vp.zoom_all(size=size)
    vp.camera_pos = (229.1, -122.8, 260.0)
    vp.camera_dir = (-0.4609, 0.7957, -0.3930)
    vp.fov = 36.0
    renderer = TachyonRenderer()
    renderer.antialiasing = True
    renderer.antialiasing_samples = 6
    renderer.ambient_occlusion = True
    renderer.ambient_occlusion_samples = 8
    renderer.shadows = True
    rendered = vp.render_image(size=size, background=(1.0, 1.0, 1.0), alpha=False,
                               renderer=renderer, stop_on_error=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    banner = 160
    canvas = Image.new("RGB", (size[0], size[1] + banner), "white")
    raw_render_path = args.output.with_name(args.output.stem + "_ovito_render.png")
    rendered.save(str(raw_render_path))
    with Image.open(raw_render_path) as raw:
        canvas.paste(raw.convert("RGB"), (0, banner))
    draw = ImageDraw.Draw(canvas)
    draw.text((40, 18), f"Ag–W V15 全模型｜TTM 末态 t={args.time_ps:g} ps", font=font(35, True), fill="#1d2b36")
    draw.text((43, 67), "四参数随机 QSGS：W 骨架 + Ag 填充；橙色仅表示 Q6/PTM/MSD 无序候选，不是已确认液相",
              font=font(22), fill="#526574")
    legend_y = 112
    legends = [((9, 12, 14), "W 骨架"), ((168, 209, 227), "Ag"), ((245, 99, 20), f"无序候选 {int(candidate.sum()):,}")]
    x = 45
    for rgb, label in legends:
        draw.ellipse((x, legend_y, x + 22, legend_y + 22), fill=rgb)
        draw.text((x + 32, legend_y - 2), label, font=font(19), fill="#354654")
        x += 235
    draw.text((800, legend_y - 2), "200 × 200 × 300 Å；x/y 周期，z 自由表面", font=font(18), fill="#526574")
    canvas.save(args.output, dpi=(180, 180))
    raw_render_path.unlink(missing_ok=True)
    start_selected_ag = (start["type"] == 2) & (start["position"][:, 2] >= 40.0) & (start["position"][:, 2] <= 200.0)
    print(
        f"saved {args.output}; atoms={len(end['id'])}; "
        f"initial Ag PTM_OTHER={int(np.count_nonzero(start_selected_ag & (start['ptm'] == int(PolyhedralTemplateMatchingModifier.Type.OTHER))))}; "
        f"Ag candidates={int(candidate.sum())}; fraction={candidate.sum()/max(int(selected_ag.sum()),1):.6f}"
    )
    scene.remove_from_scene()
    final_pipe.remove_from_scene()
    start_pipe.remove_from_scene()


if __name__ == "__main__":
    main()
