#!/usr/bin/env python3
"""Render the complete V14 atomistic box and mark exploratory Ag disorder candidates."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier, SteinhardtParametersModifier
from ovito.pipeline import Pipeline, StaticSource
from ovito.vis import ParticlesVis, TachyonRenderer, Viewport


def read_fields(path: Path):
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
    options = [
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc" if bold else "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for item in options:
        if Path(item).exists():
            return ImageFont.truetype(item, size)
    return ImageFont.load_default()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--start", type=Path, required=True)
    ap.add_argument("--snapshot", type=Path, required=True)
    ap.add_argument("--time-ps", type=float, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--q6-cutoff", type=float, default=0.40)
    ap.add_argument("--msd-cutoff-A", type=float, default=2.0)
    args = ap.parse_args()

    source_pipeline, data, end = read_fields(args.snapshot)
    _, _, start = read_fields(args.start)
    order0, order1 = np.argsort(start["id"]), np.argsort(end["id"])
    for fields in (start, end):
        ix = np.argsort(fields["id"])
        for name in fields:
            fields[name] = fields[name][ix]
    if not np.array_equal(start["id"], end["id"]):
        raise RuntimeError("Atom IDs differ between start and final snapshot")
    delta = end["position"] - start["position"]
    delta[:, :2] -= 200.0 * np.rint(delta[:, :2] / 200.0)
    displacement = np.linalg.norm(delta, axis=1)
    # Particle type 1=W, 2=Ag; select the original mapped material region.
    ag = (start["type"] == 2) & (start["position"][:, 2] >= 40.0) & (start["position"][:, 2] <= 200.0)
    other = end["ptm"] == int(PolyhedralTemplateMatchingModifier.Type.OTHER)
    low_q6 = end["q6"] < args.q6_cutoff
    candidate = ag & other & low_q6 & (displacement > args.msd_cutoff_A)
    colors = np.empty((len(end["id"]), 3), dtype=float)
    colors[end["type"] == 1] = (0.09, 0.18, 0.25)   # W scaffold
    colors[end["type"] == 2] = (0.66, 0.82, 0.89)   # Ag phase
    colors[candidate] = (0.97, 0.42, 0.10)          # exploratory candidates

    data.particles_.create_property("Color", data=colors)
    data.particles_.vis.radius = 0.90
    data.particles_.vis.shape = ParticlesVis.Shape.Sphere
    data.cell.vis.enabled = True
    data.cell.vis.render_cell = True
    data.cell.vis.line_width = 1.25
    data.cell.vis.rendering_color = (0.15, 0.18, 0.20)
    render_pipeline = Pipeline(source=StaticSource(data=data))
    render_pipeline.add_to_scene()

    panel = (1900, 1450)
    viewport = Viewport(type=Viewport.Type.Perspective)
    viewport.camera_pos = (423.6, -457.0, 425.1)
    viewport.camera_dir = (-0.4609, 0.7957, -0.3930)
    viewport.camera_up = (0.0, 0.0, 1.0)
    viewport.fov = 36.0
    # zoom_all initializes the scene's camera target; restore an explicit
    # distance/FOV afterward so the entire atomistic box fills the frame.
    viewport.zoom_all(size=panel)
    viewport.camera_pos = (423.6, -457.0, 425.1)
    viewport.camera_dir = (-0.4609, 0.7957, -0.3930)
    viewport.fov = 36.0
    renderer = TachyonRenderer()
    renderer.antialiasing = True
    renderer.antialiasing_samples = 8
    renderer.ambient_occlusion = True
    renderer.ambient_occlusion_samples = 12
    renderer.shadows = True
    rendered = viewport.render_image(size=panel, background=(1.0, 1.0, 1.0), alpha=False,
                                     renderer=renderer, stop_on_error=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    banner_h = 210
    canvas = Image.new("RGB", (panel[0], panel[1] + banner_h), "white")
    raw_render_path = args.output.with_name(args.output.stem + "_ovito_render.png")
    rendered.save(str(raw_render_path))
    with Image.open(raw_render_path) as raw_render:
        canvas.paste(raw_render.convert("RGB"), (0, banner_h))
    draw = ImageDraw.Draw(canvas)
    draw.text((52, 22), f"Ag–W V14 全模型｜{args.time_ps:g} ps", font=font(42, True), fill="#1d2b36")
    draw.text((55, 82), "高剂量探索末态；橙色仅为 Q6 + PTM + MSD 无序候选，并非已确认液相",
              font=font(26), fill="#526574")
    legend_y = 145
    legends = [
        ((23, 46, 64), "W 骨架"),
        ((168, 209, 227), "Ag"),
        ((247, 107, 26), f"Ag 无序候选 {int(candidate.sum()):,}"),
    ]
    x = 60
    for rgb, label in legends:
        draw.ellipse((x, legend_y, x + 26, legend_y + 26), fill=rgb)
        draw.text((x + 38, legend_y - 1), label, font=font(21), fill="#354654")
        x += 260
    draw.text((1120, legend_y - 1), "盒子：200 × 200 × 300 Å；x/y 周期，z 自由表面",
              font=font(19), fill="#526574")
    canvas.save(args.output, dpi=(180, 180))
    raw_render_path.unlink(missing_ok=True)
    print(f"saved {args.output}; atoms={len(end['id'])}; candidate={int(candidate.sum())}; fraction={candidate.sum()/max(int(ag.sum()), 1):.4f}")
    render_pipeline.remove_from_scene()
    source_pipeline.remove_from_scene()


if __name__ == "__main__":
    main()
