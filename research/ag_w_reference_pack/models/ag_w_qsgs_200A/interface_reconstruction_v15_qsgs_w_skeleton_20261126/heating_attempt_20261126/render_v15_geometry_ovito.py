#!/usr/bin/env python3
"""Render the full V15 QSGS atomistic geometry, colored by element."""
from pathlib import Path
import argparse
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from ovito.io import import_file
from ovito.pipeline import Pipeline, StaticSource
from ovito.vis import ParticlesVis, TachyonRenderer, Viewport


def font(size):
    for path in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--snapshot",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    src=import_file(str(a.snapshot))
    data=src.compute(); types=np.asarray(data.particles["Particle Type"],dtype=int)
    colors=np.empty((len(types),3),float)
    colors[types==1]=(0.035,0.045,0.055) # black W
    colors[types==2]=(0.66,0.82,0.89) # pale blue Ag
    data.particles_.create_property("Color",data=colors)
    data.particles_.vis.radius=0.88
    data.particles_.vis.shape=ParticlesVis.Shape.Sphere
    data.cell.vis.enabled=True; data.cell.vis.render_cell=True
    data.cell.vis.line_width=1.2; data.cell.vis.rendering_color=(0.18,0.20,0.22)
    scene=Pipeline(source=StaticSource(data=data)); scene.add_to_scene()
    size=(1700,1250); vp=Viewport(type=Viewport.Type.Perspective)
    vp.camera_pos=(229.1,-122.8,260.0); vp.camera_dir=(-0.4609,0.7957,-0.3930)
    vp.camera_up=(0.0,0.0,1.0); vp.fov=36.0; vp.zoom_all(size=size)
    vp.camera_pos=(229.1,-122.8,260.0); vp.camera_dir=(-0.4609,0.7957,-0.3930); vp.fov=36.0
    renderer=TachyonRenderer(); renderer.antialiasing=True; renderer.antialiasing_samples=6
    renderer.ambient_occlusion=True; renderer.ambient_occlusion_samples=8; renderer.shadows=True
    rendered=vp.render_image(size=size,background=(1,1,1),alpha=False,renderer=renderer,stop_on_error=True)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    tmp=a.output.with_name(a.output.stem+"_render.png"); rendered.save(str(tmp))
    canvas=Image.new("RGB",(size[0],size[1]+120),"white")
    with Image.open(tmp) as im: canvas.paste(im.convert("RGB"),(0,120))
    tmp.unlink(missing_ok=True)
    d=ImageDraw.Draw(canvas); d.text((40,18),"Ag–W V15｜四参数随机 QSGS W 骨架 + Ag 填充",font=font(34),fill="#1d2b36")
    d.ellipse((48,75,70,97),fill="#090c0e"); d.text((80,72),"W 骨架",font=font(20),fill="#354654")
    d.ellipse((250,75,272,97),fill="#a8d1e3"); d.text((282,72),"Ag",font=font(20),fill="#354654")
    d.text((510,72),"未加热结构；完整 200 × 200 × 300 Å 盒子",font=font(19),fill="#526574")
    canvas.save(a.output,dpi=(180,180))
    print(f"saved {a.output}; atoms={len(types)}; W={int(np.count_nonzero(types==1))}; Ag={int(np.count_nonzero(types==2))}")
    scene.remove_from_scene(); src.remove_from_scene()

if __name__=="__main__": main()
