#!/usr/bin/env python3
"""Render the dense rebuilt candidate in the user's reference layout."""
from pathlib import Path
import json
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from ovito.io import import_file
from ovito.pipeline import Pipeline,StaticSource
from ovito.vis import Viewport,TachyonRenderer,ParticlesVis
from ovito.nonpublic import ParticleRenderingQuality

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'qsgs_thin_1A_seed41_200x200x500.data'
REPORT=json.loads((ROOT/'qsgs_thin_1A_seed41_report.json').read_text())
OUT=ROOT/'AgTi3SiC2_更细骨架_seed41_OVITO_CN_300dpi.png'
RAW=ROOT/'AgTi3SiC2_更细骨架_seed41_raw.png'

data=import_file(str(DATA)).compute()
types=np.asarray(data.particles['Particle Type'],dtype=int)
n_all=len(types); n_ag=int(np.count_nonzero(types==4)); n_tsc=n_all-n_ag
colors=np.empty((len(types),3)); colors[types!=4]=(0.015,0.018,0.022); colors[types==4]=(0.53,0.78,0.91)
data.particles_.create_property('Color',data=colors)
# Atomic visualization radii close the apparent gaps between neighboring sites.
data.particles_.create_property('Radius',data=np.where(types==4,1.45,0.85))
data.particles.vis.shape=ParticlesVis.Shape.Sphere
data.particles.vis.periodic_image_mode=ParticlesVis.PeriodicImageMode.Off
data.particles.vis.rendering_quality=ParticleRenderingQuality.HighQuality
data.cell.vis.render_cell=True; data.cell.vis.rendering_color=(0.32,0.36,0.40); data.cell.vis.line_width=2.0
scene=Pipeline(source=StaticSource(data=data)); scene.add_to_scene()
size=(1200,1100); center=np.array([100.,100.,150.]); camera=np.array([520.,-590.,420.])
vp=Viewport(type=Viewport.Type.Perspective); vp.camera_pos=camera; vp.camera_dir=center-camera; vp.camera_up=(0.,0.,1.)
vp.zoom_all(size=size)
renderer=TachyonRenderer(); renderer.ambient_occlusion=False; renderer.antialiasing=True
renderer.antialiasing_samples=2; renderer.shadows=False
vp.render_image(filename=str(RAW),size=size,renderer=renderer,background=(1,1,1),crop=False)
scene.remove_from_scene()

font_path='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
font_bold='/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc'
raw=Image.open(RAW).convert('RGB'); RAW.unlink(missing_ok=True)
sidebar,header=560,185
canvas=Image.new('RGB',(sidebar+20+raw.width,header+raw.height),'white')
draw=ImageDraw.Draw(canvas)
title=ImageFont.truetype(font_bold,44); body=ImageFont.truetype(font_path,26)
legend=ImageFont.truetype(font_path,28); note=ImageFont.truetype(font_path,20)
draw.text((30,14),'Ag–Ti3SiC2 更细枝骨架模型',font=title,fill=(20,23,26))
draw.text((33,77),f"全新随机 QSGS｜seed {REPORT['QSGS']['seed']}",font=body,fill=(48,53,58))
draw.text((33,116),'材料 200 × 200 × 200 Å｜计算盒 200 × 200 × 500 Å',font=body,fill=(48,53,58))
draw.text((33,154),f'{n_all:,} 个原子｜Ag/Ti3SiC2 质量比 {REPORT["Ag_over_Ti3SiC2_mass_ratio"]:.3f}:1',font=body,fill=(48,53,58))
canvas.paste(raw,(sidebar+20,header))
n_skeleton=n_all-n_ag
x0,y0=32,header+165
items=[((0.015,0.018,0.022),f'Ti3SiC2 骨架  {n_skeleton:,}'),((0.53,0.78,0.91),f'Ag 银相  {n_ag:,}')]
for i,(rgb,label) in enumerate(items):
    y=y0+76*i
    draw.rectangle((x0,y,x0+36,y+28),fill=tuple(round(v*255) for v in rgb))
    draw.text((x0+53,y-5),label,font=legend,fill=(25,29,32))
draw.text((32,y0+185),'两相主连通网络均贯通 x、y、z',font=note,fill=(60,65,70))
draw.text((32,y0+219),f'顶部骨架最高占比 {REPORT["top_surface_z180_200A"]["max_20A_xy_bin"]:.1%}；无 >90% 分区',font=note,fill=(60,65,70))
draw.text((32,y0+253),f'20 Å 分区最少 {REPORT["material_20A_cube_atom_counts"]["min"]} 个原子，无空区',font=note,fill=(60,65,70))
draw.text((32,y0+287),'几何候选，尚未弛豫或 MD 验证',font=note,fill=(60,65,70))
draw.text((32,y0+335),'z 分层（Å）',font=legend,fill=(35,39,43))
draw.text((32,y0+380),'固定层：0–20（仅骨架）',font=note,fill=(60,65,70))
draw.text((32,y0+414),'恒温层：20–40（仅骨架）',font=note,fill=(60,65,70))
draw.text((32,y0+448),'牛顿层：40–200',font=note,fill=(60,65,70))
draw.text((32,y0+482),'真空区：200–500',font=note,fill=(60,65,70))
draw.text((32,y0+532),'底部固定/恒温仅作用于骨架；Ag 可动',font=note,fill=(60,65,70))
canvas.save(OUT,dpi=(300,300))
print(f'atoms={n_all} TSC={n_tsc} Ag={n_ag} output={OUT}',flush=True)
