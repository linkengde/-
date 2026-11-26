#!/usr/bin/env python3
"""Render the seed41 model after 1.8 ps of Ag-only NVT heating."""
from pathlib import Path
import json
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier
from ovito.pipeline import Pipeline, StaticSource
from ovito.vis import Viewport, TachyonRenderer, ParticlesVis
from ovito.nonpublic import ParticleRenderingQuality

root=Path(__file__).resolve().parent
full=root/'seed41_heat_final_1p8ps_all.dump'
traj=root/'seed41_heat_q6_msd_ag.dump'
report=json.loads((root/'seed41_ptm_q6_summary.json').read_text())
last=report[-1]

# PTM was computed on Ag only, matching the LAMMPS phase-analysis group.
ptm_pipe=import_file(str(traj))
ptm_pipe.modifiers.append(PolyhedralTemplateMatchingModifier())
ptm=ptm_pipe.compute(ptm_pipe.source.num_frames-1)
ag_ids=np.asarray(ptm.particles['Particle Identifier'],dtype=np.int64)
ag_struct=np.asarray(ptm.particles['Structure Type'],dtype=np.int32)
ptm_map=dict(zip(ag_ids.tolist(),ag_struct.tolist()))

data=import_file(str(full)).compute()
ids=np.asarray(data.particles['Particle Identifier'],dtype=np.int64)
types=np.asarray(data.particles['Particle Type'],dtype=np.int32)
colors=np.empty((len(ids),3),dtype=np.float32)
radii=np.empty(len(ids),dtype=np.float32)
# Black: rigid Ti3SiC2 geometry. Pale blue: residual FCC/HCP Ag.
# Orange: Ag not classified as FCC/HCP by PTM at the hot endpoint.
colors[types!=4]=(0.015,0.018,0.022); radii[types!=4]=0.85
colors[types==4]=(0.96,0.43,0.08); radii[types==4]=1.35
residual=np.fromiter((ptm_map.get(int(i),0) in (1,2) for i in ids),dtype=bool,count=len(ids)) & (types==4)
colors[residual]=(0.55,0.80,0.92)
data.particles_.create_property('Color',data=colors)
data.particles_.create_property('Radius',data=radii)
data.particles.vis.shape=ParticlesVis.Shape.Sphere
data.particles.vis.periodic_image_mode=ParticlesVis.PeriodicImageMode.Off
data.particles.vis.rendering_quality=ParticleRenderingQuality.HighQuality
data.cell.vis.render_cell=True
data.cell.vis.rendering_color=(0.30,0.34,0.38)
data.cell.vis.line_width=2.0
scene=Pipeline(source=StaticSource(data=data)); scene.add_to_scene()
raw=root/'seed41_heat_final_raw.png'
size=(1350,1200)
vp=Viewport(type=Viewport.Type.Perspective)
center=np.array([100.,100.,250.]); camera=np.array([550.,-650.,480.])
vp.camera_pos=camera; vp.camera_dir=center-camera; vp.camera_up=(0.,0.,1.)
vp.zoom_all(size=size)
renderer=TachyonRenderer(); renderer.ambient_occlusion=False; renderer.antialiasing=True
renderer.antialiasing_samples=2; renderer.shadows=False
vp.render_image(filename=str(raw),size=size,renderer=renderer,background=(1,1,1),crop=False)
scene.remove_from_scene()

font_path='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
font_bold='/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc'
base=Image.open(raw).convert('RGB'); raw.unlink(missing_ok=True)
sidebar,header=590,195
canvas=Image.new('RGB',(sidebar+base.width,header+base.height),'white')
d=ImageDraw.Draw(canvas)
title=ImageFont.truetype(font_bold,42); body=ImageFont.truetype(font_path,25)
legend=ImageFont.truetype(font_path,27); note=ImageFont.truetype(font_path,20)
d.text((28,14),'Ag–Ti₃SiC₂ seed41 加热终态',font=title,fill=(20,23,26))
d.text((31,76),'t = 1.800 ps｜Ag NVT 均匀控温',font=body,fill=(48,53,58))
d.text((31,113),'材料 200 × 200 × 200 Å｜计算盒 200 × 200 × 500 Å',font=body,fill=(48,53,58))
d.text((31,150),f"{int(last['n_Ag']):,} 个 Ag｜Ag 温度约 2407 K｜MSD {7.4063:.2f} Å²",font=body,fill=(48,53,58))
canvas.paste(base,(sidebar,header))
x0,y0=30,header+140
items=[((0.015,0.018,0.022),'Ti₃SiC₂ 固定骨架  358,428'),((0.96,0.43,0.08),f"Ag 高温失序候选  {int(last['n_Ag'])-int(last['fcc'])-int(last['hcp']):,}"),((0.55,0.80,0.92),f"Ag 残余 FCC/HCP  {int(last['fcc'])+int(last['hcp']):,}")]
for j,(rgb,label) in enumerate(items):
    y=y0+58*j
    d.rectangle((x0,y,x0+34,y+26),fill=tuple(round(v*255) for v in rgb))
    d.text((x0+48,y-5),label,font=legend,fill=(28,32,36))
d.text((30,y0+200),f"PTM FCC/HCP：{last['ordered_fcc_hcp_pct']:.3f}%",font=note,fill=(54,59,64))
d.text((30,y0+232),'橙色表示 PTM 未识别为 FCC/HCP；结合 Q6、MSD 判断高温失序。',font=note,fill=(54,59,64))
d.text((30,y0+264),'骨架固定作几何障碍；Ag–骨架为探索性排斥项，不代表真实粘附。',font=note,fill=(54,59,64))
out=root/'AgTi3SiC2_seed41_heated_1p8ps_OVITO_CN_300dpi.png'
canvas.save(out,dpi=(300,300))
print(out)
