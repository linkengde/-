#!/usr/bin/env python3
"""Render the final three-pulse seed41 state with OVITO PTM coloring."""
from pathlib import Path
import json
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier,ExpressionSelectionModifier,DeleteSelectedModifier
from ovito.pipeline import Pipeline,StaticSource
from ovito.vis import Viewport,TachyonRenderer,ParticlesVis
from ovito.nonpublic import ParticleRenderingQuality

ROOT=Path(__file__).resolve().parent
FULL=ROOT/'seed41_afterpulse3_final_14p8ps_nomsd_all.dump'
TRAJ=ROOT/'seed41_afterpulse3_NVE_ag_nomsd.dump'
SUMMARY=json.loads((ROOT/'seed41_pulse3_summary.json').read_text())
FONT='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
BOLD='/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc'

ptm_pipe=import_file(str(TRAJ))
ptm_pipe.modifiers.append(PolyhedralTemplateMatchingModifier())
ptm=ptm_pipe.compute(ptm_pipe.source.num_frames-1)
ptm_ids=np.asarray(ptm.particles['Particle Identifier'],dtype=np.int64)
ptm_struct=np.asarray(ptm.particles['Structure Type'],dtype=np.int32)
ptm_map=dict(zip(ptm_ids.tolist(),ptm_struct.tolist()))

def styled_data(data):
    ids=np.asarray(data.particles['Particle Identifier'],dtype=np.int64)
    typ=np.asarray(data.particles['Particle Type'],dtype=np.int32)
    z=np.asarray(data.particles.positions)[:,2]
    colors=np.empty((len(ids),3),dtype=np.float32)
    radii=np.empty(len(ids),dtype=np.float32)
    colors[typ!=4]=(0.012,0.016,0.020); radii[typ!=4]=0.82
    is_ag=typ==4
    colors[is_ag]=(0.96,0.43,0.08); radii[is_ag]=1.15
    ordered=np.fromiter((ptm_map.get(int(i),0) in (1,2) for i in ids),dtype=bool,count=len(ids))
    crystalline=is_ag&ordered&(z<=220)
    disordered=is_ag&(~ordered)&(z<=220)
    plume=is_ag&(z>220)
    colors[crystalline]=(0.56,0.80,0.92)
    colors[disordered]=(0.96,0.43,0.08)
    colors[plume]=(0.83,0.16,0.13)
    radii[plume]=1.30
    data.particles_.create_property('Color',data=colors)
    data.particles_.create_property('Radius',data=radii)
    data.particles.vis.shape=ParticlesVis.Shape.Sphere
    data.particles.vis.periodic_image_mode=ParticlesVis.PeriodicImageMode.Off
    data.particles.vis.rendering_quality=ParticleRenderingQuality.HighQuality
    return int(np.count_nonzero(typ==4)),int(np.count_nonzero(crystalline)),int(np.count_nonzero(disordered)),int(np.count_nonzero(plume)),int(np.count_nonzero(typ!=4))

def render(pipeline,filename,center,camera,cell,zoom=True):
    data=pipeline.compute()
    counts=styled_data(data)
    data.cell.vis.render_cell=cell
    data.cell.vis.rendering_color=(0.30,0.34,0.38); data.cell.vis.line_width=2.0
    scene=Pipeline(source=StaticSource(data=data)); scene.add_to_scene()
    raw=ROOT/(Path(filename).stem+'_raw.png'); size=(1350,1180)
    vp=Viewport(type=Viewport.Type.Perspective); vp.camera_pos=camera; vp.camera_dir=np.asarray(center)-np.asarray(camera); vp.camera_up=(0.,0.,1.)
    if zoom: vp.zoom_all(size=size)
    renderer=TachyonRenderer(); renderer.ambient_occlusion=False; renderer.antialiasing=True; renderer.antialiasing_samples=2; renderer.shadows=False
    vp.render_image(filename=str(raw),size=size,renderer=renderer,background=(1,1,1),crop=False)
    scene.remove_from_scene()
    base=Image.open(raw).convert('RGB'); raw.unlink(missing_ok=True)
    canvas=Image.new('RGB',(590+base.width,190+base.height),'white')
    d=ImageDraw.Draw(canvas); title=ImageFont.truetype(BOLD,40); body=ImageFont.truetype(FONT,24); legend=ImageFont.truetype(FONT,25); note=ImageFont.truetype(FONT,18)
    heading='Ag–Ti3SiC2 中心局部结构' if 'closeup' in Path(filename).stem else 'Ag–Ti3SiC2 三次局域脉冲后结构'
    d.text((28,12),heading,font=title,fill=(22,25,29))
    d.text((30,72),'t = 14.800 ps｜seed41｜累计输入 4,437.75 eV',font=body,fill=(48,53,58))
    scope='中心截取' if 'closeup' in Path(filename).stem else '盒子 200 × 200 × 500 Å'
    d.text((30,109),f"{scope}｜{len(data.particles):,} 个原子｜Ag {counts[0]:,}",font=body,fill=(48,53,58))
    d.text((30,146),'骨架固定；Ag 继续自由运动；无 TTM',font=body,fill=(48,53,58))
    canvas.paste(base,(590,190))
    y0=410
    items=[((0.012,0.016,0.020),f'Ti3SiC2 固定骨架  {counts[4]:,}'),((0.56,0.80,0.92),f'Ag PTM FCC/HCP  {counts[1]:,}'),((0.96,0.43,0.08),f'Ag 未识别为 FCC/HCP  {counts[2]:,}'),((0.83,0.16,0.13),f'Ag 位于 z > 220 Å  {counts[3]:,}')]
    for j,(rgb,label) in enumerate(items):
        y=y0+58*j; d.rectangle((30,y,65,y+25),fill=tuple(round(v*255) for v in rgb)); d.text((80,y-4),label,font=legend,fill=(25,29,33))
    d.text((30,y0+252),'橙色只表示 PTM 未识别为 FCC/HCP，',font=note,fill=(58,62,66))
    d.text((30,y0+280),'需结合 Q6、MSD 判断无序；不单独等同熔化。',font=note,fill=(58,62,66))
    d.text((30,y0+321),'红色为原始表面 z=200 Å 上方 Ag，',font=note,fill=(58,62,66))
    d.text((30,y0+349),'不代表这些原子都已蒸发。',font=note,fill=(58,62,66))
    d.text((30,y0+395),'Ag–骨架为排斥性交叉项；脉冲为速度注能，',font=note,fill=(58,62,66))
    d.text((30,y0+423),'结果仅作定性探索。',font=note,fill=(58,62,66))
    canvas.save(filename,dpi=(300,300))
    print(filename,counts)

full_pipe=import_file(str(FULL))
render(full_pipe,ROOT/'seed41_three_pulse_final_14p8ps_OVITO_CN_300dpi.png',center=(100,100,250),camera=(520,-650,450),cell=True)

zoom_pipe=import_file(str(FULL))
zoom_pipe.modifiers.append(ExpressionSelectionModifier(expression='Position.Z < 135 || Position.Z > 250 || Position.X < 50 || Position.X > 150 || Position.Y < 50 || Position.Y > 150'))
zoom_pipe.modifiers.append(DeleteSelectedModifier())
render(zoom_pipe,ROOT/'seed41_three_pulse_crater_closeup_14p8ps_OVITO_CN_300dpi.png',center=(100,100,190),camera=(180,-80,260),cell=False,zoom=False)
