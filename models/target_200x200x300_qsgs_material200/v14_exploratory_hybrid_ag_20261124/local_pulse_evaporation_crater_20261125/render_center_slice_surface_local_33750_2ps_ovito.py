#!/workspace/software/ovito-lavapipe/ovito-python.sh
from pathlib import Path
from ovito.io import import_file
from ovito.modifiers import ExpressionSelectionModifier, DeleteSelectedModifier, AssignColorModifier
from ovito.vis import Viewport, TachyonRenderer
import ovito

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'
OUT.mkdir(exist_ok=True)
FILES=[
    ('脉冲前（1.700 ps）', ROOT.parent/'simple_ag_fixed_skeleton_20261125/reheat_melt_assessment_20261125/trajectories/ag_reheat_q6_ptm_msd_final.dump', 'center_slice_before.png'),
    ('表层局部脉冲 33750 K（2.0 ps）', ROOT/'trajectory_surface_local_33750_2ps/final_all_atoms_surface_local_33750.dump', 'center_slice_after_surface_local_33750_2ps.png'),
]
colors={1:(0.015,0.02,0.025),2:(0.015,0.02,0.025),3:(0.015,0.02,0.025),4:(0.53,0.78,0.91)}
for title,path,name in FILES:
    pipeline=import_file(str(path))
    pipeline.modifiers.append(ExpressionSelectionModifier(expression='abs(Position.Y - 100.0) > 5.0 || Position.X < 50.0 || Position.X > 150.0 || Position.Z < 170.0 || Position.Z > 280.0'))
    pipeline.modifiers.append(DeleteSelectedModifier())
    for typ,col in colors.items():
        pipeline.modifiers.append(ExpressionSelectionModifier(expression=f'ParticleType == {typ}'))
        pipeline.modifiers.append(AssignColorModifier(color=col))
    pipeline.modifiers.append(ExpressionSelectionModifier(expression='ParticleType == 4 && Position.Z > 220.0'))
    pipeline.modifiers.append(AssignColorModifier(color=(1.0,0.25,0.06)) )
    ovito.scene.pipelines.clear()
    ovito.scene.pipelines.append(pipeline)
    vp=Viewport(type=Viewport.Type.Ortho)
    vp.camera_pos=(100.0,-500.0,225.0)
    vp.camera_dir=(0.0,1.0,0.0)
    vp.camera_up=(0.0,0.0,1.0)
    vp.fov=65.0
    vp.render_image(filename=str(OUT/name),size=(1800,1200),background=(1,1,1),renderer=TachyonRenderer(ambient_occlusion=True))
    print(title, OUT/name)

# Crop empty viewport margins and compose a labeled Chinese comparison figure.
from PIL import Image, ImageDraw, ImageFont
import numpy as np

def tight_crop(im):
    a=np.asarray(im.convert('RGB'))
    mask=np.min(a,axis=2)<242
    yy,xx=np.where(mask)
    if len(xx)==0: return im
    pad=45
    box=(max(0,int(xx.min())-pad),max(0,int(yy.min())-pad),min(im.width,int(xx.max())+pad),min(im.height,int(yy.max())+pad))
    return im.crop(box)

before=tight_crop(Image.open(OUT/'center_slice_before.png'))
after=tight_crop(Image.open(OUT/'center_slice_after_surface_local_33750_2ps.png'))
# Scale both panels to a common height while preserving the OVITO render aspect ratio.
height=780
panels=[]
for im in (before,after):
    w=round(im.width*height/im.height)
    panels.append(im.resize((w,height),Image.Resampling.LANCZOS))
gap=30
header=175
canvas=Image.new('RGB',(panels[0].width+panels[1].width+gap, height+header),'white')
canvas.paste(panels[0],(0,header)); canvas.paste(panels[1],(panels[0].width+gap,header))
d=ImageDraw.Draw(canvas)
font_path='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
f_title=ImageFont.truetype(font_path,42)
f_sub=ImageFont.truetype(font_path,27)
f_legend=ImageFont.truetype(font_path,28)
d.text((40,18),'Ag–Ti3SiC2 表层局部脉冲剂量对比',font=f_title,fill=(20,20,20))
d.text((40,78),'沿 y=100 Å、厚 10 Å 切片；黑色：Ti3SiC2 骨架；浅蓝：Ag；显示 x=50–150 Å、z=170–280 Å 区域',font=f_sub,fill=(65,65,65))
d.text((panels[0].width//2-150,128),'脉冲前（1.700 ps）',font=f_legend,fill=(20,20,20))
d.text((panels[0].width+gap+panels[1].width//2-215,128),'表层局部脉冲 33750 K（2.0 ps）',font=f_legend,fill=(20,20,20))
canvas.save(OUT/'center_slice_comparison_surface_local_33750_2ps_Chinese.png',dpi=(300,300))
