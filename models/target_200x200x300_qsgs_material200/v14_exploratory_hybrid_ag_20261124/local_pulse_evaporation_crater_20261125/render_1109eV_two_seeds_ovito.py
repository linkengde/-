#!/workspace/software/ovito-lavapipe/ovito-python.sh
from pathlib import Path
from ovito.io import import_file
from ovito.modifiers import ExpressionSelectionModifier, DeleteSelectedModifier, AssignColorModifier
from ovito.vis import Viewport, TachyonRenderer
import ovito

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
FILES=[
 ('脉冲前（1.700 ps）',ROOT.parent/'simple_ag_fixed_skeleton_20261125/reheat_melt_assessment_20261125/trajectories/ag_reheat_q6_ptm_msd_final.dump','slice_before_1109eV.png'),
 ('1.109 keV，seed 20261131（2 ps）',ROOT/'trajectory_surface_local_1109eV_seed20261131_2ps/final_all_atoms_1109eV_seed20261131.dump','slice_after_1109eV_seed20261131.png'),
 ('1.109 keV，seed 20261132（2 ps）',ROOT/'trajectory_surface_local_1109eV_seed20261132_2ps/final_all_atoms_1109eV_seed20261132.dump','slice_after_1109eV_seed20261132.png'),
]
for title,path,name in FILES:
    pipeline=import_file(str(path))
    pipeline.modifiers.append(ExpressionSelectionModifier(expression='abs(Position.Y-100.0)>5.0 || Position.X<50.0 || Position.X>150.0 || Position.Z<170.0 || Position.Z>300.0'))
    pipeline.modifiers.append(DeleteSelectedModifier())
    for typ in (1,2,3):
        pipeline.modifiers.append(ExpressionSelectionModifier(expression=f'ParticleType=={typ}'))
        pipeline.modifiers.append(AssignColorModifier(color=(0.015,0.02,0.025)))
    pipeline.modifiers.append(ExpressionSelectionModifier(expression='ParticleType==4'))
    pipeline.modifiers.append(AssignColorModifier(color=(0.53,0.78,0.91)))
    pipeline.modifiers.append(ExpressionSelectionModifier(expression='ParticleType==4 && Position.Z>220.0'))
    pipeline.modifiers.append(AssignColorModifier(color=(1.0,0.25,0.06)))
    ovito.scene.pipelines.clear(); ovito.scene.pipelines.append(pipeline)
    vp=Viewport(type=Viewport.Type.Ortho)
    vp.camera_pos=(100.0,-500.0,230.0); vp.camera_dir=(0.0,1.0,0.0); vp.camera_up=(0.0,0.0,1.0); vp.fov=65.0
    vp.render_image(filename=str(OUT/name),size=(1500,1050),background=(1,1,1),renderer=TachyonRenderer(ambient_occlusion=True))
    print(title,OUT/name)

from PIL import Image,ImageDraw,ImageFont
import numpy as np
def crop(im):
    a=np.asarray(im.convert('RGB')); yy,xx=np.where(np.min(a,axis=2)<242)
    if not len(xx): return im
    p=35; return im.crop((max(0,int(xx.min())-p),max(0,int(yy.min())-p),min(im.width,int(xx.max())+p),min(im.height,int(yy.max())+p)))
names=[x[2] for x in FILES]
panels=[crop(Image.open(OUT/n)) for n in names]
height=720
panels=[im.resize((round(im.width*height/im.height),height),Image.Resampling.LANCZOS) for im in panels]
gap=20; header=165
canvas=Image.new('RGB',(sum(im.width for im in panels)+2*gap,height+header),'white')
x=0
for im in panels: canvas.paste(im,(x,header)); x+=im.width+gap
d=ImageDraw.Draw(canvas); font='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
f1=ImageFont.truetype(font,39); f2=ImageFont.truetype(font,24); f3=ImageFont.truetype(font,24)
d.text((30,12),'Ag–Ti3SiC2：1.109 keV 局部脉冲随机种子对照',font=f1,fill=(20,20,20))
d.text((30,66),'OVITO 中心截面：y=100±5 Å，x=50–150 Å，z=170–300 Å；黑色：骨架，浅蓝：Ag，红色：z>220 Å 的 Ag',font=f2,fill=(60,60,60))
labels=['脉冲前（1.700 ps）','seed 20261131（2.0 ps）','seed 20261132（2.0 ps）']
x=0
for im,label in zip(panels,labels):
    d.text((x+im.width//2-145,122),label,font=f3,fill=(30,30,30)); x+=im.width+gap
canvas.save(OUT/'center_slice_1109eV_two_seeds_Chinese.png',dpi=(300,300))
print(OUT/'center_slice_1109eV_two_seeds_Chinese.png')
