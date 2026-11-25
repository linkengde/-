#!/workspace/software/ovito-lavapipe/ovito-python.sh
from pathlib import Path
from PIL import Image,ImageDraw,ImageFont
from ovito.io import import_file
from ovito.modifiers import ExpressionSelectionModifier, DeleteSelectedModifier, AssignColorModifier
from ovito.vis import Viewport, TachyonRenderer
import ovito

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
RUNS={
 20261133:[('脉冲后2 ps（总时刻3.7 ps）',ROOT/'trajectory_surface_local_1479eV_seed20261133_2ps/final_all_atoms_1479eV_seed20261133.dump'),
          ('脉冲后4 ps（总时刻5.7 ps）',ROOT/'trajectory_surface_local_1479eV_seed20261133_extension2ps/final_all_atoms_1479eV_seed20261133_4ps.dump')],
 20261134:[('脉冲后2 ps（总时刻3.7 ps）',ROOT/'trajectory_surface_local_1479eV_seed20261134_2ps/final_all_atoms_1479eV_seed20261134.dump'),
          ('脉冲后4 ps（总时刻5.7 ps）',ROOT/'trajectory_surface_local_1479eV_seed20261134_extension2ps/final_all_atoms_1479eV_seed20261134_4ps.dump')]
}
font='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
f1=ImageFont.truetype(font,34); f2=ImageFont.truetype(font,23)
for seed,items in RUNS.items():
    panels=[]
    for title,path in items:
        p=import_file(str(path))
        p.modifiers.append(ExpressionSelectionModifier(expression='abs(Position.Y-100.0)>5.0 || Position.X<50.0 || Position.X>150.0 || Position.Z<170.0 || Position.Z>300.0'))
        p.modifiers.append(DeleteSelectedModifier())
        for typ in (1,2,3):
            p.modifiers.append(ExpressionSelectionModifier(expression=f'ParticleType=={typ}'))
            p.modifiers.append(AssignColorModifier(color=(0.015,0.02,0.025)))
        p.modifiers.append(ExpressionSelectionModifier(expression='ParticleType==4'))
        p.modifiers.append(AssignColorModifier(color=(0.53,0.78,0.91)))
        p.modifiers.append(ExpressionSelectionModifier(expression='ParticleType==4 && Position.Z>220.0'))
        p.modifiers.append(AssignColorModifier(color=(1.0,0.25,0.06)))
        ovito.scene.pipelines.clear(); ovito.scene.pipelines.append(p)
        vp=Viewport(type=Viewport.Type.Ortho)
        vp.camera_pos=(100.0,-500.0,230.0); vp.camera_dir=(0.0,1.0,0.0); vp.camera_up=(0.0,0.0,1.0); vp.fov=65.0
        name=OUT/f'extension_seed{seed}_{len(panels)+1}.png'
        vp.render_image(filename=str(name),size=(1500,1050),background=(1,1,1),renderer=TachyonRenderer(ambient_occlusion=True))
        im=Image.open(name).convert('RGB')
        import numpy as np
        a=np.asarray(im); yy,xx=np.where(np.min(a,axis=2)<242)
        if len(xx):
            pad=25; im=im.crop((max(0,int(xx.min())-pad),max(0,int(yy.min())-pad),min(im.width,int(xx.max())+pad),min(im.height,int(yy.max())+pad)))
        panels.append((title,im))
    h=780; gap=24; head=150
    panels=[(t,im.resize((round(im.width*h/im.height),h),Image.Resampling.LANCZOS)) for t,im in panels]
    canvas=Image.new('RGB',(sum(im.width for _,im in panels)+gap,h+head),'white'); d=ImageDraw.Draw(canvas)
    d.text((24,8),f'Ag–Ti3SiC2：1.479 keV 延时对照，seed {seed}',font=f1,fill=(20,20,20))
    d.text((24,56),'OVITO中心截面 y=100±5 Å、x=50–150 Å、z=170–300 Å；黑色：骨架，浅蓝：Ag，红色：z>220 Å 的 Ag',font=f2,fill=(55,55,55))
    x=0
    for title,im in panels:
        canvas.paste(im,(x,head)); d.text((x+im.width//2-160,110),title,font=f2,fill=(35,35,35)); x+=im.width+gap
    path=OUT/f'center_slice_1479eV_extension_seed{seed}_Chinese.png'
    canvas.save(path,dpi=(300,300)); print(path)
