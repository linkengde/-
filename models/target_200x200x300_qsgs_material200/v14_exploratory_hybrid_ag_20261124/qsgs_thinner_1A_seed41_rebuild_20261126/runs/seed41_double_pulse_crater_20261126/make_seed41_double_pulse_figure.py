#!/usr/bin/env python3
"""Create Chinese diagnostic plots for the second seed41 local pulse."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from analyze_seed41_double_pulse import dump_frames

ROOT=Path(__file__).resolve().parent
font=FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
bold=FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc')

def surface_map(path,bw=4.):
    fr=dump_frames(path)[0]; idx=fr[1]; a=fr[2]
    x=a[:,idx['x']]-100.; y=a[:,idx['y']]-100.; z=a[:,idx['z']]
    keep=(np.abs(x)<50)&(np.abs(y)<50)&(z>=150)&(z<=230)
    n=int(100/bw); m=np.full((n,n),np.nan)
    ix=np.floor((x[keep]+50)/bw).astype(int); iy=np.floor((y[keep]+50)/bw).astype(int)
    for bx,by,h in zip(ix,iy,z[keep]):
        if 0<=bx<n and 0<=by<n and (np.isnan(m[by,bx]) or h>m[by,bx]): m[by,bx]=h
    return m

def read_csv(name):
    with (ROOT/name).open() as f: return list(csv.DictReader(f))

paths=[ROOT.parent/'seed41_local_pulse_crater_20261126'/'seed41_local_pulse_final_5p8ps_all.dump',
       ROOT/'seed41_pulse2_final_7p8ps_all.dump',
       ROOT/'seed41_afterpulse2_final_10p8ps_nomsd_all.dump']
maps=[surface_map(p) for p in paths]
labels=['脉冲2前 5.8 ps','脉冲2后 7.8 ps','无源松弛后 10.8 ps']
traj=read_csv('seed41_double_pulse_trajectory.csv')
qmsd=read_csv('seed41_double_pulse_q6_msd.csv')
t=np.array([float(r['time_ps']) for r in traj])
cuts=['n_z_gt_220','n_z_gt_250','n_z_gt_300','n_z_gt_400']
cutlabels=['z > 220 Å','z > 250 Å','z > 300 Å','z > 400 Å']
colors=['#5d83a5','#71819b','#e78a28','#d84c34']
fig=plt.figure(figsize=(17,9),constrained_layout=True)
gs=fig.add_gridspec(2,3)
map_axes=[fig.add_subplot(gs[0,i]) for i in range(3)]
axp=fig.add_subplot(gs[1,0]); axm=fig.add_subplot(gs[1,1]); axq=fig.add_subplot(gs[1,2])
for ax,m,label in zip(map_axes,maps,labels):
    im=ax.imshow(np.ma.masked_invalid(m),origin='lower',extent=(-50,50,-50,50),vmin=170,vmax=210,cmap='Blues',interpolation='nearest')
    ax.add_patch(plt.Circle((0,0),25,fill=False,ls='--',lw=1.2,ec='#e78a28'))
    ax.set_title(label,fontproperties=bold,fontsize=12)
    ax.set_xlabel('相对 x (Å)',fontproperties=font); ax.set_ylabel('相对 y (Å)',fontproperties=font)
cb=fig.colorbar(im,ax=map_axes,shrink=.85,pad=.02); cb.set_label('Ag 局部顶面高度 z (Å)',fontproperties=font)
for key,label,col in zip(cuts,cutlabels,colors):
    axp.plot(t,[int(r[key]) for r in traj],label=label,color=col,lw=1.8)
axp.axvline(7.8,color='#555',ls='--',lw=1); axp.axvline(10.375,color='#bd3030',ls=':',lw=1)
axp.set_title('真空区 Ag 数量（当前盒内）',fontproperties=bold,fontsize=12)
axp.set_xlabel('时间 (ps)',fontproperties=font); axp.set_ylabel('Ag 原子数',fontproperties=font)
axp.legend(prop=font,frameon=False,ncol=2)
qt=np.array([float(r['time_ps']) for r in qmsd])
axm.plot(qt,[float(r['MSD_Ag_A2']) for r in qmsd],label='全部 Ag MSD',color='#4477aa',lw=1.8)
axm.plot(qt,[float(r['MSD_initial_center_r25_A2']) for r in qmsd],label='初始中心 r<25 Å Ag MSD',color='#d56b22',lw=1.8)
axm.set_title('Ag 位移（按逐帧周期边界解包）',fontproperties=bold,fontsize=12)
axm.set_xlabel('时间 (ps)',fontproperties=font); axm.set_ylabel('MSD (Å²)',fontproperties=font); axm.legend(prop=font,frameon=False)
axq.plot(qt,[float(r['Q6mean']) for r in qmsd],label='平均 Q6',color='#6655aa',lw=1.8)
axq.plot(qt,[float(r['Q6_gt_0p30_fraction']) for r in qmsd],label='Q6 > 0.30 比例',color='#39856b',lw=1.8)
axq.set_title('Ag 局部取向序参量',fontproperties=bold,fontsize=12)
axq.set_xlabel('时间 (ps)',fontproperties=font); axq.set_ylabel('Q₆ / 原子比例',fontproperties=font); axq.legend(prop=font,frameon=False)
fig.suptitle('seed41 第二次局域能量脉冲：凹坑、喷出候选与 Ag 无序度',fontproperties=bold,fontsize=17)
fig.text(.5,-.02,'输入为 1,479 eV 局部速度脉冲，之后 NVE；中心 (100,100,180.09) Å，σxy=12 Å，衰减深度 6 Å。骨架固定、无 TTM。橙色虚圈半径 25 Å；z>500 Å 的原子从有限计算盒移除。',ha='center',fontproperties=font,fontsize=9)
out=ROOT/'seed41_double_pulse_diagnostics_CN_300dpi.png'
fig.savefig(out,dpi=300,bbox_inches='tight',facecolor='white')
print(out,out.stat().st_size)
