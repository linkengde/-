#!/usr/bin/env python3
"""Plot crater, plume, mobility and Q6 for the three-pulse sequence."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from analyze_seed41_double_pulse import dump_frames

ROOT=Path(__file__).resolve().parent
FONT='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
BOLD='/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc'
font=FontProperties(fname=FONT); bold=FontProperties(fname=BOLD)

def readcsv(p):
    with open(p) as f: return list(csv.DictReader(f))
def surface(path):
    fr=dump_frames(path)[0]; idx=fr[1]; a=fr[2]
    x=a[:,idx['x']]-100.; y=a[:,idx['y']]-100.; z=a[:,idx['z']]
    keep=(np.abs(x)<50)&(np.abs(y)<50)&(z>=150)&(z<=230)
    m=np.full((25,25),np.nan); ix=np.floor((x[keep]+50)/4).astype(int); iy=np.floor((y[keep]+50)/4).astype(int)
    for bx,by,h in zip(ix,iy,z[keep]):
        if 0<=bx<25 and 0<=by<25 and (np.isnan(m[by,bx]) or h>m[by,bx]): m[by,bx]=h
    return m

paths=[ROOT.parent/'seed41_local_pulse_crater_20261126'/'seed41_local_pulse_final_5p8ps_all.dump',
       ROOT/'seed41_pulse2_final_7p8ps_all.dump',
       ROOT/'seed41_afterpulse2_final_10p8ps_nomsd_all.dump',
       ROOT/'seed41_pulse3_final_12p8ps_all.dump',
       ROOT/'seed41_afterpulse3_final_14p8ps_nomsd_all.dump']
maps=[surface(p) for p in paths]
labels=['5.8 ps｜脉冲2前','7.8 ps｜脉冲2后','10.8 ps｜脉冲2后 NVE','12.8 ps｜脉冲3后','14.8 ps｜脉冲3后 NVE']
tr2=readcsv(ROOT/'seed41_double_pulse_trajectory.csv'); tr3=readcsv(ROOT/'seed41_pulse3_trajectory.csv')
traj={float(r['time_ps']):r for r in tr2+tr3}; ts=np.array(sorted(traj)); rows=[traj[x] for x in ts]
q2=readcsv(ROOT/'seed41_double_pulse_q6_msd.csv'); q3=readcsv(ROOT/'seed41_pulse3_q6_msd.csv')
qmap={float(r['time_ps']):r for r in q2+q3}; qt=np.array(sorted(qmap)); qr=[qmap[x] for x in qt]
endpoints=[(5.8,-10.6645),(7.8,-12.766),(10.8,-8.2345),(12.8,-9.4525),(14.8,-0.72)]

fig=plt.figure(figsize=(22,10),constrained_layout=True)
gs=fig.add_gridspec(2,5,height_ratios=[1,0.84])
axs=[fig.add_subplot(gs[0,i]) for i in range(5)]
for ax,m,label in zip(axs,maps,labels):
    im=ax.imshow(np.ma.masked_invalid(m),origin='lower',extent=(-50,50,-50,50),vmin=170,vmax=210,cmap='Blues',interpolation='nearest')
    ax.add_patch(plt.Circle((0,0),25,fill=False,ls='--',lw=1.1,ec='#e78a28'))
    ax.set_title(label,fontproperties=bold,fontsize=11)
    ax.set_xlabel('相对 x (Å)',fontproperties=font,fontsize=9)
    if ax is axs[0]: ax.set_ylabel('相对 y (Å)',fontproperties=font,fontsize=9)
    else: ax.set_yticklabels([])
cb=fig.colorbar(im,ax=axs,shrink=.82,pad=.015); cb.set_label('Ag 局部顶面高度 z (Å)',fontproperties=font)

ax=fig.add_subplot(gs[1,0:2])
for key,label,col in [('n_z_gt_220','z > 220 Å','#5d83a5'),('n_z_gt_250','z > 250 Å','#71819b'),('n_z_gt_300','z > 300 Å','#e78a28'),('n_z_gt_400','z > 400 Å','#d84c34')]:
    ax.plot(ts,[int(r[key]) for r in rows],label=label,color=col,lw=1.7)
for t in (7.8,10.8,12.8): ax.axvline(t,color='#555',ls='--',lw=.8)
ax.set_title('计算盒内高空 Ag 数',fontproperties=bold); ax.set_xlabel('时间 (ps)',fontproperties=font); ax.set_ylabel('原子数',fontproperties=font); ax.legend(prop=font,ncol=2,frameon=False)

ax=fig.add_subplot(gs[1,2])
ax.plot([float(r['time_ps']) for r in q2],[float(r['MSD_Ag_A2']) for r in q2],label='自 5.8 ps',color='#4477aa',lw=1.6)
ax.plot([float(r['time_ps']) for r in q3],[float(r['MSD_Ag_A2']) for r in q3],label='自 10.8 ps',color='#d56b22',lw=1.6)
ax.plot([float(r['time_ps']) for r in q2],[float(r['MSD_initial_center_r25_A2']) for r in q2],ls=':',color='#4477aa',lw=1.4)
ax.plot([float(r['time_ps']) for r in q3],[float(r['MSD_initial_center_r25_A2']) for r in q3],ls=':',color='#d56b22',lw=1.4)
ax.set_title('Ag MSD（实线全体，点线初始中心）',fontproperties=bold); ax.set_xlabel('时间 (ps)',fontproperties=font); ax.set_ylabel('MSD (Å²)',fontproperties=font); ax.legend(prop=font,frameon=False,fontsize=8)

ax=fig.add_subplot(gs[1,3])
ax.plot(qt,[float(r['Q6mean']) for r in qr],label='平均 Q6',color='#6655aa',lw=1.7)
ax.plot(qt,[float(r['Q6_gt_0p30_fraction']) for r in qr],label='Q6 > 0.30 比例',color='#39856b',lw=1.7)
ax.set_title('Ag 局部取向序',fontproperties=bold); ax.set_xlabel('时间 (ps)',fontproperties=font); ax.set_ylabel('Q6 / 原子比例',fontproperties=font); ax.legend(prop=font,frameon=False,fontsize=8)

ax=fig.add_subplot(gs[1,4])
ax.plot([x[0] for x in endpoints],[x[1] for x in endpoints],'-o',color='#bb5931',lw=1.8)
ax.axhline(0,color='#555',lw=.8)
ax.set_title('中心相对外环高度差',fontproperties=bold); ax.set_xlabel('时间 (ps)',fontproperties=font); ax.set_ylabel('中心 − 外环中位 z (Å)',fontproperties=font)
ax.text(.03,.05,'负值代表中心较低\n14.8 ps 时约 −0.72 Å',transform=ax.transAxes,fontproperties=font,fontsize=8)

fig.suptitle('seed41 三次局域脉冲：熔坑、喷出候选与 Ag 结构变化',fontproperties=bold,fontsize=17)
fig.text(.5,-.025,'每次加入 1,479.25 eV 局域速度脉冲，之后 NVE；脉冲中心 (100,100,180.09) Å、σxy=12 Å、衰减深度 6 Å。骨架固定、无 TTM。顶面按 4 Å 网格中最高 Ag 估算；z>500 Å 原子越界后移除。高空 Ag 不能直接视作蒸发产率。',ha='center',fontproperties=font,fontsize=9)
out=ROOT/'seed41_three_pulse_diagnostics_CN_300dpi.png'
fig.savefig(out,dpi=300,bbox_inches='tight',facecolor='white')
print(out,out.stat().st_size)
