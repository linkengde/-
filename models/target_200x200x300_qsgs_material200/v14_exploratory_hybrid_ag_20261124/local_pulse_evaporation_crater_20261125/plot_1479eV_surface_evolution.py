#!/usr/bin/env python3
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

ROOT=Path(__file__).resolve().parent; OUT=ROOT/'results'
font=FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
fig,ax=plt.subplots(figsize=(12,7.5),dpi=200)
colors={20261133:'#c23b22',20261134:'#2166ac'}
for seed in (20261133,20261134):
    old=list(csv.DictReader((OUT/'trajectory_metrics_1479eV_two_seeds.csv').open()))
    ext=list(csv.DictReader((OUT/'trajectory_metrics_1479eV_extension_2ps.csv').open()))
    rows=[r for r in old if int(r['seed'])==seed]+[r for r in ext if int(r['seed'])==seed and float(r['elapsed_since_pulse_ps'])>2.00001]
    x=[float(r['elapsed_ps'] if 'elapsed_ps' in r else r['elapsed_since_pulse_ps']) for r in rows]
    y=[float(r['center_minus_annulus_q50_A']) for r in rows]
    ax.plot(x,y,color=colors[seed],lw=2.2,marker='o',ms=3.5,label=f'随机种子 {seed}')
ax.axhline(0,color='#555555',lw=1)
ax.axhline(-5,color='#777777',lw=1.4,ls='--',label='本轮探索凹陷判读线（−5 Å）')
ax.axvline(2,color='#999999',lw=1,ls=':')
ax.text(2.03,ax.get_ylim()[1]*.90,'开始延长续算',fontproperties=font,color='#555555')
ax.set_xlim(0,4); ax.set_xlabel('脉冲后时间（ps）',fontproperties=font,fontsize=14)
ax.set_ylabel('中心表面高度 − 20–40 Å 环带高度（Å）',fontproperties=font,fontsize=14)
ax.set_title('1.479 keV 脉冲后中心表面相对高度',fontproperties=font,fontsize=17)
ax.grid(True,alpha=.25)
ax.legend(prop=font,loc='best')
fig.tight_layout()
path=OUT/'surface_height_difference_1479eV_0to4ps.png'
fig.savefig(path,dpi=300); print(path)
