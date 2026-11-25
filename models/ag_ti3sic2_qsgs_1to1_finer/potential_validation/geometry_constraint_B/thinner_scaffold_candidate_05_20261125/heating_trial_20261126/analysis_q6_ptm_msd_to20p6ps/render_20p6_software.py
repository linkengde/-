#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""使用OVITO导出的20.6 ps原子数据生成可视化图片（Matplotlib软件渲染）。"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = np.load(HERE / "frame_20p6ps_ovito_analysis_arrays.npz")
pos = DATA["positions"].astype(np.float64)
types = DATA["types"]
qbar = DATA["qbar6"]
ag = types == 4
skeleton = ~ag
norm = Normalize(vmin=0.17, vmax=0.62, clip=True)
cmap = plt.get_cmap("coolwarm_r")

def draw(ax, crop=False):
    """绘制整体模型或热源区，Ag颜色表示局部qbar6。"""
    if crop:
        r2 = (pos[:, 0]-100.0)**2 + (pos[:, 1]-100.0)**2
        mask = (r2 <= 42.0**2) & (pos[:, 2] >= 100.0) & (pos[:, 2] <= 240.0)
        ax.set_xlim(58, 142); ax.set_ylim(58, 142); ax.set_zlim(100, 240)
        ax.set_title("Heated core and ejected Ag", fontsize=15, pad=8)
    else:
        mask = np.ones(len(pos), dtype=bool)
        ax.set_xlim(0, 200); ax.set_ylim(0, 200); ax.set_zlim(0, 320)
        ax.set_title("Full model with free-surface plume", fontsize=15, pad=8)
    p, t, q = pos[mask], types[mask], qbar[mask]
    is_ag = t == 4
    ax.scatter(p[~is_ag,0],p[~is_ag,1],p[~is_ag,2],c="#173f78",alpha=0.34,
               s=0.16 if not crop else 0.48,depthshade=False,linewidths=0,rasterized=True)
    ax.scatter(p[is_ag,0],p[is_ag,1],p[is_ag,2],c=q[is_ag],cmap=cmap,norm=norm,
               s=0.28 if not crop else 1.3,depthshade=False,linewidths=0,rasterized=True)
    ax.view_init(elev=22 if not crop else 20, azim=-54)
    ax.set_box_aspect((200,200,320) if not crop else (84,84,140))
    ax.set_proj_type("persp",focal_length=1.4)
    ax.set_axis_off()

fig=plt.figure(figsize=(18,8),facecolor="white")
fig.suptitle("Ag/Ti3SiC2 candidate 05 | 20 ps heating + 0.6 ps ideal-sink cooling",fontsize=18,y=0.96)
ax1=fig.add_subplot(1,2,1,projection="3d"); ax2=fig.add_subplot(1,2,2,projection="3d")
draw(ax1,False); draw(ax2,True)
sm=matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap); sm.set_array([])
cb=fig.colorbar(sm,ax=[ax1,ax2],fraction=0.025,pad=0.02,shrink=0.74)
cb.set_label("Ag local qbar6 (relative order indicator)")
fig.legend(handles=[
 Line2D([0],[0],marker="o",color="w",markerfacecolor="#173f78",markersize=9,label="Rigid Ti/Si/C support"),
 Line2D([0],[0],marker="o",color="w",markerfacecolor="#df4b3f",markersize=9,label="Ag, lower qbar6"),
 Line2D([0],[0],marker="o",color="w",markerfacecolor="#b6d7e8",markersize=9,label="Ag, higher qbar6")],
 loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(0.48,0.02),fontsize=10)
fig.text(0.5,0.002,"OVITO data; software-rendered. qbar6 is a relative disorder indicator, not a universal phase boundary.",ha="center",fontsize=9)
fig.subplots_adjust(left=0.01,right=0.92,top=0.88,bottom=0.08,wspace=0.01)
out=HERE/"candidate05_20p6ps_OVITO_data_view.png"
fig.savefig(out,dpi=180,facecolor="white")
print(f"WROTE {out}",flush=True)
