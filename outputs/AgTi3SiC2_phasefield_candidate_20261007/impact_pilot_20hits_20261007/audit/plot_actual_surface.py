from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_crater import read_data

a=read_data('structures/cool_05_quench_0p5ps.data')
nb=96
u=np.clip((a['frac'][:,0]*nb).astype(int),0,nb-1)
v=np.clip((a['frac'][:,1]*nb).astype(int),0,nb-1)
z=a['xyz'][:,2]
top=np.full((nb,nb),-np.inf)
for i,(j,k) in enumerate(zip(v,u)):
    top[j,k]=max(top[j,k],z[i])
# Retain atoms in the outermost 2 A of populated surface columns, excluding deep internal columns.
keep=(top[v,u]>=20.0)&(z>=top[v,u]-2.0)
xyz=a['xyz'][keep]
typ=a['typ'][keep]
center=np.array([9.9573898076,24.2578814220])
d0=xyz[:,:2]-center
av=np.array([a['lx'],0.0]); bv=np.array([a['tilt'][0],a['ly']])
imgs=np.stack([d0+i*av+j*bv for i in (-1,0,1) for j in (-1,0,1)],axis=0)
idx=np.argmin(np.sum(imgs*imgs,axis=2),axis=0)
xy=imgs[idx,np.arange(len(xyz))]
fig=plt.figure(figsize=(9,7),dpi=180)
ax=fig.add_subplot(111,projection='3d')
colors=np.where(typ==4,'#d3a43c','#7650a5')
ax.scatter(xy[:,0],xy[:,1],xyz[:,2],c=colors,s=18,alpha=.95,depthshade=False,linewidths=0)
ax.set_xlim(-15,15); ax.set_ylim(-15,15); ax.set_zlim(20,55)
ax.set_xlabel('x from impact center (Å)'); ax.set_ylabel('y from impact center (Å)'); ax.set_zlabel('z (Å)')
ax.set_title('Actual top-layer atoms after 20 Ag impacts and 1 ps quench\noutermost 2 Å in occupied surface columns')
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],marker='o',color='w',label='Ti₃SiC₂',markerfacecolor='#7650a5',markersize=8),Line2D([0],[0],marker='o',color='w',label='Ag',markerfacecolor='#d3a43c',markersize=8)],loc='upper left')
ax.view_init(elev=32,azim=-52)
ax.set_box_aspect((30,30,30))
fig.tight_layout()
out=Path('audit/model_20hits_actual_surface_atoms.png')
fig.savefig(out,bbox_inches='tight')
print('atoms shown',len(xyz),'saved',out)
