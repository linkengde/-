#!/usr/bin/env python3
"""Render a 3-D P1G composite/Ag pool frame from its compressed HDF5 handoff."""
import argparse
from pathlib import Path
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('h5',type=Path)
p.add_argument('output',type=Path)
p.add_argument('--time-ps',type=float,default=None)
p.add_argument('--largest-liquid',action='store_true',help='choose the frame with greatest total liquid_fraction*Ag volume')
a=p.parse_args()
with h5py.File(a.h5,'r') as f:
    times=f['time_ps'][:]
    if a.largest_liquid:
        ag=f['fields/ag_marker'][:]
        lf=f['fields/liquid_fraction'][:]
        vol=np.sum(ag*lf,axis=(1,2,3))
        idx=int(np.argmax(vol))
    elif a.time_ps is not None:
        idx=int(np.argmin(np.abs(times-a.time_ps)))
    else:
        idx=len(times)-1
    fields={name:f['fields/'+name][idx] for name in ['ag_marker','skeleton_marker','liquid_fraction','temperature_K']}
    x=f['mesh/x_center_m'][:]*1e9; y=f['mesh/y_center_m'][:]*1e9; z=f['mesh/z_center_m'][:]*1e9
    dx=float(f.attrs['cell_size_m'])*1e9
    t=float(times[idx]); mode=bool(f.attrs.get('marangoni_enabled',False))

ag,sc,lf,T=(fields[k] for k in ['ag_marker','skeleton_marker','liquid_fraction','temperature_K'])
# Background phases are sampled every second cell; liquid and mushy Ag remain full resolution.
ix,iy,iz=np.indices(ag.shape)
coarse=(ix%2==0)&(iy%2==0)&(iz%2==0)
skel=(sc>0.5)&coarse
solid=(ag>0.5)&(lf<=1e-4)&coarse
mushy=(ag>0.5)&(lf>1e-4)&(lf<0.5)
liquid=(ag>0.5)&(lf>=0.5)

def points(mask):
    ijk=np.argwhere(mask)
    return (x[ijk[:,0]],y[ijk[:,1]],z[ijk[:,2]]) if len(ijk) else ([],[],[])

fig=plt.figure(figsize=(10,8),dpi=180)
ax=fig.add_subplot(111,projection='3d')
for mask,color,label,size,alpha in [(skel,'#17232c','Ti3SiC2 scaffold',3.2,0.96),(solid,'#8bb9d0','solid Ag',2.8,0.38),(mushy,'#ffe21a','partly molten Ag',8,0.98),(liquid,'#f28b21','liquid Ag',8,0.48)]:
    xx,yy,zz=points(mask)
    if len(xx): ax.scatter(xx,yy,zz,c=color,s=size,alpha=alpha,depthshade=False,label=f'{label} ({len(xx):,})',edgecolors='none')
ax.set(xlabel='X (nm)',ylabel='Y (nm)',zlabel='Z (nm)')
ax.set_xlim(float(x.min()),float(x.max())); ax.set_ylim(float(y.min()),float(y.max())); ax.set_zlim(float(z.min()),float(z.max()))
ax.set_box_aspect((len(x),len(y),len(z)))
ax.view_init(elev=27,azim=-55)
ax.grid(False)
ax.set_title(f'Ag–Ti3SiC2 3-D phases | t={t:.6f} ps | Marangoni {"on" if mode else "off"}\nPeak T={float(np.max(T)):.1f} K | liquid Ag={float(np.sum(ag*lf)*dx**3):.2f} nm³')
ax.legend(loc='upper left',bbox_to_anchor=(0.0,1.0),frameon=False,fontsize=9)
fig.tight_layout()
a.output.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(a.output,bbox_inches='tight')
print(f'frame={idx} time_ps={t:.9g} output={a.output} bytes={a.output.stat().st_size} liquid_cells={int(liquid.sum())}')
