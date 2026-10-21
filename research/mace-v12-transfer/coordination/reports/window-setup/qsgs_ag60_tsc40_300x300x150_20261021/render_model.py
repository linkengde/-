from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes
from ase.io import read
ROOT=Path(__file__).resolve().parent;M=ROOT/'model'
audit=json.loads((M/'GEOMETRY_AUDIT.json').read_text())
arrival=np.load(M/'qsgs_growth.npz')['arrival'];threshold=audit['composition_trials'][-1]['threshold'];mask=arrival<=threshold
blue='#2158cf';gray='#bfc1c3';cmap=ListedColormap([gray,blue])
fig=plt.figure(figsize=(13,6),dpi=180)
ax=fig.add_subplot(121,projection='3d')
# Three exterior faces of the initial QSGS morphology, no extrusion.
for axis,idx in [(0,149),(1,0),(2,74)]:
 if axis==0:
  yy,zz=np.meshgrid(np.arange(150)*2+1,np.arange(75)*2+1,indexing='ij');xx=np.full_like(yy,300);colors=cmap(mask[idx,:,:].astype(int))
 elif axis==1:
  xx,zz=np.meshgrid(np.arange(150)*2+1,np.arange(75)*2+1,indexing='ij');yy=np.zeros_like(xx);colors=cmap(mask[:,idx,:].astype(int))
 else:
  xx,yy=np.meshgrid(np.arange(150)*2+1,np.arange(150)*2+1,indexing='ij');zz=np.full_like(xx,150);colors=cmap(mask[:,:,idx].astype(int))
 ax.plot_surface(xx,yy,zz,facecolors=colors,rstride=1,cstride=1,shade=False,linewidth=0,antialiased=False)
ax.set_title('Two-phase morphology: 300 x 300 x 150 A')
ax2=fig.add_subplot(122,projection='3d')
v,f,_,_=marching_cubes(np.pad(mask[::2,::2,::2].astype(np.float32),1),level=.5,spacing=(4,4,4))
v-=4
poly=Poly3DCollection(v[f],facecolor=blue,edgecolor='none',alpha=1)
ax2.add_collection3d(poly);ax2.set_title('Ti3SiC2 skeleton (Ag hidden)')
for aa in [ax,ax2]:
 aa.set(xlim=(0,300),ylim=(0,300),zlim=(0,250),xlabel='x / A',ylabel='y / A',zlabel='z / A');aa.set_box_aspect((300,300,250));aa.view_init(elev=24,azim=-55)
 # box top including 100 A vacuum
 for x in (0,300):
  for y in (0,300):aa.plot([x,x],[y,y],[0,250],color='#666666',lw=.6)
 aa.text(0,300,215,'100 A vacuum',fontsize=8)
fig.suptitle('QSGS geometry preview | Blue: Ti3SiC2; Gray: Ag | Not relaxed')
fig.tight_layout();fig.savefig(M/'MODEL_3D.png');plt.close(fig)
# Atom-based thin section verifies the retained structure rather than just mask.
a=read(M/'AG_TSC_300x300x150.extxyz');p=a.positions;ph=a.arrays['phase_id'];sel=(p[:,2]>=74)&(p[:,2]<76)
fig,ax=plt.subplots(figsize=(7,7),dpi=180)
ax.set_facecolor('white');order=np.flatnonzero(sel)
ax.scatter(p[order,0],p[order,1],c=np.where(ph[order]==1,blue,gray),s=2,linewidths=0)
ax.set(xlim=(0,300),ylim=(0,300),xlabel='x / A',ylabel='y / A',title='Retained atomic slice: z=74-76 A\nBlue: Ti3SiC2; Gray: Ag; White: construction gaps')
ax.set_aspect('equal');fig.tight_layout();fig.savefig(M/'ATOMIC_XY_SLICE.png');plt.close(fig)
print('Rendered MODEL_3D.png and ATOMIC_XY_SLICE.png')
