from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

here=Path(__file__).resolve().parent
with np.load(here/'nano_geometry_ag607wt_ti3sic2_393wt.npz',allow_pickle=False) as data:
 skeleton=np.asarray(data['skeleton'],dtype=bool)
 size=np.asarray(data['size_A'],dtype=float)
a=~skeleton
nx,ny,nz=skeleton.shape
dx,dy,dz=size/np.array(skeleton.shape)
cmap=ListedColormap(['#b1843b','#547b9a']) # Ti3SiC2, Ag
fig=plt.figure(figsize=(15,12),layout='constrained')
gs=fig.add_gridspec(3,2,width_ratios=[1.1,1])
ax3=fig.add_subplot(gs[:,0],projection='3d')
step=2
ix,iy,iz=np.mgrid[0:nx:step,0:ny:step,0:nz:step]
subsk=skeleton[::step,::step,::step]
for phase,color,label in [(False,'#547b9a','Ag'),(True,'#b1843b','Ti₃SiC₂')]:
 m=subsk==phase
 ax3.scatter((ix[m]+.5*step)*dx,(iy[m]+.5*step)*dy,(iz[m]+.5*step)*dz,
             c=color,s=3.0,alpha=.62,linewidths=0,rasterized=True)
ax3.set(xlim=(0,size[0]),ylim=(0,size[1]),zlim=(0,size[2]),
        xlabel='x (Å)',ylabel='y (Å)',zlabel='z (Å)')
ax3.set_box_aspect(size.copy())
ax3.view_init(elev=23,azim=-58)
ax3.set_title('3D phase point cloud (every second voxel)')
ax3.legend(handles=[Patch(color='#547b9a',label='Ag'),Patch(color='#b1843b',label='Ti₃SiC₂')],loc='upper left')

axes=[fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,1]),fig.add_subplot(gs[2,1])]
sections=[(skeleton[:,:,nz//2].T, (0,size[0],0,size[1]),'XY midplane (z = {:.1f} Å)'.format(size[2]/2),'x (Å)','y (Å)'),
          (skeleton[:,ny//2,:].T, (0,size[0],0,size[2]),'XZ midplane (y = {:.1f} Å)'.format(size[1]/2),'x (Å)','z (Å)'),
          (skeleton[nx//2,:,:].T, (0,size[1],0,size[2]),'YZ midplane (x = {:.1f} Å)'.format(size[0]/2),'y (Å)','z (Å)')]
for ax,(sec,extent,title,xlab,ylab) in zip(axes,sections):
 ax.imshow(sec.astype(int),origin='lower',extent=extent,interpolation='nearest',aspect='equal',cmap=cmap,vmin=0,vmax=1)
 ax.set(title=title,xlabel=xlab,ylabel=ylab)
fig.suptitle('Ag 60.7 wt% / Ti₃SiC₂ 39.3 wt% | pore-free target Ag 39.9586 vol%',fontsize=16)
fig.savefig(here/'geometry_preview.png',dpi=180)
