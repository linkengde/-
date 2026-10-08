from pathlib import Path
import importlib.util
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('audit_crater', ROOT/'audit/audit_crater.py')
ac = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ac)
CENTER = ac.CENTER
IMPACT = ROOT
CONTROL = ROOT.parent/'noimpact_control_1ps_after105ps_20261008'
START = ac.read_data(IMPACT/'checkpoints/after_hit_471.data')
HIT = ac.read_data(IMPACT/'checkpoints/after_hit_476.data')
CTRL = ac.read_data(CONTROL/'checkpoints/after_control_window_476.data')


def top_grid(atoms, n=96):
    u = np.mod(atoms['frac'][:,0],1.0)
    v = np.mod(atoms['frac'][:,1],1.0)
    z = atoms['xyz'][:,2]
    grid = np.full((n,n),np.nan)
    # Exclude sputtered/incoming atoms in vacuum; reconstruct solid top surface.
    valid = (z > 3.0) & (z < 55.0)
    for i,j,zz in zip((u[valid]*n).astype(int),(v[valid]*n).astype(int),z[valid]):
        if not np.isfinite(grid[j,i]) or zz > grid[j,i]: grid[j,i]=zz
    good = np.isfinite(grid)
    num = gaussian_filter(np.where(good,grid,0.0),3.0,mode='wrap')
    den = gaussian_filter(good.astype(float),3.0,mode='wrap')
    out = np.full_like(grid,np.nan)
    ok=den>0.08
    out[ok]=num[ok]/den[ok]
    return out


def local_mesh(atoms, grid):
    n=grid.shape[0]
    f=(np.arange(n)+0.5)/n
    vv,uu=np.meshgrid(f,f,indexing='ij')
    x=atoms['bounds']['x'][0]+uu*atoms['lx']+vv*atoms['tilt'][0]
    y=atoms['bounds']['y'][0]+vv*atoms['ly']
    d=np.stack([x-CENTER[0]+i*atoms['lx']+j*atoms['tilt'][0] for i in (-2,-1,0,1,2) for j in (-2,-1,0,1,2)],axis=0)
    dy=np.stack([y-CENTER[1]+j*atoms['ly'] for i in (-2,-1,0,1,2) for j in (-2,-1,0,1,2)],axis=0)
    d2=d*d+dy*dy
    idx=np.argmin(d2,axis=0)
    ii,jj=np.indices((n,n))
    X=d[idx,ii,jj]; Y=dy[idx,ii,jj]
    R=np.sqrt(X*X+Y*Y)
    return X,Y,R


def radial(atoms,grid):
    X,Y,R=local_mesh(atoms,grid)
    edges=np.arange(0,20.0001,.5)
    centers=(edges[:-1]+edges[1:])/2
    h=[]
    for a,b in zip(edges[:-1],edges[1:]):
        m=(R>=a)&(R<b)&np.isfinite(grid)
        h.append(float(np.median(grid[m])) if m.sum()>=5 else np.nan)
    return centers,np.array(h)


g0,g1,gc=top_grid(START),top_grid(HIT),top_grid(CTRL)
X,Y,R=local_mesh(HIT,g1)
r0,h0=radial(START,g0)
r1,h1=radial(HIT,g1)
rc,hc=radial(CTRL,gc)

def depth(atoms,grid):
    X,Y,R=local_mesh(atoms,grid)
    center=(R<2.0)&np.isfinite(grid)
    rim=(R>=8.0)&(R<12.0)&np.isfinite(grid)
    return float(np.median(grid[rim])-np.median(grid[center])),float(np.median(grid[center])),float(np.median(grid[rim]))

metrics={name:depth(at,g) for name,at,g in [('105ps impact start',START,g0),('106ps impact',HIT,g1),('106ps no-impact',CTRL,gc)]}

# Actual near-surface atoms, positioned in the nearest periodic image around the impact.
xyz=HIT['xyz']; frac=HIT['frac']
u=frac[:,0]; v=frac[:,1]
xx=HIT['bounds']['x'][0]+u*HIT['lx']+v*HIT['tilt'][0]
yy=HIT['bounds']['y'][0]+v*HIT['ly']
dx=xx-CENTER[0]; dy=yy-CENTER[1]
cands=[]
for i in range(-2,3):
  for j in range(-2,3): cands.append(np.column_stack((dx+i*HIT['lx']+j*HIT['tilt'][0],dy+j*HIT['ly'])))
cands=np.stack(cands,axis=0); k=np.argmin(np.sum(cands*cands,axis=2),axis=0)
xy=cands[k,np.arange(len(xyz))]; rho=np.linalg.norm(xy,axis=1)
# Atoms within 4 A below/above the local reconstructed surface; show both phases.
ii=np.minimum((u*g1.shape[1]).astype(int),g1.shape[1]-1)
jj=np.minimum((v*g1.shape[0]).astype(int),g1.shape[0]-1)
local_top=g1[jj,ii]
shell=(rho<19.0)&(xyz[:,2]>3.0)&(xyz[:,2]<55.0)&np.isfinite(local_top)&(xyz[:,2]>=local_top-4.5)&(xyz[:,2]<=local_top+1.5)
colors={'TSC':'#6840a0','Ag':'#d49b16'}
fig=plt.figure(figsize=(15,7.5),dpi=190,constrained_layout=True)
ax=fig.add_subplot(121,projection='3d')
surf=ax.plot_surface(X,Y,np.ma.masked_invalid(g1),cmap='inferno',vmin=10,vmax=40,linewidth=0,antialiased=True,alpha=.46,rcount=100,ccount=100,shade=True)
for phase,mask,color in [('Ti₃SiC₂',HIT['typ']!=4,colors['TSC']),('Ag',HIT['typ']==4,colors['Ag'])]:
    m=shell&mask
    if m.any(): ax.scatter(xy[m,0],xy[m,1],xyz[m,2],s=16,c=color,alpha=.98,edgecolors='none',depthshade=False,label=phase)
ax.scatter([0],[0],[np.nanmedian(g1[(R<2)&np.isfinite(g1)])],marker='x',c='cyan',s=45,linewidths=2,label='Impact centre')
ax.set_xlim(-19,19);ax.set_ylim(-19,19);ax.set_zlim(5,45)
ax.set_xlabel('x from impact (Å)');ax.set_ylabel('y from impact (Å)');ax.set_zlabel('z (Å)')
ax.set_title('106 ps: 3D crater; purple = Ti₃SiC₂, gold = Ag')
ax.view_init(elev=30,azim=-48);ax.set_box_aspect((38,38,31));ax.legend(loc='upper left',fontsize=8,ncol=2)
cb=fig.colorbar(surf,ax=ax,shrink=.68,pad=.025);cb.set_label('smoothed top-surface height z (Å)')

ax2=fig.add_subplot(122)
# Vertical atomic section through impact center, with phase colors.
section=(np.abs(xy[:,1])<2.0)&(np.abs(xy[:,0])<20)&(xyz[:,2]>5)&(xyz[:,2]<40)
for phase,mask,color in [('Ti₃SiC₂',HIT['typ']!=4,colors['TSC']),('Ag',HIT['typ']==4,colors['Ag'])]:
    m=section&mask
    ax2.scatter(xy[m,0],xyz[m,2],s=16,c=color,alpha=.82,edgecolors='none',label=phase)
ax2.plot(r0,h0,color='#555555',lw=1.5,ls='--',label='105 ps start profile')
ax2.plot(-r0,h0,color='#555555',lw=1.5,ls='--')
ax2.plot(r1,h1,color='#bb352d',lw=2,label='106 ps impact profile')
ax2.plot(-r1,h1,color='#bb352d',lw=2)
ax2.plot(rc,hc,color='#177e89',lw=1.7,ls=':',label='106 ps no-impact profile')
ax2.plot(-rc,hc,color='#177e89',lw=1.7,ls=':')
ax2.set_xlim(-20,20);ax2.set_ylim(5,40);ax2.set_xlabel('distance along section from impact centre (Å)');ax2.set_ylabel('z (Å)')
ax2.set_title('Atomic section and median top-surface profiles')
ax2.grid(alpha=.25);ax2.legend(fontsize=8,ncol=2)
fig.suptitle('Ag–Ti₃SiC₂ impact-pit geometry at 106 ps: two phase colors',fontsize=15,fontweight='bold')
out=ROOT/'images/crater_3d_two_phase_colors_106ps.png';out.parent.mkdir(exist_ok=True);fig.savefig(out,bbox_inches='tight',facecolor='white');plt.close(fig)
print('IMAGE',out,out.stat().st_size)
print('METRICS depth/center-med/rim-med Å')
for name,val in metrics.items():print(name,*(round(x,3) for x in val))
