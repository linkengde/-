#!/usr/bin/env python3
"""Compare atomistic deformation against the unbombarded reference structure."""
from pathlib import Path
import importlib.util, json
import numpy as np
from scipy.ndimage import gaussian_filter, gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
PILOT=ROOT.parents[1]
TO100=ROOT.parent/'to100ps_200fs_20261008'
spec=importlib.util.spec_from_file_location('audit_crater',TO100/'audit/audit_crater.py')
ac=importlib.util.module_from_spec(spec);spec.loader.exec_module(ac)
CENTER=np.array([9.9573898076,24.2578814220])
BASE=ac.read_data(PILOT/'structures/base_active_local_relaxed.data')
STATES=[('80 Ag impacts',TO100/'checkpoints/after_hit_080.data'),
        ('230 Ag impacts',TO100/'checkpoints/after_hit_230.data'),
        ('476 Ag impacts',ROOT/'checkpoints/after_hit_476.data')]
X0,Y0=-20.,-20.; Z0=-5.; DX=DY=DZ=2.0
NX=NY=20; NZ=25


def local_xy(atoms):
    H=atoms['H']; origin=np.array([atoms['bounds']['x'][0],atoms['bounds']['y'][0]])
    center_frac=np.linalg.solve(H,CENTER-origin)
    df=atoms['frac']-center_frac
    cand=np.stack([(df+np.array([i,j]))@H.T for i in (-2,-1,0,1,2) for j in (-2,-1,0,1,2)])
    ix=np.argmin(np.sum(cand*cand,axis=2),axis=0)
    return cand[ix,np.arange(len(df))]


def displacement_field(state):
    common,ib,is_=np.intersect1d(BASE['ids'],state['ids'],assume_unique=True,return_indices=True)
    if len(common)<0.99*len(BASE['ids']):raise RuntimeError('many original atoms are missing')
    H=BASE['H'];df=state['frac'][is_]-BASE['frac'][ib]
    candidates=np.stack([(df+np.array([i,j]))@H.T for i in (-2,-1,0,1,2) for j in (-2,-1,0,1,2)])
    im=np.argmin(np.sum(candidates*candidates,axis=2),axis=0)
    dxy=candidates[im,np.arange(len(common))]
    dz=state['xyz'][is_,2]-BASE['xyz'][ib,2]
    u=np.column_stack((dxy,dz))
    # Bottom atoms act as a fixed reference; remove any residual common translation.
    deep=BASE['xyz'][ib,2]<5.0
    drift=np.median(u[deep],axis=0) if np.any(deep) else np.zeros(3)
    u-=drift
    magnitude=np.linalg.norm(u,axis=1)
    # Spread displacement magnitude into a local 3D field for legibility.
    pos=local_xy(state)[is_]
    xyz=np.column_stack((pos,state['xyz'][is_,2]))
    ix=np.floor((xyz[:,0]-X0)/DX).astype(int)
    iy=np.floor((xyz[:,1]-Y0)/DY).astype(int)
    iz=np.floor((xyz[:,2]-Z0)/DZ).astype(int)
    keep=(ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)&(iz>=0)&(iz<NZ)
    linear=(iz[keep]*NY+iy[keep])*NX+ix[keep]
    vals=magnitude[keep]
    smag=np.bincount(linear,weights=vals,minlength=NX*NY*NZ).reshape(NZ,NY,NX)
    count=np.bincount(linear,minlength=NX*NY*NZ).reshape(NZ,NY,NX).astype(float)
    smag=gaussian_filter(smag,1.0,mode='constant',cval=0.0)
    count=gaussian_filter(count,1.0,mode='constant',cval=0.0)
    field=np.full_like(smag,np.nan)
    valid=count>=0.03
    field[valid]=smag[valid]/count[valid]
    # Map the local field back to atom positions for atom-cloud rendering.
    fmag=np.full(len(state['ids']),np.nan)
    jx=np.floor((local_xy(state)[:,0]-X0)/DX).astype(int)
    jy=np.floor((local_xy(state)[:,1]-Y0)/DY).astype(int)
    jz=np.floor((state['xyz'][:,2]-Z0)/DZ).astype(int)
    ok=(jx>=0)&(jx<NX)&(jy>=0)&(jy<NY)&(jz>=0)&(jz<NZ)
    fmag[ok]=field[jz[ok],jy[ok],jx[ok]]
    return common,ib,is_,u,magnitude,pos,field,fmag


def axes_box(ax,xlim=(-20,20),ylim=(-20,20),zlim=(-5,45)):
    xs=[xlim[0],xlim[1]];ys=[ylim[0],ylim[1]];zs=[zlim[0],zlim[1]]
    for x in xs:
      for y in ys:
       ax.plot([x,x],[y,y],zs,color='#555555',lw=.65,alpha=.45)
    for x in xs:
      for z in zs:
       ax.plot([x,x],ys,[z,z],color='#555555',lw=.65,alpha=.45)
    for y in ys:
      for z in zs:
       ax.plot(xs,[y,y],[z,z],color='#555555',lw=.65,alpha=.45)


def radial_surface_profile(atoms,nbin=96):
    u=np.mod(atoms['frac'][:,0],1.0);v=np.mod(atoms['frac'][:,1],1.0);z=atoms['xyz'][:,2]
    grid=np.full((nbin,nbin),np.nan)
    keep=(z>-5.0)&(z<55.0)
    for i,j,zz in zip((u[keep]*nbin).astype(int),(v[keep]*nbin).astype(int),z[keep]):
        if not np.isfinite(grid[j,i]) or zz>grid[j,i]:grid[j,i]=zz
    valid=np.isfinite(grid)
    num=gaussian_filter(np.where(valid,grid,0.0),3.0,mode='wrap')
    den=gaussian_filter(valid.astype(float),3.0,mode='wrap')
    surface=np.full_like(grid,np.nan);good=den>.08;surface[good]=num[good]/den[good]
    f=(np.arange(nbin)+.5)/nbin;vv,uu=np.meshgrid(f,f,indexing='ij')
    x=atoms['bounds']['x'][0]+uu*atoms['lx']+vv*atoms['tilt'][0]
    y=atoms['bounds']['y'][0]+vv*atoms['ly']
    H=atoms['H'];origin=np.array([atoms['bounds']['x'][0],atoms['bounds']['y'][0]])
    center_frac=np.linalg.solve(H,CENTER-origin)
    frac=np.stack((uu.ravel(),vv.ravel()),axis=1)
    df=frac-center_frac
    cand=np.stack([(df+np.array([i,j]))@H.T for i in (-1,0,1) for j in (-1,0,1)])
    ix=np.argmin(np.sum(cand*cand,axis=2),axis=0)
    closest=cand[ix,np.arange(len(df))]
    radius=np.linalg.norm(closest,axis=1).reshape(nbin,nbin)
    edges=np.arange(0,20.5001,.5);centers=(edges[:-1]+edges[1:])/2;med=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        m=(radius>=lo)&(radius<hi)&np.isfinite(surface)
        med.append(float(np.median(surface[m])) if m.sum()>=3 else np.nan)
    med=np.asarray(med);g=np.isfinite(med)
    if g.sum()>2:
        idx=np.arange(len(med));med[~g]=np.interp(idx[~g],idx[g],med[g]);med=gaussian_filter1d(med,1.0,mode='nearest')
    return centers,med


plt.rcParams.update({'font.size':9,'axes.titleweight':'bold'})
fig=plt.figure(figsize=(15.3,9.2),dpi=190)
gs=fig.add_gridspec(2,3,height_ratios=[1.15,1.0],left=.045,right=.985,top=.91,bottom=.15,hspace=.34,wspace=.22)
base_r,base_surface=radial_surface_profile(BASE)
reports=[]
for col,(label,path) in enumerate(STATES):
    state=ac.read_data(path)
    common,ib,is_,u,mag,pos,field,fmag=displacement_field(state)
    typ=state['typ'];xyz=state['xyz'];new=state['ids']>len(BASE['ids'])
    loc=local_xy(state)
    view=(np.abs(loc[:,0])<20)&(np.abs(loc[:,1])<20)&(xyz[:,2]>=-5)&(xyz[:,2]<45)
    old=view&~new&np.isfinite(fmag)
    newag=view&new&(typ==4)
    ax=fig.add_subplot(gs[0,col],projection='3d')
    sc=ax.scatter(loc[old,0],loc[old,1],xyz[old,2],c=np.clip(fmag[old],0,15),cmap='turbo',vmin=0,vmax=15,s=5,alpha=.87,linewidths=0,depthshade=False,label='original atoms: displacement')
    ax.scatter(loc[newag,0],loc[newag,1],xyz[newag,2],c='#e1ad27',s=8,alpha=.85,linewidths=0,depthshade=False,label='Ag added by impacts')
    axes_box(ax)
    ax.set(xlim=(-20,20),ylim=(-20,20),zlim=(-5,45),xlabel='x (Å)',ylabel='y (Å)',zlabel='z (Å)')
    ax.set_title(label,fontsize=11,pad=2)
    ax.tick_params(labelsize=7,pad=0)
    ax.xaxis.labelpad=-1;ax.yaxis.labelpad=-1;ax.zaxis.labelpad=-1
    ax.view_init(elev=22,azim=-50);ax.set_box_aspect((40,40,50))
    if col==0:ax.legend(loc='upper left',fontsize=7)

    ax2=fig.add_subplot(gs[1,col])
    section=view&(np.abs(loc[:,1])<2.0)
    so=section&~new&np.isfinite(fmag)
    sn=section&newag
    ax2.scatter(loc[so,0],xyz[so,2],c=np.clip(fmag[so],0,15),cmap='turbo',vmin=0,vmax=15,s=8,alpha=.88,linewidths=0)
    ax2.scatter(loc[sn,0],xyz[sn,2],c='#e1ad27',s=10,alpha=.88,linewidths=0)
    ax2.set(xlim=(-20,20),ylim=(-5,45),xlabel='x at y≈0 (Å)',ylabel='z (Å)')
    sr,sh=radial_surface_profile(state)
    ax2.plot(base_r,base_surface,color='#666666',lw=1.1,ls='--',label='unbombarded radial median' if col==0 else None)
    ax2.plot(-base_r,base_surface,color='#666666',lw=1.1,ls='--')
    ax2.plot(sr,sh,color='#111111',lw=1.5,label='current radial median' if col==0 else None)
    ax2.plot(-sr,sh,color='#111111',lw=1.5)
    ax2.set_title('Atom slice + radial median surface',fontsize=10,pad=3)
    if col==0:ax2.legend(fontsize=7,loc='upper right',framealpha=.85)
    ax2.grid(alpha=.18)
    local=mag[np.linalg.norm(pos,axis=1)<8]
    far=mag[np.linalg.norm(pos,axis=1)>15]
    center_h=float(np.nanmedian(sh[sr<2.0]));rim_h=float(np.nanmedian(sh[(sr>=8.0)&(sr<12.0)]))
    reports.append({'state':label,'path':str(path),'total_atoms':len(state['ids']),'original_atoms_matched':len(common),'new_Ag_shown':int(newag.sum()),'local_r_lt_8_displacement_median_p90_A':[float(x) for x in np.percentile(local,[50,90])],'far_r_gt_15_displacement_median_p90_A':[float(x) for x in np.percentile(far,[50,90])],'radial_surface_center_median_A':center_h,'radial_surface_outer_ring_median_A':rim_h,'radial_surface_depression_A':rim_h-center_h,'bottom_drift_removed_A':[float(x) for x in np.median(u[BASE['xyz'][ib,2]<5.0],axis=0)]})

cax=fig.add_axes([.22,.065,.56,.025])
cb=fig.colorbar(sc,cax=cax,orientation='horizontal',aspect=45)
cb.set_label('Local-smoothed displacement magnitude |Δr| from unbombarded structure (Å; capped at 15 Å)')
fig.suptitle('Ag–Ti₃SiC₂: progressive deformation from actual atom coordinates',fontsize=15,fontweight='bold',y=.975)
out=ROOT/'images/deformation_history_80_230_476impacts.png'
out.parent.mkdir(exist_ok=True);fig.savefig(out,bbox_inches='tight',facecolor='white');plt.close(fig)
summary={'reference':str(PILOT/'structures/base_active_local_relaxed.data'),'reference_atoms':len(BASE['ids']),'method':'Match baseline atoms by ID; unwrap x/y periodic displacement; subtract median bottom-layer drift; show actual snapshot coordinates. Color for original atoms is the 3D Gaussian-smoothed local mean displacement magnitude. Newly added Ag is gold. This is deformation relative to the unbombarded reference, not temperature.','states':reports}
(ROOT/'analysis_deformation_history.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
print('IMAGE',out,out.stat().st_size)
print(json.dumps(summary,ensure_ascii=False,indent=2))
