from pathlib import Path
import re
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter, gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path('.')
OUT = ROOT / 'audit'
CENTER = np.array([9.9573898076, 24.2578814220])


def read_data(path):
    text = Path(path).read_text(errors='replace').splitlines()
    bounds = {}
    tilt = (0.0, 0.0, 0.0)
    atoms_start = None
    for i, line in enumerate(text):
        s = line.split()
        if len(s) >= 4 and s[-2:] == ['xlo', 'xhi']:
            bounds['x'] = (float(s[0]), float(s[1]))
        elif len(s) >= 4 and s[-2:] == ['ylo', 'yhi']:
            bounds['y'] = (float(s[0]), float(s[1]))
        elif len(s) >= 4 and s[-2:] == ['zlo', 'zhi']:
            bounds['z'] = (float(s[0]), float(s[1]))
        elif len(s) >= 6 and s[-3:] == ['xy', 'xz', 'yz']:
            tilt = tuple(map(float, s[:3]))
        if line.strip().startswith('Atoms'):
            atoms_start = i + 1
            break
    if atoms_start is None:
        raise ValueError(f'Atoms section missing: {path}')
    rows = []
    for line in text[atoms_start:]:
        s = line.split('#', 1)[0].split()
        if not s:
            if rows:
                break
            continue
        if len(s) < 5:
            if rows:
                break
            continue
        try:
            rows.append((int(s[0]), int(s[1]), *map(float, s[2:5])))
        except ValueError:
            if rows:
                break
    a = np.asarray(rows, float)
    ids, typ = a[:, 0].astype(int), a[:, 1].astype(int)
    xyz = a[:, 2:5]
    lx = bounds['x'][1] - bounds['x'][0]
    ly = bounds['y'][1] - bounds['y'][0]
    xy = tilt[0]
    # LAMMPS restricted triclinic coordinates: r = f_a*a + f_b*b
    fb = (xyz[:, 1] - bounds['y'][0]) / ly
    fa = (xyz[:, 0] - bounds['x'][0] - xy * fb) / lx
    frac = np.column_stack((fa, fb))
    frac %= 1.0
    H = np.array([[lx, xy], [0.0, ly]])
    return {'path': str(path), 'ids': ids, 'typ': typ, 'xyz': xyz, 'frac': frac,
            'bounds': bounds, 'tilt': tilt, 'H': H, 'lx': lx, 'ly': ly}


def min_image_xy(atoms, center=CENTER):
    H = atoms['H']
    f0 = np.linalg.solve(H, center - np.array([atoms['bounds']['x'][0], atoms['bounds']['y'][0]]))
    df = atoms['frac'] - f0
    cand = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            q = (df + np.array([i, j])) @ H.T
            cand.append(np.sum(q*q, axis=1))
    return np.sqrt(np.min(cand, axis=0))


def summary(atoms, name):
    r = min_image_xy(atoms)
    z = atoms['xyz'][:, 2]
    print(f'[{name}] N={len(z)} type_counts=' + ','.join(f'{t}:{np.sum(atoms["typ"]==t)}' for t in sorted(set(atoms['typ']))))
    for lo, hi, label in [(0,2,'r<=2'),(0,3,'r<=3'),(2,4,'2-4'),(4,8,'4-8'),(8,12,'8-12'),(12,18,'12-18')]:
        m=(r>=lo)&(r<hi)
        if np.any(m):
            print(f'  {label}: n={m.sum()} zmax={np.max(z[m]):.3f} zmed={np.median(z[m]):.3f}')
    center = r<2
    print(f'  r<2 atoms above 35 A: {np.sum(center & (z>35))}; zmax={np.max(z[center]) if center.any() else float("nan"):.3f}')


def pair_audit(atoms, radius=2.5):
    ids, typ, xyz = atoms['ids'], atoms['typ'], atoms['xyz']
    H, lx, ly = atoms['H'], atoms['lx'], atoms['ly']
    xy = atoms['tilt'][0]
    a = np.array([lx, 0.0, 0.0]); b = np.array([xy, ly, 0.0])
    shifts = [(i*a+j*b) for i in (-1,0,1) for j in (-1,0,1)]
    points = np.concatenate([xyz+s for s in shifts], axis=0)
    labels = np.tile(ids, len(shifts))
    image_ix = np.repeat(np.arange(len(shifts)), len(ids))
    tree = cKDTree(points)
    best = {}
    for i, p in enumerate(xyz):
        for k in tree.query_ball_point(p, radius):
            j = int(labels[k])
            if j == ids[i]:
                continue
            lo, hi = sorted((int(ids[i]), j))
            if lo == hi:
                continue
            d = float(np.linalg.norm(p-points[k]))
            old = best.get((lo,hi))
            if old is None or d < old[0]:
                t1 = int(typ[np.where(ids==lo)[0][0]])
                t2 = int(typ[np.where(ids==hi)[0][0]])
                best[(lo,hi)] = (d,t1,t2)
    classes = {'Ag-Ag': lambda t1,t2: t1==4 and t2==4,
               'Ag-TSC': lambda t1,t2: (t1==4) != (t2==4),
               'TSC-TSC': lambda t1,t2: t1!=4 and t2!=4,
               'C-C': lambda t1,t2: t1==3 and t2==3}
    for label, pred in classes.items():
        arr=[(d,i,j) for (i,j),(d,t1,t2) in best.items() if pred(t1,t2)]
        if not arr:
            print(f'  {label}: no pairs <{radius} A')
            continue
        arr.sort()
        d,i,j=arr[0]
        under15=sum(x[0]<1.5 for x in arr)
        under20=sum(x[0]<2.0 for x in arr)
        print(f'  {label}: min={d:.4f} A IDs={i},{j}; pairs<1.5={under15}, <2.0={under20}')


def read_last_dump(path):
    lines=Path(path).read_text().splitlines()
    starts=[i for i,s in enumerate(lines) if s.startswith('ITEM: TIMESTEP')]
    if not starts: raise ValueError(f'no dump frames in {path}')
    i=starts[-1]
    timestep=int(lines[i+1])
    n=int(lines[i+3])
    j=i+4
    while not lines[j].startswith('ITEM: ATOMS'):
        j+=1
    cols=lines[j].split()[2:]
    data=np.array([[float(x) for x in row.split()] for row in lines[j+1:j+1+n]])
    return timestep, cols, data


def make_plot(base, final, output, label):
    def grid_surface(atoms, nbin=72):
        u=np.clip((atoms['frac'][:,0]*nbin).astype(int),0,nbin-1)
        v=np.clip((atoms['frac'][:,1]*nbin).astype(int),0,nbin-1)
        z=atoms['xyz'][:,2]
        grid=np.full((nbin,nbin),np.nan)
        for k,(i,j) in enumerate(zip(u,v)):
            if np.isnan(grid[j,i]) or z[k]>grid[j,i]: grid[j,i]=z[k]
        valid=np.isfinite(grid)&(grid>18.0)
        num=gaussian_filter(np.where(valid,grid,0.0),1.25,mode='wrap')
        den=gaussian_filter(valid.astype(float),1.25,mode='wrap')
        smooth=np.full_like(grid,np.nan)
        good=den>0.12
        smooth[good]=num[good]/den[good]
        return smooth

    nbin=72
    grids={id(base):grid_surface(base,nbin),id(final):grid_surface(final,nbin)}
    def radial_profile(atoms):
        grid=grids[id(atoms)]
        u=(np.arange(nbin)+0.5)/nbin
        vv,uu=np.meshgrid(u,u,indexing='ij')
        x=atoms['bounds']['x'][0]+uu*atoms['lx']+vv*atoms['tilt'][0]
        y=atoms['bounds']['y'][0]+vv*atoms['ly']
        pts=np.column_stack((x.ravel()-CENTER[0],y.ravel()-CENTER[1]))
        avec=np.array([atoms['lx'],0.0]); bvec=np.array([atoms['tilt'][0],atoms['ly']])
        images=[]
        for i in (-1,0,1):
            for j in (-1,0,1): images.append(pts+i*avec+j*bvec)
        images=np.stack(images,axis=0)
        im=np.argmin(np.sum(images*images,axis=2),axis=0)
        closest=images[im,np.arange(len(pts))]
        radius=np.sqrt(np.sum(closest*closest,axis=1)).reshape((nbin,nbin))
        edges=np.arange(0,16.5,.75); centers=(edges[:-1]+edges[1:])/2
        med=[]; q25=[]; q75=[]
        for lo,hi in zip(edges[:-1],edges[1:]):
            m=(radius>=lo)&(radius<hi)&np.isfinite(grid)
            vals=grid[m]
            if len(vals)<3: med.append(np.nan); q25.append(np.nan); q75.append(np.nan)
            else:
                med.append(np.median(vals)); q25.append(np.percentile(vals,25)); q75.append(np.percentile(vals,75))
        med=np.array(med); q25=np.array(q25); q75=np.array(q75)
        good=np.isfinite(med)
        if good.sum()>2:
            ix=np.arange(len(med)); med[~good]=np.interp(ix[~good],ix[good],med[good])
            med=gaussian_filter1d(med,1.0,mode='nearest')
        return centers,med,q25,q75
    rb,mb,_,_=radial_profile(base)
    rf,mf,q25,q75=radial_profile(final)

    fig=plt.figure(figsize=(14,6),dpi=170,constrained_layout=True)
    ax=fig.add_subplot(121,projection='3d')
    theta=np.linspace(0,2*np.pi,100)
    radius=np.linspace(0,15,90)
    RR,TT=np.meshgrid(radius,theta)
    ZZ=np.interp(RR,rf,mf,left=mf[0],right=mf[-1])
    XX=RR*np.cos(TT); YY=RR*np.sin(TT)
    surf=ax.plot_surface(XX,YY,ZZ,cmap='inferno',vmin=20,vmax=52,linewidth=0,antialiased=True,shade=False)
    ax.scatter([0],[0],[mf[0]],marker='x',s=65,c='cyan',label='impact center')
    ax.set_title('Radially averaged crater surface')
    ax.set_xlabel('x from impact (Å)'); ax.set_ylabel('y from impact (Å)'); ax.set_zlabel('median top z (Å)')
    ax.set_xlim(-15,15); ax.set_ylim(-15,15); ax.set_zlim(20,55)
    ax.view_init(elev=30,azim=-52); ax.set_box_aspect((30,30,30)); ax.legend(loc='upper left',fontsize=8)
    cb=fig.colorbar(surf,ax=ax,shrink=.68,pad=.02); cb.set_label('surface height z (Å)')

    ax2=fig.add_subplot(122)
    xyz=final['xyz']; d0=xyz[:,:2]-CENTER
    avec=np.array([final['lx'],0.0]); bvec=np.array([final['tilt'][0],final['ly']])
    images=[]
    for i in (-1,0,1):
        for j in (-1,0,1): images.append(d0+i*avec+j*bvec)
    images=np.stack(images,axis=0)
    score=np.abs(images[:,:,1])+1.0e-7*np.abs(images[:,:,0])
    im=np.argmin(score,axis=0)
    chosen=images[im,np.arange(len(xyz))]
    xs,ys=chosen[:,0],chosen[:,1]
    zmin=final['bounds']['z'][0]
    mask=(np.abs(ys)<1.2)&(np.abs(xs)<15.5)&(xyz[:,2]>max(zmin+1,15))&(xyz[:,2]<56)
    shown=set()
    for typ in sorted(set(final['typ'])):
        m=mask&(final['typ']==typ)
        if np.any(m):
            phase='Ag' if typ==4 else 'Ti₃SiC₂'
            if phase not in shown:
                ax2.scatter(xs[m],xyz[m,2],s=12,c=('#d3a43c' if typ==4 else '#7650a5'),label=phase,alpha=.82,linewidths=0)
                shown.add(phase)
            else:
                ax2.scatter(xs[m],xyz[m,2],s=12,c=('#d3a43c' if typ==4 else '#7650a5'),alpha=.82,linewidths=0)
    ax2.plot(rf,mf,color='black',lw=2.1,label='post-hit radial median')
    ax2.plot(-rf,mf,color='black',lw=2.1)
    ax2.plot(rb,mb,color='#777777',lw=1.6,ls='--',label='initial radial median')
    ax2.plot(-rb,mb,color='#777777',lw=1.6,ls='--')
    ax2.set_xlim(-15,15); ax2.set_ylim(15,55)
    ax2.set_xlabel('distance along section from impact (Å)'); ax2.set_ylabel('z (Å)')
    ax2.set_title('Atomic section: Ag and Ti₃SiC₂ skeleton')
    ax2.grid(alpha=.22); ax2.legend(fontsize=8,ncol=2)
    fig.suptitle(f'Ag–Ti₃SiC₂ pilot: {label}',fontsize=14)
    fig.savefig(output,bbox_inches='tight')
    plt.close(fig)

def max_image_r(atoms):
    return min_image_xy(atoms)

if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--final', default='structures/cool_02_0p5ps.data')
    parser.add_argument('--label', default='8 impacts + 0.5 ps cooling')
    parser.add_argument('--plot', default='audit/crater_after_8impacts_cooling.png')
    parser.add_argument('--dump', default='dumps/cool_02_0p5ps.dump')
    args=parser.parse_args()
    base=read_data('structures/base_active_local_relaxed.data')
    final=read_data(args.final)
    summary(base,'initial')
    print()
    summary(final,args.label)
    pair_audit(final)
    print()
    if args.dump and Path(args.dump).exists():
        step,cols,frame=read_last_dump(args.dump)
        force_cols=[cols.index(c) for c in ('fx','fy','fz') if c in cols]
        if len(force_cols)==3:
            force=np.linalg.norm(frame[:,force_cols],axis=1)
            i=int(np.argmax(force))
            xyz=frame[i,[cols.index('x'),cols.index('y'),cols.index('z')]]
            print(f'[final dump] step={step} N={len(frame)} max_force={force[i]:.6g} eV/A atom_id={int(frame[i,cols.index("id")])} type={int(frame[i,cols.index("type")])} xyz={xyz}')
    if args.plot:
        make_plot(base,final,args.plot,args.label)
        print('plot:',args.plot)
