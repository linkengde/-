from pathlib import Path
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'generation'))
import map_atomic_candidate as ma
LX,LY,XY=37.28597486102,32.290601434511,-18.64298743051
CX,CY=9.9573898076,24.2578814220
SH=ma.cell_translations(LX,LY,XY)
def read(p):
    l=Path(p).read_text().splitlines();j=next(i for i,s in enumerate(l) if s.strip().startswith('Atoms'))+1;r=[]
    for s in l[j:]:
        a=s.split('#',1)[0].split()
        if len(a)<5:
            if r:break
            continue
        try:r.append([int(a[0]),int(a[1]),*map(float,a[2:5])])
        except ValueError:
            if r:break
    a=np.array(r);a=a[np.argsort(a[:,0])];return a[:,0].astype(int),a[:,1].astype(int),a[:,2:5]
def radial(x):
    v=(x[:,1]/LY)%1;u=((x[:,0]-XY*(x[:,1]/LY))/LX)%1;cv=CY/LY;cu=((CX-XY*cv)/LX)%1;du=(u-cu+.5)%1-.5;dv=(v-cv+.5)%1-.5
    return np.sqrt((du*LX+dv*XY)**2+(dv*LY)**2)
source=read(ROOT/'structures/phasefield_clearance_2p25.data')
relaxed=read(ROOT/'structures/active_local_relaxed.data')
ids,t,x=relaxed
ids0,t0,x0=source
lookup={int(v):i for i,v in enumerate(ids0)}
indices=np.array([lookup[int(v)] for v in ids])
d=x-x0[indices]
dv=d[:,1]/LY;du=(d[:,0]-XY*dv)/LX;du-=np.rint(du);dv-=np.rint(dv);d[:,0]=du*LX+dv*XY;d[:,1]=dv*LY
dmag=np.linalg.norm(d,axis=1);r=radial(x)
active=radial(x0)<=8.0

def stats(v):return {'n':int(v.size),'rms_A':float(np.sqrt(np.mean(v*v))),'p50_A':float(np.percentile(v,50)),'p95_A':float(np.percentile(v,95)),'max_A':float(np.max(v))}
force_lines=Path(ROOT/'audit/active_local_relaxed_forces.dump').read_text().splitlines()
j=next(i for i,s in enumerate(force_lines) if s.startswith('ITEM: ATOMS'))+1
force=np.array([[float(q) for q in s.split()] for s in force_lines[j:]])
summary={
 'input_clearance_variant': 'Ag/TSC >= 2.25 A geometric mapping threshold before local FIRE',
 'impact_center_A':[CX,CY],
 'active_region': 'radius <= 8 A around actual impact center; outside constrained',
 'atoms':int(len(ids)),
 'type_counts':{str(int(k)):int(np.sum(t==k)) for k in np.unique(t)},
 'max_force_active_eV_per_A':float(np.max(force[radial(force[:,2:5])<=8,8])),
 'local_fire':{'iterations':3000,'stopping':'maximum iteration count; not converged to force tolerance'},
 'displacement_from_clearance_mapped_geometry':{'all':stats(dmag),'within_r_le8_A':stats(dmag[active]),'outside_r_gt8_max_A':float(np.max(dmag[~active]))},
 'connectivity':{
  'TSC_before':ma.component_summary(x0[t0!=4],3.35,SH),
  'TSC_after':ma.component_summary(x[t!=4],3.35,SH),
  'Ag_before':ma.component_summary(x0[t0==4],3.2,SH),
  'Ag_after':ma.component_summary(x[t==4],3.2,SH),
 },
 'cross_minima_after_local_fire':ma.all_pair_stats(x[t!=4],t[t!=4],x[t==4],t[t==4],SH,3.2),
}
(ROOT/'audit/centered_local_relax_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
