#!/usr/bin/env python3
"""Analyze the 100-step Ag-only relaxation with a rigid Ti3SiC2 scaffold."""
import gzip, json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
HERE=Path(__file__).resolve().parent
BOX=np.asarray(json.loads((HERE/'geometry_report.json').read_text())['simulation_box_A'],float)
ELEMENT={1:'Ti',2:'Si',3:'C',4:'Ag'}
def load(p):
 with gzip.open(p,'rt') as f:return np.loadtxt(f,skiprows=9)
def stats(x):
 return {'n':int(len(x)),'p50_eV_A':float(np.percentile(x,50)),'p95_eV_A':float(np.percentile(x,95)),'p99_eV_A':float(np.percentile(x,99)),'max_eV_A':float(np.max(x)),'n_gt_2':int(np.sum(x>2))}
def topology(xyz,ids,types,cutoff,box):
 pairs=cKDTree(xyz,boxsize=box).query_pairs(cutoff,output_type='ndarray')
 g=coo_matrix((np.ones(len(pairs)*2,dtype=np.uint8),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(xyz),len(xyz))).tocsr()
 n,lab=connected_components(g,directed=False);sz=np.bincount(lab);large=int(sz.max());small=np.flatnonzero(sz<large);out=np.flatnonzero(np.isin(lab,small))
 near,ni=cKDTree(xyz,boxsize=box).query(xyz[out],k=2) if len(out) else (np.empty((0,2)),np.empty((0,2),int))
 return {'components':int(n),'largest_component_atoms':large,'outside_largest':int(len(xyz)-large),'cutoff_A':cutoff,'detached_ids':[int(ids[i]) for i in out]}
ini=load(HERE/'initial_forces.dump.gz'); fin=load(HERE/'ag_only_rigid_skeleton_100.dump.gz')
assert np.array_equal(ini[:,0],fin[:,0]) and np.array_equal(ini[:,1],fin[:,1])
ids=ini[:,0].astype(int); t=ini[:,1].astype(int); ip=ini[:,2:5]; fp=fin[:,2:5]
d=fp-ip;d[:,0]-=BOX[0]*np.rint(d[:,0]/BOX[0]);d[:,1]-=BOX[1]*np.rint(d[:,1]/BOX[1]);disp=np.linalg.norm(d,axis=1);f=np.linalg.norm(fin[:,5:8],axis=1)
mobile=ip[:,2]>=20
xyz=fp.copy();xyz[:,0]%=BOX[0];xyz[:,1]%=BOX[1];pbox=[BOX[0],BOX[1],1000.]
mobile_stats={}
for ty,name in ELEMENT.items():
 sel=mobile&(t==ty);v=f[sel];q=disp[sel]
 mobile_stats[name]={'displacement_A':{'p50':float(np.percentile(q,50)),'p95':float(np.percentile(q,95)),'max':float(np.max(q))},'force':stats(v)}
minpair={}
for a,b in [('Ag','Ti'),('Ag','Si'),('Ag','C'),('Ti','C')]:
 ta=next(k for k,v in ELEMENT.items() if v==a);tb=next(k for k,v in ELEMENT.items() if v==b)
 tree=cKDTree(xyz[t==tb],boxsize=pbox); minpair[a+'-'+b]=float(np.min(tree.query(xyz[t==ta],k=1)[0]))
sel_sk=t!=4;sel_ag=t==4
res={'status':'STATIC_100_STEP_AG_ONLY_RIGID_SCAFFOLD_SCREEN_NOT_A_COLD_BASELINE','atoms':int(len(ids)),'max_mobile_Ag_force_eV_A':float(np.max(f[mobile&sel_ag])),'max_mobile_skeleton_force_after_releasing_skeleton_constraint_eV_A':float(np.max(f[mobile&sel_sk])),'max_Ag_displacement_A':float(np.max(disp[mobile&sel_ag])),'by_mobile_element':mobile_stats,'minimum_pair_distances_A':minpair,'connectivity_after_relaxation':{'Ti3SiC2':topology(xyz[sel_sk],ids[sel_sk],np.array([ELEMENT[int(v)] for v in t[sel_sk]]),3.25,pbox),'Ag':topology(xyz[sel_ag],ids[sel_ag],np.array(['Ag']*np.sum(sel_ag)),3.15,pbox)},'limitations':['Ag was minimized with the entire Ti3SiC2 scaffold fixed; scaffold forces are only evaluated after releasing that constraint.','100 FIRE iterations; no finite-temperature dynamics or TTM.','This is a geometric-confinement sensitivity screen under calibration-v0, not a validated material baseline.']}
(HERE/'rigid_scaffold_ag_minimize_metrics.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
