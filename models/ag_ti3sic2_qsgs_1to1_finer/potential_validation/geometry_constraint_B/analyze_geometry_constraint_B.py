#!/usr/bin/env python3
"""Check geometry, forces, and WCA activity for exploratory Model B."""
import gzip,json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
HERE=Path(__file__).resolve().parent
MODEL=HERE.parents[1]
GEOM=json.loads((MODEL/'potential_validation/relaxed_backbone_candidate/geometry_report.json').read_text())
BOX=np.asarray(GEOM['simulation_box_A'],float); TREEBOX=[BOX[0],BOX[1],1000.]
MASS={1:47.867,2:28.085,3:12.011,4:107.8682};EL={1:'Ti',2:'Si',3:'C',4:'Ag'}
RC={1:2.33708864,2:2.40407837,3:2.15255924};EPS=0.5

def load(name):
 with gzip.open(HERE/name,'rt') as f:return np.loadtxt(f,skiprows=9)
def force_stats(f):return {'p50':float(np.percentile(f,50)),'p95':float(np.percentile(f,95)),'p99':float(np.percentile(f,99)),'max':float(np.max(f)),'n_gt_2':int(np.sum(f>2.0))}
def topology(pos,cut):
 pairs=cKDTree(pos,boxsize=TREEBOX).query_pairs(cut,output_type='ndarray')
 graph=sparse.coo_matrix((np.ones(len(pairs)*2,dtype=np.uint8),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(pos),len(pos))).tocsr()
 n,labels=connected_components(graph,directed=False);sizes=np.bincount(labels);big=int(np.argmax(sizes))
 return {'components':int(n),'largest_component_atoms':int(sizes[big]),'atoms_outside_largest':int(len(pos)-sizes[big]),'cutoff_A':cut}

ini=load('initial_forces.dump.gz');fin=load('ag_minimized_100.dump.gz')
if not np.array_equal(ini[:,:2],fin[:,:2]):raise RuntimeError('Atom ID/type mismatch')
if np.max(np.abs(ini[:,4]-fin[:,4]))>2:raise RuntimeError('Unexpected z displacement')
typeid=fin[:,1].astype(int);pos=fin[:,2:5].copy();pos[:,0]%=BOX[0];pos[:,1]%=BOX[1]
f=np.linalg.norm(fin[:,5:8],axis=1)
delta=fin[:,2:5]-ini[:,2:5];delta[:,0]-=BOX[0]*np.rint(delta[:,0]/BOX[0]);delta[:,1]-=BOX[1]*np.rint(delta[:,1]/BOX[1]);disp=np.linalg.norm(delta,axis=1)
mobile=(typeid==4)&(ini[:,4]>=20.0);mobile_sk=(typeid!=4)&(ini[:,4]>=20.0)
minpairs={}
for a,b in [('Ag','Ti'),('Ag','Si'),('Ag','C'),('Ti','C'),('Ag','Ag')]:
 ta=next(t for t,e in EL.items() if e==a);tb=next(t for t,e in EL.items() if e==b)
 tree=cKDTree(pos[typeid==tb],boxsize=TREEBOX)
 distances=tree.query(pos[typeid==ta],k=2 if a==b else 1)[0]
 minpairs[a+'-'+b]=float(distances[:,1].min() if a==b else distances.min())
pair_idx=cKDTree(pos,boxsize=TREEBOX).query_pairs(max(RC.values()),output_type='ndarray')
wca_by={EL[t]:{'contacts':0,'energy_eV':0.0,'max_pair_force_eV_A':0.0,'minimum_distance_A':None} for t in RC}
for i,j in pair_idx:
 ti,tj=int(typeid[i]),int(typeid[j])
 if ti==4 and tj in RC:cer=tj
 elif tj==4 and ti in RC:cer=ti
 else:continue
 dr=pos[i]-pos[j];dr[0]-=BOX[0]*np.rint(dr[0]/BOX[0]);dr[1]-=BOX[1]*np.rint(dr[1]/BOX[1]);r=float(np.linalg.norm(dr));rc=RC[cer]
 if r>=rc:continue
 sigma=rc/(2**(1/6));x=sigma/r;x6=x**6;x12=x6*x6
 u=4*EPS*(x12-x6)+EPS;force=24*EPS/r*(2*x12-x6)
 rec=wca_by[EL[cer]];rec['contacts']+=1;rec['energy_eV']+=u;rec['max_pair_force_eV_A']=max(rec['max_pair_force_eV_A'],force);rec['minimum_distance_A']=r if rec['minimum_distance_A'] is None else min(rec['minimum_distance_A'],r)
report={'status':'GEOMETRY_CONSTRAINED_WCA_EXPLORATORY_STATIC_SCREEN','model_atoms':int(len(fin)),'counts':{EL[t]:int(np.sum(typeid==t)) for t in EL},'mass_ratio_Ag_over_Ti3SiC2':float(MASS[4]*np.sum(typeid==4)/sum(MASS[t]*np.sum(typeid==t) for t in (1,2,3))),'fixed_skeleton':True,'bottom_test_constraint_z_0_to_20_A':True,'wca_parameterization':{'epsilon_eV':EPS,'zero_force_cutoffs_A':{EL[t]:RC[t] for t in RC},'status':'EXPLORATORY GEOMETRIC WALL, NOT MATERIAL CALIBRATION'},'mobile_Ag_force_eV_A':force_stats(f[mobile]),'mobile_skeleton_force_after_releasing_scaffold_eV_A':force_stats(f[mobile_sk]),'mobile_Ag_displacement_A':{'p50':float(np.percentile(disp[mobile],50)),'p95':float(np.percentile(disp[mobile],95)),'max':float(np.max(disp[mobile]))},'minimum_pair_distances_A':minpairs,'phase_connectivity':{'Ti3SiC2':topology(pos[typeid!=4],3.25),'Ag':topology(pos[typeid==4],3.15)},'active_WCA_pairs':wca_by,'limitations':['Only 100 FIRE iterations; the force tolerance was not reached.','Rigid ceramic has no thermal response or heat transport.','WCA contact distances reuse Morse r0 and epsilon=0.5 eV is an exploratory choice; this is not an Ag/Ti3SiC2 interface potential.','No NVT or TTM was run by this static report.']}
(HERE/'static_screen_metrics.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
