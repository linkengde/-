#!/usr/bin/env python3
"""Summarize the 25 fs no-source NVT screen for exploratory geometry Model B."""
import gzip,json,re
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
HERE=Path(__file__).resolve().parent;MODEL=HERE.parents[1]
GEOM=json.loads((MODEL/'potential_validation/relaxed_backbone_candidate/geometry_report.json').read_text())
BOX=np.asarray(GEOM['simulation_box_A'],float); TREEBOX=[BOX[0],BOX[1],1000.]

def load(name):
 with gzip.open(HERE/name,'rt') as f:return np.loadtxt(f,skiprows=9)
def stats(x):return {'p50':float(np.percentile(x,50)),'p95':float(np.percentile(x,95)),'p99':float(np.percentile(x,99)),'max':float(np.max(x))}
def topology(xyz,cut):
 p=cKDTree(xyz,boxsize=TREEBOX).query_pairs(cut,output_type='ndarray')
 g=sparse.coo_matrix((np.ones(len(p)*2,dtype=np.uint8),(np.r_[p[:,0],p[:,1]],np.r_[p[:,1],p[:,0]])),shape=(len(xyz),len(xyz))).tocsr()
 n,l=connected_components(g,directed=False);s=np.bincount(l);k=int(s.argmax())
 return {'components':int(n),'largest_component_atoms':int(s[k]),'atoms_outside_largest':int(len(xyz)-s[k]),'cutoff_A':cut}

def thermo_rows():
 pat=re.compile(r'^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$')
 out=[]
 for line in (HERE/'log.ag_nvt_0p025ps').read_text().splitlines():
  m=pat.match(line)
  if m:out.append(list(map(float,m.groups())))
 return np.asarray(out)

start=load('ag_minimized_100.dump.gz');end=load('ag_nvt_0p025ps_final.dump.gz')
if not np.all(end[:,1]==4):raise RuntimeError('NVT endpoint must contain Ag')
lookup=np.searchsorted(start[:,0],end[:,0])
if not np.array_equal(start[lookup,0],end[:,0]):raise RuntimeError('Ag IDs mismatch')
delta=end[:,2:5]-start[lookup,2:5];delta[:,0]-=BOX[0]*np.rint(delta[:,0]/BOX[0]);delta[:,1]-=BOX[1]*np.rint(delta[:,1]/BOX[1]);disp=np.linalg.norm(delta,axis=1)
f=np.linalg.norm(end[:,8:11],axis=1)
rows=thermo_rows();first,last=rows[0],rows[-1]
# combine fixed bottom Ag positions with the mobile Ag NVT endpoint
bottom=(start[:,1]==4)&(start[:,4]<20.0);all_ag=np.vstack([start[bottom,:5],end[:,:5]])
_,unique=np.unique(all_ag[:,0],return_index=True);all_ag=all_ag[np.sort(unique)]
pos=all_ag[:,2:5].copy();pos[:,0]%=BOX[0];pos[:,1]%=BOX[1]
report={'status':'GEOMETRY_CONSTRAINED_WCA_NO_SOURCE_NVT_25FS_EXPLORATORY','duration_fs':float((last[0]-first[0])*0.5),'timestep_fs':0.5,'thermo_initial':{'Ag_temperature_K':float(first[1]),'potential_energy_eV':float(first[2]),'kinetic_energy_eV':float(first[3]),'total_energy_eV':float(first[4]),'max_Ag_force_eV_A':float(first[5])},'thermo_final':{'Ag_temperature_K':float(last[1]),'potential_energy_eV':float(last[2]),'kinetic_energy_eV':float(last[3]),'total_energy_eV':float(last[4]),'max_Ag_force_eV_A':float(last[5])},'energy_change_eV':{'potential':float(last[2]-first[2]),'kinetic':float(last[3]-first[3]),'total':float(last[4]-first[4])},'mobile_Ag_displacement_A':stats(disp),'mobile_Ag_force_eV_A':stats(f),'Ag_connectivity_after_test':topology(pos,3.15),'limitations':['No external heat source or TTM; skeleton and bottom slab were fixed.','Only 25 fs; not sufficient to establish equilibrium or melting.','The WCA interaction is an exploratory geometric wall, not calibrated interfacial physics.']}
(HERE/'nvt_screen_metrics.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
