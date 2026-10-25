#!/usr/bin/env python3
"""Summarize cumulative 50 fs no-source NVT for exploratory geometry Model B."""
import gzip,json,re
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
HERE=Path(__file__).resolve().parent;MODEL=HERE.parents[1]
GEOM=json.loads((MODEL/'potential_validation/relaxed_backbone_candidate/geometry_report.json').read_text())
BOX=np.asarray(GEOM['simulation_box_A'],float);TREEBOX=[BOX[0],BOX[1],1000.]
def load(name):
 with gzip.open(HERE/name,'rt') as f:return np.loadtxt(f,skiprows=9)
def stats(x):return {'p50':float(np.percentile(x,50)),'p95':float(np.percentile(x,95)),'p99':float(np.percentile(x,99)),'max':float(np.max(x))}
def rows(name):
 pat=re.compile(r'^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$');out=[]
 for line in (HERE/name).read_text().splitlines():
  m=pat.match(line)
  if m:out.append(list(map(float,m.groups())))
 return np.asarray(out)
def topology(xyz,cut):
 p=cKDTree(xyz,boxsize=TREEBOX).query_pairs(cut,output_type='ndarray');g=sparse.coo_matrix((np.ones(len(p)*2,dtype=np.uint8),(np.r_[p[:,0],p[:,1]],np.r_[p[:,1],p[:,0]])),shape=(len(xyz),len(xyz))).tocsr();n,l=connected_components(g,directed=False);s=np.bincount(l);k=int(s.argmax());return {'components':int(n),'largest_component_atoms':int(s[k]),'atoms_outside_largest':int(len(xyz)-s[k]),'cutoff_A':cut}
start=load('ag_minimized_100.dump.gz');end=load('ag_nvt_0p05ps_final.dump.gz')
lookup=np.searchsorted(start[:,0],end[:,0]);assert np.array_equal(start[lookup,0],end[:,0])
d=end[:,2:5]-start[lookup,2:5];d[:,0]-=BOX[0]*np.rint(d[:,0]/BOX[0]);d[:,1]-=BOX[1]*np.rint(d[:,1]/BOX[1]);disp=np.linalg.norm(d,axis=1);force=np.linalg.norm(end[:,8:11],axis=1)
thermo=np.vstack([rows('log.ag_nvt_0p025ps'),rows('log.ag_nvt_continue_to_0p05ps')]);first=thermo[0];last=thermo[-1]
bottom=(start[:,1]==4)&(start[:,4]<20.0);allag=np.vstack([start[bottom,:5],end[:,:5]]);_,ix=np.unique(allag[:,0],return_index=True);allag=allag[np.sort(ix)];xyz=allag[:,2:5].copy();xyz[:,0]%=BOX[0];xyz[:,1]%=BOX[1]
report={'status':'GEOMETRY_CONSTRAINED_WCA_NO_SOURCE_NVT_50FS_EXPLORATORY','duration_fs':25.0,'cumulative_duration_fs':50.0,'timestep_fs':0.5,'thermo_initial':{'Ag_temperature_K':float(first[1]),'potential_energy_eV':float(first[2]),'kinetic_energy_eV':float(first[3]),'total_energy_eV':float(first[4]),'max_Ag_force_eV_A':float(first[5])},'thermo_final':{'Ag_temperature_K':float(last[1]),'potential_energy_eV':float(last[2]),'kinetic_energy_eV':float(last[3]),'total_energy_eV':float(last[4]),'max_Ag_force_eV_A':float(last[5])},'energy_change_eV':{'potential':float(last[2]-first[2]),'kinetic':float(last[3]-first[3]),'total':float(last[4]-first[4])},'mobile_Ag_displacement_from_minimized_state_A':stats(disp),'mobile_Ag_force_eV_A':stats(force),'Ag_connectivity_after_test':topology(xyz,3.15),'limitations':['No external heat or TTM; skeleton and bottom slab fixed.','50 fs is only a short cold-start screen, not an equilibrium or melt-pool validation.','WCA contact wall is exploratory, not calibrated interface physics.']}
(HERE/'nvt_0p05ps_metrics.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
