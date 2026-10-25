#!/usr/bin/env python3
"""Analyze the bounded continuous 0.1 ps B-route cold-start test."""
import gzip,json,re
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
HERE=Path(__file__).resolve().parent;MODEL=HERE.parents[1]
GEOM=json.loads((MODEL/'potential_validation/relaxed_backbone_candidate/geometry_report.json').read_text())
BOX=np.asarray(GEOM['simulation_box_A'],float);TREEBOX=[BOX[0],BOX[1],1000.]
def load(n):
 with gzip.open(HERE/n,'rt') as f:return np.loadtxt(f,skiprows=9)
def stats(x):return {'p50':float(np.percentile(x,50)),'p95':float(np.percentile(x,95)),'p99':float(np.percentile(x,99)),'max':float(np.max(x))}
def thermo():
 pat=re.compile(r'^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$');out=[]
 for line in (HERE/'log.ag_nvt_0p1ps_continuous').read_text().splitlines():
  m=pat.match(line)
  if m:out.append(list(map(float,m.groups())))
 return np.asarray(out)
def topology(xyz,cut):
 p=cKDTree(xyz,boxsize=TREEBOX).query_pairs(cut,output_type='ndarray');g=sparse.coo_matrix((np.ones(len(p)*2,dtype=np.uint8),(np.r_[p[:,0],p[:,1]],np.r_[p[:,1],p[:,0]])),shape=(len(xyz),len(xyz))).tocsr();n,l=connected_components(g,directed=False);s=np.bincount(l);k=int(s.argmax());return {'components':int(n),'largest_component_atoms':int(s[k]),'atoms_outside_largest':int(len(xyz)-s[k]),'cutoff_A':cut}
start=load('ag_minimized_100.dump.gz');end=load('ag_nvt_0p1ps_continuous_final.dump.gz');assert np.all(end[:,1]==4)
ix=np.searchsorted(start[:,0],end[:,0]);assert np.array_equal(start[ix,0],end[:,0]);d=end[:,2:5]-start[ix,2:5];d[:,0]-=BOX[0]*np.rint(d[:,0]/BOX[0]);d[:,1]-=BOX[1]*np.rint(d[:,1]/BOX[1]);disp=np.linalg.norm(d,axis=1);force=np.linalg.norm(end[:,8:11],axis=1)
r=thermo();first,last=r[0],r[-1]
bottom=(start[:,1]==4)&(start[:,4]<20);ag=np.vstack([start[bottom,:5],end[:,:5]]);_,u=np.unique(ag[:,0],return_index=True);ag=ag[np.sort(u)];xyz=ag[:,2:5].copy();xyz[:,0]%=BOX[0];xyz[:,1]%=BOX[1]
report={'status':'GEOMETRY_CONSTRAINED_WCA_NO_SOURCE_NVT_100FS_EXPLORATORY_NOT_EQUILIBRATED','duration_fs':float((last[0]-first[0])*0.5),'timestep_fs':0.5,'thermo_initial':{'Ag_temperature_K':float(first[1]),'potential_energy_eV':float(first[2]),'kinetic_energy_eV':float(first[3]),'total_energy_eV':float(first[4]),'max_Ag_force_eV_A':float(first[5])},'thermo_final':{'Ag_temperature_K':float(last[1]),'potential_energy_eV':float(last[2]),'kinetic_energy_eV':float(last[3]),'total_energy_eV':float(last[4]),'max_Ag_force_eV_A':float(last[5])},'energy_change_eV':{'potential':float(last[2]-first[2]),'kinetic':float(last[3]-first[3]),'total':float(last[4]-first[4])},'mobile_Ag_displacement_A':stats(disp),'mobile_Ag_force_eV_A':stats(force),'Ag_connectivity_after_test':topology(xyz,3.15),'limitations':['The WCA wall is an exploratory geometry-only interaction, not calibrated interface physics.','The Ti3SiC2 skeleton and 0-20 A bottom slab are fixed; their thermal response/heat transport is absent.','The 100 fs trajectory avoids catastrophic force/temperature runaway but does not equilibrate at the 300 K thermostat target.','No TTM, external heating, or melting was tested.']}
(HERE/'nvt_0p1ps_metrics.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
