#!/usr/bin/env python3
"""Summarize the 25 fs no-source Ag NVT test after the 300-step static screen."""
import gzip
import json
import re
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
BOX = np.asarray(json.loads((HERE / 'geometry_report.json').read_text())['simulation_box_A'], float)
TREE_BOX = [BOX[0], BOX[1], 1000.]

def load_dump(name):
    with gzip.open(HERE / name, 'rt') as f:
        return np.loadtxt(f, skiprows=9)

def load_log():
    rows=[]
    pattern=re.compile(r'^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$')
    for line in (HERE/'log.rigid_skeleton_Ag_NVT_after300_0p025ps').read_text().splitlines():
        m=pattern.match(line)
        if m:
            v=list(map(float,m.groups())); rows.append(v)
    if not rows: raise RuntimeError('No thermo rows parsed')
    return np.asarray(rows)

def stats(v):
    return {'p50':float(np.percentile(v,50)),'p95':float(np.percentile(v,95)),'p99':float(np.percentile(v,99)),'max':float(np.max(v))}

start=load_dump('ag_only_rigid_skeleton_continue_200.dump.gz')
end=load_dump('rigid_skeleton_Ag_NVT_after300_0p025ps_final.dump.gz')
if not np.all(end[:,1]==4): raise RuntimeError('Endpoint dump should contain only Ag')
lookup=np.searchsorted(start[:,0],end[:,0])
if not np.array_equal(start[lookup,0],end[:,0]): raise RuntimeError('Ag IDs not found in starting snapshot')
start_ag=start[lookup]; delta=end[:,2:5]-start_ag[:,2:5]
delta[:,0]-=BOX[0]*np.rint(delta[:,0]/BOX[0]); delta[:,1]-=BOX[1]*np.rint(delta[:,1]/BOX[1])
disp=np.linalg.norm(delta,axis=1); f=np.linalg.norm(end[:,8:11],axis=1)
thermo=load_log(); first,last=thermo[0],thermo[-1]
dPE=last[2]-first[2]; dKE=last[3]-first[3]; dE=last[4]-first[4]
# Rebuild total Ag geometry by combining mobile endpoint atoms with the frozen bottom Ag.
mobile_ids=set(end[:,0].astype(int).tolist())
bottom_ag=(start[:,1]==4)&(start[:,4]<20.)
fixed_bottom=start[bottom_ag]
all_ag=np.vstack([fixed_bottom[:,0:5],end[:,0:5]])
# Remove accidental ID duplicates defensively (the groups should be disjoint).
_, unique=np.unique(all_ag[:,0],return_index=True); all_ag=all_ag[np.sort(unique)]
xyz=all_ag[:,2:5].copy(); xyz[:,0]%=BOX[0]; xyz[:,1]%=BOX[1]
pairs=cKDTree(xyz,boxsize=TREE_BOX).query_pairs(3.15,output_type='ndarray')
graph=sparse.coo_matrix((np.ones(len(pairs)*2,dtype=np.uint8),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(xyz),len(xyz))).tocsr()
components,labels=connected_components(graph,directed=False); sizes=np.bincount(labels); largest=int(sizes.argmax())
report={
 'status':'NO_SOURCE_300K_NVT_25FS_TEST_NOT_A_COLD_BASELINE',
 'timesteps':{'initial':int(first[0]),'final':int(last[0]),'dt_fs':0.5,'segment_duration_fs':float((last[0]-first[0])*0.5)},
 'thermo_initial':{'Ag_temperature_K':float(first[1]),'potential_energy_eV':float(first[2]),'kinetic_energy_eV':float(first[3]),'total_energy_eV':float(first[4]),'max_mobile_Ag_force_eV_A':float(first[5])},
 'thermo_final':{'Ag_temperature_K':float(last[1]),'potential_energy_eV':float(last[2]),'kinetic_energy_eV':float(last[3]),'total_energy_eV':float(last[4]),'max_mobile_Ag_force_eV_A':float(last[5])},
 'energy_change_eV':{'potential':float(dPE),'kinetic':float(dKE),'total':float(dE),'fraction_of_potential_drop_appearing_as_kinetic_increase':float(dKE/(-dPE))},
 'mobile_Ag_displacement_A':stats(disp),
 'mobile_Ag_force_eV_A':stats(f),
 'Ag_connectivity_after_test':{'components':int(components),'largest_component_atoms':int(sizes[largest]),'atoms_outside_largest':int(len(xyz)-sizes[largest]),'cutoff_A':3.15},
 'limitations':['No external heat source or TTM; fixed skeleton and fixed bottom layer.','The short test measures continued structural energy release under the current calibration-v0 model; it does not establish liquid Ag.']
}
(HERE/'rigid_scaffold_ag_nvt_after300_metrics.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
