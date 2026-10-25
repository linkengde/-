#!/usr/bin/env python3
"""Resolve selected atoms' Morse cross-force vectors by neighboring species."""
from __future__ import annotations
import gzip,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
GEOM=json.loads((HERE/'geometry_report.json').read_text())
BOX=np.asarray(GEOM['simulation_box_A'],float)
TREE_BOX=[BOX[0],BOX[1],1000.]
ELEMENT={1:'Ti',2:'Si',3:'C',4:'Ag'}
PARAM={1:(0.94814194,2.0,2.33708864),2:(0.93944948,2.0,2.40407837),3:(1.35077230,2.0,2.15255924)}

def load(name):
    with gzip.open(HERE/name,'rt') as f:return np.loadtxt(f,skiprows=9)
full=load('forces_after300_FULL.dump.gz'); nocross=load('forces_after300_NO_CROSS.dump.gz')
if not np.array_equal(full[:,:2],nocross[:,:2]):raise RuntimeError('ID/type mismatch')
ids=full[:,0].astype(int); types=full[:,1].astype(int); pos=full[:,2:5].copy(); pos[:,0]%=BOX[0];pos[:,1]%=BOX[1]
full_force=full[:,5:8]; no_force=nocross[:,5:8]; cross=full_force-no_force
sk_idx=np.flatnonzero(types!=4); ag_idx=np.flatnonzero(types==4)
sk_tree=cKDTree(pos[sk_idx],boxsize=TREE_BOX);ag_tree=cKDTree(pos[ag_idx],boxsize=TREE_BOX)
loc=json.loads((HERE/'high_force_localization_after300.json').read_text())

results=[]
for center_name in ('Ag','Ti','Si','C'):
    for entry in loc['by_element'][center_name]['top_10_full_force_atoms'][:5]:
        atom_id=entry['id']; i=int(np.searchsorted(ids,atom_id)); ctype=int(types[i]);
        tree=sk_tree if ctype==4 else ag_tree
        pool=sk_idx if ctype==4 else ag_idx
        neighbors=tree.query_ball_point(pos[i],r=8.0)
        records=[]; by_species={}
        for local in neighbors:
            j=int(pool[local]); jtype=int(types[j]); D,alpha,r0=PARAM[ctype if ctype!=4 else jtype]
            dr=pos[i]-pos[j];dr[0]-=BOX[0]*np.rint(dr[0]/BOX[0]);dr[1]-=BOX[1]*np.rint(dr[1]/BOX[1])
            r=float(np.linalg.norm(dr))
            if r==0:raise RuntimeError('Overlapping atoms')
            ex=np.exp(-alpha*(r-r0)); dudr=2*D*alpha*ex*(1.0-ex)
            force=(-dudr/r)*dr
            species=ELEMENT[jtype]
            by_species.setdefault(species,np.zeros(3));by_species[species]+=force
            records.append({'neighbor_id':int(ids[j]),'neighbor_element':species,'distance_A':r,'morse_force_on_center_eV_A':float(np.linalg.norm(force)),'force_vector_eV_A':[float(x) for x in force]})
        predicted=np.sum(np.asarray([x['force_vector_eV_A'] for x in records]),axis=0) if records else np.zeros(3)
        cross_vec=cross[i]
        records.sort(key=lambda x:x['morse_force_on_center_eV_A'],reverse=True)
        results.append({'center_id':atom_id,'center_element':center_name,'position_A':[float(x) for x in full[i,2:5]],'full_force_eV_A':float(np.linalg.norm(full_force[i])),'no_cross_force_eV_A':float(np.linalg.norm(no_force[i])),'observed_cross_force_vector_eV_A':[float(x) for x in cross_vec],'observed_cross_force_eV_A':float(np.linalg.norm(cross_vec)),'summed_pairwise_morse_vector_eV_A':[float(x) for x in predicted],'pair_sum_residual_eV_A':float(np.linalg.norm(predicted-cross_vec)),'cross_force_by_neighbor_element_eV_A':{k:{'vector':[float(x) for x in v],'magnitude':float(np.linalg.norm(v))} for k,v in by_species.items()},'largest_individual_neighbor_contributions':records[:8]})
report={'status':'TOP_ATOM_PAIR_FORCE_DECOMPOSITION','method':'Analytic pair Morse force at the same dumped coordinates; exact same calibration-v0 D0/alpha/r0 values, summed over opposite-phase neighbors within the 8 A Morse cutoff. Compared with full-minus-zero-cross LAMMPS force vectors.','pair_sum_residual_max_eV_A':max(x['pair_sum_residual_eV_A'] for x in results),'atoms':results}
(HERE/'top_atom_morse_force_decomposition.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
