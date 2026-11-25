#!/usr/bin/env python3
"""Track high-z Ag atoms and Ag connectivity in the local-pulse pilot."""
from pathlib import Path
import csv, json
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT=Path(__file__).resolve().parent
FILES={
 'initial_1p7ps':ROOT.parent/'simple_ag_fixed_skeleton_20261125/reheat_melt_assessment_20261125/trajectories/ag_reheat_q6_ptm_msd_final.dump',
 'pulse_3000K_1ps':ROOT/'trajectory_1ps/final_all_atoms_1ps.dump',
 'pulse_9000K_1ps':ROOT/'trajectory_9000K/final_all_atoms_9000K.dump',
 'pulse_22500K_1ps':ROOT/'trajectory_22500K/final_all_atoms_22500K.dump',
 'pulse_surface_local_2ps':ROOT/'trajectory_surface_local_2ps/final_all_atoms_surface_local.dump',
}

def read_dump(path):
    with path.open() as f:
        while f.readline().strip()!='ITEM: NUMBER OF ATOMS': pass
        n=int(f.readline()); f.readline()
        bounds=np.array([[float(v) for v in f.readline().split()[:2]] for _ in range(3)])
        names=f.readline().split()[2:]
        arr=np.loadtxt(f,max_rows=n)
    return {name:i for i,name in enumerate(names)},arr,bounds

results={}; atom_states={}
for label,path in FILES.items():
    idx,all_atoms,bounds=read_dump(path)
    ag=all_atoms[all_atoms[:,idx['type']]==4]
    ids=ag[:,idx['id']].astype(int)
    xyz=ag[:,[idx['x'],idx['y'],idx['z']]].copy()
    xy=xyz[:,:2].copy()%200.0
    z=xyz[:,2].copy()
    r=np.hypot((xy[:,0]-100+100)%200-100,(xy[:,1]-100+100)%200-100)
    tree=cKDTree(np.column_stack([xy,z]),boxsize=[200.0,200.0,1000.0])
    pairs=tree.query_pairs(3.8,output_type='ndarray')
    adjacency=coo_matrix((np.ones(len(pairs)*2,dtype=np.uint8),
        (np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(ag),len(ag))).tocsr()
    ncomp,labels=connected_components(adjacency,directed=False)
    sizes=np.bincount(labels)
    high=z>250.0
    result={
      'N_Ag':int(len(ag)),
      'n_connected_components_at_3p8A':int(ncomp),
      'largest_Ag_component':int(sizes.max()),
      'Ag_outside_largest_component':int(len(ag)-sizes.max()),
      'N_Ag_z_gt_220A':int((z>220).sum()),
      'N_Ag_z_gt_250A':int(high.sum()),
      'N_Ag_z_gt_220A_inside_source_r25A':int(((z>220)&(r<25)).sum()),
      'N_Ag_z_gt_250A_inside_source_r25A':int((high&(r<25)).sum()),
      'Ag_ids_z_gt_250A':ids[high].tolist(),
    }
    results[label]=result
    atom_states[label]={int(atom_id):row for atom_id,row in zip(ids,ag)}

initial=atom_states['initial_1p7ps']
idx_final,final_arr,_=read_dump(FILES['pulse_9000K_1ps'])
final_ag=final_arr[final_arr[:,idx_final['type']]==4]
tracks=[]
for atom_id,row in atom_states['pulse_9000K_1ps'].items():
    if row[idx_final['z']]<=250: continue
    prev=initial.get(atom_id)
    tracks.append({
      'id':atom_id,
      'initial_x_A':float(prev[2]) if prev is not None else None,
      'initial_y_A':float(prev[3]) if prev is not None else None,
      'initial_z_A':float(prev[4]) if prev is not None else None,
      'final_x_A':float(row[2]),'final_y_A':float(row[3]),'final_z_A':float(row[4]),
      'distance_from_pulse_xy_A':float(np.hypot((row[2]-100+100)%200-100,(row[3]-100+100)%200-100)),
      'same_atom_already_z_gt_220_initial':bool(prev is not None and prev[4]>220),
    })
    # row layout is id,type,x,y,z,vx,vy,vz
(ROOT/'results/ejection_cluster_summary.json').write_text(json.dumps(results,indent=2)+'\n')
with (ROOT/'results/high_z_ag_id_tracks_9000K.csv').open('w',newline='') as f:
    fields=['id','initial_x_A','initial_y_A','initial_z_A','final_x_A','final_y_A','final_z_A','distance_from_pulse_xy_A','same_atom_already_z_gt_220_initial']
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(tracks)
print(json.dumps(results,indent=2))
print(json.dumps(tracks,indent=2))
