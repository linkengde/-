#!/usr/bin/env python3
"""Summarize high-force atoms and nearest-neighbor distances for the LAMMPS run-0 audit."""
import collections
import gzip
import json
import pathlib
import statistics
import numpy as np
from scipy.spatial import cKDTree

DATA = pathlib.Path('/workspace/-/models/ag_ti3sic2_qsgs_1to1/ag_ti3sic2_qsgs.data.gz')
DUMP = pathlib.Path(__file__).with_name('force_hotspots.lammpstrj')
OUT = pathlib.Path(__file__).with_name('force_audit_summary.json')

# Read the LAMMPS atomic data section.
coords=[]; types=[]; atom_ids=[]
with gzip.open(DATA, 'rt') as f:
    in_atoms=False
    for line in f:
        if line.startswith('Atoms #'):
            in_atoms=True
            continue
        if not in_atoms:
            continue
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        p=line.split()
        if len(p) < 5:
            if atom_ids:
                break
            continue
        try:
            atom_ids.append(int(p[0])); types.append(int(p[1])); coords.append([float(p[2]),float(p[3]),float(p[4])])
        except ValueError:
            if atom_ids: break
            continue
xyz=np.asarray(coords,dtype=np.float64)
types=np.asarray(types,dtype=np.int8)
ids=np.asarray(atom_ids,dtype=np.int64)

# Parse the thresholded force dump (LAMMPS screened v_fmag > 0.5 eV/A).
force_rows=[]
with DUMP.open() as f:
    lines=f.read().splitlines()
for line in lines[9:]:
    p=line.split()
    if len(p) >= 9:
        force_rows.append((int(p[0]),int(p[1]),float(p[2]),float(p[3]),float(p[4]),float(p[5])))
hot_ids=np.asarray([r[0] for r in force_rows],dtype=np.int64)
hot_types=np.asarray([r[1] for r in force_rows],dtype=np.int8)
hot_force=np.asarray([r[5] for r in force_rows],dtype=np.float64)
hot_z=np.asarray([r[4] for r in force_rows],dtype=np.float64)

# Enforce the LAMMPS x/y periodic boundaries and keep z effectively non-periodic.
box=np.array([299.635,301.857484616085,258.6115],dtype=np.float64)
query_xyz=xyz.copy(); query_xyz[:,2]+=box[2]
periodic_box=[box[0],box[1],3*box[2]]
tree=cKDTree(query_xyz,boxsize=periodic_box)
neighbor_d,neighbor_i=tree.query(query_xyz,k=2,workers=4)
nearest=neighbor_d[:,1]
# Data file atom IDs map to rows by explicit dictionary for robustness.
id_to_row={int(aid):i for i,aid in enumerate(ids)}
hot_idx=np.fromiter((id_to_row[int(a)] for a in hot_ids),dtype=np.int64,count=len(hot_ids))
hot_near=nearest[hot_idx]
hot_neighbor_types=types[neighbor_i[hot_idx,1]]

names={1:'Ti',2:'Si',3:'C',4:'Ag'}
summary={
  'atoms':int(len(ids)),
  'total_atoms_by_type':{names[int(k)]:int(v) for k,v in collections.Counter(types.tolist()).items()},
  'force_dump_threshold_eV_A':0.5,
  'dumped_atoms_above_0p5_eV_A':int(len(force_rows)),
  'force_threshold_counts_from_lammps_thermo':{'gt_0p5_eV_A':330452,'gt_1_eV_A':203536,'gt_5_eV_A':2143},
  'max_force_eV_A':float(hot_force.max()),
  'hot_force_counts_by_type':{},
  'nearest_neighbor_distance_all_A':{
     'min':float(nearest.min()),
     'p01':float(np.quantile(nearest,.01)),
     'p05':float(np.quantile(nearest,.05)),
     'median':float(np.median(nearest)),
     'count_lt_1p5_A':int(np.sum(nearest<1.5)),
     'count_lt_2p0_A':int(np.sum(nearest<2.0)),
  },
  'nearest_neighbor_distance_for_dumped_hot_atoms_A':{
     'min':float(hot_near.min()),
     'p01':float(np.quantile(hot_near,.01)),
     'p05':float(np.quantile(hot_near,.05)),
     'median':float(np.median(hot_near)),
     'count_lt_1p5_A':int(np.sum(hot_near<1.5)),
     'count_lt_2p0_A':int(np.sum(hot_near<2.0)),
  },
  'hot_atom_z_bins_10':{},
  'top_20_force_atoms':[],
  'caveat':'Nearest distances use x/y periodic boundaries and treat z as non-periodic; forces are model forces on the unrelaxed QSGS geometry, not fitting errors.'
}
for typ,name in names.items():
    mask=hot_types==typ
    total=int(np.sum(types==typ))
    summary['hot_force_counts_by_type'][name]={
      'total_atoms':total,
      'gt_0p5_eV_A':int(mask.sum()),
      'fraction_gt_0p5':float(mask.sum()/total) if total else None,
      'gt_1_eV_A':int(np.sum(mask&(hot_force>1.0))),
      'gt_5_eV_A':int(np.sum(mask&(hot_force>5.0))),
      'max_force_eV_A':float(hot_force[mask].max()) if mask.any() else None,
    }
for b in range(10):
    mask=(hot_z>=b*box[2]/10)&(hot_z<(b+1)*box[2]/10)
    summary['hot_atom_z_bins_10'][str(b)]={'z_interval_A':[b*box[2]/10,(b+1)*box[2]/10],'gt_0p5':int(mask.sum()),'gt_1':int(np.sum(mask&(hot_force>1.0)))}
for i in np.argsort(hot_force)[-20:][::-1]:
    row=hot_idx[i]; nearrow=neighbor_i[row,1]
    summary['top_20_force_atoms'].append({
      'id':int(hot_ids[i]),'type':names[int(hot_types[i])],
      'xyz_A':[float(x) for x in xyz[row]],'force_eV_A':float(hot_force[i]),
      'nearest_neighbor_distance_A':float(hot_near[i]),
      'nearest_neighbor_type':names[int(types[nearrow])],
    })
OUT.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
