#!/usr/bin/env python3
from __future__ import annotations
import csv, importlib.util, json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parent; CONTROL=ROOT.parent/'noimpact_control_1ps_after105ps_20261008'
spec=importlib.util.spec_from_file_location('al',ROOT/'analyze_extension.py'); al=importlib.util.module_from_spec(spec); spec.loader.exec_module(al)
def initial_bonds(frame):
 ids=al.column(frame,'id').astype(np.int64);typ=al.column(frame,'type').astype(np.int64);xyz=frame['rows'][:,[frame['cols'].index(c) for c in ('x','y','z')]]
 ag=np.flatnonzero(typ==4);rho=np.linalg.norm(al.min_image_xy(xyz[ag,:2],frame),axis=1);centers=ag[(rho<8)&(xyz[ag,2]>=12)&(xyz[ag,2]<55)]
 avec=np.array([frame['a'][0],frame['a'][1],0]);bvec=np.array([frame['b'][0],frame['b'][1],0]);shifts=np.array([i*avec+j*bvec for i in (-1,0,1) for j in (-1,0,1)])
 points=np.concatenate([xyz[ag]+s for s in shifts]);pointids=np.tile(ids[ag],len(shifts));tree=cKDTree(points);pairs=set()
 for ix in centers:
  for k in tree.query_ball_point(xyz[ix],3.5):
   other=int(pointids[k]);center=int(ids[ix])
   if other!=center:pairs.add(tuple(sorted((center,other))))
 return np.asarray(sorted(pairs),dtype=np.int64),int(len(centers))
def retained(frame,pairs):
 ids=al.column(frame,'id').astype(np.int64);xyz=frame['rows'][:,[frame['cols'].index(c) for c in ('x','y','z')]];order=np.argsort(ids);ids=ids[order];xyz=xyz[order];ix=np.searchsorted(ids,pairs)
 if np.any(ix>=len(ids)) or not np.array_equal(ids[ix],pairs):return float('nan')
 dv=xyz[ix[:,0]]-xyz[ix[:,1]];shifts=np.array([i*frame['a']+j*frame['b'] for i in (-1,0,1) for j in (-1,0,1)]);cand=dv[:,None,:2]+shifts[None,:,:];k=np.argmin(np.sum(cand*cand,axis=2),axis=1);dxy=cand[np.arange(len(dv)),k];return float(np.mean(np.sum(dxy*dxy,axis=1)+dv[:,2]**2<=3.5**2))
def trajectory(branch,prefix):
 out=[]
 for hit in range(472,477):
  frames=al.read_dump(branch/'dumps'/f'{prefix}_{hit}_62eV.dump.gz')
  if len(frames)!=41:raise RuntimeError('unexpected segment frame count')
  out.extend(frames if hit==472 else frames[1:])
 if len(out)!=201:raise RuntimeError(f'expected 201 frames, got {len(out)}')
 return out
impact=trajectory(ROOT,'impact');control=trajectory(CONTROL,'control');pairs,ncenter=initial_bonds(impact[0]);p2,n2=initial_bonds(control[0])
if ncenter!=n2 or not np.array_equal(pairs,p2):raise RuntimeError('initial cohorts differ')
series={k:[retained(f,pairs) for f in frames] for k,frames in [('impact',impact),('noimpact',control)]}
ref_path=ROOT.parent/'ag_eam_phase_references_20261008/analysis/ag_reference_comparison.json';refs=json.loads(ref_path.read_text())['references']
summary={'initial_local_Ag_centers':ncenter,'initial_unique_Ag_Ag_pairs':len(pairs),'neighbor_cutoff_A':3.5,'retained_fraction':{k:{f'{lag:g}_ps':float(v[round(lag/.005)]) for lag in (.1,.5,1.)} for k,v in series.items()},'pure_Ag_reference':{k:refs[k]['neighbor_pair_retention'] for k in ('solid_300K','hot_1800K')}}
(ROOT/'analysis_neighbor_survival_extension.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
with (ROOT/'analysis_neighbor_survival_extension.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['time_ps','lag_from_105ps','impact_retained_fraction','noimpact_retained_fraction'])
 for i,(a,b) in enumerate(zip(series['impact'],series['noimpact'])):w.writerow([105+i*.005,i*.005,a,b])
print(json.dumps(summary,indent=2,ensure_ascii=False))
