#!/workspace/venvs/dft-labels/bin/python
"""Extract C-C precursor at step 18 from the v12m cold replay."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import write
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'cold_v12m_cc_hotspot_trace/C_C_local_by_step.dump'
TARGET_IDS=[128985,134363]
STEP=18
RADIUS_A=5.5
CELL_A=22.0
lines=SOURCE.read_text().splitlines(); frame=None; i=0
while i<len(lines):
    if lines[i]!='ITEM: TIMESTEP': i+=1; continue
    step=int(lines[i+1]); n=int(lines[i+3]); h=i+8; cols=lines[h].split()[2:]
    rows=np.loadtxt(lines[h+1:h+1+n],dtype=float)
    if step==STEP: frame=(cols,rows); break
    i=h+1+n
if frame is None: raise RuntimeError(f'Missing timestep {STEP}')
cols,rows=frame
ix={k:cols.index(k) for k in ('id','type','x','y','z','fx','fy','fz')}
ids=rows[:,ix['id']].astype(int); types=rows[:,ix['type']].astype(int)
xyz=rows[:,[ix[k] for k in ('x','y','z')]]; forces=rows[:,[ix[k] for k in ('fx','fy','fz')]]
lookup={int(atom_id):j for j,atom_id in enumerate(ids)}
if any(a not in lookup for a in TARGET_IDS): raise RuntimeError('C-C targets missing from local dump')
if [types[lookup[a]] for a in TARGET_IDS] != [2,2]: raise RuntimeError('Target atoms are not both C')
anchor=xyz[lookup[TARGET_IDS[0]]]
unwrapped=[]
for atom_id in TARGET_IDS:
    delta=xyz[lookup[atom_id]]-anchor
    delta[:2]-=200*np.rint(delta[:2]/200)
    unwrapped.append(anchor+delta)
unwrapped=np.asarray(unwrapped); center=unwrapped.mean(axis=0)
relative=xyz-center; relative[:,:2]-=200*np.rint(relative[:,:2]/200)
keep=np.linalg.norm(relative,axis=1)<=RADIUS_A
kept_ids=ids[keep]; kept_types=types[keep]
if len(kept_ids)>64 or not np.isin(TARGET_IDS,kept_ids).all(): raise RuntimeError('Bad crop size/target inclusion')
target=np.isin(kept_ids,TARGET_IDS).astype(np.int8)
atoms=Atoms(symbols=[{1:'Ag',2:'C',3:'Si',4:'Ti'}[int(t)] for t in kept_types],positions=relative[keep]+CELL_A/2,cell=[CELL_A]*3,pbc=True)
atoms.new_array('source_atom_id',kept_ids.astype(np.int32)); atoms.new_array('POD_v12m_forces',forces[keep]); atoms.new_array('is_DFT_force_target',target)
dist=float(np.linalg.norm(unwrapped[0]-unwrapped[1]))
atoms.info.update({'config_type':'full200_v12m_POD_cold_C-C_collapse_precursor','source_timestep':STEP,'cumulative_source_time_fs':10+STEP,'center_atom_ids':json.dumps(TARGET_IDS),'elements':json.dumps(['C','C']),'center_pair_distance_A':dist,'crop_radius_A':RADIUS_A,'crop_cell_A':CELL_A,'geometry_source':str(SOURCE),'source_pbc':'x/y periodic, z open; local crop embedded in cubic vacuum box','DFT_target_rows':2,'label_role':'v12n training candidate; Si01 center 104359 remains held out'})
out=HERE/'candidates/cc_pair_v12m_t28fs_step18.extxyz'; out.parent.mkdir(exist_ok=True); write(out,atoms,format='extxyz')
print(json.dumps({'status':'EXTRACTED','output':str(out),'natoms':len(atoms),'target_ids':TARGET_IDS,'target_rows':int(target.sum()),'distance_A':dist,'types':[int(types[lookup[a]]) for a in TARGET_IDS]},indent=2))
