#!/workspace/venvs/dft-labels/bin/python
"""Extract representative Ag-skeleton short contacts from SNAP v12o cold 10 fs."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import write

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'cold_snap_v12o_initial_458k_1fs'/'end_0p01ps.dump'
BOX=200.0
RADIUS=5.5
CELL=22.0
SYMBOL={1:'Ag',2:'C',3:'Si',4:'Ti'}
PAIRS={
 'ag_si_430045_250182':([430045,250182],['Ag','Si']),
 'ag_si_455657_339883':([455657,339883],['Ag','Si']),
 'ag_c_454089_336261':([454089,336261],['Ag','C']),
 'ag_ti_412351_186486':([412351,186486],['Ag','Ti']),
}

def main():
    lines=SOURCE.read_text().splitlines()
    h=next(i for i,s in enumerate(lines) if s.startswith('ITEM: ATOMS '))
    n=int(lines[3]); cols=lines[h].split()[2:]; ix={name:cols.index(name) for name in ('id','type','x','y','z','fx','fy','fz')}
    a=np.loadtxt(lines[h+1:h+1+n],ndmin=2)
    ids=a[:,ix['id']].astype(np.int64); typ=a[:,ix['type']].astype(np.int32)
    xyz=a[:,[ix['x'],ix['y'],ix['z']]]; f=a[:,[ix['fx'],ix['fy'],ix['fz']]]
    by_id={int(v):i for i,v in enumerate(ids)}
    if len(by_id)!=n or not np.isfinite(xyz).all() or not np.isfinite(f).all(): raise RuntimeError('Bad dump IDs or values')
    records=[]
    for name,(target_ids,elements) in PAIRS.items():
        if any(i not in by_id for i in target_ids): raise RuntimeError(f'{name}: target ID absent')
        ix_target=[by_id[i] for i in target_ids]
        actual=[SYMBOL[int(typ[i])] for i in ix_target]
        if actual!=elements: raise RuntimeError(f'{name}: expected {elements}, found {actual}')
        anchor=xyz[ix_target[0]].copy(); target_xyz=[]
        for i in ix_target:
            d=xyz[i]-anchor; d[:2]-=BOX*np.rint(d[:2]/BOX); target_xyz.append(anchor+d)
        target_xyz=np.asarray(target_xyz); midpoint=target_xyz.mean(axis=0)
        rel=xyz-midpoint; rel[:,:2]-=BOX*np.rint(rel[:,:2]/BOX)
        keep=np.linalg.norm(rel,axis=1)<=RADIUS
        kept_ids=ids[keep]
        if not np.isin(target_ids,kept_ids).all() or int(keep.sum())>96: raise RuntimeError(f'{name}: crop size/target mapping error {keep.sum()}')
        d=float(np.linalg.norm(target_xyz[0]-target_xyz[1])); mask=np.isin(kept_ids,target_ids).astype(np.int8)
        atoms=Atoms(symbols=[SYMBOL[int(t)] for t in typ[keep]],positions=rel[keep]+CELL/2,cell=[CELL]*3,pbc=True)
        atoms.new_array('source_atom_id',kept_ids.astype(np.int32)); atoms.new_array('SNAP_v12o_forces',f[keep]); atoms.new_array('is_DFT_force_target',mask)
        atoms.info.update({'config_type':f'SNAP_v12o_cold_10fs_{name}','source_timestep':10,'source_time_fs':10,'center_atom_ids':json.dumps(target_ids),'elements':json.dumps(elements),'center_pair_distance_A':d,'crop_radius_A':RADIUS,'crop_cell_A':CELL,'geometry_source':str(SOURCE),'source_pbc':'x/y periodic, z open; local crop embedded in cubic vacuum box','source_model':'SNAP v12o','DFT_target_rows':len(target_ids),'label_role':'v12p training candidate; not validation or holdout'})
        out=HERE/'candidates'/f'{name}_t10fs.extxyz'; write(out,atoms,format='extxyz')
        model_force={str(i):f[by_id[i]].tolist() for i in target_ids}
        records.append({'name':name,'target_ids':target_ids,'elements':elements,'pair_distance_A':d,'natoms':len(atoms),'target_model_forces_eV_A':model_force,'max_crop_model_force_eV_A':float(np.linalg.norm(f[keep],axis=1).max()),'output':str(out)})
    out={'status':'EXTRACTED','source_model':'SNAP v12o','source_time_fs':10,'radius_A':RADIUS,'cell_A':CELL,'pairs':records}
    (HERE/'extraction_manifest.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
if __name__=='__main__':main()
