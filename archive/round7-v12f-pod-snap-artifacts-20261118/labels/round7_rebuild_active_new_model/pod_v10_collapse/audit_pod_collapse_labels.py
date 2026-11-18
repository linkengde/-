#!/workspace/venvs/dft-labels/bin/python
"""Audit the two PW350/Gamma labels for POD-v10 FIRE collapse environments."""
from __future__ import annotations
import hashlib
import json
import re
from pathlib import Path
import numpy as np
from ase.io import read

HERE=Path(__file__).resolve().parent
LABELS=HERE/'labels'
SOURCE=HERE/'candidate_geometries.extxyz'
CASES=(0,2)

def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def audit_one(idx:int):
    stem=f'pod_collapse_local_{idx:02d}_pw350_gamma'
    out=LABELS/f'{stem}.extxyz'; summary=LABELS/f'{stem}.json'; log=LABELS/f'{stem}.txt'
    missing=[str(p) for p in (out,summary,log) if not p.is_file()]
    if missing:return {'frame_index':idx,'status':'pending','missing':missing}
    src=read(SOURCE,index=idx); atoms=read(out); record=json.loads(summary.read_text())
    issues=[]
    if not np.array_equal(src.numbers,atoms.numbers):issues.append('element order/composition changed')
    if not np.array_equal(src.positions,atoms.positions):issues.append('coordinates changed')
    if not np.array_equal(src.cell.array,atoms.cell.array):issues.append('cell changed')
    if not np.array_equal(src.pbc,atoms.pbc):issues.append('PBC changed')
    forces=np.asarray(atoms.arrays.get('REF_forces',[]),dtype=float)
    if forces.shape!=(len(atoms),3) or not np.isfinite(forces).all():issues.append('DFT force array invalid')
    ids=np.asarray(atoms.arrays.get('source_atom_id',[]),dtype=int)
    center=int(atoms.info.get('center_atom_id',-1)); partner=int(atoms.info.get('partner_atom_id',-1))
    for label,atom_id in (('center',center),('partner',partner)):
        if np.count_nonzero(ids==atom_id)!=1:issues.append(f'{label} ID not uniquely mapped')
    lines=[line for line in log.read_text(errors='replace').splitlines() if '|iter:' in line]
    markers=re.findall(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?c(?=\s*[|])',lines[-1]) if lines else []
    scf_converged=len(markers)>=3
    if not scf_converged:issues.append('last GPAW SCF row does not show three convergence markers')
    if int(record.get('natoms',-1))!=len(atoms):issues.append('JSON atom count disagrees')
    if int(record.get('center_atom_id',-1))!=center or int(record.get('partner_atom_id',-1))!=partner:issues.append('JSON center/partner IDs disagree')
    if abs(float(record.get('energy_eV',np.nan))-float(atoms.info.get('REF_energy',np.nan)))>1e-10:issues.append('JSON energy disagrees')
    cim=np.flatnonzero(ids==center); pim=np.flatnonzero(ids==partner)
    projection=None
    if len(cim)==1 and len(pim)==1 and np.isfinite(forces).all():
        ci=int(cim[0]); pi=int(pim[0])
        vec=atoms.positions[pi]-atoms.positions[ci]
        vec-=atoms.cell.lengths()*np.rint(vec/atoms.cell.lengths())
        unit=vec/np.linalg.norm(vec)
        projection={'center_toward_partner_eV_A':float(forces[ci]@unit),
            'partner_toward_center_eV_A':float(forces[pi]@(-unit))}
    return {'frame_index':idx,'status':'PASS' if not issues else 'FAIL','issues':issues,
        'center_atom_id':center,'partner_atom_id':partner,'pair_distance_A':float(atoms.info['center_partner_distance_A']),
        'natoms':len(atoms),'SCF_iterations':len(lines),'SCF_converged':scf_converged,
        'center_force_norm_eV_A':float(np.linalg.norm(forces[np.flatnonzero(ids==center)[0]])),
        'partner_force_norm_eV_A':float(np.linalg.norm(forces[np.flatnonzero(ids==partner)[0]])),
        'DFT_pair_force_projection_eV_A':projection,
        'max_force_norm_eV_A':float(np.linalg.norm(forces,axis=1).max()),
        'sha256':{'extxyz':sha256(out),'json':sha256(summary),'log':sha256(log)}}

def main():
    results=[audit_one(i) for i in CASES]
    complete=all(r['status']=='PASS' for r in results)
    report={'status':'PASS' if complete else 'INCOMPLETE',
        'method':'GPAW PBE/PW350 Gamma; Fermi-Dirac sigma 0.10 eV','expected_labels':len(CASES),
        'passed_labels':sum(r['status']=='PASS' for r in results),'frames':results,
        'role':'v11 training candidate only; finite-crop energy has zero fit weight; not blind test'}
    out=HERE/'pod_collapse_label_audit.json'
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(report,indent=2,ensure_ascii=False))
    return 0 if complete else 2

if __name__=='__main__':
    raise SystemExit(main())
