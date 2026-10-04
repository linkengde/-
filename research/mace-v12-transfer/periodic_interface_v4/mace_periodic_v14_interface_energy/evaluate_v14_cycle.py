#!/usr/bin/env python3
"""Evaluate frozen test frames and strictly held-out interface geometries."""
from hashlib import sha256
import json
from pathlib import Path
import sys
import csv
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parent
if len(sys.argv)!=2 or sys.argv[1] not in ('v13','v14'):
    raise SystemExit('usage: python evaluate_v14_cycle.py v13|v14')
mode=sys.argv[1]
VERSION_DIR=PROJECT/('mace_periodic_v13_interface_energy' if mode=='v13' else 'mace_periodic_v14_interface_energy')
MODELS={n:PROJECT/f'mace_periodic_{n}_{suffix}'/path for n,suffix,path in [
 ('v9','forceonly','models/Ag_Ti_Si_C_MACE_periodic_v9_forceonly.model'),
 ('v10','energycal','checkpoints/MACE_periodic_v10_energycal_run-43.model'),
 ('v11','energyfocus','checkpoints/MACE_periodic_v11_energyfocus_run-44.model'),
 ('v12','interface_energy','checkpoints/MACE_periodic_v12_interface_energy_run-45.model'),
 ('v13','interface_energy','checkpoints/MACE_periodic_v13_interface_energy_run-45.model'),
]}
if mode=='v14': MODELS['v14']=PROJECT/'mace_periodic_v14_interface_energy/checkpoints/MACE_periodic_v14_interface_energy_run-45.model'
if any(not p.is_file() for p in MODELS.values()): raise SystemExit('A required checkpoint is missing.')

def source_frames():
    records=[]
    if mode=='v13':
        for folder,labels in [
          (PROJECT/'pbe_interface_v13_holdouts',['AgC_lateral_registry_holdout_v13','AgSi_lateral_registry_holdout_v13']),
          (PROJECT/'pbe_interface_v14_parallel_acquisition',['AgTi_registry_probe_v13_01','AgC_registry_strain_acq_v14_01','AgSi_registry_strain_acq_v14_01'])]:
            manifest=json.loads((folder/'input_manifest.json').read_text())
            for label in labels:
                meta=next(r for r in manifest['records'] if r['label']==label)
                result=folder/'calculations'/label
                verification=json.loads((result/'verification.json').read_text())
                if verification.get('status')!='PASS': raise SystemExit('Unverified v13 geometry: '+label)
                ext=next(result.glob('*_PW_PBE.extxyz'))
                if sha256(ext.read_bytes()).hexdigest()!=verification['sha256'][ext.name]: raise SystemExit('Hash mismatch: '+label)
                records.append((label,meta['role'],ext,verification['sha256'][ext.name]))
    else:
        for folder,label in [
          (PROJECT/'pbe_interface_v14_main_holdouts','AgC_registry_holdout_v14_01'),
          (PROJECT/'pbe_interface_v14_main_holdouts','AgSi_registry_holdout_v14_01'),
          (PROJECT/'pbe_interface_v14_parallel_acquisition','AgTi_registry_holdout_v14_01')]:
            result=folder/'calculations'/label
            verification=json.loads((result/'verification.json').read_text())
            if verification.get('status')!='PASS' or verification.get('role')!='v14_holdout': raise SystemExit('Not a reserved, verified v14 holdout: '+label)
            ext=next(result.glob('*_PW_PBE.extxyz'))
            if sha256(ext.read_bytes()).hexdigest()!=verification['sha256'][ext.name]: raise SystemExit('Hash mismatch: '+label)
            records.append((label,'v14_holdout',ext,verification['sha256'][ext.name]))
    return records

def reference(a):
    energy=next((a.info[k] for k in ('REF_energy','PW_PBE_energy_eV','energy') if k in a.info),None)
    forces=next((a.arrays[k] for k in ('REF_forces','PW_PBE_forces','forces') if k in a.arrays),None)
    if energy is None or forces is None:
        results=getattr(a.calc,'results',{}) if a.calc is not None else {}
        energy=results.get('energy',energy);forces=results.get('forces',forces)
    forces=np.asarray(forces,dtype=float)
    if not np.isfinite(float(energy)) or forces.shape!=(len(a),3) or not np.isfinite(forces).all(): raise SystemExit('Invalid reference labels.')
    return float(energy),forces

torch.set_default_dtype(torch.float64)
calculators={}
for name,path in MODELS.items():
    model=torch.load(path,map_location='cpu',weights_only=False).double().eval()
    calculators[name]=MACECalculator(models=[model],device='cpu',default_dtype='float64')

rows=[]
for a in read(VERSION_DIR/'data/test.extxyz',index=':'):
    rows.append((a.info.get('config_type',''), 'frozen_test', a, *reference(a)))
holdouts=source_frames()
for label,role,path,_hash in holdouts:
    a=read(path);rows.append((label,role,a,*reference(a)))
output=[]
for label,role,source,eref,fref in rows:
    a=source.copy();symbols=a.get_chemical_symbols()
    central=np.flatnonzero(np.asarray(a.arrays.get('central_pair',np.zeros(len(a))),dtype=bool))
    radial_ref=radial_error=None
    if len(central)==2:
        i,j=map(int,central);v=a.positions[j]-a.positions[i];d=float(np.linalg.norm(v))
        radial_ref=float(np.dot(fref[j]-fref[i],v/d))
    row={'label':label,'role':role,'atoms':len(a),'formula':a.get_chemical_formula(),'energy_DFT_eV_cell':eref,
         'geometry_sha256':sha256(Path(source.info['source_path']).read_bytes()).hexdigest() if 'source_path' in source.info else None}
    for name,calc in calculators.items():
        a.calc=calc;ep=float(a.get_potential_energy());fp=np.asarray(a.get_forces());delta=fp-fref
        row[name+'_energy_abs_error_meV_atom']=abs(eref-ep)*1000/len(a)
        row[name+'_force_vector_RMSE_eV_A']=float(np.sqrt(np.mean(np.sum(delta**2,axis=1))))
        row[name+'_force_vector_max_error_eV_A']=float(np.linalg.norm(delta,axis=1).max())
        if len(central)==2:
            radial=float(np.dot(fp[j]-fp[i],v/d));row[name+'_separating_force_eV_A']=radial
            row[name+'_separating_force_abs_error_eV_A']=abs(radial-radial_ref)
    output.append(row)

threshold={'energy_MAE_meV_atom':10.0,'force_vector_RMSE_eV_A':0.05,'separating_force_abs_error_eV_A':0.10}
types={}
for kind in ('AgC','AgSi','AgTi'):
    subset=[r for r in output if r['role']!='frozen_test' and r['label'].startswith(kind)]
    if not subset: continue
    for name in MODELS:
        energy=float(np.mean([r[name+'_energy_abs_error_meV_atom'] for r in subset]))
        fr=np.sqrt(sum(r['atoms']*r[name+'_force_vector_RMSE_eV_A']**2 for r in subset)/sum(r['atoms'] for r in subset))
        radial=max(r[name+'_separating_force_abs_error_eV_A'] for r in subset if name+'_separating_force_abs_error_eV_A' in r)
        types.setdefault(kind,{})[name]={'structures':len(subset),'energy_MAE_meV_atom':energy,'force_vector_RMSE_eV_A':float(fr),
                                       'separating_force_max_abs_error_eV_A':radial,
                                       'screen_pass':bool(energy<=threshold['energy_MAE_meV_atom'] and fr<=threshold['force_vector_RMSE_eV_A'] and radial<=threshold['separating_force_abs_error_eV_A'])}

out_dir=VERSION_DIR/'results';out_dir.mkdir(exist_ok=True)
basename='v13_independent_holdout_comparison' if mode=='v13' else 'v14_independent_holdout_comparison'
summary={'mode':mode,'models':{n:{'path':str(p.relative_to(PROJECT)),'sha256':sha256(p.read_bytes()).hexdigest()} for n,p in MODELS.items()},
         'thresholds':threshold,'per_structure':output,'by_interface':types,
         'independent_geometry_holdouts':{r[0]:{'role':r[1],'path':str(r[2].relative_to(PROJECT)),'sha256':r[3]} for r in holdouts},
         'scope':'Single geometries from existing cluster motifs are narrow checks, not broad transfer validation for an extended periodic interface or thermal/liquid configurations.',
         'overall_status':'SCREEN_FAIL' if any(not info[name]['screen_pass'] for info in types.values() for name in [mode]) else 'UNDETERMINED'}
(out_dir/(basename+'.json')).write_text(json.dumps(summary,indent=2)+'\n')
with (out_dir/(basename+'.csv')).open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in output))),lineterminator='\n');writer.writeheader();writer.writerows(output)
print(json.dumps(summary,indent=2))
