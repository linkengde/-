#!/usr/bin/env python3
"""Score five fresh acquisition labels with frozen v9-v13 models before v14 use."""
from hashlib import sha256
import json
from pathlib import Path
import csv
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parent
MODELS={
 'v9':PROJECT/'mace_periodic_v9_forceonly/models/Ag_Ti_Si_C_MACE_periodic_v9_forceonly.model',
 'v10':PROJECT/'mace_periodic_v10_energycal/checkpoints/MACE_periodic_v10_energycal_run-43.model',
 'v11':PROJECT/'mace_periodic_v11_energyfocus/checkpoints/MACE_periodic_v11_energyfocus_run-44.model',
 'v12':PROJECT/'mace_periodic_v12_interface_energy/checkpoints/MACE_periodic_v12_interface_energy_run-45.model',
 'v13':PROJECT/'mace_periodic_v13_interface_energy/checkpoints/MACE_periodic_v13_interface_energy_run-45.model',
}
FROZEN=PROJECT/'mace_periodic_v13_interface_energy/data/test.extxyz'
SOURCES=[
 (PROJECT/'pbe_interface_v13_holdouts','AgC_lateral_registry_holdout_v13'),
 (PROJECT/'pbe_interface_v13_holdouts','AgSi_lateral_registry_holdout_v13'),
 (PROJECT/'pbe_interface_v14_parallel_acquisition','AgTi_registry_probe_v13_01'),
 (PROJECT/'pbe_interface_v14_parallel_acquisition','AgC_registry_strain_acq_v14_01'),
 (PROJECT/'pbe_interface_v14_parallel_acquisition','AgSi_registry_strain_acq_v14_01'),
]
if any(not p.is_file() for p in list(MODELS.values())+[FROZEN]): raise SystemExit('Missing model or frozen test data.')
torch.set_default_dtype(torch.float64)
calculators={}
for name,path in MODELS.items():
    model=torch.load(path,map_location='cpu',weights_only=False).double().eval()
    calculators[name]=MACECalculator(models=[model],device='cpu',default_dtype='float64')

def ref(atoms):
    energy=atoms.info.get('PW_PBE_energy_eV',atoms.info.get('REF_energy'))
    force=atoms.arrays.get('PW_PBE_forces',atoms.arrays.get('REF_forces'))
    force=np.asarray(force,dtype=float)
    if energy is None or force.shape!=(len(atoms),3) or not np.isfinite(float(energy)) or not np.isfinite(force).all():
        raise SystemExit('Missing or invalid DFT labels.')
    return float(energy),force

sources=[]
for folder,label in SOURCES:
    manifest=json.loads((folder/'input_manifest.json').read_text())
    meta=next((r for r in manifest['records'] if r['label']==label),None)
    if meta is None: raise SystemExit('Assigned input record missing: '+label)
    result=folder/'calculations'/label
    verification=json.loads((result/'verification.json').read_text())
    if verification.get('status')!='PASS': raise SystemExit('Archive has no PASS verification: '+label)
    ext=next(result.glob('*_PW_PBE.extxyz'))
    if sha256(ext.read_bytes()).hexdigest()!=verification['sha256'][ext.name]: raise SystemExit('Reference hash mismatch: '+label)
    sources.append((label,meta.get('role','v13_independent'),ext,verification['sha256'][ext.name]))

frames=[('frozen_test',a.info.get('config_type',''),'frozen_test',a,*ref(a)) for a in read(FROZEN,index=':')]
for label,role,path,hashval in sources:
    a=read(path); e,f=ref(a);frames.append(('acquisition',label,role,a,e,f))
rows=[]
for split,label,role,source,eref,fref in frames:
    a=source.copy(); n=len(a); central=np.flatnonzero(np.asarray(a.arrays.get('central_pair',np.zeros(n)),dtype=bool))
    radial_ref=None
    if len(central)==2:
        i,j=map(int,central);vec=a.positions[j]-a.positions[i];dist=float(np.linalg.norm(vec))
        radial_ref=float(np.dot(fref[j]-fref[i],vec/dist))
    row={'split':split,'label':label,'role':role,'atoms':n,'formula':a.get_chemical_formula()}
    if radial_ref is not None: row['DFT_separating_force_eV_A']=radial_ref
    for name,calc in calculators.items():
        a.calc=calc;ep=float(a.get_potential_energy());fp=np.asarray(a.get_forces());delta=fp-fref
        row[name+'_energy_abs_error_meV_atom']=abs(eref-ep)*1000/n
        row[name+'_force_vector_RMSE_eV_A']=float(np.sqrt(np.mean(np.sum(delta**2,axis=1))))
        row[name+'_force_vector_max_error_eV_A']=float(np.linalg.norm(delta,axis=1).max())
        if len(central)==2:
            row[name+'_separating_force_eV_A']=float(np.dot(fp[j]-fp[i],vec/dist))
            row[name+'_separating_force_abs_error_eV_A']=abs(row[name+'_separating_force_eV_A']-radial_ref)
    rows.append(row)

threshold={'energy_abs_error_meV_atom':10.0,'force_vector_RMSE_eV_A':0.05,'separating_force_abs_error_eV_A':0.10}
by_interface={}
for kind in ('AgC','AgSi','AgTi'):
    subset=[r for r in rows if r['split']=='acquisition' and r['label'].startswith(kind)]
    if not subset: raise SystemExit('No acquisition geometry for '+kind)
    by_interface[kind]={}
    for model in MODELS:
        e=max(r[model+'_energy_abs_error_meV_atom'] for r in subset)
        f=np.sqrt(sum(r['atoms']*r[model+'_force_vector_RMSE_eV_A']**2 for r in subset)/sum(r['atoms'] for r in subset))
        radial=max(r.get(model+'_separating_force_abs_error_eV_A',0) for r in subset)
        by_interface[kind][model]={'structures':len(subset),'max_energy_abs_error_meV_atom':e,
             'force_vector_RMSE_eV_A':float(f),'max_separating_force_abs_error_eV_A':radial,
             'screen_pass':bool(e<=threshold['energy_abs_error_meV_atom'] and f<=threshold['force_vector_RMSE_eV_A'] and radial<=threshold['separating_force_abs_error_eV_A'])}
out={'models':{n:{'path':str(p.relative_to(PROJECT)),'sha256':sha256(p.read_bytes()).hexdigest()} for n,p in MODELS.items()},
     'frozen_test_sha256':sha256(FROZEN.read_bytes()).hexdigest(),'thresholds':threshold,
     'per_structure':rows,'by_interface':by_interface,
     'holdouts':{label:{'role':role,'path':str(path.relative_to(PROJECT)),'sha256':h} for label,role,path,h in sources},
     'v13_screen':'FAIL' if any(not by_interface[k]['v13']['screen_pass'] for k in by_interface) else 'LOCAL_SCREEN_PASS',
     'overall_status':'UNDETERMINED' if all(by_interface[k]['v13']['screen_pass'] for k in by_interface) else 'FAIL_SCREENING',
     'scope':'Five new registry/local-perturbation labels are related to existing cluster motifs. These fixed geometries do not establish transfer to an extended interface, thermal disorder, liquid Ag or pressure.'}
outdir=PROJECT/'mace_periodic_v13_interface_energy/results';outdir.mkdir(exist_ok=True)
(outdir/'v13_independent_holdout_comparison.json').write_text(json.dumps(out,indent=2)+'\n')
with (outdir/'v13_independent_holdout_comparison.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in rows))),lineterminator='\n');writer.writeheader();writer.writerows(rows)
print(json.dumps(out,indent=2))
