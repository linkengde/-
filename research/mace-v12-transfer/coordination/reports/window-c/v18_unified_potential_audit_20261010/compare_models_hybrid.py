#!/usr/bin/env python3
"""Score two frozen MACE models on the exact six public legacy-hybrid diagnostic DFT structures."""
import csv, hashlib, importlib.util, json, sys, os, tempfile
from pathlib import Path
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT=Path(os.environ.get('CLOUD_C_REPO_ROOT',Path(__file__).resolve().parents[6])).resolve()
BASE=ROOT/'research/mace-v12-transfer'
OUT=Path(os.environ.get('CLOUD_C_AUDIT_OUTPUT',Path(tempfile.gettempdir())/'cloud-c-v18-audit')).resolve(); OUT.mkdir(parents=True,exist_ok=True)
torch.set_num_threads(4); torch.set_num_interop_threads(1)
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod
helper=load('hybrid_helpers',BASE/'coordination/reports/window-b/v19_existing_AgSi_candidate_force_error_screen_install_and_infer_window_b/infer_and_evaluate.py')
MODEL=ROOT/helper.MODEL_PATH
FOUNDATION=BASE/'periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model'
ROOT_CALC=BASE/'periodic_interface_v4/pbe_interface_v19_residual_acquisition/calculations'
labels=['AgC_residual_mode_m_06_v19_acq','AgC_residual_mode_p_06_v19_acq','AgSi_residual_mode_m_03_v19_acq','AgSi_residual_mode_p_03_v19_acq','AgTi_residual_mode_m_03_v19_acq','AgTi_residual_mode_p_03_v19_acq']
old={r['label']:r for r in csv.DictReader((BASE/'coordination/reports/window-a/candidate02_cross_potential_dft_diagnosis/metrics.csv').open())}
assert len(old)==16
records=[]
for label in labels:
    folder=ROOT_CALC/label
    summary=json.loads((folder/'summary.json').read_text()); verification=json.loads((folder/'verification.json').read_text())
    assert summary.get('label')==label and summary.get('scf_converged') is True
    assert verification.get('status')=='PASS' and all(verification.get('checks',{}).values())
    checks=[]
    for line in (folder/'SHA256SUMS.txt').read_text().splitlines():
        if not line.strip(): continue
        expected,name=line.split(maxsplit=1); name=name.lstrip('*'); member=(folder/name).resolve()
        assert folder.resolve() in member.parents and member.is_file()
        actual=hashlib.sha256(member.read_bytes()).hexdigest(); assert actual==expected
        checks.append({'file':name,'sha256':actual})
    assert checks
    out=folder/f'{label}_PW_PBE.extxyz'; assert verification.get('sha256',{}).get(out.name)==hashlib.sha256(out.read_bytes()).hexdigest()
    atoms=read(out,index=0)
    src_text=str(summary['source']).replace('\\','/'); ix=src_text.find('research/mace-v12-transfer/'); assert ix>=0
    source=(ROOT/src_text[ix:]).resolve(); assert ROOT.resolve() in source.parents
    assert hashlib.sha256(source.read_bytes()).hexdigest()==summary['source_sha256']
    source_atoms=read(source,index=0)
    for a,b in [('numbers',atoms.numbers),('positions',atoms.positions),('cell',atoms.cell.array),('pbc',atoms.pbc)]:
        x=b; y={'numbers':source_atoms.numbers,'positions':source_atoms.positions,'cell':source_atoms.cell.array,'pbc':source_atoms.pbc}[a]
        assert np.array_equal(x,y),(label,a)
    assert np.array_equal(helper.atom_ids(atoms),helper.atom_ids(source_atoms))
    force=np.asarray(atoms.arrays['PW_PBE_forces'],float); assert force.shape==(len(atoms),3) and np.isfinite(force).all()
    records.append((label,folder,atoms,force,summary,checks))
models={}
for name,path in [('foundation',FOUNDATION),('V18',MODEL)]:
    models[name]=MACECalculator(model_paths=str(path),device='cpu',default_dtype='float64')
rows=[]
for label,folder,atoms,reference,summary,checks in records:
    for name,calc in models.items():
        a=atoms.copy(); a.calc=calc; pred=np.asarray(a.get_forces(),float); assert np.isfinite(pred).all()
        m=helper.force_metrics(pred,reference,a)
        row={'label':label,'contact':old[label]['contact'],'role':old[label]['role'],'atoms':len(a),'DFT_source_sha256':summary['source_sha256'],'DFT_result_sha256':hashlib.sha256((folder/f'{label}_PW_PBE.extxyz').read_bytes()).hexdigest(),'model':name,'model_force_vector_RMSE_eV_A':m['force_vector_RMSE_eV_A'],'model_marked_pair_signed_error_eV_A':m['marked_pair_separating_force_eV_A']['signed_error_model_minus_dft'],'model_marked_pair_abs_error_eV_A':m['marked_pair_separating_force_eV_A']['absolute_error'],'old_candidate02_hybrid_force_vector_RMSE_eV_A':float(old[label]['force_vector_rmse_eV_A']),'old_candidate02_hybrid_marked_pair_signed_error_eV_A':float(old[label]['separation_error_eV_A']),'old_candidate02_hybrid_abs_pair_error_eV_A':abs(float(old[label]['separation_error_eV_A'])),'old_hybrid_DFT_energies_eV_cell_unscored':json.dumps([float(old[label]['model_energy_eV']),float(old[label]['native_DFT_eV']),float(old[label]['free_DFT_eV'])]),'DFT_converged_hashes_checked':len(checks),'independent_validation':False}
        for element,values in m['per_species'].items(): row[f'{element}_force_vector_RMSE_eV_A']=values['force_vector_RMSE_eV_A']
        rows.append(row)
with (OUT/'same_geometry_old_potential_comparison.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={'scope':'Frozen CPU float64 MACE inference and comparison with already published, unchanged candidate02 EAM+Tersoff+Morse run0 metrics on the same six nonsealed, converged DFT structures. No LAMMPS/DFT/MD/training launched.','candidate02_identity':'This is the published candidate02 heating calibration-v0 package. Published records explicitly do not establish that its exact files/configuration were Stage69; therefore this comparison cannot be attributed to an authenticated Stage69 run.','structures':len(labels),'DFT_archives_verified':True,'metrics_csv_hash':hashlib.sha256((BASE/'coordination/reports/window-a/candidate02_cross_potential_dft_diagnosis/metrics.csv').read_bytes()).hexdigest(),'potential_file_hashes':json.loads((BASE/'coordination/reports/window-b/current_heating_potential_handoff_20261009/manifest.json').read_text())['files'],'models':{name:{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()} for name,path in [('V18',MODEL),('foundation',FOUNDATION)]},'records':len(rows),'output':'same_geometry_old_potential_comparison.csv','constraints':['Six local points reuse earlier training/development parent families and are not independent validation.','Old hybrid metrics were precomputed and verified in candidate02 report; this script did not run LAMMPS or re-evaluate hybrid.','Different energy zeros mean absolute energies are not scored; same-parent relative response and force errors are separate.']}
(OUT/'same_geometry_old_potential_comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Compared',len(labels),'verified structures with two frozen MACE models against the existing candidate02 hybrid results.')
