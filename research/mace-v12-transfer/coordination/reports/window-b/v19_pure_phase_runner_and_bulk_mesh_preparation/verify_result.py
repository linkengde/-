"""Immutable archive inventory and complete pilot-label verification."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('label');ap.add_argument('result_dir',type=Path);args=ap.parse_args();from entry import load
manifest,r=load(args.label);expected=r['method'];folder=args.result_dir;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
if (folder/'SHA256SUMS.txt').exists():
 inventory={}
 for line in (folder/'SHA256SUMS.txt').read_text().splitlines():
  h,n=line.split('  ',1);p=folder/n;assert p.resolve().is_relative_to(folder.resolve()) and p.is_file() and sha(p)==h,'Archive inventory mismatch: '+n
  assert n not in inventory;inventory[n]=h
 actual={p.name for p in folder.iterdir() if p.is_file() and p.name!='SHA256SUMS.txt' and p.suffix!='.gpw'};assert actual==set(inventory),'Incomplete/extraneous archive inventory'
src=ROOT/r['input'];a=read(src);p=folder/(args.label+'_PW_PBE.extxyz');b=read(p);s=json.loads((folder/'summary.json').read_text());progress=json.loads((folder/'progress.json').read_text());log=(folder/'gpaw.log').read_text();assert (folder/'SHA256SUMS.txt').exists(), 'Complete immutable inventory required'
f=b.arrays['PW_PBE_forces'];expected=r['method']
checks={'input_hash':sha(src)==r['input_sha256']==s['source_sha256'],'label':s['label']==args.label==progress['label'],'count_formula':s['atoms']==len(a)==len(b) and s['formula']==a.get_chemical_formula(),'symbols':a.get_chemical_symbols()==b.get_chemical_symbols(),'geometry':np.array_equal(a.positions,b.positions) and np.array_equal(a.cell,b.cell) and np.array_equal(a.pbc,b.pbc),'ids':np.array_equal(a.arrays['local_new_id'],b.arrays['local_new_id']),'convergence':s['scf_converged'] is True and progress['scf_converged'] is True and progress['status']=='complete' and f"Converged in {s['scf_iterations']} steps" in log,'progress_consistency':progress['iteration']==s['scf_iterations'] and progress['source_sha256']==r['input_sha256'],'mpi4':s['mpi_ranks']==progress['mpi_ranks']==4,'energy_exact':b.info['PW_PBE_energy_eV']==s['energy_eV_cell'],'finite':f.shape==(len(a),3) and np.isfinite(f).all() and np.isfinite(s['energy_eV_cell']) and np.isfinite(s['free_energy_eV_cell']),'method':s['method']==expected,'pbc_metadata':s['method']['pbc']==b.pbc.tolist(),'poisson':s['method']['poissonsolver']==r['poissonsolver'],'energy_convention':s['energy_convention']==manifest['energy_convention'],'role':s['dataset_role']==r['role'],'force_summaries':np.isclose(s['fmax_eV_A'],np.linalg.norm(f,axis=1).max(),rtol=0,atol=1e-12) and np.isclose(s['rms_force_eV_A'],np.sqrt(np.mean(np.sum(f*f,axis=1))),rtol=0,atol=1e-12)}
checks={k:bool(v) for k,v in checks.items()}
if not all(checks.values()):raise SystemExit(json.dumps(checks))
print(json.dumps({'status':'PASS','checks':checks,'owner_task':r['owner_task'],'input_sha256':r['input_sha256'],'sha256':{p.name:sha(p) for p in folder.iterdir() if p.is_file() and p.name not in ['verification.json','SHA256SUMS.txt'] and p.suffix!='.gpw'},'large_gpw_archived':False},indent=2))
