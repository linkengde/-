#!/usr/bin/env python3
"""Assemble the pure Ti-Si-C train/dev files after every registered DFT label passes."""
from __future__ import annotations
import hashlib, io, json, subprocess
from pathlib import Path
import numpy as np
from ase.io import read, write

ROOT=Path(__file__).resolve().parent
REPO=Path('/workspace/jeitschko-main-run')
MANIFEST=json.loads((ROOT/'input_manifest.json').read_text())
SOURCE_COMMIT=MANIFEST['source']['source_commit']

def blob(path):
    return subprocess.check_output(['git','-C',str(REPO),'show',f'{SOURCE_COMMIT}:{path}'])
def sha(data): return hashlib.sha256(data).hexdigest()
def normalize(atoms,label,split,source_record):
    a=atoms.copy(); a.calc=None
    free=float(a.info.get('PW_PBE_free_energy_eV',a.info.get('REF_energy',float('nan'))))
    forces=np.asarray(a.arrays.get('PW_PBE_forces',a.arrays.get('REF_forces')),dtype=float)
    if set(a.get_chemical_symbols())!={'Ti','Si','C'} or len(a)!=12:
        raise ValueError(f'{label}: wrong species or atom count')
    if not np.isfinite(free) or forces.shape!=(len(a),3) or not np.isfinite(forces).all():
        raise ValueError(f'{label}: missing/nonfinite free-energy or force targets')
    a.info.update({'config_type':label,'REF_energy':free,'dataset_split':split,
        'dataset_role':'skeleton_train_candidate_pending_audit' if split=='train' else 'skeleton_development_candidate',
        'energy_target_convention':'force-consistent finite-width free energy at Fermi-Dirac sigma 0.10 eV',
        'source_record_sha256':source_record['source_result_sha256']})
    a.arrays['REF_forces']=forces.copy()
    return a

status_path=ROOT/'sequence_status.json'
status=json.loads(status_path.read_text())
if status.get('status')!='all_skeleton_archives_verified':
    raise SystemExit(f"refusing dataset assembly; sequence status is {status.get('status')}")
verified={r['label'] for r in status.get('completed_labels',[]) if r.get('archive_verification')=='PASS'}
expected={r['label'] for r in MANIFEST['preflight_jobs']}
if verified!=expected: raise SystemExit(f'missing verified labels: {sorted(expected-verified)}')

source_labels=[]
for row in MANIFEST['previously_verified_training_candidates']:
    raw=blob(row['source_result_path'])
    a=read(io.StringIO(raw.decode()),format='extxyz')
    source_labels.append((row['label'],row['split'],normalize(a,row['label'],row['split'],{
        'source_result_sha256':sha(raw)}),sha(raw),row))

new_labels=[]
for row in MANIFEST['preflight_jobs']:
    folder=ROOT/'calculations'/row['label']
    file=folder/f"{row['label']}_PW_PBE.extxyz"
    raw=file.read_bytes()
    a=read(io.StringIO(raw.decode()),format='extxyz')
    new_labels.append((row['label'],row['split'],normalize(a,row['label'],row['split'],{
        'source_result_sha256':sha(raw)}),sha(raw),row))

train=[row[2] for row in source_labels+new_labels if row[1]=='train']
dev=[row[2] for row in source_labels+new_labels if row[1]=='development']
if len(train)!=19 or len(dev)!=4 or not train or not dev:
    raise SystemExit(f'expected 19 train/4 dev; received {len(train)}/{len(dev)}')

def geometry_key(a):
    return (tuple(a.get_chemical_symbols()),tuple(np.round(a.cell.array.ravel(),8)),
            tuple(np.round(a.positions.ravel(),8)),tuple(a.pbc.tolist()))
keys=[geometry_key(a) for a in train]
if len(set(keys))!=len(keys): raise SystemExit('exact duplicate geometries in train split')
dev_keys=[geometry_key(a) for a in dev]
if len(set(dev_keys))!=len(dev_keys): raise SystemExit('exact duplicate geometries in dev split')
if set(keys)&set(dev_keys): raise SystemExit('train/dev exact geometry overlap')

dataset=ROOT/'dataset'; dataset.mkdir(exist_ok=True)
write(dataset/'train.extxyz',train,format='extxyz')
write(dataset/'valid.extxyz',dev,format='extxyz')
source_rows=[]
for label,split,a,result_sha,row in source_labels+new_labels:
    source_rows.append({'label':label,'split':split,'source_result_sha256':result_sha,
        'input_sha256':row.get('input_sha256'),'reused_source_path':row.get('source_result_path'),
        'scf_iterations':row.get('scf_iterations')})
dataset_manifest={
  'schema':'tisic_skeleton_dataset_v25_v1','status':'NON_OFFICIAL_SKELETON_CANDIDATE_PENDING_AUDIT',
  'train_frames':len(train),'development_frames':len(dev),'test_frames':0,
  'species':['Ti','Si','C'],'train_sha256':sha((dataset/'train.extxyz').read_bytes()),
  'development_sha256':sha((dataset/'valid.extxyz').read_bytes()),
  'energy_target':'finite-width free energy at Fermi-Dirac sigma 0.10 eV, force-consistent',
  'force_target':'GPAW/ASE forces at Fermi-Dirac sigma 0.10 eV',
  'source_family':'Jeitschko_Nowotny_1967_single_crystal','source_cif_sha256':MANIFEST['source']['source_cif_sha256'],
  'acceptance_standard':MANIFEST['acceptance_standard_fixed_before_fit'],
  'source_rows':source_rows,
  'limitations':MANIFEST['limitations'],
  'official_data_changed':False,'sealed_test_opened':False,
}
(dataset/'manifest.json').write_text(json.dumps(dataset_manifest,indent=2,ensure_ascii=False)+'\n')
files=[dataset/'train.extxyz',dataset/'valid.extxyz',dataset/'manifest.json']
(dataset/'SHA256SUMS.txt').write_text(''.join(f'{sha(f.read_bytes())}  {f.name}\n' for f in files))
print(json.dumps({'status':dataset_manifest['status'],'train_frames':len(train),'development_frames':len(dev),
                  'train_sha256':dataset_manifest['train_sha256'],'development_sha256':dataset_manifest['development_sha256']},indent=2))
