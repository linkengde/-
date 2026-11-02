#!/usr/bin/env python3
"""Verify the assembled Ti-Si-C candidate and pinned MACE foundation before fitting."""
import hashlib,json
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parent
DATA=ROOT/'dataset'
FOUNDATION=Path('/workspace/-/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model')
EXPECTED_FOUNDATION_SHA256='2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
if not (DATA/'manifest.json').is_file(): raise SystemExit('dataset not assembled')
m=json.loads((DATA/'manifest.json').read_text())
if m.get('status')!='NON_OFFICIAL_SKELETON_CANDIDATE_PENDING_AUDIT': raise SystemExit('unexpected dataset status')
if m.get('acceptance_standard',{}).get('per_configuration_force_vector_RMSE_max_eV_A')!=0.07: raise SystemExit('skeleton 0.07 eV/A acceptance standard mismatch')
if sha(FOUNDATION)!=EXPECTED_FOUNDATION_SHA256: raise SystemExit('MACE foundation hash mismatch')
if m.get('species')!=['Ti','Si','C'] or m.get('train_frames')!=19 or m.get('development_frames')!=4: raise SystemExit('unexpected dataset scope/counts')
rows={line.split('  ',1)[1]:line.split('  ',1)[0] for line in (DATA/'SHA256SUMS.txt').read_text().splitlines()}
for name in ('train.extxyz','valid.extxyz','manifest.json'):
 if rows.get(name)!=sha(DATA/name): raise SystemExit(f'dataset checksum mismatch: {name}')
train=read(DATA/'train.extxyz',index=':'); dev=read(DATA/'valid.extxyz',index=':')
for a in train+dev:
 if set(a.get_chemical_symbols())!={'Ti','Si','C'} or not np.isfinite(float(a.info['REF_energy'])) or a.arrays['REF_forces'].shape!=(len(a),3) or not np.isfinite(a.arrays['REF_forces']).all(): raise SystemExit('invalid skeleton energy/force row')
print(json.dumps({'status':'PASS','train_frames':len(train),'development_frames':len(dev),'foundation_model_sha256':sha(FOUNDATION)},indent=2))
