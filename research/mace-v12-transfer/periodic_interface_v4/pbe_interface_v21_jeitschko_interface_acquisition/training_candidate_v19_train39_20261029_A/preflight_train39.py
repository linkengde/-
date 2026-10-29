#!/usr/bin/env python3
"""Read-only preflight for the user-authorized non-official V19 train39 experiment."""
from __future__ import annotations
import hashlib, importlib.metadata as metadata, json, shutil, subprocess, sys
from pathlib import Path
import numpy as np
import torch
from ase.io import read

HERE = Path(__file__).resolve().parent
TASK = HERE.parent
REPO = Path('/workspace/jeitschko-main-run')
DATA = HERE / 'dataset'
FOUNDATION = Path('/workspace/-/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model')
EXPECTED_FOUNDATION = '2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'
LABELS = [
    'AgSi_JN1967_registryB_k6x6x1_sigma0p10_preflight_v21',
    'AgC_JN1967_registryB_k6x6x1_sigma0p10_acquisition_v21',
    'AgTi_JN1967_registryB_k6x6x1_sigma0p10_acquisition_v21',
]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
manifest = json.loads((DATA / 'manifest.json').read_text())
if manifest['status'] != 'NON_OFFICIAL_V19_TRAIN39_EXPERIMENT_CANDIDATE':
    raise SystemExit('Unexpected dataset status')
if manifest['roles'] != {'train_candidate':39,'development_validation_frozen':3,'test':0}:
    raise SystemExit('Unexpected dataset roles or frame counts')
if manifest['official_train36_and_dev3_mutated'] or manifest['withheld_labels_opened']:
    raise SystemExit('Official inputs changed or withheld labels opened')
if sha(DATA/'train.extxyz') != manifest['train_sha256'] or sha(DATA/'valid.extxyz') != manifest['dev_sha256']:
    raise SystemExit('Dataset file digest mismatch')
if sha(FOUNDATION) != EXPECTED_FOUNDATION or EXPECTED_FOUNDATION != manifest['foundation_model_sha256']:
    raise SystemExit('Foundation model digest mismatch')
train = read(DATA/'train.extxyz', index=':')
dev = read(DATA/'valid.extxyz', index=':')
if len(train) != 39 or len(dev) != 3:
    raise SystemExit(f'Frame count mismatch: train={len(train)} dev={len(dev)}')
for atoms in train + dev:
    e = float(atoms.info['REF_energy'])
    f = np.asarray(atoms.arrays['REF_forces'], dtype=float)
    if not np.isfinite(e) or f.shape != (len(atoms),3) or not np.isfinite(f).all():
        raise SystemExit('Nonfinite or malformed energy/force labels')
for label, row, atoms in zip(LABELS, manifest['candidate_rows'], train[-3:]):
    folder = TASK/'calculations'/label
    result = subprocess.check_output([
        '/workspace/.venvs/gpaw-mpi/bin/python', str(TASK/'verify_preflight.py'), label
    ], text=True)
    verification = json.loads(result)
    if verification['status'] != 'PASS' or row['config_type'] != label:
        raise SystemExit(f'Archived result not verified: {label}')
    if row['source_family'] != 'Jeitschko_Nowotny_1967_single_crystal':
        raise SystemExit(f'Unexpected source family: {label}')
    summary = json.loads((folder/'summary.json').read_text())
    if atoms.info['REF_energy'] != summary['native_extrapolated_energy_eV_cell']:
        raise SystemExit(f'Energy target mismatch: {label}')
    if not np.array_equal(atoms.arrays['REF_forces'], atoms.arrays['PW_PBE_forces']):
        raise SystemExit(f'Force target mismatch: {label}')
    if not np.isclose(float(atoms.info['free_energy_eV']), summary['free_energy_eV_cell'], rtol=0, atol=1e-12):
        raise SystemExit(f'Free-energy diagnostic mismatch: {label}')
if (DATA/'test.extxyz').exists():
    raise SystemExit('Unexpected test data file in experiment dataset')
for name, version in [('torch','2.5.1+cpu'),('mace-torch','0.3.16'),('ase','3.29.0')]:
    if metadata.version(name) != version:
        raise SystemExit(f'Runtime version mismatch for {name}: {metadata.version(name)}')
if torch.cuda.is_available():
    raise SystemExit('Expected CPU-only training environment')
free = shutil.disk_usage(HERE).free
if free < 2_950_000_000:
    raise SystemExit(f'Insufficient free disk for retained epoch checkpoint verification: {free}')
run = HERE/'run'
if run.exists() and any(run.iterdir()):
    raise SystemExit('Refusing to overwrite a nonempty run directory')
if subprocess.run(['pgrep','-af','[m]ace_run_train'],capture_output=True,text=True).stdout.strip():
    raise SystemExit('Another MACE training process appears active')
print(json.dumps({
    'status':'V19_TRAIN39_PREFLIGHT_PASS_EXPERIMENT_ONLY',
    'train_frames':len(train),'dev_frames':len(dev),'test_frames':0,
    'train_sha256':manifest['train_sha256'],'dev_sha256':manifest['dev_sha256'],
    'foundation_model_sha256':EXPECTED_FOUNDATION,
    'torch':torch.__version__,'mace_torch':metadata.version('mace-torch'),'ase':metadata.version('ase'),
    'free_disk_bytes':free,'official_split_mutated':False,'withheld_labels_opened':False,
},indent=2))
