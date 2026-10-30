#!/usr/bin/env python3
"""Read-only provenance, dataset and runtime preflight for experiment-only train49."""
from __future__ import annotations
import hashlib, importlib.metadata as metadata, io, json, shutil, subprocess
from pathlib import Path
import numpy as np
import torch
from ase.io import read

HERE = Path(__file__).resolve().parent
TASK = HERE.parent
REPO = Path('/workspace/jeitschko-main-run')
DATA = HERE / 'dataset'
BASE = TASK / 'training_candidate_v19_train43_20261029_A/dataset'
FOUNDATION = Path('/workspace/-/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model')
FOUNDATION_SHA256 = '2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'
SOURCE_COMMIT = '5d1e734c16d17dc9b460e3184eab4670ccc2739c'
SOURCE_TASK = 'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v22_force_response_acquisition'
LABELS = [
    'AgTi_JN1967_registryB_gap_m0p20_sigma0p10_k6x6x1_v22',
    'AgTi_JN1967_registryB_gap_p0p20_sigma0p10_k6x6x1_v22',
    'AgSi_JN1967_registryB_gap_m0p20_sigma0p10_k6x6x1_v22',
    'AgSi_JN1967_registryB_gap_p0p20_sigma0p10_k6x6x1_v22',
    'AgC_JN1967_registryB_gap_m0p20_sigma0p10_k6x6x1_v22',
    'AgC_JN1967_registryB_gap_p0p20_sigma0p10_k6x6x1_v22',
]

def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def blob(path: str) -> bytes:
    return subprocess.check_output(['git', '-C', str(REPO), 'show', f'{SOURCE_COMMIT}:{path}'])

def main() -> None:
    manifest = json.loads((DATA / 'manifest.json').read_text())
    if manifest.get('status') != 'NON_OFFICIAL_V19_TRAIN49_EXPERIMENT_CANDIDATE':
        raise SystemExit('Unexpected candidate status')
    if manifest.get('roles') != {'train_candidate': 49, 'development_validation_frozen': 3, 'test': 0}:
        raise SystemExit('Unexpected frame roles or counts')
    if manifest.get('official_train36_and_dev3_mutated') or manifest.get('withheld_labels_opened'):
        raise SystemExit('Official split was changed or withheld data was opened')
    if digest((DATA / 'train.extxyz').read_bytes()) != manifest['train_sha256']:
        raise SystemExit('train49 hash mismatch')
    if digest((DATA / 'valid.extxyz').read_bytes()) != manifest['dev_sha256']:
        raise SystemExit('Frozen dev3 hash mismatch')
    base_train = (BASE / 'train.extxyz').read_bytes()
    if digest(base_train) != manifest['base_train43_sha256']:
        raise SystemExit('train43 base hash mismatch')
    if digest((BASE / 'valid.extxyz').read_bytes()) != manifest['base_dev3_sha256']:
        raise SystemExit('Frozen dev3 base hash mismatch')
    if digest((DATA / 'train.extxyz').read_bytes()[:len(base_train)]) != digest(base_train):
        raise SystemExit('train43 byte prefix changed')
    if digest(FOUNDATION.read_bytes()) != FOUNDATION_SHA256 or manifest['foundation_model_sha256'] != FOUNDATION_SHA256:
        raise SystemExit('Foundation model hash mismatch')
    if subprocess.run(['git', '-C', str(REPO), 'merge-base', '--is-ancestor', SOURCE_COMMIT, 'HEAD']).returncode != 0:
        raise SystemExit('Pinned v22 source commit is not an ancestor of current HEAD')
    sums = {}
    for line in (DATA / 'SHA256SUMS.txt').read_text().splitlines():
        sha, name = line.split('  ', 1)
        sums[name] = sha
    for name in ['train.extxyz', 'valid.extxyz', 'manifest.json']:
        if sums.get(name) != digest((DATA / name).read_bytes()):
            raise SystemExit(f'Dataset checksum inventory mismatch: {name}')
    train = read(DATA / 'train.extxyz', index=':')
    dev = read(DATA / 'valid.extxyz', index=':')
    if len(train) != 49 or len(dev) != 3 or len(manifest['added_candidate_rows']) != 6:
        raise SystemExit(f'Frame count mismatch: train={len(train)} dev={len(dev)}')
    if [r['config_type'] for r in manifest['added_candidate_rows']] != LABELS:
        raise SystemExit('Candidate label order mismatch')
    for atoms in train + dev:
        energy = float(atoms.info['REF_energy'])
        forces = np.asarray(atoms.arrays['REF_forces'], dtype=float)
        if not np.isfinite(energy) or forces.shape != (len(atoms), 3) or not np.isfinite(forces).all():
            raise SystemExit('Nonfinite or malformed energy/force labels')

    for label, row, atoms in zip(LABELS, manifest['added_candidate_rows'], train[-6:]):
        folder = f'{SOURCE_TASK}/calculations/{label}'
        result_name = f'{label}_PW_PBE.extxyz'
        raw = blob(f'{folder}/{result_name}')
        inventory_raw = blob(f'{folder}/SHA256SUMS.txt')
        inventory = {line.split('  ', 1)[1]: line.split('  ', 1)[0]
                     for line in inventory_raw.decode().splitlines()}
        for name, expected in inventory.items():
            if digest(blob(f'{folder}/{name}')) != expected:
                raise SystemExit(f'Archived result checksum mismatch: {label}/{name}')
        summary = json.loads(blob(f'{folder}/summary.json'))
        verification_raw = blob(f'{folder}/verification.json')
        verification = json.loads(verification_raw)
        input_raw = blob(f"{SOURCE_TASK}/{summary['input']}")
        if (digest(raw) != inventory.get(result_name) or digest(raw) != row['source_label_sha256']
                or digest(inventory_raw) != row['source_result_archive_sha256']
                or digest(verification_raw) != row['source_verification_sha256']):
            raise SystemExit(f'Archived source hashes differ from manifest: {label}')
        if (verification.get('status') != 'PASS' or summary.get('scf_converged') is not True
                or summary.get('label') != label or digest(input_raw) != summary.get('source_sha256')):
            raise SystemExit(f'Archived source verification/input mismatch: {label}')
        source = read(io.StringIO(raw.decode()), format='extxyz')
        if (len(source) != 52 or summary['source_family'] != manifest['candidate_source_family']
                or atoms.info.get('source_label_sha256') != digest(raw)
                or atoms.info.get('source_input_sha256') != summary['source_sha256']
                or atoms.info.get('REF_energy') != summary['native_extrapolated_energy_eV_cell']
                or atoms.info.get('free_energy_eV') != summary['free_energy_eV_cell']):
            raise SystemExit(f'Label identity or energy mapping mismatch: {label}')
        if (not np.array_equal(atoms.positions, source.positions)
                or not np.array_equal(atoms.cell.array, source.cell.array)
                or not np.array_equal(atoms.pbc, source.pbc)
                or not np.array_equal(atoms.arrays['REF_forces'], source.arrays['PW_PBE_forces'])):
            raise SystemExit(f'Geometry or force mapping mismatch: {label}')
    if (DATA / 'test.extxyz').exists():
        raise SystemExit('Unexpected test/sealed dataset file present')
    for name, version in [('torch', '2.5.1+cpu'), ('mace-torch', '0.3.16'), ('ase', '3.29.0')]:
        if metadata.version(name) != version:
            raise SystemExit(f'Runtime version mismatch for {name}: {metadata.version(name)}')
    if torch.cuda.is_available():
        raise SystemExit('Expected CPU-only training environment')
    free = shutil.disk_usage(HERE).free
    if free < 4_000_000_000:
        raise SystemExit(f'Insufficient free disk for epoch checkpoints: {free}')
    run = HERE / 'run'
    if run.exists() and any(run.iterdir()):
        raise SystemExit('Refusing to overwrite a nonempty train49 run directory')
    active = subprocess.run(['pgrep', '-af', '[m]ace_run_train'], capture_output=True, text=True).stdout.strip()
    if active:
        raise SystemExit(f'Another MACE training process appears active: {active}')
    print(json.dumps({'status': 'V19_TRAIN49_PREFLIGHT_PASS_EXPERIMENT_ONLY',
                      'train_frames': len(train), 'dev_frames': len(dev), 'test_frames': 0,
                      'train_sha256': manifest['train_sha256'], 'dev_sha256': manifest['dev_sha256'],
                      'foundation_model_sha256': FOUNDATION_SHA256, 'source_commit': SOURCE_COMMIT,
                      'source_label_count': len(LABELS), 'torch': torch.__version__,
                      'mace_torch': metadata.version('mace-torch'), 'ase': metadata.version('ase'),
                      'free_disk_bytes': free, 'official_split_mutated': False,
                      'withheld_labels_opened': False}, indent=2))

if __name__ == '__main__':
    main()
