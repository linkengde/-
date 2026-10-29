#!/usr/bin/env python3
"""Read-only preflight for the experiment-only train43 refit."""
from __future__ import annotations
import hashlib, importlib.metadata as metadata, json, shutil, subprocess
from pathlib import Path
import numpy as np
import torch
from ase.io import read

HERE = Path(__file__).resolve().parent
TASK = HERE.parent
REPO = Path('/workspace/jeitschko-main-run')
DATA = HERE / 'dataset'
TRAIN39 = TASK / 'training_candidate_v19_train39_20261029_A/dataset'
FOUNDATION = Path('/workspace/-/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model')
EXPECTED_FOUNDATION = '2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'
SOURCE_COMMIT = '847676691c89faa78750802d2217bea9d8702f85'
SOURCE_TASK = 'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_bulk_training_acquisition'
LABELS = [
    'Ti3SiC2_JN1967_primitive_Ti4f_id0003_z_m002A_k10x10x2_v21',
    'Ti3SiC2_JN1967_primitive_Ti4f_id0003_z_p002A_k10x10x2_v21',
    'Ti3SiC2_JN1967_primitive_C4f_id0009_z_m002A_k10x10x2_v21',
    'Ti3SiC2_JN1967_primitive_C4f_id0009_z_p002A_k10x10x2_v21',
]

def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def blob(path: str) -> bytes:
    return subprocess.check_output(['git', '-C', str(REPO), 'show', f'{SOURCE_COMMIT}:{path}'])

def main() -> None:
    (DATA / 'SHA256SUMS.txt').read_text()
    manifest = json.loads((DATA / 'manifest.json').read_text())
    if manifest['status'] != 'NON_OFFICIAL_V19_TRAIN43_EXPERIMENT_CANDIDATE':
        raise SystemExit('Unexpected train43 data status')
    if manifest['roles'] != {'train_candidate': 43, 'development_validation_frozen': 3, 'test': 0}:
        raise SystemExit('Unexpected roles or frame counts')
    if manifest['official_train36_and_dev3_mutated'] or manifest['withheld_labels_opened']:
        raise SystemExit('Official inputs mutated or withheld labels opened')
    if digest((DATA / 'train.extxyz').read_bytes()) != manifest['train_sha256']:
        raise SystemExit('train43 file digest mismatch')
    if digest((DATA / 'valid.extxyz').read_bytes()) != manifest['dev_sha256']:
        raise SystemExit('Frozen dev3 file digest mismatch')
    if digest((TRAIN39 / 'train.extxyz').read_bytes()) != manifest['base_train39_sha256']:
        raise SystemExit('train39 base digest mismatch')
    if digest((DATA / 'train.extxyz').read_bytes()[:(TRAIN39 / 'train.extxyz').stat().st_size]) != manifest['base_train39_sha256']:
        raise SystemExit('train39 exact prefix not preserved')
    if digest(FOUNDATION.read_bytes()) != EXPECTED_FOUNDATION or EXPECTED_FOUNDATION != manifest['foundation_model_sha256']:
        raise SystemExit('Foundation-model digest mismatch')
    if subprocess.run(['git', '-C', str(REPO), 'merge-base', '--is-ancestor', SOURCE_COMMIT, 'HEAD']).returncode != 0:
        raise SystemExit('Pinned source commit is not an ancestor of current main worktree HEAD')
    train = read(DATA / 'train.extxyz', index=':')
    dev = read(DATA / 'valid.extxyz', index=':')
    if len(train) != 43 or len(dev) != 3:
        raise SystemExit(f'Frame count mismatch: train={len(train)} dev={len(dev)}')
    for atoms in train + dev:
        energy = float(atoms.info['REF_energy'])
        forces = np.asarray(atoms.arrays['REF_forces'], dtype=float)
        if not np.isfinite(energy) or forces.shape != (len(atoms), 3) or not np.isfinite(forces).all():
            raise SystemExit('Nonfinite or malformed energy/force labels')
    rows = manifest['added_candidate_rows']
    if [r['config_type'] for r in rows] != LABELS:
        raise SystemExit('Added candidate order or labels differ from manifest')
    for label, row, atoms in zip(LABELS, rows, train[-4:]):
        folder = f'{SOURCE_TASK}/calculations/{label}'
        source_name = f'{label}_PW_PBE.extxyz'
        raw = blob(f'{folder}/{source_name}')
        sums_text = blob(f'{folder}/SHA256SUMS.txt').decode()
        sums = {line.split('  ', 1)[1]: line.split('  ', 1)[0] for line in sums_text.splitlines()}
        summary = json.loads(blob(f'{folder}/summary.json'))
        verification_bytes = blob(f'{folder}/verification.json')
        verification = json.loads(verification_bytes)
        input_bytes = blob(f"{SOURCE_TASK}/{summary['input']}")
        if digest(raw) != sums.get(source_name) or digest(raw) != row['source_label_sha256']:
            raise SystemExit(f'Raw result archive hash mismatch: {label}')
        if digest(verification_bytes) != row['source_verification_sha256'] or verification.get('status') != 'PASS':
            raise SystemExit(f'Archived source verification missing or invalid: {label}')
        if digest(blob(f'{folder}/SHA256SUMS.txt')) != row['source_result_archive_sha256']:
            raise SystemExit(f'Source checksum inventory changed: {label}')
        if digest(input_bytes) != summary['input_sha256'] or summary['source_family'] != manifest['candidate_source_family']:
            raise SystemExit(f'Input provenance or source family mismatch: {label}')
        if row['source_repository_commit'] != SOURCE_COMMIT:
            raise SystemExit(f'Unexpected source commit: {label}')
        if atoms.info['REF_energy'] != summary['native_extrapolated_energy_eV_cell']:
            raise SystemExit(f'Energy target mismatch: {label}')
        if atoms.info['source_label_sha256'] != digest(raw):
            raise SystemExit(f'Embedded source label hash mismatch: {label}')
        if not np.array_equal(atoms.arrays['REF_forces'], atoms.arrays['PW_PBE_forces']):
            raise SystemExit(f'Force target mismatch: {label}')
        if not np.isclose(float(atoms.info['free_energy_eV']), summary['free_energy_eV_cell'], rtol=0, atol=1e-12):
            raise SystemExit(f'Free-energy diagnostic mismatch: {label}')
    if (DATA / 'test.extxyz').exists():
        raise SystemExit('Unexpected test/sealed dataset file present')
    if digest(FOUNDATION.read_bytes()) != EXPECTED_FOUNDATION:
        raise SystemExit('Foundation model changed')
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
        raise SystemExit('Refusing to overwrite a nonempty run directory')
    active = subprocess.run(['pgrep', '-af', '[m]ace_run_train'], capture_output=True, text=True).stdout.strip()
    if active:
        raise SystemExit(f'Another MACE training process appears active: {active}')
    print(json.dumps({
        'status': 'V19_TRAIN43_PREFLIGHT_PASS_EXPERIMENT_ONLY',
        'train_frames': len(train), 'dev_frames': len(dev), 'test_frames': 0,
        'train_sha256': manifest['train_sha256'], 'dev_sha256': manifest['dev_sha256'],
        'foundation_model_sha256': EXPECTED_FOUNDATION, 'source_commit': SOURCE_COMMIT,
        'source_label_count': len(rows), 'torch': torch.__version__,
        'mace_torch': metadata.version('mace-torch'), 'ase': metadata.version('ase'),
        'free_disk_bytes': free, 'official_split_mutated': False, 'withheld_labels_opened': False,
    }, indent=2))

if __name__ == '__main__':
    main()
