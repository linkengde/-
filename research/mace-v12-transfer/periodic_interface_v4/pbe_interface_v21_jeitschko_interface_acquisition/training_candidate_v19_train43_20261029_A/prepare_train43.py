#!/usr/bin/env python3
"""Build an experiment-only train43 set by appending four verified Ti/C labels to train39."""
from __future__ import annotations
import hashlib, io, json, subprocess
from pathlib import Path
import numpy as np
from ase.io import read, write

HERE = Path(__file__).resolve().parent
TASK = HERE.parent
REPO = Path('/workspace/jeitschko-main-run')
TRAIN39 = TASK / 'training_candidate_v19_train39_20261029_A/dataset'
DEST = HERE / 'dataset'
SOURCE_COMMIT = '847676691c89faa78750802d2217bea9d8702f85'
SOURCE_TASK = 'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_bulk_training_acquisition'
FOUNDATION_SHA256 = '2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'
TRAIN39_SHA256 = 'd2336a650abadf46ee0c4f74496bff6264e0c406af2a5f5d2495cd1e67af7434'
DEV3_SHA256 = 'e9e86e3a353c50c8173cefede370e2653ec2e81506f4b8122b01b7a8efe76e4c'
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

def parsed_sums(data: bytes) -> dict[str, str]:
    result = {}
    for line in data.decode().splitlines():
        sha, name = line.split('  ', 1)
        result[name] = sha
    return result

def main() -> None:
    if DEST.exists():
        raise SystemExit(f'Refusing to overwrite existing candidate dataset: {DEST}')
    if digest((TRAIN39 / 'train.extxyz').read_bytes()) != TRAIN39_SHA256:
        raise SystemExit('Frozen train39 file hash mismatch')
    if digest((TRAIN39 / 'valid.extxyz').read_bytes()) != DEV3_SHA256:
        raise SystemExit('Frozen dev3 file hash mismatch')
    base_bytes = (TRAIN39 / 'train.extxyz').read_bytes()
    dev_bytes = (TRAIN39 / 'valid.extxyz').read_bytes()
    train39_manifest = json.loads((TRAIN39 / 'manifest.json').read_text())
    if train39_manifest.get('status') != 'NON_OFFICIAL_V19_TRAIN39_EXPERIMENT_CANDIDATE':
        raise SystemExit('Unexpected train39 base status')
    if subprocess.run(['git', '-C', str(REPO), 'merge-base', '--is-ancestor', SOURCE_COMMIT, 'HEAD']).returncode != 0:
        raise SystemExit('Pinned bulk-label source commit is not an ancestor of current main worktree HEAD')

    blocks, rows = [], []
    for label in LABELS:
        folder = f'{SOURCE_TASK}/calculations/{label}'
        result_name = f'{label}_PW_PBE.extxyz'
        raw = blob(f'{folder}/{result_name}')
        sums = parsed_sums(blob(f'{folder}/SHA256SUMS.txt'))
        if digest(raw) != sums.get(result_name):
            raise SystemExit(f'Per-result archive checksum mismatch: {label}')
        summary = json.loads(blob(f'{folder}/summary.json'))
        verification = json.loads(blob(f'{folder}/verification.json'))
        if verification.get('status') != 'PASS' or summary.get('label') != label:
            raise SystemExit(f'Archived label verifier did not pass: {label}')
        if summary.get('source_family') != 'Jeitschko_Nowotny_1967_single_crystal':
            raise SystemExit(f'Unexpected source family: {label}')
        input_raw = blob(f"{SOURCE_TASK}/{summary['input']}")
        if digest(input_raw) != summary.get('input_sha256'):
            raise SystemExit(f'Registered calculation input hash mismatch: {label}')
        atoms = read(io.StringIO(raw.decode()), format='extxyz')
        if len(atoms) != 12 or atoms.get_chemical_formula() != 'C4Si2Ti6':
            raise SystemExit(f'Unexpected formula or atom count: {label}')
        native = float(summary['native_extrapolated_energy_eV_cell'])
        free = float(summary['free_energy_eV_cell'])
        forces = np.asarray(atoms.arrays['PW_PBE_forces'], dtype=float)
        if atoms.info.get('PW_PBE_energy_eV') != native or not np.isfinite(forces).all():
            raise SystemExit(f'Nonfinite or mismatched energy/force label: {label}')
        atoms.info['REF_energy'] = native
        atoms.info['free_energy_eV'] = free
        atoms.info['dataset_role'] = 'v19_experiment_train43_candidate'
        atoms.info['source_family'] = summary['source_family']
        atoms.info['source_doi'] = summary['source_doi']
        atoms.info['source_label_sha256'] = digest(raw)
        atoms.info['source_input_sha256'] = summary['input_sha256']
        atoms.info['source_repository_commit'] = SOURCE_COMMIT
        atoms.info['energy_convention'] = 'native GPAW extrapolated energy (force_consistent=False); free energy retained separately'
        atoms.info['source_method'] = 'GPAW 26.7.0 PBE/PW 500 eV, 10x10x2, Fermi-Dirac sigma 0.10 eV, TTT'
        atoms.arrays['REF_forces'] = forces.copy()
        tmp = HERE / f'.{label}.tmp.extxyz'
        write(tmp, atoms, format='extxyz')
        block_bytes = tmp.read_bytes()
        tmp.unlink()
        if not block_bytes.endswith(b'\n'):
            block_bytes += b'\n'
        blocks.append(block_bytes)
        rows.append({
            'config_type': label,
            'source_family': summary['source_family'],
            'source_doi': summary['source_doi'],
            'source_repository_commit': SOURCE_COMMIT,
            'source_repository_path': f'{folder}/{result_name}',
            'source_input_sha256': summary['input_sha256'],
            'source_label_sha256': digest(raw),
            'source_result_archive_sha256': digest(blob(f'{folder}/SHA256SUMS.txt')),
            'source_verification_sha256': digest(blob(f'{folder}/verification.json')),
            'native_extrapolated_energy_eV_cell': native,
            'free_energy_eV_cell': free,
            'ref_energy_key': 'REF_energy (native extrapolated energy)',
            'ref_force_key': 'REF_forces (GPAW/ASE forces at sigma 0.10 eV)',
            'geometry_label': summary.get('geometry_label', summary['label']),
            'method': summary['method'],
            'scf_iterations': summary['scf_iterations'],
            'role': 'experiment-only training candidate; owner-side external role review waived by user; not an official train assignment',
        })

    DEST.mkdir(parents=True)
    train_bytes = base_bytes + b''.join(blocks)
    (DEST / 'train.extxyz').write_bytes(train_bytes)
    (DEST / 'valid.extxyz').write_bytes(dev_bytes)
    train = read(DEST / 'train.extxyz', index=':')
    dev = read(DEST / 'valid.extxyz', index=':')
    if len(train) != 43 or len(dev) != 3:
        raise SystemExit(f'Unexpected frame counts: train={len(train)} dev={len(dev)}')
    if (DEST / 'train.extxyz').read_bytes()[:len(base_bytes)] != base_bytes:
        raise SystemExit('Train39 byte prefix changed')
    for atoms in train + dev:
        energy = float(atoms.info['REF_energy'])
        forces = np.asarray(atoms.arrays['REF_forces'], dtype=float)
        if not np.isfinite(energy) or forces.shape != (len(atoms), 3) or not np.isfinite(forces).all():
            raise SystemExit('Nonfinite or malformed training labels')
    elements = ['Ag', 'C', 'Si', 'Ti']
    matrix = np.asarray([[a.get_chemical_symbols().count(e) for e in elements] for a in train], dtype=float)
    matrix39 = matrix[:39]
    sv43 = np.linalg.svd(matrix, compute_uv=False)
    sv39 = np.linalg.svd(matrix39, compute_uv=False)
    manifest = {
        'status': 'NON_OFFICIAL_V19_TRAIN43_EXPERIMENT_CANDIDATE',
        'base_train39_sha256': TRAIN39_SHA256,
        'base_dev3_sha256': DEV3_SHA256,
        'base_candidate_manifest_sha256': digest((TRAIN39 / 'manifest.json').read_bytes()),
        'foundation_model_sha256': FOUNDATION_SHA256,
        'roles': {'train_candidate': 43, 'development_validation_frozen': 3, 'test': 0},
        'withheld_labels_opened': False,
        'official_train36_and_dev3_mutated': False,
        'train39_byte_prefix_preserved': True,
        'energy_target': 'native extrapolated PW_PBE_energy_eV copied to REF_energy; free_energy_eV retained separately',
        'force_target': 'PW_PBE_forces copied to REF_forces; GPAW/ASE forces at sigma 0.10 eV',
        'added_candidate_rows': rows,
        'candidate_source_group_count': 1,
        'candidate_source_family': 'Jeitschko_Nowotny_1967_single_crystal',
        'source_task_commit': SOURCE_COMMIT,
        'train_sha256': digest(train_bytes),
        'dev_sha256': digest(dev_bytes),
        'composition_matrix': {
            'elements': elements,
            'rank': int(np.linalg.matrix_rank(matrix)),
            'train39_singular_values_descending': sv39.tolist(),
            'train39_condition_number': float(sv39[0] / sv39[-1]),
            'train43_singular_values_descending': sv43.tolist(),
            'train43_condition_number': float(sv43[0] / sv43[-1]),
        },
        'limitations': [
            'The four added displacement labels share one Jeitschko–Nowotny 1967 source family and geometry reconstruction.',
            'The original paper discusses possible nonstoichiometry; these calculations assume ideal full occupancy.',
            'The four bulk labels use 10x10x2/TTT; train39 contains Gamma/TTT and 6x6x1/TTF labels.',
            'Native extrapolated energy and finite-smearing forces follow the existing dataset convention; exact energy-force gradient consistency is not claimed.',
            'Direct geometry/record review is provenance evidence, not a universal supercell-equivalence audit.',
            'This versioned experiment is not an official split or production model approval.',
        ],
    }
    (DEST / 'manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + '\n')
    entries = []
    for name in ['train.extxyz', 'valid.extxyz', 'manifest.json']:
        entries.append(f'{digest((DEST/name).read_bytes())}  {name}\n')
    (DEST / 'SHA256SUMS.txt').write_text(''.join(entries))
    print(json.dumps({
        'status': manifest['status'], 'train_frames': len(train), 'dev_frames': len(dev),
        'train_sha256': manifest['train_sha256'], 'dev_sha256': manifest['dev_sha256'],
        'train39_condition_number': manifest['composition_matrix']['train39_condition_number'],
        'train43_condition_number': manifest['composition_matrix']['train43_condition_number'],
        'official_split_mutated': False, 'withheld_labels_opened': False,
    }, indent=2))

if __name__ == '__main__':
    main()
