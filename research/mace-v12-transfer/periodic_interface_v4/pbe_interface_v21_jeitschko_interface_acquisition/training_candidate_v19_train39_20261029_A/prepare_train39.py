#!/usr/bin/env python3
"""Assemble a non-official V19 train39 experiment from frozen train36 plus three audited labels."""
from __future__ import annotations
import hashlib, json, subprocess
from pathlib import Path
import numpy as np
from ase.io import read, write

HERE = Path(__file__).resolve().parent
TASK = HERE.parent
REPO = Path('/workspace/jeitschko-main-run')
BASE_ENTRY = Path('/workspace/-/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v19_residual_response')
BASE_DATA = BASE_ENTRY / 'data'
DEST = HERE / 'dataset'
LABELS = [
    'AgSi_JN1967_registryB_k6x6x1_sigma0p10_preflight_v21',
    'AgC_JN1967_registryB_k6x6x1_sigma0p10_acquisition_v21',
    'AgTi_JN1967_registryB_k6x6x1_sigma0p10_acquisition_v21',
]
EXPECTED_BASE_TRAIN = 'a247f7619e613297a1c11d7dc97abec888aa0bbb8cc6889bcbbe05825ade0c9c'
EXPECTED_BASE_DEV = 'e9e86e3a353c50c8173cefede370e2653ec2e81506f4b8122b01b7a8efe76e4c'
EXPECTED_FOUNDATION = '2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'

def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()
def sha(path: Path) -> str:
    return sha_bytes(path.read_bytes())
def json_write(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + '\n')
def main() -> None:
    if not BASE_DATA.is_dir():
        raise SystemExit(f'Base data not found: {BASE_DATA}')
    if sha(BASE_DATA / 'train.extxyz') != EXPECTED_BASE_TRAIN:
        raise SystemExit('Frozen V19 train36 hash mismatch')
    if sha(BASE_DATA / 'valid.extxyz') != EXPECTED_BASE_DEV:
        raise SystemExit('Frozen V19 dev3 hash mismatch')
    if sha(BASE_ENTRY.parent / 'mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model') != EXPECTED_FOUNDATION:
        raise SystemExit('Foundation model digest mismatch')
    head = subprocess.check_output(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip()
    if subprocess.run(['git', '-C', str(REPO), 'merge-base', '--is-ancestor', '3562256258aa16998bea3a48a16dd59cff62e4fa', head]).returncode != 0:
        raise SystemExit(f'Pinned main-worktree baseline is not an ancestor: {head}')
    if DEST.exists():
        raise SystemExit(f'Refusing to overwrite existing candidate dataset: {DEST}')
    DEST.mkdir(parents=True)
    base_bytes = (BASE_DATA / 'train.extxyz').read_bytes()
    if not base_bytes.endswith(b'\n'):
        raise SystemExit('Frozen train36 lacks final newline; preserve the source bytes exactly')
    blocks = []
    new_rows = []
    for label in LABELS:
        subprocess.run([
            '/workspace/.venvs/gpaw-mpi/bin/python', str(TASK / 'verify_preflight.py'), label
        ], check=True, stdout=subprocess.DEVNULL)
        folder = TASK / 'calculations' / label
        summary = json.loads((folder / 'summary.json').read_text())
        source = folder / f'{label}_PW_PBE.extxyz'
        inventory = {name: digest for digest, name in (line.split('  ', 1) for line in (folder / 'SHA256SUMS.txt').read_text().splitlines())}
        if sha(source) != inventory[source.name]:
            raise SystemExit(f'Archive checksum mismatch: {label}')
        atoms = read(source)
        native = float(summary['native_extrapolated_energy_eV_cell'])
        free = float(summary['free_energy_eV_cell'])
        forces = np.asarray(atoms.arrays['PW_PBE_forces'], dtype=float)
        if atoms.info.get('PW_PBE_energy_eV') != native or not np.isfinite(forces).all():
            raise SystemExit(f'Label values do not match summary: {label}')
        atoms.info['REF_energy'] = native
        atoms.info['free_energy_eV'] = free
        atoms.info['dataset_role'] = 'v19_experiment_train39_candidate'
        atoms.info['source_family'] = summary['source_family']
        atoms.info['source_doi'] = summary['source_doi']
        atoms.info['source_label_sha256'] = sha(source)
        atoms.info['source_input_sha256'] = summary['source_sha256']
        atoms.info['energy_convention'] = 'native GPAW extrapolated energy (force_consistent=False); free energy retained separately'
        atoms.info['source_method'] = 'GPAW 26.7.0 PBE/PW 500 eV, 6x6x1, Fermi-Dirac sigma 0.10 eV, xy dipole, TTF'
        atoms.arrays['REF_forces'] = forces.copy()
        temp = HERE / f'.{label}.extxyz.tmp'
        write(temp, atoms, format='extxyz')
        block = temp.read_bytes()
        temp.unlink()
        if not block.endswith(b'\n'):
            block += b'\n'
        blocks.append(block)
        new_rows.append({
            'config_type': label,
            'source_family': summary['source_family'],
            'source_doi': summary['source_doi'],
            'source_input_sha256': summary['source_sha256'],
            'source_label_sha256': sha(source),
            'result_archive_sha256': sha(folder / 'SHA256SUMS.txt'),
            'native_extrapolated_energy_eV_cell': native,
            'free_energy_eV_cell': free,
            'ref_energy_key': 'REF_energy (native extrapolated energy)',
            'ref_force_key': 'REF_forces (GPAW/ASE forces at sigma 0.10 eV)',
            'geometry_label': summary['geometry_label'],
            'method': summary['method'],
            'scf_iterations': summary['scf_iterations'],
        })
    train_bytes = base_bytes + b''.join(blocks)
    (DEST / 'train.extxyz').write_bytes(train_bytes)
    dev_bytes = (BASE_DATA / 'valid.extxyz').read_bytes()
    (DEST / 'valid.extxyz').write_bytes(dev_bytes)
    train = read(DEST / 'train.extxyz', index=':')
    dev = read(DEST / 'valid.extxyz', index=':')
    if len(train) != 39 or len(dev) != 3:
        raise SystemExit(f'Unexpected frame counts train={len(train)} dev={len(dev)}')
    if (DEST / 'train.extxyz').read_bytes()[:len(base_bytes)] != base_bytes:
        raise SystemExit('Original train36 bytes were not preserved as the exact file prefix')
    for atoms in train + dev:
        if not np.isfinite(float(atoms.info['REF_energy'])):
            raise SystemExit('Nonfinite REF_energy found')
        f = np.asarray(atoms.arrays['REF_forces'], dtype=float)
        if f.shape != (len(atoms), 3) or not np.isfinite(f).all():
            raise SystemExit('Invalid REF_forces found')
    elements = ['Ag', 'C', 'Si', 'Ti']
    matrix = np.asarray([[atoms.get_chemical_symbols().count(e) for e in elements] for atoms in train], dtype=float)
    base_matrix = matrix[:36]
    singular = np.linalg.svd(matrix, compute_uv=False)
    base_singular = np.linalg.svd(base_matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix))
    base_condition = float(base_singular[0] / base_singular[-1])
    condition = float(singular[0] / singular[-1])
    manifest = {
        'status': 'NON_OFFICIAL_V19_TRAIN39_EXPERIMENT_CANDIDATE',
        'base_main_commit': '3562256258aa16998bea3a48a16dd59cff62e4fa',
        'base_data_source_checkout': str(BASE_ENTRY),
        'base_train_sha256': EXPECTED_BASE_TRAIN,
        'base_dev3_sha256': EXPECTED_BASE_DEV,
        'foundation_model_sha256': EXPECTED_FOUNDATION,
        'roles': {'train_candidate': 39, 'development_validation_frozen': 3, 'test': 0},
        'withheld_labels_opened': False,
        'official_train36_and_dev3_mutated': False,
        'base_train36_byte_prefix_preserved': True,
        'energy_target': 'native extrapolated PW_PBE_energy_eV copied to REF_energy; free_energy_eV retained as a separate field',
        'force_target': 'PW_PBE_forces copied to REF_forces; GPAW/ASE forces at sigma 0.10 eV',
        'candidate_source_group_count': 1,
        'candidate_source_family': 'Jeitschko_Nowotny_1967_single_crystal',
        'candidate_source_doi': '10.1007/BF00899949',
        'candidate_rows': new_rows,
        'train_sha256': sha_bytes(train_bytes),
        'dev_sha256': sha_bytes(dev_bytes),
        'composition_matrix': {'elements': elements, 'rank': rank, 'base36_singular_values_descending': base_singular.tolist(), 'base36_condition_number': base_condition, 'singular_values_descending': singular.tolist(), 'condition_number': condition, 'condition_number_reduction_factor': base_condition / condition},
        'limitations': [
            'All three added frames are one crystallographic source family.',
            'Direct geometry screening found no exact/near hits against frozen train36 or dev3 in the audited scope; no universal supercell-equivalence search was performed.',
            'Candidate settings match dev3 slab settings and differ from train36 Gamma/TTT settings.',
            'Native extrapolated energy matches the existing V19 target convention; exact finite-smearing energy-force gradient consistency is not claimed.',
            'This is an experiment-only training set, not a versioned official split or production approval.'
        ]
    }
    json_write(DEST / 'manifest.json', manifest)
    inventory = []
    for name in ['train.extxyz', 'valid.extxyz', 'manifest.json']:
        inventory.append(f'{sha(DEST / name)}  {name}\n')
    (DEST / 'SHA256SUMS.txt').write_text(''.join(inventory))
    print(json.dumps({'status': manifest['status'], 'train_frames': len(train), 'dev_frames': len(dev), 'train_sha256': manifest['train_sha256'], 'dev_sha256': manifest['dev_sha256'], 'composition_rank': rank, 'condition_number': manifest['composition_matrix']['condition_number'], 'base36_condition_number': manifest['composition_matrix']['base36_condition_number']}, indent=2))
if __name__ == '__main__':
    main()
