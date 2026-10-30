#!/usr/bin/env python3
"""Build a separate train49 experiment from train43 plus six verified v22 labels."""
from __future__ import annotations
import hashlib, io, json, subprocess
from pathlib import Path
import numpy as np
from ase.io import read, write

HERE = Path(__file__).resolve().parent
TASK = HERE.parent
REPO = Path('/workspace/jeitschko-main-run')
BASE_DATA = TASK / 'training_candidate_v19_train43_20261029_A/dataset'
DEST = HERE / 'dataset'
SOURCE_COMMIT = '5d1e734c16d17dc9b460e3184eab4670ccc2739c'
SOURCE_TASK = 'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v22_force_response_acquisition'
FOUNDATION_SHA256 = '2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'
BASE_TRAIN_SHA256 = '5895c130b9ff1350e1da88dbbab1ca6b671b42a4fdc0a349870472a1f1c1b97b'
DEV3_SHA256 = 'e9e86e3a353c50c8173cefede370e2653ec2e81506f4b8122b01b7a8efe76e4c'
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

def parse_sums(raw: bytes) -> dict[str, str]:
    return {line.split('  ', 1)[1]: line.split('  ', 1)[0]
            for line in raw.decode().splitlines()}

def main() -> None:
    if DEST.exists():
        raise SystemExit(f'Refusing to overwrite existing train49 dataset: {DEST}')
    if digest((BASE_DATA / 'train.extxyz').read_bytes()) != BASE_TRAIN_SHA256:
        raise SystemExit('Frozen train43 file hash mismatch')
    if digest((BASE_DATA / 'valid.extxyz').read_bytes()) != DEV3_SHA256:
        raise SystemExit('Frozen dev3 file hash mismatch')
    base_bytes = (BASE_DATA / 'train.extxyz').read_bytes()
    dev_bytes = (BASE_DATA / 'valid.extxyz').read_bytes()
    base_manifest_raw = (BASE_DATA / 'manifest.json').read_bytes()
    base_manifest = json.loads(base_manifest_raw)
    if base_manifest.get('status') != 'NON_OFFICIAL_V19_TRAIN43_EXPERIMENT_CANDIDATE':
        raise SystemExit('Unexpected base candidate status')
    if subprocess.run(['git', '-C', str(REPO), 'merge-base', '--is-ancestor', SOURCE_COMMIT, 'HEAD']).returncode != 0:
        raise SystemExit('Pinned v22 source commit is not an ancestor of current worktree HEAD')

    base_frames = read(io.StringIO(base_bytes.decode()), index=':', format='extxyz')
    dev_frames = read(io.StringIO(dev_bytes.decode()), index=':', format='extxyz')
    blocks, rows = [], []
    for label in LABELS:
        folder = f'{SOURCE_TASK}/calculations/{label}'
        result_name = f'{label}_PW_PBE.extxyz'
        raw = blob(f'{folder}/{result_name}')
        inventory_raw = blob(f'{folder}/SHA256SUMS.txt')
        inventory = parse_sums(inventory_raw)
        if result_name not in inventory or digest(raw) != inventory[result_name]:
            raise SystemExit(f'Result checksum mismatch: {label}')
        for name, expected in inventory.items():
            if digest(blob(f'{folder}/{name}')) != expected:
                raise SystemExit(f'Archived file checksum mismatch: {label}/{name}')
        summary = json.loads(blob(f'{folder}/summary.json'))
        verification = json.loads(blob(f'{folder}/verification.json'))
        if verification.get('status') != 'PASS' or summary.get('label') != label or summary.get('scf_converged') is not True:
            raise SystemExit(f'Archived DFT verification did not pass: {label}')
        if summary.get('source_family') != 'Jeitschko_Nowotny_1967_single_crystal':
            raise SystemExit(f'Unexpected source family: {label}')
        input_raw = blob(f"{SOURCE_TASK}/{summary['input']}")
        if digest(input_raw) != summary.get('source_sha256'):
            raise SystemExit(f'Registered input hash mismatch: {label}')
        atoms = read(io.StringIO(raw.decode()), format='extxyz')
        if len(atoms) != 52 or atoms.get_chemical_formula() != 'C16Ag4Si8Ti24':
            raise SystemExit(f'Unexpected geometry/formula: {label}')
        forces = np.asarray(atoms.arrays['PW_PBE_forces'], dtype=float)
        native = float(summary['native_extrapolated_energy_eV_cell'])
        free = float(summary['free_energy_eV_cell'])
        if not np.isfinite(forces).all() or atoms.info.get('PW_PBE_energy_eV') != native:
            raise SystemExit(f'Invalid force or energy label: {label}')
        for previous in base_frames + dev_frames:
            if (len(previous) == len(atoms) and previous.get_chemical_symbols() == atoms.get_chemical_symbols()
                    and np.array_equal(previous.positions, atoms.positions)
                    and np.array_equal(previous.cell.array, atoms.cell.array)
                    and np.array_equal(previous.pbc, atoms.pbc)):
                raise SystemExit(f'Exact duplicate geometry found in frozen base: {label}')

        atoms.info['REF_energy'] = native
        atoms.info['free_energy_eV'] = free
        atoms.info['dataset_role'] = 'v19_experiment_train49_candidate'
        atoms.info['source_family'] = summary['source_family']
        atoms.info['source_doi'] = summary['source_doi']
        atoms.info['source_label_sha256'] = digest(raw)
        atoms.info['source_input_sha256'] = summary['source_sha256']
        atoms.info['source_repository_commit'] = SOURCE_COMMIT
        atoms.info['source_result_archive_sha256'] = digest(inventory_raw)
        atoms.info['energy_convention'] = 'native GPAW extrapolated energy (force_consistent=False); free energy retained separately'
        atoms.info['source_method'] = 'GPAW 26.7.0 PBE/PW 500 eV, 6x6x1, Fermi-Dirac sigma 0.10 eV, xy dipole correction'
        atoms.arrays['REF_forces'] = forces.copy()
        block = io.StringIO()
        write(block, atoms, format='extxyz')
        block_bytes = block.getvalue().encode()
        if not block_bytes.endswith(b'\n'):
            block_bytes += b'\n'
        blocks.append(block_bytes)
        rows.append({
            'config_type': label,
            'source_family': summary['source_family'],
            'source_doi': summary['source_doi'],
            'source_repository_commit': SOURCE_COMMIT,
            'source_repository_path': f'{folder}/{result_name}',
            'source_input_sha256': summary['source_sha256'],
            'source_label_sha256': digest(raw),
            'source_result_archive_sha256': digest(inventory_raw),
            'source_verification_sha256': digest(blob(f'{folder}/verification.json')),
            'native_extrapolated_energy_eV_cell': native,
            'free_energy_eV_cell': free,
            'ref_energy_key': 'REF_energy (native extrapolated energy)',
            'ref_force_key': 'REF_forces (GPAW/ASE forces at sigma 0.10 eV)',
            'geometry_label': summary.get('geometry_label', label),
            'method': summary['method'],
            'scf_iterations': summary['scf_iterations'],
            'role': 'experiment-only acquisition candidate; not an official train assignment',
        })

    DEST.mkdir(parents=True)
    train_bytes = base_bytes + b''.join(blocks)
    (DEST / 'train.extxyz').write_bytes(train_bytes)
    (DEST / 'valid.extxyz').write_bytes(dev_bytes)
    train = read(DEST / 'train.extxyz', index=':')
    dev = read(DEST / 'valid.extxyz', index=':')
    if len(train) != 49 or len(dev) != 3 or train_bytes[:len(base_bytes)] != base_bytes:
        raise SystemExit(f'Unexpected frame counts or changed train43 prefix: train={len(train)} dev={len(dev)}')
    for atoms in train + dev:
        energy = float(atoms.info['REF_energy'])
        forces = np.asarray(atoms.arrays['REF_forces'], dtype=float)
        if not np.isfinite(energy) or forces.shape != (len(atoms), 3) or not np.isfinite(forces).all():
            raise SystemExit('Nonfinite or malformed training labels')
    elements = ['Ag', 'C', 'Si', 'Ti']
    matrix = np.asarray([[a.get_chemical_symbols().count(e) for e in elements] for a in train], dtype=float)
    singular = np.linalg.svd(matrix, compute_uv=False)
    manifest = {
        'status': 'NON_OFFICIAL_V19_TRAIN49_EXPERIMENT_CANDIDATE',
        'base_train43_sha256': BASE_TRAIN_SHA256,
        'base_dev3_sha256': DEV3_SHA256,
        'base_train43_manifest_sha256': digest(base_manifest_raw),
        'foundation_model_sha256': FOUNDATION_SHA256,
        'roles': {'train_candidate': 49, 'development_validation_frozen': 3, 'test': 0},
        'withheld_labels_opened': False,
        'official_train36_and_dev3_mutated': False,
        'train43_byte_prefix_preserved': True,
        'energy_target': 'native extrapolated PW_PBE_energy_eV copied to REF_energy; finite-width free_energy_eV retained separately',
        'force_target': 'PW_PBE_forces copied to REF_forces; GPAW/ASE forces at sigma 0.10 eV',
        'added_candidate_rows': rows,
        'candidate_source_group_count': 1,
        'candidate_source_family': 'Jeitschko_Nowotny_1967_single_crystal',
        'source_task_commit': SOURCE_COMMIT,
        'train_sha256': digest(train_bytes),
        'dev_sha256': digest(dev_bytes),
        'composition_matrix': {'elements': elements, 'rank': int(np.linalg.matrix_rank(matrix)),
                               'singular_values_descending': singular.tolist(),
                               'condition_number': float(singular[0] / singular[-1])},
        'exact_geometry_duplicate_screen_pass': True,
        'limitations': [
            'The six added interface labels share the Jeitschko–Nowotny 1967 source family and registry-B slab construction.',
            'All six labels sample rigid Ag displacement at fixed geometry; they are not independent relaxed interfaces.',
            'Native extrapolated energy and finite-smearing forces follow the existing convention; exact energy-force gradient consistency is not claimed.',
            'The candidate mixes Gamma, 6x6x1, and 10x10x2 labels and differing boundary conditions.',
            'The three frozen development labels share the COD 9009647 source family.',
            'This experiment is not an official split or production model approval; slab convergence is not fully established.',
        ],
    }
    (DEST / 'manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + '\n')
    entries = [f'{digest((DEST/name).read_bytes())}  {name}\n'
               for name in ['train.extxyz', 'valid.extxyz', 'manifest.json']]
    (DEST / 'SHA256SUMS.txt').write_text(''.join(entries))
    print(json.dumps({'status': manifest['status'], 'train_frames': len(train), 'dev_frames': len(dev),
                      'train_sha256': manifest['train_sha256'], 'dev_sha256': manifest['dev_sha256'],
                      'condition_number': manifest['composition_matrix']['condition_number'],
                      'exact_geometry_duplicate_screen_pass': True,
                      'official_split_mutated': False, 'withheld_labels_opened': False}, indent=2))

if __name__ == '__main__':
    main()
