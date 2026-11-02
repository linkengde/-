#!/usr/bin/env python3
"""Build a pinned Ti3SiC2-only acquisition set from the archived 1967 geometry."""
from __future__ import annotations

import hashlib
import io
import json
import subprocess
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
REPO = Path('/workspace/jeitschko-main-run')
SOURCE_COMMIT = '847676691c89faa78750802d2217bea9d8702f85'
SOURCE_ROOT = 'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_bulk_training_acquisition'
SOURCE_MANIFEST_PATH = f'{SOURCE_ROOT}/input_manifest.json'
SOURCE_CIF_SHA256 = '26960fc12d2de584c6a13ba5fe0d15713c07f697f455dbe29791b747f07359a1'
SOURCE_PAPER_SHA256 = '6765f99b3cc83daf3e48e4674dbeba109e5fffef8fa66cd53e266733ac4a4a63'


def git_blob(commit_path: str) -> bytes:
    return subprocess.check_output(['git', '-C', str(REPO), 'show', commit_path])


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_input(label: str, atoms, role: str, split: str, perturbation: dict) -> dict:
    atoms = atoms.copy()
    atoms.calc = None
    atoms.info.clear()
    atoms.info.update({
        'config_type': label,
        'candidate_role': 'skeleton_train_candidate_pending_audit',
        'dataset_split': split,
        'source_family': 'Jeitschko_Nowotny_1967_single_crystal',
        'source_doi': '10.1007/BF00899949',
        'source_cif_sha256': SOURCE_CIF_SHA256,
        'geometry_parent_sha256': BASELINE_INPUT_SHA256,
        'perturbation': json.dumps(perturbation, sort_keys=True, separators=(',', ':')),
    })
    path = ROOT / 'inputs' / f'{label}.extxyz'
    write(path, atoms, format='extxyz')
    reread = read(path)
    assert len(reread) == 12 and reread.get_chemical_formula() == 'C4Si2Ti6'
    assert reread.pbc.tolist() == [True, True, True]
    assert not any(k in reread.info for k in ('energy', 'REF_energy', 'PW_PBE_energy_eV'))
    assert not any(k in reread.arrays for k in ('forces', 'REF_forces', 'PW_PBE_forces'))
    return {
        'label': label,
        'input': f'inputs/{path.name}',
        'input_sha256': sha(path.read_bytes()),
        'atom_count': len(reread),
        'formula': reread.get_chemical_formula(),
        'role': role,
        'split': split,
        'perturbation': perturbation,
        'launch_enabled': True,
        'status': 'queued',
    }


source_manifest = json.loads(git_blob(f'{SOURCE_COMMIT}:{SOURCE_MANIFEST_PATH}'))
baseline_label = 'Ti3SiC2_JN1967_primitive_baseline_k10x10x2_mesh_control_v21'
baseline_record = next(r for r in source_manifest['records'] if r['label'] == baseline_label)
baseline_blob_path = f"{SOURCE_COMMIT}:{SOURCE_ROOT}/{baseline_record['input']}"
baseline_blob = git_blob(baseline_blob_path)
BASELINE_INPUT_SHA256 = sha(baseline_blob)
baseline = read(io.StringIO(baseline_blob.decode()), format='extxyz')
assert len(baseline) == 12 and baseline.get_chemical_formula() == 'C4Si2Ti6'
assert np.allclose(baseline.cell.lengths(), [3.068, 3.068, 17.669], atol=1e-8)
assert baseline.pbc.tolist() == [True, True, True]

input_dir = ROOT / 'inputs'
input_dir.mkdir(parents=True, exist_ok=True)
records = []

# Paired local displacements cover the sites/directions absent from the first four labels.
site_jobs = [
    ('Ti2a_z_m0p03_v25', 0, 2, -0.03), ('Ti2a_z_p0p03_v25', 0, 2, 0.03),
    ('Si2b_z_m0p03_v25', 6, 2, -0.03), ('Si2b_z_p0p03_v25', 6, 2, 0.03),
    ('Ti4f_x_m0p03_v25', 3, 0, -0.03), ('Ti4f_x_p0p03_v25', 3, 0, 0.03),
    ('C4f_x_m0p03_v25', 9, 0, -0.03), ('C4f_x_p0p03_v25', 9, 0, 0.03),
]
for label, atom_index, axis, delta in site_jobs:
    atoms = baseline.copy()
    atoms.positions[atom_index, axis] += delta
    records.append(write_input(label, atoms, 'train_candidate_pending_audit', 'train', {
        'kind': 'single_atom_cartesian_displacement', 'atom_index_zero_based': atom_index,
        'element': atoms[atom_index].symbol, 'axis': ['x','y','z'][axis], 'delta_A': delta,
    }))

# Cell strains are kept as training controls in both signs; stress labels will be archived.
for axis_name, axis_index in [('basal_ab', 0), ('c_axis', 2)]:
    for sign, factor, suffix in [(-1, 0.99, 'm1p0'), (1, 1.01, 'p1p0')]:
        atoms = baseline.copy()
        new_cell = atoms.cell.array.copy()
        if axis_name == 'basal_ab':
            new_cell[0:2, :] *= factor
        else:
            new_cell[axis_index, :] *= factor
        atoms.set_cell(new_cell, scale_atoms=True)
        label = f'strain_{axis_name}_{suffix}pct_v25'
        records.append(write_input(label, atoms, 'train_candidate_pending_audit', 'train', {
            'kind': 'homogeneous_cell_strain', 'strain_axis': axis_name,
            'linear_scale': factor, 'nominal_percent': sign * 1.0,
        }))

def random_snapshot(label: str, seed: int, target_rms_A: float, split: str) -> dict:
    atoms = baseline.copy()
    rng = np.random.default_rng(seed)
    delta = rng.normal(size=(len(atoms), 3))
    delta -= delta.mean(axis=0, keepdims=True)
    rms = float(np.sqrt(np.mean(np.sum(delta**2, axis=1))))
    delta *= target_rms_A / rms
    atoms.positions += delta
    return write_input(label, atoms, 'train_candidate_pending_audit' if split == 'train' else 'development_candidate', split, {
        'kind': 'zero_mean_all_atom_random_displacement', 'seed': seed,
        'target_rms_displacement_A': target_rms_A,
        'realized_rms_displacement_A': float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
    })

records.append(random_snapshot('bulk_random_train_seed4101_rms0p02_v25', 4101, 0.02, 'train'))
records.append(random_snapshot('bulk_random_train_seed4102_rms0p04_v25', 4102, 0.04, 'train'))
for seed, amplitude in [(8201, 0.02), (8202, 0.03), (8203, 0.045), (8204, 0.06)]:
    records.append(random_snapshot(f'bulk_random_dev_seed{seed}_rms{str(amplitude).replace(".","p")}_v25', seed, amplitude, 'development'))

manifest = {
    'schema': 'tisic_skeleton_acquisition_v25_v1',
    'task_id': 'Ti3SiC2_TiSiC_skeleton_potential_candidate_v25',
    'status': 'inputs_prepared',
    'launch_enabled': True,
    'fit_enabled': False,
    'data_role': 'non_official_skeleton_train_and_development_candidate',
    'source': {
        'doi': '10.1007/BF00899949',
        'source_family': 'Jeitschko_Nowotny_1967_single_crystal',
        'source_commit': SOURCE_COMMIT,
        'source_manifest_path': SOURCE_MANIFEST_PATH,
        'source_cif_sha256': SOURCE_CIF_SHA256,
        'source_paper_sha256': SOURCE_PAPER_SHA256,
        'baseline_input_path_at_source_commit': baseline_record['input'],
        'baseline_input_sha256': BASELINE_INPUT_SHA256,
        'occupancy_assumption': 'full occupancy; inherited from reconstructed source model',
    },
    'method': {
        'xc': 'PBE', 'basis': 'plane wave', 'cutoff_eV': 500,
        'kpts': [10, 10, 2], 'smearing_eV': 0.10,
        'energy_target': 'finite-width free energy, force-consistent with ASE forces at sigma 0.10 eV',
        'also_record': ['native_extrapolated_energy_eV_cell', 'free_energy_eV_cell', 'forces_eV_A', 'stress_eV_A3'],
        'pbc': [True, True, True], 'poissonsolver': {},
        'mixer': {'type': 'Pulay', 'beta': 0.05, 'nmaxold': 8, 'weight': 100},
        'convergence': {'energy': 1e-5, 'density': 1e-5, 'eigenstates': 1e-5},
        'convergence_note': 'GPAW native criterion keys; each threshold is 1e-5 in its GPAW-defined units.',
        'maxiter': 160, 'GPAW_version': '26.7.0', 'ASE_version': '3.29.0',
        'MPI_ranks': 4, 'threads_per_rank': 1,
        'PAW_data_version': '1.2.1',
        'PAW_setup_sha256': source_manifest['method']['PAW_setup_sha256'],
        'shared_dft_lock': '/workspace/.mace-v19-dft.lock',
        'cgroup_limit_bytes': 34359738368, 'memory_stop_guard_bytes': 32212254720,
        'checkpoint_interval_scf_steps': 20,
    },
    'working_gate_fixed_before_fit': {
        'per_configuration_force_vector_RMSE_max_eV_A': 0.07,
        'note': 'Working fit screen only; the development set is not a blind test set.'
    },
    'previously_verified_training_candidates': [
        {'label': baseline_label, 'source_result_path': f'{SOURCE_ROOT}/calculations/{baseline_label}/{baseline_label}_PW_PBE.extxyz',
         'source_result_commit': SOURCE_COMMIT, 'source_result_role': '10x10x2 fixed-geometry mesh control, re-used as a candidate baseline only in this separate skeleton experiment',
         'split': 'train', 'target_energy_field': 'PW_PBE_free_energy_eV'},
        *[{'label': n, 'source_result_path': f'{SOURCE_ROOT}/calculations/{n}/{n}_PW_PBE.extxyz',
           'source_result_commit': SOURCE_COMMIT, 'source_result_role': 'previous Ti/C displacement candidate; pending source-parent audit',
           'split': 'train', 'target_energy_field': 'PW_PBE_free_energy_eV'} for n in [
            'Ti3SiC2_JN1967_primitive_Ti4f_id0003_z_m002A_k10x10x2_v21',
            'Ti3SiC2_JN1967_primitive_Ti4f_id0003_z_p002A_k10x10x2_v21',
            'Ti3SiC2_JN1967_primitive_C4f_id0009_z_m002A_k10x10x2_v21',
            'Ti3SiC2_JN1967_primitive_C4f_id0009_z_p002A_k10x10x2_v21',
        ]],
    ],
    'preflight_jobs': records,
    'splits': {'train': 19, 'development': 4, 'test': 0,
               'note': 'Counts include the 5 reused historical labels. Random development seeds are held out as whole configurations.'},
    'limitations': [
        'This initial panel is for near-structure Ti3SiC2 bulk response, not defects, surfaces, other phases, or high-temperature MD.',
        'The 0.07 eV/A screen is a per-configuration force-vector RMSE criterion and does not bound each atom maximum error.',
        'No sealed test set exists yet; a passing development screen does not establish general production validity.',
        'Stress is archived for future variable-cell training/audit but the initial MACE recipe may use energy and forces only.',
    ],
}
(ROOT / 'input_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n')
print(json.dumps({'status': manifest['status'], 'registered_new_jobs': len(records),
                  'train_frames_including_reused': manifest['splits']['train'],
                  'development_frames': manifest['splits']['development'],
                  'labels': [r['label'] for r in records]}, indent=2))
