#!/usr/bin/env python3
"""Build six C-orbit response inputs from registered training-interface parents."""
from __future__ import annotations

import hashlib
import io
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from ase.io import read, write

from geometry_utils import compare_geometry, extxyz_geometry_frames

ROOT = Path(__file__).resolve().parent
REPO = Path('/workspace/jeitschko-main-run')
V22_COMMIT = '5d1e734c16d17dc9b460e3184eab4670ccc2739c'
TRAIN55_COMMIT = '6b8fe5520679808cd98506eb08b8ea72cd71f61a'
V22_TASK = 'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v22_force_response_acquisition'
TRAIN55_TASK = 'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_interface_acquisition/training_candidate_v19_train55_20261101_A'
TRAIN55_DATA = REPO / TRAIN55_TASK / 'dataset'
INTERFACES = ['AgTi', 'AgSi', 'AgC']
SIGNS = [('m', -0.02), ('p', 0.02)]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_blob(commit: str, path: str) -> bytes:
    return subprocess.check_output(['git', '-C', str(REPO), 'show', f'{commit}:{path}'])


def label_for(interface: str, sign: str) -> str:
    return f'{interface}_JN1967_registryB_surfaceC_dz_{sign}0p02_sigma0p10_k6x6x1_v24'


def main() -> None:
    inputs = ROOT / 'inputs'
    if inputs.exists() and any(inputs.iterdir()):
        raise SystemExit(f'Refusing to overwrite existing input files: {inputs}')
    train_manifest = json.loads((TRAIN55_DATA / 'manifest.json').read_text())
    train_raw = (TRAIN55_DATA / 'train.extxyz').read_bytes()
    dev_raw = (TRAIN55_DATA / 'valid.extxyz').read_bytes()
    if digest(train_raw) != train_manifest['train_sha256'] or digest(dev_raw) != train_manifest['dev_sha256']:
        raise SystemExit('Pinned train55/development dataset hash mismatch')
    if subprocess.run(['git', '-C', str(REPO), 'merge-base', '--is-ancestor', TRAIN55_COMMIT, 'HEAD']).returncode != 0:
        raise SystemExit('Pinned train55 commit is not an ancestor of current HEAD')

    v22_manifest_raw = git_blob(V22_COMMIT, f'{V22_TASK}/input_manifest.json')
    v22_manifest = json.loads(v22_manifest_raw)
    train_frames = extxyz_geometry_frames(TRAIN55_DATA / 'train.extxyz')
    dev_frames = extxyz_geometry_frames(TRAIN55_DATA / 'valid.extxyz')
    if len(train_frames) != 55 or len(dev_frames) != 3:
        raise SystemExit(f'Geometry-only screen expected 55+3 frames; got {len(train_frames)}+{len(dev_frames)}')

    inputs.mkdir(parents=True, exist_ok=True)
    records, jobs = [], []
    for interface_index, interface in enumerate(INTERFACES):
        parent_row = next(row for row in v22_manifest['records'] if row['interface'] == interface)
        parent_path = parent_row['parent_input']
        parent_raw = git_blob(V22_COMMIT, parent_path)
        if digest(parent_raw) != parent_row['parent_input_sha256']:
            raise SystemExit(f'Pinned parent input hash mismatch: {interface}')
        parent = read(io.StringIO(parent_raw.decode()), format='extxyz')
        if len(parent) != 52 or parent.get_chemical_formula() != 'C16Ag4Si8Ti24':
            raise SystemExit(f'Unexpected parent geometry: {interface}')
        if parent.pbc.tolist() != [True, True, False]:
            raise SystemExit(f'Unexpected PBC: {interface}')
        if any(k in parent.info for k in ('energy', 'REF_energy', 'PW_PBE_energy_eV')):
            raise SystemExit(f'Parent geometry unexpectedly contains energy labels: {interface}')
        if any(k in parent.arrays for k in ('forces', 'REF_forces', 'PW_PBE_forces')):
            raise SystemExit(f'Parent geometry unexpectedly contains force labels: {interface}')

        symbols = np.asarray(parent.get_chemical_symbols())
        c_indices = np.flatnonzero(symbols == 'C')
        center = np.flatnonzero(np.asarray(parent.arrays.get('central_pair', np.zeros(len(parent))), dtype=bool))
        if len(center) != 2:
            raise SystemExit(f'Missing registered central contact pair: {interface}')
        z_levels = sorted(set(round(float(parent.positions[i, 2]), 4) for i in c_indices))
        c_orbits = []
        for z in z_levels:
            orbit = c_indices[np.isclose(parent.positions[c_indices, 2], z, rtol=0, atol=1e-4)]
            if len(orbit) != 4:
                continue
            dists = [[parent.get_distance(int(i), int(j), mic=True) for j in center] for i in orbit]
            score = float(np.mean(np.min(np.asarray(dists), axis=1)))
            c_orbits.append((score, float(z), orbit))
        if not c_orbits:
            raise SystemExit(f'No four-atom C orbit found: {interface}')
        score, zref, target = min(c_orbits, key=lambda row: (row[0], row[1]))
        target_ids = [int(x) for x in parent.arrays.get('lammps_id', np.arange(len(parent)))[target]]

        parent_hits = [idx for idx, frame in enumerate(train_frames[:49])
                       if (compare_geometry(parent, frame) or {}).get('exact_duplicate')]
        if len(parent_hits) != 1:
            raise SystemExit(f'Expected one registered parent in train55 train49-prefix, got {parent_hits} for {interface}')
        parent_frame = parent_hits[0]
        expected_response_frames = {49 + 2 * interface_index, 50 + 2 * interface_index}
        allowed_near = {parent_frame} | expected_response_frames

        panel_geometry = []
        for sign, displacement in SIGNS:
            label = label_for(interface, sign)
            atoms = parent.copy()
            atoms.positions[target, 2] += displacement
            atoms.info.update({
                'config_type': label,
                'geometry_label': f'{interface}_registryB_contact_C_orbit_z{displacement:+.2f}A',
                'dataset_role': 'train_acquisition_candidate_pending_audit',
                'parent_geometry_sha256': digest(parent_raw),
                'perturbed_atom_ids': target_ids,
                'perturbation': 'four C atoms in the contact-nearest symmetry orbit displaced together along surface normal',
            })
            output = inputs / f'{label}.extxyz'
            write(output, atoms, format='extxyz')
            input_raw = output.read_bytes()
            frames = extxyz_geometry_frames(output)
            if len(frames) != 1:
                raise SystemExit(f'Expected one geometry frame: {label}')

            train_comparisons = []
            for idx, candidate in enumerate(train_frames):
                comp = compare_geometry(frames[0], candidate)
                if comp is not None:
                    train_comparisons.append({'frame_index': idx, **comp})
            exact_train = [x for x in train_comparisons if x['exact_duplicate']]
            near_train = [x for x in train_comparisons if x['near_duplicate']]
            unexpected_near = [x for x in near_train if x['frame_index'] not in allowed_near]
            if exact_train or unexpected_near:
                raise SystemExit(f'Unexpected train55 hit for {label}: exact={exact_train}, near={unexpected_near}')

            dev_comparisons = []
            for idx, candidate in enumerate(dev_frames):
                comp = compare_geometry(frames[0], candidate)
                if comp is not None:
                    dev_comparisons.append({'frame_index': idx, **comp})
            exact_dev = [x for x in dev_comparisons if x['exact_duplicate']]
            near_dev = [x for x in dev_comparisons if x['near_duplicate']]
            if exact_dev or near_dev:
                raise SystemExit(f'Exact/near geometry hit against frozen dev3 for {label}')

            distances = frames[0].get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            minimum = float(distances.min())
            if minimum < 1.8:
                raise SystemExit(f'Unphysically short pair distance {minimum:.4f} Å: {label}')
            contact = float(atoms.get_distance(int(center[0]), int(center[1]), mic=True))
            panel_geometry.append(frames[0])
            record = {
                'label': label, 'interface': interface, 'input': f'inputs/{output.name}',
                'geometry_label': atoms.info['geometry_label'],
                'input_sha256': digest(input_raw), 'parent_input': parent_path,
                'parent_input_sha256': digest(parent_raw), 'parent_source_commit': V22_COMMIT,
                'parent_source_manifest_sha256': digest(v22_manifest_raw),
                'atoms': len(frames[0]), 'formula': frames[0].get_chemical_formula(),
                'cell_A': frames[0].cell.array.tolist(), 'pbc': frames[0].pbc.tolist(),
                'target_element': 'C', 'target_atom_ids': target_ids,
                'target_original_z_A': zref, 'target_contact_nearest_mean_distance_A': score,
                'target_displacement_z_A': displacement,
                'intentional_train55_parent_frame_index': parent_frame,
                'intentional_train55_v23_response_frame_indices': sorted(expected_response_frames),
                'intentional_train55_near_hit_frame_indices': [x['frame_index'] for x in near_train],
                'marked_contact_pair_distance_A': contact,
                'minimum_periodic_pair_distance_A': minimum,
                'role': 'train_acquisition_candidate_pending_audit', 'launch_enabled': True,
                'geometry_screens': {
                    'train55_exact_hits': len(exact_train), 'train55_near_hits': len(near_train),
                    'train55_unexpected_near_hits': len(unexpected_near),
                    'frozen_dev3_exact_hits': len(exact_dev), 'frozen_dev3_near_hits': len(near_dev),
                    'train55_comparisons': train_comparisons, 'frozen_dev3_comparisons': dev_comparisons,
                },
            }
            records.append(record)
            jobs.append({
                'label': label, 'geometry_label': record['geometry_label'],
                'input': record['input'], 'input_sha256': record['input_sha256'],
                'kpts': [6, 6, 1], 'smearing_eV': 0.10,
                'role': record['role'], 'launch_enabled': True,
            })
        paired = compare_geometry(panel_geometry[0], panel_geometry[1])
        if paired is None or not paired['near_duplicate'] or paired['exact_duplicate']:
            raise SystemExit(f'Paired ± displacement geometries did not screen as intended: {interface}')

    jobs.sort(key=lambda row: (INTERFACES.index(row['label'].split('_', 1)[0]),
                               0 if '_m0p02_' in row['label'] else 1))
    manifest = {
        'task_id': 'pbe_interface_v24_carbon_force_response_acquisition',
        'status': 'V24_GEOMETRY_SCREENED_ACQUISITION_CANDIDATES',
        'source': {
            'source_family': 'Jeitschko_Nowotny_1967_single_crystal',
            'source_doi': '10.1007/BF00899949', 'parent_source_commit': V22_COMMIT,
            'train55_screen_commit': TRAIN55_COMMIT,
            'train55_manifest_sha256': digest((TRAIN55_DATA / 'manifest.json').read_bytes()),
        },
        'reference_datasets': {
            'train55_frames': 55, 'train55_train_sha256': train_manifest['train_sha256'],
            'train55_manifest_sha256': digest((TRAIN55_DATA / 'manifest.json').read_bytes()),
            'frozen_dev3_frames': 3, 'frozen_dev3_sha256': train_manifest['dev_sha256'],
        },
        'method_candidate': v22_manifest['method_candidate'],
        'runtime_guard': {'cgroup_limit_bytes': 34359738368, 'stop_threshold_bytes': 32212254720,
                          'prelaunch_memory_limit_bytes': 3221225472, 'minimum_free_disk_bytes': 5368709120,
                          'maximum_concurrent_jobs': 1, 'automatic_retries': 0},
        'role_policy': {'official_train36_dev3_mutated': False, 'test_or_sealed_data_opened': False,
                        'dataset_integration_enabled': False},
        'geometry_perturbation_policy': 'On each registered AgTi, AgSi and AgC training-interface parent, displace the four C atoms in the symmetry orbit nearest the registered central contact pair by ±0.02 Å along the surface normal. Geometry generation reads only parent geometry; frozen dev3 is used only for geometry-only duplicate screening.',
        'geometry_screens': {
            'screened_against': 'train55 and frozen dev3 using geometry-only parser',
            'train55_expected_near_hit_types': ['registered parent', 'same-interface v23 topmost-Ti response pair'],
            'train55_unexpected_near_hits': sum(r['geometry_screens']['train55_unexpected_near_hits'] for r in records),
            'train55_exact_hits': sum(r['geometry_screens']['train55_exact_hits'] for r in records),
            'frozen_dev3_exact_hits': sum(r['geometry_screens']['frozen_dev3_exact_hits'] for r in records),
            'frozen_dev3_near_hits': sum(r['geometry_screens']['frozen_dev3_near_hits'] for r in records),
            'within_panel_paired_displacements_are_intentional': True,
        },
        'records': records, 'preflight_jobs': jobs,
        'launch_enabled': True, 'training_enabled': False,
        'dataset_integration_enabled': False, 'test_or_sealed_data_opened': False,
        'labels_or_forces_read': False,
        'limitations': [
            'All six configurations share one crystallographic source family and registry-B slab construction.',
            'The panel samples one contact-nearest C orbit at fixed geometry and does not establish broad interface coverage.',
            'Native extrapolated energy and finite-smearing forces follow the existing convention; exact energy-force gradient consistency is not claimed.',
            'Slab k-point, smearing and vacuum convergence are not fully established.',
        ],
    }
    (ROOT / 'input_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + '\n')
    audit = {
        'status': 'PASS', 'source_v22_commit': V22_COMMIT,
        'train55_screen_commit': TRAIN55_COMMIT,
        'train55_train_sha256': train_manifest['train_sha256'],
        'frozen_dev3_sha256': train_manifest['dev_sha256'],
        'screen_parser': 'geometry_utils.py; symbols, positions, cell, PBC and persistent IDs only',
        'candidate_count': len(records),
        'candidates': [{
            'label': r['label'], 'input_sha256': r['input_sha256'],
            'target_atom_ids': r['target_atom_ids'], 'displacement_z_A': r['target_displacement_z_A'],
            'minimum_distance_A': r['minimum_periodic_pair_distance_A'],
            'marked_contact_pair_distance_A': r['marked_contact_pair_distance_A'],
            'geometry_screens': r['geometry_screens'],
        } for r in records],
        'test_or_sealed_data_opened': False,
    }
    (ROOT / 'geometry_audit.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False, sort_keys=True) + '\n')
    sequence = {'status': 'not_started', 'created_utc': datetime.now(timezone.utc).isoformat(),
                'completed_labels': [], 'current_label': None, 'automatic_retries': 0}
    (ROOT / 'sequence_status.json').write_text(json.dumps(sequence, indent=2) + '\n')
    print(json.dumps({'status': manifest['status'], 'labels': [j['label'] for j in jobs],
                      'train55_unexpected_near_hits': manifest['geometry_screens']['train55_unexpected_near_hits'],
                      'frozen_dev3_exact_hits': manifest['geometry_screens']['frozen_dev3_exact_hits'],
                      'frozen_dev3_near_hits': manifest['geometry_screens']['frozen_dev3_near_hits'],
                      'official_split_mutated': False, 'withheld_data_opened': False}, indent=2))


if __name__ == '__main__':
    main()
