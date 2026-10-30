#!/usr/bin/env python3
"""Build six symmetry-preserving Ti surface displacement inputs from registered training geometry."""
from __future__ import annotations
import hashlib, io, json, subprocess
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from ase.io import read, write
from geometry_utils import compare_geometry, extxyz_geometry_frames

ROOT = Path(__file__).resolve().parent
REPO = Path('/workspace/jeitschko-main-run')
V22_COMMIT = '5d1e734c16d17dc9b460e3184eab4670ccc2739c'
TRAIN49_COMMIT = 'e1a6c134b240a8afa9b44febc0bb646a4015ca3b'
V22_TASK = 'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v22_force_response_acquisition'
TRAIN49_TASK = 'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_interface_acquisition/training_candidate_v19_train49_20261030_A'
TRAIN49_DATA = REPO / TRAIN49_TASK / 'dataset'
INTERFACES = ['AgTi', 'AgSi', 'AgC']
SIGNS = [('m', -0.02), ('p', 0.02)]

def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def git_blob(commit: str, path: str) -> bytes:
    return subprocess.check_output(['git', '-C', str(REPO), 'show', f'{commit}:{path}'])

def label_for(interface: str, sign: str) -> str:
    return f'{interface}_JN1967_registryB_surfaceTi_dz_{sign}0p02_sigma0p10_k6x6x1_v23'

def main() -> None:
    inputs = ROOT / 'inputs'
    if inputs.exists() and any(inputs.iterdir()):
        raise SystemExit(f'Refusing to overwrite existing input files: {inputs}')
    train_manifest = json.loads((TRAIN49_DATA / 'manifest.json').read_text())
    train_raw = (TRAIN49_DATA / 'train.extxyz').read_bytes()
    dev_raw = (TRAIN49_DATA / 'valid.extxyz').read_bytes()
    if digest(train_raw) != train_manifest['train_sha256'] or digest(dev_raw) != train_manifest['dev_sha256']:
        raise SystemExit('Pinned train49/development dataset hash mismatch')
    if subprocess.run(['git', '-C', str(REPO), 'merge-base', '--is-ancestor', TRAIN49_COMMIT, 'HEAD']).returncode != 0:
        raise SystemExit('Pinned train49 data commit is not an ancestor of current HEAD')
    v22_manifest_raw = git_blob(V22_COMMIT, f'{V22_TASK}/input_manifest.json')
    v22_manifest = json.loads(v22_manifest_raw)
    if v22_manifest.get('status') != 'V22_GEOMETRY_SCREENED_ACQUISITION_CANDIDATES':
        if not v22_manifest.get('preflight_jobs'):
            raise SystemExit('Pinned v22 parent manifest is malformed')
    train_frames = extxyz_geometry_frames(TRAIN49_DATA / 'train.extxyz')
    dev_frames = extxyz_geometry_frames(TRAIN49_DATA / 'valid.extxyz')
    if len(train_frames) != 49 or len(dev_frames) != 3:
        raise SystemExit(f'Geometry-only screen expected 49+3 frames; got {len(train_frames)}+{len(dev_frames)}')
    screened_frames = train_frames + dev_frames

    inputs.mkdir(parents=True, exist_ok=True)
    records, jobs = [], []
    for interface in INTERFACES:
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
        ti = np.flatnonzero(symbols == 'Ti')
        zmax = float(parent.positions[ti, 2].max())
        target = ti[np.isclose(parent.positions[ti, 2], zmax, rtol=0, atol=1e-6)]
        if len(target) != 4:
            raise SystemExit(f'Expected four symmetry-equivalent topmost Ti atoms, got {len(target)} for {interface}')
        parent_matches = []
        for idx, frame in enumerate(train_frames):
            comp = compare_geometry(parent, frame)
            if comp is not None and comp['exact_duplicate']:
                parent_matches.append(idx)
        if len(parent_matches) != 1:
            raise SystemExit(f'Expected one exact geometry-only parent frame in train49, got {parent_matches} for {interface}')
        parent_train_frame = parent_matches[0]
        target_ids = [int(x) for x in parent.arrays.get('lammps_id', np.arange(len(parent)))[target]]

        for sign, displacement in SIGNS:
            label = label_for(interface, sign)
            atoms = parent.copy()
            atoms.positions[target, 2] += displacement
            atoms.info.update({
                'config_type': label,
                'geometry_label': f'{interface}_registryB_topmost_Ti_orbit_z{displacement:+.2f}A',
                'dataset_role': 'train_acquisition_candidate_pending_audit',
                'parent_geometry_sha256': digest(parent_raw),
                'perturbed_atom_ids': target_ids,
                'perturbation': 'all four symmetry-equivalent topmost Ti atoms displaced together along surface normal',
            })
            output = inputs / f'{label}.extxyz'
            write(output, atoms, format='extxyz')
            input_raw = output.read_bytes()
            input_hash = digest(input_raw)
            frames = extxyz_geometry_frames(output)
            if len(frames) != 1:
                raise SystemExit(f'Expected one geometry frame: {label}')
            comparisons = []
            for category, dataset in [('train49_geometry_only', train_frames), ('frozen_dev3_geometry_only', dev_frames)]:
                category_comparisons = []
                for idx, candidate in enumerate(dataset):
                    comp = compare_geometry(frames[0], candidate)
                    if comp is not None:
                        category_comparisons.append({'frame_index': idx, **comp})
                exact = [x for x in category_comparisons if x['exact_duplicate']]
                near = [x for x in category_comparisons if x['near_duplicate']]
                intentional_parent_near = category == 'train49_geometry_only' and [x['frame_index'] for x in near] == [parent_train_frame] and not exact
                if category == 'train49_geometry_only':
                    if not intentional_parent_near:
                        raise SystemExit(f'Unexpected train49 geometry hits: {label}, exact={exact}, near={near}')
                elif exact or near:
                    raise SystemExit(f'Exact/near geometry hit against {category}: {label}')
                comparisons.append({'dataset': category, 'frames_read_geometry_only': len(dataset),
                                    'same_composition_normalized_comparisons': len(category_comparisons),
                                    'exact_hits': len(exact), 'near_hits': len(near),
                                    'intentional_parent_near_hit_frame_index': parent_train_frame if intentional_parent_near else None,
                                    'comparisons': category_comparisons})
            distances = frames[0].get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            minimum = float(distances.min())
            if minimum < 1.8:
                raise SystemExit(f'Unphysically short pair distance {minimum:.4f} Å: {label}')
            # The geometry-only parser intentionally drops non-geometric arrays.
            # Read the registered contact marker from the original ASE candidate,
            # while keeping all train/dev duplicate checks on geometry-only frames.
            central = np.flatnonzero(np.asarray(atoms.arrays.get('central_pair', np.zeros(len(atoms))), dtype=bool))
            if len(central) != 2:
                raise SystemExit(f'Missing registered central contact pair: {label}')
            contact = float(atoms.get_distance(int(central[0]), int(central[1]), mic=True))
            record = {
                'label': label, 'interface': interface, 'input': f'inputs/{output.name}',
                'input_sha256': input_hash, 'parent_input': parent_path,
                'parent_input_sha256': digest(parent_raw), 'parent_source_commit': V22_COMMIT,
                'parent_source_manifest_sha256': digest(v22_manifest_raw),
                'atoms': len(frames[0]), 'formula': frames[0].get_chemical_formula(),
                'cell_A': frames[0].cell.array.tolist(), 'pbc': frames[0].pbc.tolist(),
                'target_element': 'Ti', 'target_atom_ids': target_ids,
                'target_original_z_A': zmax, 'target_displacement_z_A': displacement,
                'intentional_parent_train49_geometry_frame_index': parent_train_frame,
                'marked_contact_pair_distance_A': contact,
                'minimum_periodic_pair_distance_A': minimum,
                'role': 'train_acquisition_candidate_pending_audit',
                'launch_enabled': True, 'geometry_screens': comparisons,
            }
            records.append(record)
            jobs.append({'label': label, 'geometry_label': record['geometry_label'] if 'geometry_label' in record else f'{interface}_surfaceTi_{displacement:+.2f}',
                         'input': record['input'], 'input_sha256': input_hash,
                         'kpts': [6, 6, 1], 'smearing_eV': 0.10,
                         'role': record['role'], 'launch_enabled': True})

    jobs.sort(key=lambda x: (INTERFACES.index(x['label'].split('_', 1)[0]),
                             0 if '_m0p02_' in x['label'] else 1))
    parent_manifest = json.loads(git_blob(V22_COMMIT, f'{V22_TASK}/input_manifest.json'))
    manifest = {
        'task_id': 'pbe_interface_v23_ti_substrate_force_response_acquisition',
        'status': 'V23_GEOMETRY_SCREENED_ACQUISITION_CANDIDATES',
        'source': {'source_family': 'Jeitschko_Nowotny_1967_single_crystal',
                   'source_doi': '10.1007/BF00899949',
                   'parent_inputs_from': [{'interface': interface,
                                           'path': next(row for row in v22_manifest['records'] if row['interface'] == interface)['parent_input'],
                                           'sha256': next(row for row in v22_manifest['records'] if row['interface'] == interface)['parent_input_sha256']}
                                          for interface in INTERFACES],
                   'v22_source_commit': V22_COMMIT,
                   'train49_geometry_screen_commit': TRAIN49_COMMIT},
        'reference_datasets': {
            'train49_manifest_sha256': digest((TRAIN49_DATA / 'manifest.json').read_bytes()),
            'train49_train_sha256': train_manifest['train_sha256'],
            'frozen_dev3_sha256': train_manifest['dev_sha256'],
            'train49_frames': 49, 'frozen_dev3_frames': 3,
        },
        'method_candidate': parent_manifest['method_candidate'],
        'runtime_guard': {'cgroup_limit_bytes': 34359738368, 'stop_threshold_bytes': 32212254720,
                          'prelaunch_memory_limit_bytes': 3221225472, 'minimum_free_disk_bytes': 5368709120,
                          'maximum_concurrent_jobs': 1, 'automatic_retries': 0},
        'role_policy': {'official_train36_dev3_mutated': False, 'test_or_sealed_data_opened': False,
                        'owner_side_external_role_review': 'waived by user for experiment-only work',
                        'dataset_integration_enabled': False},
        'geometry_perturbation_policy': 'Use three registered neutral-gap registry-B training geometries and displace the four symmetry-equivalent topmost Ti atoms together by ±0.02 Å along surface normal. Geometry generation and screening read only atomic identities, positions, cell and PBC; no dev3 coordinates or labels were copied into inputs.',
        'geometry_screens': {'screened_against': 'train49 and frozen dev3 using geometry-only parser',
                             'train49_exact_hits': 0, 'train49_intentional_parent_near_hits': 6,
                             'train49_unexpected_near_hits': 0, 'frozen_dev3_exact_hits': 0,
                             'frozen_dev3_near_hits': 0,
                             'within_panel_paired_displacements_are_intentional': True},
        'records': records, 'preflight_jobs': jobs,
        'launch_enabled': True, 'training_enabled': False,
        'dataset_integration_enabled': False, 'test_or_sealed_data_opened': False,
        'labels_or_forces_read': False,
        'limitations': [
            'All six candidate configurations derive from one 1967 source-family geometry and one registry-B slab family.',
            'Paired local Ti displacements sample a narrow force-response mode; they do not establish broad interface coverage.',
            'Native extrapolated energy and finite-smearing forces follow the existing convention; exact energy-force gradient consistency is not claimed.',
            'Slab k-point, smearing and vacuum convergence are not fully established.',
        ],
    }
    (ROOT / 'input_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + '\n')
    audit = {
        'status': 'PASS', 'source_v22_commit': V22_COMMIT,
        'train49_data_commit': TRAIN49_COMMIT,
        'train49_train_sha256': train_manifest['train_sha256'],
        'frozen_dev3_sha256': train_manifest['dev_sha256'],
        'screen_parser': 'geometry_utils.py; symbols, positions, cell, PBC and persistent IDs only',
        'candidate_count': len(records),
        'candidates': [{'label': r['label'], 'input_sha256': r['input_sha256'],
                        'target_atom_ids': r['target_atom_ids'],
                        'displacement_z_A': r['target_displacement_z_A'],
                        'minimum_distance_A': r['minimum_periodic_pair_distance_A'],
                        'marked_contact_pair_distance_A': r['marked_contact_pair_distance_A'],
                        'geometry_screens': r['geometry_screens']} for r in records],
        'test_or_sealed_data_opened': False,
    }
    (ROOT / 'geometry_audit.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False, sort_keys=True) + '\n')
    sequence = {'status': 'not_started', 'created_utc': datetime.now(timezone.utc).isoformat(),
                'completed_labels': [], 'current_label': None, 'automatic_retries': 0}
    (ROOT / 'sequence_status.json').write_text(json.dumps(sequence, indent=2) + '\n')
    print(json.dumps({'status': manifest['status'], 'candidate_count': len(records),
                      'labels': [r['label'] for r in records], 'no_exact_hits_or_unexpected_near_hits': True,
                      'intentional_parent_near_hits': 6, 'frozen_dev3_near_hits': 0,
                      'train49_sha256': train_manifest['train_sha256'],
                      'frozen_dev3_sha256': train_manifest['dev_sha256']}, indent=2))

if __name__ == '__main__':
    main()
