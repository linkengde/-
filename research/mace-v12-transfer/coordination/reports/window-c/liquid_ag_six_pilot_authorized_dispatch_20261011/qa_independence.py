#!/usr/bin/env python3
"""Compare only six explicit liquid-pilot candidate manifests and geometries.

Unique RNG seeds support trajectory sampling independence, never independent
validation-family eligibility. No directories are enumerated for structures.
No calculators, DFT labels, dynamics, or model inference are used.
"""
from __future__ import annotations
import argparse
import csv
import itertools
import json
from pathlib import Path
import sys
from qa_sampling import Hold, sha256_file, canonical_sha, geometry_identity, save_json, save_csv


def periodic_assignment(first, second, edge, np, linear_sum_assignment):
    """Translation and permutation invariant MIC assignment RMS distance.

    Try all atom alignments to atom zero. This detects exact/near duplicates;
    it is a conservative screen, not a complete search over rotations or a
    rigorous optimal periodic registration. Rotational equivalence needs C.
    """
    best = float('inf')
    best_maximum = None
    for anchor in second:
        shift = anchor - first[0]
        vectors = first[:, None, :] + shift - second[None, :, :]
        vectors -= np.rint(vectors / edge) * edge
        squared = (vectors * vectors).sum(axis=2)
        rows, columns = linear_sum_assignment(squared)
        matched = np.sqrt(squared[rows, columns])
        rms = float(np.sqrt((matched * matched).mean()))
        if rms < best:
            best, best_maximum = rms, float(matched.max())
    return best, best_maximum


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate-manifests', type=Path, nargs=6, required=True,
                        help='Exactly six explicit, public, newly generated pilot manifests')
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True, help='A new comparison folder')
    args = parser.parse_args(argv)
    if args.out.exists():
        parser.error('Output folder exists; never overwrite an earlier review')
    args.out.mkdir(parents=True, exist_ok=False)
    summary = {'status': 'HOLD_FOR_C_REVIEW', 'INDEPENDENCE_SCREEN_PASS': False,
               'manualC_REQUIRED': True, 'STRUCTURE_PASS': False,
               'independent_test_eligible': False, 'DFT_LABELS': 'ABSENT',
               'reasons': [], 'warnings': ['Seeds are not independent test families',
                    'Registration checks periodic translation and permutation, not all rotations',
                    'All common source ancestors remain grouped for future data splitting']}
    try:
        if not args.config.is_file():
            raise Hold('RELEASE_CONFIG_MISSING')
        config = json.loads(args.config.read_text())
        if config.get('source_release') is not True or config.get('density_reference_pass') is not True:
            raise Hold('SOURCE_OR_DENSITY_RELEASE_HOLD')
        tasks_path = args.config.parent / config.get('tasks_file', 'AUTHORIZED_SIX_TASKS.csv')
        if not tasks_path.is_file() or sha256_file(tasks_path) != config.get('tasks_sha256'):
            raise Hold('TASK_MATRIX_MISSING_OR_SHA_MISMATCH')
        with tasks_path.open(newline='') as stream:
            tasks = list(csv.DictReader(stream))
        if len(tasks) != 6 or len({task['task_id'] for task in tasks}) != 6:
            raise Hold('EXPECTED_EXACTLY_SIX_UNIQUE_PILOT_TASKS')
        task_map = {task['task_id']: task for task in tasks}
        summary['config_sha256'] = sha256_file(args.config)
        import numpy as np
        from scipy.optimize import linear_sum_assignment
        from ase.io import read
        candidates = []
        for path in args.candidate_manifests:
            if not path.is_file():
                raise Hold('MISSING_ACTUAL_CANDIDATE_MANIFEST:' + str(path))
            manifest = json.loads(path.read_text())
            task_id = manifest.get('task_id')
            if task_id not in task_map:
                raise Hold('CANDIDATE_NOT_IN_ORIGINAL_SIX_TASKS')
            task = task_map[task_id]
            if manifest.get('owner') != task['owner'] or manifest.get('RNG_seed') != int(task['RNG_seed']):
                raise Hold('OWNER_OR_SEED_CHANGED:' + task_id)
            if manifest.get('trajectory_parent') != task['trajectory_parent']:
                raise Hold('PARENT_CHANGED:' + task_id)
            expected_root = config.get('lineage', {}).get('source_root_by_task', {}).get(task_id, task.get('source_root'))
            if not expected_root or 'TO_BE_FROZEN' in expected_root or manifest.get('source_root') != expected_root:
                raise Hold('SOURCE_ROOT_UNFROZEN_OR_CHANGED:' + task_id)
            if manifest.get('AUTO_DIAGNOSTIC_PASS') is not True or manifest.get('manualC_REQUIRED') is not True:
                raise Hold('CANDIDATE_HAS_NO_VALID_LOCAL_DIAGNOSTIC_REVIEW:' + task_id)
            relative = Path(manifest.get('candidate_relative_path', ''))
            if relative.is_absolute() or '..' in relative.parts or not relative.name:
                raise Hold('INVALID_CANDIDATE_RELATIVE_PATH')
            candidate_path = (path.parent / relative).resolve()
            if path.parent.resolve() not in candidate_path.parents or not candidate_path.is_file():
                raise Hold('CANDIDATE_PATH_ESCAPES_QA_FOLDER_OR_MISSING')
            if sha256_file(candidate_path) != manifest.get('candidate_file_sha256'):
                raise Hold('CANDIDATE_FILE_SHA_MISMATCH:' + task_id)
            atoms = read(candidate_path)
            if atoms.calc is not None:
                raise Hold('CANDIDATE_CONTAINS_CALCULATOR_LABELS:' + task_id)
            if (len(atoms) != 64 or not np.all(atoms.numbers == 47) or not np.all(atoms.pbc) or
                    not np.allclose(atoms.cell.array, np.eye(3) * float(task['box_edge_A']), atol=1e-7, rtol=0)):
                raise Hold('CANDIDATE_ELEMENT_CELL_OR_PBC_MISMATCH:' + task_id)
            identity = geometry_identity(atoms, np)
            if identity['atom_ids'] != list(range(1, 65)):
                raise Hold('ATOM_ORDER_CHANGED:' + task_id)
            if canonical_sha(identity) != manifest.get('geometry_sha256'):
                raise Hold('SERIALIZED_GEOMETRY_SHA_MISMATCH:' + task_id)
            for key in ('source_trajectory_sha256', 'source_manifest_sha256', 'source_config_sha256'):
                value = manifest.get(key, '')
                if len(value) != 64 or any(character not in '0123456789abcdef' for character in value):
                    raise Hold('SOURCE_SHA_MISSING:' + key)
            qa_report = path.parent / 'qa_report.json'
            if not qa_report.is_file():
                raise Hold('QA_REPORT_MISSING')
            qa = json.loads(qa_report.read_text())
            if (qa.get('AUTO_DIAGNOSTIC_PASS') is not True or qa.get('task_id') != task_id or
                    qa.get('input_sha256', {}).get('trajectory.traj') != manifest['source_trajectory_sha256'] or
                    qa.get('candidate_file_sha256') != manifest['candidate_file_sha256']):
                raise Hold('CANDIDATE_AND_QA_REPORT_MISMATCH:' + task_id)
            candidates.append((manifest, atoms))
        if len({manifest['task_id'] for manifest, _ in candidates}) != 6:
            raise Hold('TASK_ID_DUPLICATED_OR_MISSING')
        if len({manifest['RNG_seed'] for manifest, _ in candidates}) != 6:
            raise Hold('RNG_SEED_DUPLICATED')
        if len({manifest['trajectory_parent'] for manifest, _ in candidates}) != 6:
            raise Hold('INDEPENDENT_TRAJECTORY_PARENT_DUPLICATED')
        if len({manifest['source_trajectory_sha256'] for manifest, _ in candidates}) != 6:
            raise Hold('FULL_TRAJECTORY_SHA_DUPLICATED; possible shared trajectory')
        qa_config = config.get('qa', {})
        limit = float(qa_config.get('near_duplicate_RMS_A', 0.10))
        if not np.isfinite(limit) or not 0 < limit < 0.5:
            raise Hold('INVALID_NEAR_DUPLICATE_THRESHOLD')
        pairs = []
        for (first, atoms_first), (second, atoms_second) in itertools.combinations(candidates, 2):
            edge_first = float(atoms_first.cell[0, 0])
            edge_second = float(atoms_second.cell[0, 0])
            # Fractional comparison detects a reused configuration simply scaled
            # to another density. Report equivalent A-distance at mean box size.
            fractional_first = atoms_first.get_scaled_positions(wrap=True)
            fractional_second = atoms_second.get_scaled_positions(wrap=True)
            fractional_rms, fractional_maximum = periodic_assignment(fractional_first, fractional_second,
                                                                      1.0, np, linear_sum_assignment)
            mean_edge = (edge_first + edge_second) / 2
            scaled_rms = fractional_rms * mean_edge
            same_density = abs(edge_first - edge_second) < 1e-7
            exact = first['geometry_sha256'] == second['geometry_sha256']
            duplicate = exact or scaled_rms < limit
            pairs.append({'task_id_1': first['task_id'], 'task_id_2': second['task_id'],
                          'same_density': same_density, 'exact_geometry_sha_equal': exact,
                          'fractional_assignment_RMS': fractional_rms,
                          'equivalent_assignment_RMS_A': scaled_rms,
                          'equivalent_assignment_max_A': fractional_maximum * mean_edge,
                          'near_duplicate_threshold_A': limit, 'duplicate_screen': duplicate,
                          'same_source_root': first['source_root'] == second['source_root']})
        save_csv(args.out / 'pairwise_duplicate_screen.csv', pairs, list(pairs[0]))
        summary['pair_count'] = len(pairs)
        summary['candidate_manifest_sha256'] = {str(path): sha256_file(path) for path in args.candidate_manifests}
        summary['family_groups'] = {}
        for manifest, _ in candidates:
            summary['family_groups'].setdefault(manifest['source_root'], []).append(manifest['task_id'])
        summary['near_duplicate_RMS_A'] = limit
        summary['threshold_source'] = 'C_config' if 'near_duplicate_RMS_A' in qa_config else 'provisional_screen_requires_C'
        duplicates = [pair for pair in pairs if pair['duplicate_screen']]
        if duplicates:
            summary['reasons'].append('EXACT_OR_NEAR_DUPLICATE_CANDIDATE_PAIR; C review required')
        elif qa_config.get('thresholds_frozen_by_C') is not True:
            summary['reasons'].append('INDEPENDENCE_THRESHOLDS_NOT_FROZEN_BY_C')
        else:
            summary['INDEPENDENCE_SCREEN_PASS'] = True
            summary['status'] = 'INDEPENDENCE_SCREEN_PASS_C_SIGNOFF_REQUIRED'
    except Exception as error:
        summary['INDEPENDENCE_SCREEN_PASS'] = False
        summary['reasons'].append(type(error).__name__ + ':' + str(error))
    save_json(args.out / 'independence_report.json', summary)
    files = sorted(path for path in args.out.iterdir() if path.is_file())
    with (args.out / 'SHA256SUMS.txt').open('x') as stream:
        for path in files:
            stream.write(f'{sha256_file(path)}  {path.name}\n')
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0 if summary['INDEPENDENCE_SCREEN_PASS'] else 2


if __name__ == '__main__':
    sys.exit(main())
