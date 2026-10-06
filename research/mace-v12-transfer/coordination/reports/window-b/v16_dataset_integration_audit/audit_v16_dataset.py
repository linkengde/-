#!/usr/bin/env python3
"""Read-only audit of V16 training data; never opens fresh holdout reference files."""
import argparse
import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path

import numpy as np
from ase.io import read


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(v) for v in value]
    return value


def frame_identity(a, b):
    require(np.array_equal(a.numbers, b.numbers), 'Atomic species/order changed')
    require(np.array_equal(a.pbc, b.pbc), 'PBC changed')
    require(np.array_equal(a.cell.array, b.cell.array), 'Cell changed')
    require(np.array_equal(a.positions, b.positions), 'Positions changed')
    require(set(a.arrays) == set(b.arrays), 'Atom array fields changed')
    for key in a.arrays:
        require(np.array_equal(a.arrays[key], b.arrays[key]), f'Atom array changed: {key}')
    require(jsonable(a.info) == jsonable(b.info), 'Frame metadata changed')


def finite_labels(frames, split):
    for i, atoms in enumerate(frames):
        require('REF_energy' in atoms.info, f'{split}[{i}] missing REF_energy')
        forces = np.asarray(atoms.arrays.get('REF_forces'), dtype=float)
        require(forces.shape == (len(atoms), 3), f'{split}[{i}] invalid REF_forces shape')
        require(np.isfinite(float(atoms.info['REF_energy'])) and np.isfinite(forces).all(),
                f'{split}[{i}] has nonfinite labels')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--output-json', type=Path, required=True)
    args = ap.parse_args()
    repo = args.repo_root.resolve()
    project = repo / 'research/mace-v12-transfer/periodic_interface_v4'
    v16 = project / 'mace_periodic_v16_interface_energy'
    v15 = project / 'mace_periodic_v15_interface_energy'
    v14 = project / 'mace_periodic_v14_interface_energy'
    data = v16 / 'data'
    manifest = json.loads((data / 'dataset_manifest.json').read_text())
    parent = json.loads((v15 / 'data/dataset_manifest.json').read_text())
    pins = json.loads((v16 / 'input_pins.json').read_text())['sha256']
    for relative, expected in pins.items():
        require(sha(repo / relative) == expected, f'Pinned input hash mismatch: {relative}')

    expected_counts = {'train': 35, 'valid': 2, 'test': 3}
    require(manifest['split_sizes'] == expected_counts, 'V16 split counts differ')
    frames = {}
    for split, expected in expected_counts.items():
        path = data / f'{split}.extxyz'
        require(sha(path) == manifest['output_sha256'][split], f'{split} manifest hash mismatch')
        frames[split] = read(path, index=':')
        require(len(frames[split]) == expected, f'{split} frame count differs')
        finite_labels(frames[split], split)
        require(manifest['energy_label_counts'][split] == expected, f'{split} energy count differs')
        require(manifest['force_label_counts'][split] == expected, f'{split} force count differs')
    for line in (data / 'SHA256SUMS.txt').read_text().splitlines():
        digest, filename = line.split('  ', 1)
        require(sha(data / filename) == digest, f'Dataset SHA256SUMS mismatch: {filename}')

    composition = np.array([
        [atoms.get_chemical_symbols().count(el) for el in ['Ag', 'C', 'Si', 'Ti']]
        for atoms in frames['train']
    ])
    rank = int(np.linalg.matrix_rank(composition))
    require(rank == 4 and rank == manifest['energy_composition_matrix']['rank'],
            'Training composition rank differs')
    require(composition.tolist() == manifest['energy_composition_matrix']['matrix'],
            'Training composition matrix differs')

    parent_train = read(v15 / 'data/train.extxyz', index=':')
    require(len(parent_train) == 29 and len(frames['train']) == len(parent_train) + 6,
            'Expected V15 retained 29 plus six V16 additions')
    require(manifest['parent_V15_manifest_sha256'] == sha(v15 / 'data/dataset_manifest.json'),
            'V15 manifest source hash mismatch')
    inherited_method_fields = {
        'source_method_present': sum(bool(a.info.get('source_method')) for a in parent_train),
        'source_method_missing': sum(not bool(a.info.get('source_method')) for a in parent_train),
        'dft_method_present': sum(bool(a.info.get('dft_method')) for a in parent_train),
        'dft_method_missing': sum(not bool(a.info.get('dft_method')) for a in parent_train),
        'either_method_field_present': sum(bool(a.info.get('source_method') or a.info.get('dft_method'))
                                           for a in parent_train),
        'both_method_fields_present': sum(bool(a.info.get('source_method') and a.info.get('dft_method'))
                                          for a in parent_train),
        'neither_method_field_present': sum(not bool(a.info.get('source_method') or a.info.get('dft_method'))
                                            for a in parent_train),
        'distinct_source_method_values': sorted({str(a.info.get('source_method')) for a in parent_train}),
        'distinct_dft_method_values': sorted({str(a.info.get('dft_method')) for a in parent_train}),
    }
    for i, (old, new) in enumerate(zip(parent_train, frames['train'][:29])):
        try:
            frame_identity(old, new)
        except AssertionError as exc:
            raise AssertionError(f'Retained V15 training frame {i}: {exc}') from exc

    historical_hashes = {}
    for split in ('valid', 'test'):
        hashes = {
            'V16': sha(data / f'{split}.extxyz'),
            'V15': sha(v15 / f'data/{split}.extxyz'),
            'V14': sha(v14 / f'data/{split}.extxyz'),
        }
        require(len(set(hashes.values())) == 1, f'{split} bytes differ across V14/V15/V16')
        historical_hashes[split] = hashes

    require(manifest['per_frame_provenance'] == parent['per_frame_provenance'],
            'Inherited V15 provenance ledger changed')
    require(manifest['excluded_frame_ledger'] == parent['excluded_frame_ledger'],
            'V15 excluded-row ledger changed')
    require(len(manifest['per_frame_provenance']) == 34, 'Expected 29/2/3 inherited provenance rows')
    require(len(manifest['excluded_frame_ledger']) == 5, 'Expected four LCAO and one unresolved exclusion')
    excluded_types = {x['config_type'] for x in manifest['excluded_frame_ledger']}
    require(excluded_types == {
        'AgC_primary_2.0A_LCAO', 'AgC_primary_2.4A_LCAO',
        'AgTi_2p1A_LCAO', 'AgTi_2p7A_LCAO',
        'AgTi_2p3871A_periodic_PW_PBE_force_only',
    }, 'Unexpected excluded-row identity')

    # Verify all six additions against their separate, archived PW-PBE sources.
    additions = {x['label']: x for x in manifest['V16_additions']}
    require(len(additions) == 6, 'Expected six V16 additions')
    source_dirs = ['pbe_interface_v16_main_acquisition', 'pbe_interface_v16_parallel_acquisition']
    source_records = []
    source_methods = set()
    energy_conventions = set()
    source_hashes = manifest['V16_sources_sha256']
    for dirname in source_dirs:
        folder = project / dirname
        im = json.loads((folder / 'input_manifest.json').read_text())
        archive_manifest = json.loads((folder / 'archive_manifest.json').read_text())
        archive_records = {r['label']: r for r in archive_manifest['records']}
        require(archive_manifest.get('complete') is True and archive_manifest.get('missing_labels') == [],
                f'{dirname} archive manifest incomplete')
        require(source_hashes.get(str((folder / 'input_manifest.json').relative_to(repo))) ==
                sha(folder / 'input_manifest.json'), f'{dirname} manifest hash mismatch')
        for rec in im['records']:
            label = rec['label']
            require(label in additions, f'Unlisted V16 addition source: {label}')
            result = folder / 'calculations' / label
            verification = json.loads((result / 'verification.json').read_text())
            summary = json.loads((result / 'summary.json').read_text())
            progress = json.loads((result / 'progress.json').read_text())
            archived_record = archive_records.get(label)
            require(archived_record is not None and archived_record['status'] == 'PASS' and
                    all(archived_record['checks'].values()), f'{label} not in a passing archive manifest')
            source_input = folder / rec['input']
            result_path = result / f'{label}_PW_PBE.extxyz'
            require(sha(source_input) == rec['input_sha256'] == summary['source_sha256'],
                    f'{label} input identity mismatch')
            input_relative = str(source_input.relative_to(repo))
            require(source_hashes.get(input_relative) == sha(source_input),
                    f'{label} input hash not included in V16 source ledger')
            require(verification['status'] == 'PASS' and all(verification['checks'].values()),
                    f'{label} source verification failed')
            require(summary['scf_converged'] is True and progress['status'] == 'complete' and
                    progress['mpi_ranks'] == 4 and summary['mpi_ranks'] == 4,
                    f'{label} convergence/rank evidence failed')
            for filename, digest in verification['sha256'].items():
                member = result / filename
                require(sha(member) == digest, f'{label} archived hash mismatch: {filename}')
                relative = str(member.relative_to(repo))
                require(source_hashes.get(relative) == digest,
                        f'{label} archive hash not included in V16 source ledger: {filename}')
            require(sha(result_path) == additions[label]['source_sha256'],
                    f'{label} addition source hash mismatch')
            raw_input = read(source_input)
            archived = read(result_path)
            integrated = frames['train'][29 + len(source_records)]
            for left, right, what in [(raw_input, archived, 'input/archive'),
                                      (archived, integrated, 'archive/train')]:
                require(np.array_equal(left.numbers, right.numbers), f'{label} {what} species/order mismatch')
                require(np.array_equal(left.pbc, right.pbc) and np.array_equal(left.cell.array, right.cell.array),
                        f'{label} {what} boundary/cell mismatch')
                require(np.allclose(left.positions, right.positions, atol=1e-12, rtol=0),
                        f'{label} {what} positions mismatch')
                if 'lammps_id' in left.arrays:
                    require(np.array_equal(left.arrays['lammps_id'], right.arrays['lammps_id']),
                            f'{label} {what} atom IDs mismatch')
                if 'central_pair' in left.arrays:
                    require(np.array_equal(left.arrays['central_pair'], right.arrays['central_pair']),
                            f'{label} {what} marked-pair identity mismatch')
            require(integrated.info.get('config_type') == label + '_PW_PBE_v16_train' and
                    integrated.info.get('dataset_role') == 'training_acquisition',
                    f'{label} V16 role metadata mismatch')
            require(integrated.info.get('source_method') == archived.info.get('source_method'),
                    f'{label} source_method metadata mismatch')
            require(integrated.info.get('REF_energy') == archived.info.get('PW_PBE_energy_eV') and
                    np.array_equal(integrated.arrays['REF_forces'], archived.arrays['PW_PBE_forces']),
                    f'{label} integrated labels do not match archive')
            method = summary['method']
            require(method.get('xc') == 'PBE' and method.get('basis') == 'plane wave' and
                    method.get('cutoff_eV') == 500 and method.get('kpts') == [1, 1, 1] and
                    method.get('smearing_eV') == 0.1 and method.get('GPAW_version') == '26.7.0',
                    f'{label} source method metadata differs')
            require(additions[label]['energy_convention'] == summary['energy_convention'],
                    f'{label} energy convention differs')
            source_methods.add(json.dumps(method, sort_keys=True))
            energy_conventions.add(summary['energy_convention'])
            source_records.append({
                'label': label,
                'source_archive_sha256': sha(result_path),
                'input_sha256': sha(source_input),
                'provenance_status': additions[label]['provenance_status'],
                'mpi_ranks': summary['mpi_ranks'],
                'method': method,
                'energy_convention': summary['energy_convention'],
            })
    require(set(additions) == {x['label'] for x in source_records}, 'V16 addition/source identity sets differ')
    blind_source_entries = [name for name in source_hashes if 'pbe_interface_v16_blind_holdouts/' in name]
    require(all(('/inputs/' in name or name.endswith('/pbe_interface_v16_blind_holdouts/input_manifest.json'))
                and '/calculations/' not in name for name in blind_source_entries),
            'V16 source ledger reaches a blind calculation output')

    # Read only the three frozen INPUT geometries. Never open holdout calculations/.
    blind = project / 'pbe_interface_v16_blind_holdouts'
    blind_manifest = json.loads((blind / 'input_manifest.json').read_text())
    comparison_path = repo / 'research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision/build_v15_draft.py'
    spec = importlib.util.spec_from_file_location('v15_geometry_compare', comparison_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    heldout = []
    for rec in blind_manifest['records']:
        inp = (blind / rec['input']).resolve()
        require(inp.parent == (blind / 'inputs').resolve(), 'Attempted path outside frozen input folder')
        require(sha(inp) == rec['input_sha256'], f'{rec["label"]} frozen input hash mismatch')
        atoms = read(inp)
        require(atoms.calc is None and 'REF_energy' not in atoms.info and 'REF_forces' not in atoms.arrays and
                'PW_PBE_energy_eV' not in atoms.info and 'PW_PBE_forces' not in atoms.arrays,
                f'{rec["label"]} input unexpectedly contains labels')
        ids = atoms.arrays.get('lammps_id')
        require(ids is not None and len(set(map(int, ids))) == len(atoms),
                f'{rec["label"]} input IDs are absent or duplicated')
        collisions = []
        pbc_invariance = True
        flipped = atoms.copy()
        flipped.pbc = ~np.asarray(flipped.pbc, dtype=bool)
        for split, other_frames in frames.items():
            for index, other in enumerate(other_frames):
                result = module.comparison(other, atoms)
                changed_pbc = module.comparison(other, flipped)
                classify = None if result is None else (result['exact_geometry'], result['near_geometry'])
                flipped_classify = None if changed_pbc is None else \
                    (changed_pbc['exact_geometry'], changed_pbc['near_geometry'])
                if classify != flipped_classify:
                    pbc_invariance = False
                if result is not None and (result['exact_geometry'] or result['near_geometry']):
                    collisions.append({'split': split, 'frame': index,
                                       'exact': result['exact_geometry'], 'near': result['near_geometry']})
        require(pbc_invariance, f'{rec["label"]} geometry classification depends on PBC')
        require(not collisions, f'{rec["label"]} frozen input overlaps a committed split')
        heldout.append({'label': rec['label'], 'input_sha256': sha(inp),
                        'atoms': len(atoms), 'pbc': atoms.pbc.tolist(),
                        'split_collisions': collisions, 'geometry_screen_pbc_independent': pbc_invariance})

    script = (v16 / 'run_v16_training.sh').read_text()
    required_tokens = [
        '--foundation_model', '--seed 45', '--device cpu', '--batch_size 1',
        '--max_num_epochs 80', '--lr 0.0001', '--weight_decay 5e-7',
        '--energy_weight 100', '--forces_weight 1000', '--valid_file', '--test_file',
    ]
    require(all(token in script for token in required_tokens), 'A controlled training setting is absent')
    require('--check-blind-archives' in script and 'blind' in script.lower(),
            'Training does not gate on blind archive integrity')
    require('blind_holdouts/inputs' not in script and 'blind_holdouts/calculations' not in script,
            'Training script directly references blind files')
    wrapper_mismatches = {
        'actual_script': 'run_v16_training.sh',
        'usage_text': 'run_v15_training.sh' in script,
        'output_error_mentions_v15': 'V15 run directory' in script,
    }
    evaluator = (v16 / 'evaluate_v16.py').read_text()
    evaluator_message_mismatch = 'Output must be V15 run child' in evaluator
    preflight_text = (v16 / 'preflight.py').read_text()
    status_counts = manifest.get('provenance_status_counts_train', {})
    evidence_counts = manifest.get('evidence_tier_counts_train', {})
    provenance_summary = {
        'reported_provenance_status_counts': status_counts,
        'reported_provenance_status_total': sum(status_counts.values()),
        'reported_evidence_tier_counts': evidence_counts,
        'reported_evidence_tier_total': sum(evidence_counts.values()),
        'current_train_frames': len(frames['train']),
        'new_verified_archive_additions': len(source_records),
        'count_summaries_match_current_training_count':
            sum(status_counts.values()) == len(frames['train']) and
            sum(evidence_counts.values()) == len(frames['train']),
        'preflight_checks_provenance_count_summaries':
            'provenance_status_counts_train' in preflight_text and
            'evidence_tier_counts_train' in preflight_text,
    }
    unpinned_entry_files = [name for name in ['run_v16_training.sh', 'preflight.py', 'evaluate_v16.py']
                            if str((v16 / name).relative_to(repo)) not in pins]

    report = {
        'audit': 'V16 dataset and training-entry integration; read-only',
        'blind_reference_content_read': False,
        'counts': expected_counts,
        'energy_and_force_label_counts': {
            split: {'energy': len(items), 'forces': len(items)} for split, items in frames.items()
        },
        'all_labels_finite': True,
        'training_composition_rank': rank,
        'retained_v15_training_frames_semantically_identical': len(parent_train),
        'retained_v15_method_metadata_coverage': inherited_method_fields,
        'historical_valid_test_byte_identity': historical_hashes,
        'inherited_provenance_rows_preserved': len(manifest['per_frame_provenance']),
        'inherited_excluded_rows_preserved': len(manifest['excluded_frame_ledger']),
        'excluded_config_types': sorted(excluded_types),
        'inherited_physical_method_consistency': manifest.get('physical_method_consistency'),
        'inherited_provenance_scope': manifest.get('inherited_provenance_ledger_scope'),
        'provenance_summary_audit': provenance_summary,
        'v16_addition_count': len(source_records),
        'input_pin_count_verified': len(pins),
        'v16_additions': source_records,
        'distinct_new_source_methods': len(source_methods),
        'distinct_new_energy_conventions': sorted(energy_conventions),
        'frozen_holdout_input_isolation': heldout,
        'training_entry_settings_present': required_tokens,
        'training_has_direct_blind_file_argument': False,
        'training_archive_preflight_gate': True,
        'training_wrapper_text_mismatch': wrapper_mismatches,
        'evaluator_v15_wording_mismatch': evaluator_message_mismatch,
        'unpinned_entry_files': unpinned_entry_files,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({
        'counts': expected_counts,
        'rank': rank,
        'retained_v15_frames': len(parent_train),
        'verified_v16_additions': len(source_records),
        'fresh_holdout_input_collisions': sum(len(x['split_collisions']) for x in heldout),
        'blind_reference_content_read': False,
        'training_wrapper_text_mismatch': wrapper_mismatches,
        'unpinned_entry_files': unpinned_entry_files,
    }, indent=2))


if __name__ == '__main__':
    main()
