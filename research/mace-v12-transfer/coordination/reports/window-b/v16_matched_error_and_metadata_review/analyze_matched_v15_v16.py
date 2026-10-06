#!/usr/bin/env python3
"""Compare frozen V15/V16 models on their six already-scored registry probes."""
import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np

PROJECT = Path('research/mace-v12-transfer/periodic_interface_v4')
REPORT = Path('research/mace-v12-transfer/coordination/reports/window-b/v16_matched_error_and_metadata_review')
TOL_FORCE = 1e-6
TOL_ENERGY = 1e-6
CUTOFF = 4.5
KINDS = ['AgC', 'AgSi', 'AgTi']


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text())


def dump_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def region_metrics(mask, delta, atom_ids, symbols):
    values = np.asarray(mask, dtype=bool)
    errors = np.linalg.norm(delta, axis=1)
    sse = errors ** 2
    return {
        'atoms': int(values.sum()),
        'force_vector_RMSE_eV_A': float(np.sqrt(np.mean(sse[values]))) if values.any() else None,
        'force_error_SSE_share': float(sse[values].sum() / sse.sum()) if sse.sum() else 0.0,
        'max_atom_error_eV_A': float(errors[values].max()) if values.any() else None,
    }


def run(repo, out):
    from ase.io import read
    import torch
    from mace.calculators import MACECalculator

    repo = repo.resolve()
    project = repo / PROJECT
    draft = project / 'mace_periodic_v16_interface_energy'
    v15 = project / 'mace_periodic_v15_interface_energy'
    v15_archive = project / 'pbe_interface_v15_registry_holdouts'
    v16_archive = project / 'pbe_interface_v16_blind_holdouts'
    accessed = {}

    def record(path):
        p = Path(path).resolve()
        require(p.is_relative_to(repo), 'Evidence path escapes repository')
        accessed[str(p.relative_to(repo))] = sha(p)
        return p

    def read_json(path):
        return json.loads(record(path).read_text())

    # Recheck the limited metadata repair contract before any model inference.
    manifest = read_json(draft / 'data/dataset_manifest.json')
    pins = read_json(draft / 'input_pins.json')
    before = read_json(draft / 'metadata_revision_before.json')
    require(sum(manifest['provenance_status_counts_train'].values()) == 35, 'V16 provenance total mismatch')
    require(sum(manifest['evidence_tier_counts_train'].values()) == 35, 'V16 evidence-tier total mismatch')
    require(manifest['split_sizes'] == {'train': 35, 'valid': 2, 'test': 3}, 'V16 split size mismatch')
    for relative, expected in pins['sha256'].items():
        require(sha(record(repo / relative)) == expected, 'V16 pinned input mismatch: ' + relative)
    old_manifest = json.loads(__import__('subprocess').check_output(
        ['git', 'show', 'b6ac255:' + str(draft.relative_to(repo) / 'data/dataset_manifest.json')],
        cwd=repo, text=True))
    old_pins = json.loads(__import__('subprocess').check_output(
        ['git', 'show', 'b6ac255:' + str(draft.relative_to(repo) / 'input_pins.json')],
        cwd=repo, text=True))
    require(before['manifest'] == old_manifest and before['input_pins'] == old_pins,
            'Stored pre-repair snapshot does not match b6ac255')
    require(before['manifest'] != manifest and before['input_pins'] != pins,
            'Metadata repair snapshot is not distinct from current metadata')
    split_hashes = {}
    for split in ('train', 'valid', 'test'):
        path = draft / 'data' / (split + '.extxyz')
        current = sha(record(path))
        old_bytes = __import__('subprocess').check_output(
            ['git', 'show', 'b6ac255:' + str(path.relative_to(repo))], cwd=repo)
        require(current == hashlib.sha256(old_bytes).hexdigest(), 'V16 split changed since b6ac255: ' + split)
        require(current == manifest['output_sha256'][split], 'V16 split differs from manifest: ' + split)
        split_hashes[split] = current
    metadata_check = {
        'repair_verified': True,
        'provenance_status_total': sum(manifest['provenance_status_counts_train'].values()),
        'evidence_tier_total': sum(manifest['evidence_tier_counts_train'].values()),
        'split_sizes': manifest['split_sizes'],
        'pre_repair_snapshot_matches_commit_b6ac255': True,
        'train_valid_test_bytes_unchanged_from_b6ac255': True,
        'split_sha256': split_hashes,
        'relevant_pin_sha256': {k: v for k, v in pins['sha256'].items()
                                if k.endswith(('run_v16_training.sh', 'preflight.py', 'evaluate_v16.py'))},
    }

    import subprocess
    # Capture the reviewed frozen model and score records before loading either model.
    versions = {}
    for version, root in [('v15', v15), ('v16', draft)]:
        selection_path = root / 'selected_model_record.json'
        eval_path = root / 'evaluation_01/evaluation.json'
        selection = read_json(selection_path)
        evaluation = read_json(eval_path)
        model = record(root / 'training_01/checkpoints' / f'MACE_periodic_{version}_interface_energy_run-45.model')
        require(sha(model) == selection['model_sha256'] == evaluation['model_sha256'], version + ' model hash mismatch')
        require(evaluation['selection_record_sha256'] == accessed[str(selection_path.relative_to(repo))],
                version + ' selected-model record hash mismatch')
        require(selection['reviewed_by'] == 'window-a' and selection['completed_epochs'] == 80
                and selection['blind_labels_used_for_selection'] is False, version + ' selection record not frozen/reviewed')
        if version == 'v15':
            require(selection['training_exit_code'] == 0, 'V15 selected model lacks successful training exit')
        else:
            require(selection['training_exit_code'] in (0, None), 'V16 selected model has invalid training exit')
        for evidence in (selection['training_stdout'], selection['epoch_completion_evidence']):
            require(sha(record(evidence['path'])) == evidence['sha256'], version + ' training completion evidence changed')
        for name, expected in [
                ('evaluation.json', None), ('per_structure.csv', None)]:
            file = root / 'evaluation_01' / name
            recorded = next(line.split('  ', 1) for line in record(root / 'evaluation_01/SHA256SUMS.txt').read_text().splitlines()
                            if line.endswith('  ' + name))
            require(sha(record(file)) == recorded[0], version + ' published evaluation checksum mismatch: ' + name)
        model_bytes = model
        versions[version] = {
            'root': root,
            'selection': selection,
            'selection_sha256': accessed[str(selection_path.relative_to(repo))],
            'evaluation': evaluation,
            'model': model_bytes,
            'model_sha256': sha(model_bytes),
            'rows': {r['label']: r for r in evaluation['per_structure']},
        }

    # Validate all six archived structures and their immutable inputs/hashes.
    structures = {}
    reference_hashes = {}
    for version, archive, expected_role in [
            ('v15', v15_archive, 'v15_blind_registry_holdout'),
            ('v16', v16_archive, 'v16_blind_local_registry_holdout')]:
        archive_manifest = read_json(archive / 'archive_manifest.json')
        input_manifest = read_json(archive / 'input_manifest.json')
        require(archive_manifest['complete'] is True and len(archive_manifest['records']) == 3,
                version + ' archive manifest incomplete')
        inputs = {r['label']: r for r in input_manifest['records']}
        require(set(inputs) == {f'{kind}_registry_holdout_{version}_01' for kind in KINDS},
                version + ' input labels unexpected')
        for row in archive_manifest['records']:
            label = row['label']
            require(row['status'] == 'PASS' and row['role'] == expected_role
                    and all(value is True for value in row['checks'].values()),
                    version + ' archive not verified: ' + label)
            item = inputs[label]
            folder = archive / 'calculations' / label
            reference = folder / (label + '_PW_PBE.extxyz')
            input_path = archive / item['input']
            record(reference)
            ref_hash = sha(reference)
            require(ref_hash == row['sha256'][reference.name], 'Archive reference hash mismatch: ' + label)
            require(ref_hash == versions[version]['rows'][label]['reference_sha256'],
                    'Reference differs from frozen evaluation: ' + label)
            require(sha(record(input_path)) == item['input_sha256'] == row['input_sha256'],
                    'Input hash chain mismatch: ' + label)
            verification = read_json(folder / 'verification.json')
            summary = read_json(folder / 'summary.json')
            progress = read_json(folder / 'progress.json')
            require(verification['status'] == 'PASS' and all(x is True for x in verification['checks'].values()),
                    'Archive verification changed: ' + label)
            require(summary['scf_converged'] is True and progress['scf_converged'] is True,
                    'Reference SCF not converged: ' + label)
            for member, digest in row['sha256'].items():
                require(sha(record(folder / member)) == digest, 'Archive member changed: ' + label + '/' + member)
            ref_atoms = read(reference)
            input_atoms = read(input_path)
            require(np.array_equal(ref_atoms.numbers, input_atoms.numbers)
                    and np.array_equal(ref_atoms.pbc, input_atoms.pbc)
                    and np.array_equal(ref_atoms.arrays['lammps_id'], input_atoms.arrays['lammps_id'])
                    and np.array_equal(ref_atoms.arrays['central_pair'], input_atoms.arrays['central_pair'])
                    and np.allclose(ref_atoms.positions, input_atoms.positions, atol=1e-12, rtol=0)
                    and np.allclose(ref_atoms.cell.array, input_atoms.cell.array, atol=1e-12, rtol=0),
                    'Reference geometry differs from immutable input: ' + label)
            require('PW_PBE_energy_eV' in ref_atoms.info and 'PW_PBE_forces' in ref_atoms.arrays,
                    'Reference lacks archived labels: ' + label)
            structures[label] = {
                'version': version,
                'atoms': ref_atoms,
                'reference_path': reference,
                'input_path': input_path,
                'input_record': item,
                'reference_sha256': ref_hash,
            }
            reference_hashes[label] = ref_hash

    # Keep inference deterministic, CPU-only and single-threaded.
    torch.set_default_dtype(torch.float64)
    torch.set_num_threads(1)
    calculators = {}
    runtime = {'torch': torch.__version__, 'mace-torch': importlib.metadata.version('mace-torch'),
               'ase': importlib.metadata.version('ase'), 'torch_threads': torch.get_num_threads(),
               'device': 'cpu', 'dtype': 'float64'}
    for version, info in versions.items():
        loaded = torch.load(info['model'], map_location='cpu', weights_only=False).double().eval()
        calculators[version] = MACECalculator(models=[loaded], device='cpu', default_dtype='float64')

    summaries = []
    atom_rows = []
    for label, item in structures.items():
        atoms = item['atoms']
        symbols = np.asarray(atoms.get_chemical_symbols())
        atom_ids = np.asarray(atoms.arrays['lammps_id'], dtype=int)
        central = np.flatnonzero(np.asarray(atoms.arrays['central_pair'], dtype=bool))
        require(len(central) == 2 and len(set(atom_ids)) == len(atoms), 'Invalid atom IDs/central pair: ' + label)
        ag = [int(i) for i in central if symbols[i] == 'Ag']
        contact = [int(i) for i in central if symbols[i] != 'Ag']
        require(len(ag) == len(contact) == 1, 'Marked pair is not one Ag-X pair: ' + label)
        ia, ix = ag[0], contact[0]
        pair_ids = {'Ag_id': int(atom_ids[ia]), 'X_id': int(atom_ids[ix]), 'X_element': str(symbols[ix])}
        pair_meta = item['input_record']['central_pair']
        require(pair_ids['Ag_id'] == pair_meta['Ag_id'] and pair_ids['X_id'] == pair_meta.get('contact_id', pair_meta.get('X_id')),
                'Marked pair IDs differ from input manifest: ' + label)
        fref = np.asarray(atoms.arrays['PW_PBE_forces'], dtype=float)
        eref = float(atoms.info['PW_PBE_energy_eV'])
        distance = atoms.get_all_distances(mic=True)
        local_distance = np.min(distance[:, central], axis=1)
        local = local_distance <= CUTOFF
        u = atoms.positions[ix] - atoms.positions[ia]
        u /= np.linalg.norm(u)
        base = {'label': label, 'source_holdout': item['version'], 'atoms': len(atoms),
                'formula': atoms.get_chemical_formula(), 'reference_sha256': item['reference_sha256'],
                'input_sha256': item['input_record']['input_sha256'], 'marked_pair': pair_ids,
                'local_cutoff_A': CUTOFF}
        for version, calc in calculators.items():
            predicted = atoms.copy()
            predicted.calc = calc
            ep = float(predicted.get_potential_energy())
            fp = np.asarray(predicted.get_forces(), dtype=float)
            require(np.isfinite(ep) and fp.shape == fref.shape and np.isfinite(fp).all(),
                    'Invalid prediction: ' + version + '/' + label)
            delta = fp - fref
            norms = np.linalg.norm(delta, axis=1)
            sse = norms ** 2
            radial_ref = float(np.dot(fref[ix] - fref[ia], u))
            radial_pred = float(np.dot(fp[ix] - fp[ia], u))
            signed_sep_error = radial_pred - radial_ref
            ag_projection = float(np.dot(delta[ia], u))
            x_projection = float(np.dot(delta[ix], u))
            ag_transverse = delta[ia] - ag_projection * u
            x_transverse = delta[ix] - x_projection * u
            groups = {
                'all': np.ones(len(atoms), dtype=bool),
                'within_4p5A': local,
                'outside_4p5A': ~local,
                'Ag_within_4p5A': (symbols == 'Ag') & local,
                'framework_within_4p5A': (symbols != 'Ag') & local,
                'Ag_outside_4p5A': (symbols == 'Ag') & (~local),
                'framework_outside_4p5A': (symbols != 'Ag') & (~local),
            }
            region = {name: region_metrics(mask, delta, atom_ids, symbols) for name, mask in groups.items()}
            species = {element: region_metrics(symbols == element, delta, atom_ids, symbols)
                       for element in sorted(set(symbols))}
            max_indices = np.argsort(norms)[::-1][:10]
            top = []
            for index in max_indices:
                neighbors = [int(j) for j in np.argsort(distance[index]) if j != index][:4]
                top.append({'id': int(atom_ids[index]), 'element': str(symbols[index]),
                            'error_norm_eV_A': float(norms[index]),
                            'error_vector_xyz_eV_A': delta[index].tolist(),
                            'distance_to_marked_pair_A': float(local_distance[index]),
                            'local_within_4p5A': bool(local[index]),
                            'nearest_4_MIC_neighbors': [
                                {'id': int(atom_ids[j]), 'element': str(symbols[j]),
                                 'distance_A': float(distance[index, j])} for j in neighbors]})
            metrics = {
                'force_vector_RMSE_eV_A': float(np.sqrt(np.mean(sse))),
                'force_vector_max_error_eV_A': float(norms.max()),
                'energy_abs_error_meV_atom': abs(ep - eref) * 1000.0 / len(atoms),
                'energy_prediction_eV_cell': ep,
                'energy_reference_eV_cell': eref,
                'separating_force_reference_eV_A': radial_ref,
                'separating_force_prediction_eV_A': radial_pred,
                'separating_force_signed_error_pred_minus_ref_eV_A': signed_sep_error,
                'separating_force_abs_error_eV_A': abs(signed_sep_error),
                'marked_pair_error_decomposition': {
                    'Ag_signed_longitudinal_contribution_eV_A': -ag_projection,
                    'X_signed_longitudinal_contribution_eV_A': x_projection,
                    'sum_signed_separation_error_eV_A': -ag_projection + x_projection,
                    'Ag_transverse_error_norm_eV_A': float(np.linalg.norm(ag_transverse)),
                    'X_transverse_error_norm_eV_A': float(np.linalg.norm(x_transverse)),
                    'Ag_full_vector_error_norm_eV_A': float(norms[ia]),
                    'X_full_vector_error_norm_eV_A': float(norms[ix]),
                },
                'regions': region,
                'by_element': species,
                'top_10_atoms': top,
            }
            require(abs((-ag_projection + x_projection) - signed_sep_error) < 1e-10,
                    'Marked-pair projection decomposition failed: ' + version + '/' + label)
            # Reproduce the original frozen score when this model is on its own three probes.
            if version == item['version']:
                prior = versions[version]['rows'][label]
                checks = {
                    'force_vector_RMSE_eV_A': abs(metrics['force_vector_RMSE_eV_A'] - prior['force_vector_RMSE_eV_A']),
                    'force_vector_max_error_eV_A': abs(metrics['force_vector_max_error_eV_A'] - prior['force_vector_max_error_eV_A']),
                    'energy_prediction_eV_cell': abs(ep - prior['energy_prediction_eV_cell']),
                    'energy_abs_error_meV_atom': abs(metrics['energy_abs_error_meV_atom'] - prior['energy_abs_error_meV_atom']),
                    'separating_force_abs_error_eV_A': abs(metrics['separating_force_abs_error_eV_A'] - prior['separating_force_abs_error_eV_A']),
                }
                require(max(checks[k] for k in checks if k != 'energy_prediction_eV_cell') <= TOL_FORCE
                        and checks['energy_prediction_eV_cell'] <= TOL_ENERGY,
                        'Frozen score reproduction failed: ' + version + '/' + label + ' ' + json.dumps(checks))
                metrics['published_score_reproduction_abs_deltas'] = checks
            summaries.append({**base, 'model': version, **metrics})
            for i in range(len(atoms)):
                atom_rows.append({
                    'label': label, 'source_holdout': item['version'], 'model': version,
                    'atom_index': i, 'atom_id': int(atom_ids[i]), 'element': str(symbols[i]),
                    'distance_to_marked_pair_A': float(local_distance[i]), 'within_4p5A': bool(local[i]),
                    'DFT_Fx_eV_A': float(fref[i, 0]), 'DFT_Fy_eV_A': float(fref[i, 1]), 'DFT_Fz_eV_A': float(fref[i, 2]),
                    'MACE_Fx_eV_A': float(fp[i, 0]), 'MACE_Fy_eV_A': float(fp[i, 1]), 'MACE_Fz_eV_A': float(fp[i, 2]),
                    'error_Fx_eV_A': float(delta[i, 0]), 'error_Fy_eV_A': float(delta[i, 1]), 'error_Fz_eV_A': float(delta[i, 2]),
                    'error_norm_eV_A': float(norms[i]),
                })

    # Pairwise matched comparison across the exact same six references.
    pairs = []
    aggregate = []
    for label in structures:
        a = next(x for x in summaries if x['label'] == label and x['model'] == 'v15')
        b = next(x for x in summaries if x['label'] == label and x['model'] == 'v16')
        pairs.append({
            'label': label, 'source_holdout': structures[label]['version'], 'atoms': a['atoms'],
            'reference_sha256': a['reference_sha256'],
            'v15_force_RMSE_eV_A': a['force_vector_RMSE_eV_A'],
            'v16_force_RMSE_eV_A': b['force_vector_RMSE_eV_A'],
            'delta_v16_minus_v15_force_RMSE_eV_A': b['force_vector_RMSE_eV_A'] - a['force_vector_RMSE_eV_A'],
            'v15_energy_abs_error_meV_atom': a['energy_abs_error_meV_atom'],
            'v16_energy_abs_error_meV_atom': b['energy_abs_error_meV_atom'],
            'v15_separation_abs_error_eV_A': a['separating_force_abs_error_eV_A'],
            'v16_separation_abs_error_eV_A': b['separating_force_abs_error_eV_A'],
            'v15_regions': a['regions'], 'v16_regions': b['regions'],
            'v15_marked_pair': a['marked_pair'], 'v16_marked_pair': b['marked_pair'],
        })
    for scope, selected in [('all_six', list(structures)),
                            *[(kind, [label for label in structures if label.startswith(kind + '_')]) for kind in KINDS]]:
        require(selected, 'Empty comparison scope: ' + scope)
        row = {'scope': scope, 'structures': len(selected), 'atoms': sum(len(structures[x]['atoms']) for x in selected)}
        for version in ('v15', 'v16'):
            chosen = [x for x in summaries if x['model'] == version and x['label'] in selected]
            row[version] = {
                'pooled_force_vector_RMSE_eV_A': float(np.sqrt(sum(x['atoms'] * x['force_vector_RMSE_eV_A'] ** 2 for x in chosen) / sum(x['atoms'] for x in chosen))),
                'energy_MAE_meV_atom': float(np.mean([x['energy_abs_error_meV_atom'] for x in chosen])),
                'max_separation_abs_error_eV_A': float(max(x['separating_force_abs_error_eV_A'] for x in chosen)),
                'max_force_vector_error_eV_A': float(max(x['force_vector_max_error_eV_A'] for x in chosen)),
            }
        row['delta_v16_minus_v15_pooled_force_RMSE_eV_A'] = row['v16']['pooled_force_vector_RMSE_eV_A'] - row['v15']['pooled_force_vector_RMSE_eV_A']
        aggregate.append(row)

    # Include hashes for the existing targeted proposal design used by recommendations.
    design = repo / 'research/mace-v12-transfer/coordination/reports/window-b/targeted_transverse_framework_design/candidate_set_final'
    design_report = design / 'proposal_report.md'
    design_manifest = design / 'proposal_manifest.json'
    design_proposals = read_json(design_manifest)
    v17_manifest_path = project / 'pbe_interface_v17_parallel_targeted_acquisition/input_manifest.json'
    v17_manifest = read_json(v17_manifest_path)
    recommendations = {
        'evidence_scope': 'Retrospective comparison of frozen V15/V16 models on the six already-scored V15/V16 registry references; not a new validation set.',
        'interpretation': [
            'Use the already-approved paired V17 framework/transverse acquisitions as a small diagnostic training batch; compare paired signs and preserve the existing training recipe for the first data-only comparison.',
            'Prioritize both AgSi arms because vector errors can be transverse even when the scalar separation error passes; the matched per-atom/contact projections below test whether this persists in V16.',
            'Keep AgC and AgTi framework-contact pairs as the focused coverage checks. If errors remain outside the 4.5 A contact shell, propose a separate framework-shell pair around the identified high-error atom only after A reviews the matched IDs/neighbors.',
            'Do not use these six retrospective probes as fresh blind validation or select a model by them. Freeze a new independent geometry set and reference labels before comparing V17 models.',
        ],
        'design_proposal_sha256': {'proposal_report.md': sha(record(design_report)), 'proposal_manifest.json': sha(record(design_manifest))},
        'already_registered_v17_inputs': [
            {'label': x['label'], 'input_sha256': x['input_sha256'], 'mode': x['design']['mode'], 'sign': x['design']['sign'],
             'changed_atoms': x['design']['changed_atoms']}
            for x in v17_manifest['records']
        ],
        'design_candidate_examples': [
            {'label': x['label'], 'mode': x['mode'], 'sign': x['sign'],
             'changed_atoms': x['candidate']['changed_atoms']}
            for group in design_proposals['priority_subset'] for x in group['variants']
        ],
        'causal_test_needed': {
            'coverage_arm': 'Add only the approved paired-sign acquisition labels, retrain with unchanged data recipe, and compare per-atom/local/outside metrics on the same frozen diagnostic references plus a separate fresh independent set.',
            'loss_weight_arm': 'On the same fixed dataset and seed schedule, change one loss weight at a time; keep model size, optimizer, training steps, and split fixed. Compare vector, longitudinal, and transverse errors on the same independent labels.',
            'decision_rule': 'Coverage is supported only if the data-only addition improves relevant errors on independent structures. A loss-weight effect requires a controlled fixed-data one-factor comparison; the current model-to-model change is confounded by data and selected epoch.'
        }
    }

    for version, info in versions.items():
        require(sha(info['model']) == info['model_sha256'], version + ' model changed during inference')
    report = {
        'task': 'v16_matched_error_and_metadata_review',
        'scope': 'Both window-A-reviewed frozen selected models evaluated on the identical union of three archived V15 registry probes and three archived V16 registry probes. Retrospective diagnosis only; no selection, training, new DFT, MD or TTM.',
        'metadata_repair_check': metadata_check,
        'models': {v: {'model_sha256': x['model_sha256'], 'selected_model_record_sha256': x['selection_sha256'],
                       'selected_epoch': x['selection']['selected_epoch'],
                       'selection_record_reviewed_by': x['selection']['reviewed_by'],
                       'blind_labels_used_for_selection': x['selection']['blind_labels_used_for_selection']}
                   for v, x in versions.items()},
        'references': reference_hashes,
        'runtime': runtime,
        'geometry_rule': {'local_cutoff_A': CUTOFF, 'distance': 'ASE minimum-image distance to either marked Ag/X atom',
                          'longitudinal_axis': 'direct marked Ag-to-X coordinate vector, matching A evaluation convention',
                          'framework': 'all non-Ag atoms'},
        'reproduction_tolerances': {'force_metrics_eV_A': TOL_FORCE, 'energy_prediction_eV_cell': TOL_ENERGY},
        'matched_pairs': pairs,
        'aggregate_comparison': aggregate,
        'per_model_per_geometry': summaries,
        'recommendations': recommendations,
        'input_sha256': accessed,
    }
    dump_json(out / 'matched_comparison.json', report)
    with (out / 'per_atom_errors.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(atom_rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(atom_rows)
    with (out / 'matched_metrics.csv').open('w', newline='') as f:
        columns = ['label', 'source_holdout', 'atoms', 'reference_sha256',
                   'v15_force_RMSE_eV_A', 'v16_force_RMSE_eV_A', 'delta_v16_minus_v15_force_RMSE_eV_A',
                   'v15_energy_abs_error_meV_atom', 'v16_energy_abs_error_meV_atom',
                   'v15_separation_abs_error_eV_A', 'v16_separation_abs_error_eV_A']
        writer = csv.DictWriter(f, fieldnames=columns, lineterminator='\n')
        writer.writeheader()
        for row in pairs:
            writer.writerow({key: row[key] for key in columns})
    dump_json(out / 'recommendations.json', recommendations)
    print(json.dumps({'aggregate_comparison': aggregate, 'matched_metrics': [
        {k: x[k] for k in ['label', 'source_holdout', 'v15_force_RMSE_eV_A', 'v16_force_RMSE_eV_A',
                           'delta_v16_minus_v15_force_RMSE_eV_A', 'v15_separation_abs_error_eV_A',
                           'v16_separation_abs_error_eV_A']} for x in pairs]}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    out = args.output_dir.resolve()
    require(out.is_relative_to(repo / REPORT) and out != repo / REPORT, 'Output must be a child of report directory')
    require(not out.exists() or (out.is_dir() and not any(out.iterdir())), 'Nonempty output refused')
    out.mkdir(parents=True, exist_ok=True)
    try:
        run(repo, out)
    except Exception as exc:
        dump_json(out / 'diagnostics_failure.json', {'error': type(exc).__name__, 'message': str(exc)})
        raise


if __name__ == '__main__':
    main()
