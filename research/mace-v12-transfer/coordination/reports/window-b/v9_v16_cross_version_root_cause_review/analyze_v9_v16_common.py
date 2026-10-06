#!/usr/bin/env python3
"""Matched retrospective inference on the six already-scored V15/V16 probes."""
import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np

PROJECT = Path('research/mace-v12-transfer/periodic_interface_v4')
REPORT = Path('research/mace-v12-transfer/coordination/reports/window-b/v9_v16_cross_version_root_cause_review')
PREVIOUS = Path('research/mace-v12-transfer/coordination/reports/window-b/v16_matched_error_and_metadata_review/analysis_04')
V14_PREVIOUS = Path('research/mace-v12-transfer/coordination/reports/window-b/force_separation_driver_diagnosis/analysis_final')
TOL = 1e-6
CUTOFF = 4.5
KINDS = ('AgC', 'AgSi', 'AgTi')

MODEL_FILES = {
    'v9': PROJECT / 'mace_periodic_v9_forceonly/models/Ag_Ti_Si_C_MACE_periodic_v9_forceonly.model',
    'v10': PROJECT / 'mace_periodic_v10_energycal/checkpoints/MACE_periodic_v10_energycal_run-43.model',
    'v11': PROJECT / 'mace_periodic_v11_energyfocus/checkpoints/MACE_periodic_v11_energyfocus_run-44.model',
    'v12': PROJECT / 'mace_periodic_v12_interface_energy/checkpoints/MACE_periodic_v12_interface_energy_run-45.model',
    'v13': PROJECT / 'mace_periodic_v13_interface_energy/checkpoints/MACE_periodic_v13_interface_energy_run-45.model',
    'v14': PROJECT / 'mace_periodic_v14_interface_energy/checkpoints/MACE_periodic_v14_interface_energy_run-45.model',
    'v15': PROJECT / 'mace_periodic_v15_interface_energy/training_01/checkpoints/MACE_periodic_v15_interface_energy_run-45.model',
    'v16': PROJECT / 'mace_periodic_v16_interface_energy/training_01/checkpoints/MACE_periodic_v16_interface_energy_run-45.model',
}

DATA_NOTES = {
    'v9': 'Initial four-element force-focused baseline; original dataset and complete training command are absent from the handoff.',
    'v10': 'Added four PBE elemental anchors, an Ag-Ti 2.55 A energy/force frame and an energy label for the 2.3871 A force-only frame; kept 2.70 A same-motif holdout.',
    'v11': 'Energy-focused model; packaged split/data manifest is the V10 23/2/3 dataset. Exact isolated V10-to-V11 recipe change is undocumented.',
    'v12': 'Replaced four Ag-Si/Ag-C force-only LCAO frames with same-geometry PW-PBE energy/force labels.',
    'v13': 'Added three verified PW-PBE distance labels: AgC 2.30 A, AgSi 2.60 A, AgTi 2.60 A.',
    'v14': 'Added five PW-PBE training labels: two lateral-registry holdouts, AgTi probe, AgC/AgSi strain acquisitions.',
    'v15': 'Excluded four provenance-unknown force-only LCAO frames plus one unresolved AgTi force-only frame; added three verified residual-shell PW-PBE labels.',
    'v16': 'Added three wide-span registry and three wide-span distance PW-PBE labels; corrected metadata counts without changing train/valid/test bytes.',
}

ROLE_NOTES = {
    'v9': 'Retrospective scores inherit the old test roles recorded in the V12 comparison; original V9 split membership is not fully reconstructed.',
    'v10': 'Retrospective scores on later registry references are not V10-era tests; elemental reference meshes differ from interface meshes.',
    'v11': 'Same retrospective scope as V10; split files copied from V10 and exact V11 training change is unknown.',
    'v12': 'The three frozen regression frames and five independent v13 geometries have different roles; midpoint checks share distance-scan motifs.',
    'v13': 'Independent registry/acquisition screen is narrow and related to existing cluster motifs; earlier V12 validations were promoted into training.',
    'v14': 'Three V14 holdouts; AgTi holdout duplicates an AgTi training geometry, so it is not independent geometric validation.',
    'v15': 'Three post-freeze registry references; each is a small perturbation of an existing parent motif, not independent morphology.',
    'v16': 'Three frozen local-registry references were screened against split inputs; same six-probe inference is retrospective, not a new blind test.',
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def dump_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def rmse_rows(error):
    norms = np.linalg.norm(error, axis=1)
    return float(np.sqrt(np.mean(norms ** 2))), float(norms.max()), norms


def summarize_atoms(atoms, fref, fpred, energy_ref, energy_pred):
    symbols = np.asarray(atoms.get_chemical_symbols())
    ids = np.asarray(atoms.arrays['lammps_id'], dtype=int)
    central = np.flatnonzero(np.asarray(atoms.arrays['central_pair'], dtype=bool))
    require(len(central) == 2, 'Expected exactly two marked pair atoms')
    ia = next(int(i) for i in central if symbols[i] == 'Ag')
    ix = next(int(i) for i in central if symbols[i] != 'Ag')
    direction = atoms.positions[ix] - atoms.positions[ia]
    direction /= np.linalg.norm(direction)
    dist = atoms.get_all_distances(mic=True)
    local_distance = np.min(dist[:, central], axis=1)
    local = local_distance <= CUTOFF
    error = np.asarray(fpred) - np.asarray(fref)
    force_rmse, force_max, norms = rmse_rows(error)
    ref_sep = float(np.dot(fref[ix] - fref[ia], direction))
    pred_sep = float(np.dot(fpred[ix] - fpred[ia], direction))
    agproj = float(np.dot(error[ia], direction))
    xproj = float(np.dot(error[ix], direction))

    def group(mask):
        mask = np.asarray(mask, dtype=bool)
        return {'atoms': int(mask.sum()),
                'force_vector_RMSE_eV_A': float(np.sqrt(np.mean(norms[mask] ** 2))) if mask.any() else None,
                'force_error_SSE_share': float(np.sum(norms[mask] ** 2) / np.sum(norms ** 2)) if np.sum(norms ** 2) else 0.0}

    top = []
    for i in np.argsort(norms)[::-1][:10]:
        near = [int(j) for j in np.argsort(dist[i]) if j != i][:4]
        top.append({'atom_id': int(ids[i]), 'element': str(symbols[i]), 'error_norm_eV_A': float(norms[i]),
                    'error_xyz_eV_A': error[i].tolist(), 'within_4p5A': bool(local[i]),
                    'nearest_4_MIC_neighbors': [{'atom_id': int(ids[j]), 'element': str(symbols[j]),
                                                  'distance_A': float(dist[i, j])} for j in near]})
    return {
        'energy_abs_error_meV_atom': abs(float(energy_pred) - float(energy_ref)) * 1000 / len(atoms)
        if energy_pred is not None else None,
        'energy_prediction_eV_cell': float(energy_pred) if energy_pred is not None else None,
        'force_vector_RMSE_eV_A': force_rmse,
        'force_vector_max_error_eV_A': force_max,
        'separating_force_reference_eV_A': ref_sep,
        'separating_force_prediction_eV_A': pred_sep,
        'separating_force_signed_error_eV_A': pred_sep - ref_sep,
        'separating_force_abs_error_eV_A': abs(pred_sep - ref_sep),
        'marked_pair': {'Ag_id': int(ids[ia]), 'X_id': int(ids[ix]), 'X_element': str(symbols[ix]),
                        'Ag_longitudinal_error_contribution_eV_A': -agproj,
                        'X_longitudinal_error_contribution_eV_A': xproj,
                        'sum_longitudinal_error_eV_A': -agproj + xproj,
                        'Ag_transverse_error_norm_eV_A': float(np.linalg.norm(error[ia] - agproj * direction)),
                        'X_transverse_error_norm_eV_A': float(np.linalg.norm(error[ix] - xproj * direction)),
                        'Ag_vector_error_norm_eV_A': float(norms[ia]), 'X_vector_error_norm_eV_A': float(norms[ix])},
        'regions': {
            'all': group(np.ones(len(atoms), dtype=bool)), 'within_4p5A': group(local),
            'outside_4p5A': group(~local), 'Ag_within_4p5A': group((symbols == 'Ag') & local),
            'framework_within_4p5A': group((symbols != 'Ag') & local),
            'Ag_outside_4p5A': group((symbols == 'Ag') & (~local)),
            'framework_outside_4p5A': group((symbols != 'Ag') & (~local)),
        },
        'by_element': {el: group(symbols == el) for el in sorted(set(symbols))},
        'top_10_atoms': top,
    }


def build_ledger(repo, evidence):
    import re
    root = repo / PROJECT
    v13_compare_path = root / 'mace_periodic_v13_interface_energy/results/v9_v10_v11_v12_v13_comparison.json'
    v14_compare_path = root / 'mace_periodic_v14_interface_energy/results/v14_independent_holdout_comparison.json'
    v14_assessment_path = root / 'mace_periodic_v14_interface_energy/results/v14_validation_assessment.json'
    v15_selection_path = root / 'mace_periodic_v15_interface_energy/selected_model_record.json'
    v16_selection_path = root / 'mace_periodic_v16_interface_energy/selected_model_record.json'
    v13_compare = json.loads(evidence(repo / v13_compare_path).read_text())
    v14_compare = json.loads(evidence(repo / v14_compare_path).read_text())
    v14_assessment = json.loads(evidence(repo / v14_assessment_path).read_text())
    v15_selection = json.loads(evidence(repo / v15_selection_path).read_text())
    v16_selection = json.loads(evidence(repo / v16_selection_path).read_text())
    split_info = {}
    for version, directory in [('v11', 'mace_periodic_v11_energyfocus'), ('v12', 'mace_periodic_v12_interface_energy'),
                               ('v13', 'mace_periodic_v13_interface_energy'), ('v14', 'mace_periodic_v14_interface_energy'),
                               ('v15', 'mace_periodic_v15_interface_energy'), ('v16', 'mace_periodic_v16_interface_energy')]:
        path = root / directory / 'data/dataset_manifest.json'
        data = json.loads(evidence(path).read_text())
        split_info[version] = {'manifest_path': str(path.relative_to(repo)), 'manifest_sha256': sha(path),
                               'splits': data.get('split_sizes', data.get('splits')),
                               'energy_labels': data.get('energy_label_counts', data.get('energy_frame_counts')),
                               'force_labels': data.get('force_label_counts', data.get('force_frame_counts')),
                               'energy_composition_rank': data.get('energy_composition_matrix', data.get('training_energy_composition_matrix', {})).get('rank')}
    train_paths = {
        12: (PROJECT / 'mace_periodic_v12_interface_energy/run_v12_training.sh', PROJECT / 'mace_periodic_v12_interface_energy/train.stdout'),
        13: (PROJECT / 'mace_periodic_v13_interface_energy/run_v13_training.sh', PROJECT / 'mace_periodic_v13_interface_energy/train.stdout'),
        14: (PROJECT / 'mace_periodic_v14_interface_energy/run_v14_training.sh', PROJECT / 'mace_periodic_v14_interface_energy/train.stdout'),
        15: (PROJECT / 'mace_periodic_v15_interface_energy/run_v15_training.sh', PROJECT / 'mace_periodic_v15_interface_energy/training_01/train.stdout'),
        16: (PROJECT / 'mace_periodic_v16_interface_energy/run_v16_training.sh', PROJECT / 'mace_periodic_v16_interface_energy/training_01/train.stdout'),
    }
    rows = []
    settings = 'unknown: training command/log not preserved'
    for v in range(9, 17):
        version = f'v{v}'
        model_path = MODEL_FILES[version]
        hash_source = 'v13 comparison model manifest' if v <= 13 else (
            'v14 comparison model manifest' if v == 14 else 'window-A selected model record')
        if v <= 13:
            expected_hash = v13_compare['models'][version]['sha256']
            hash_path = root / 'mace_periodic_v13_interface_energy/results/v9_v10_v11_v12_v13_comparison.json'
        elif v == 14:
            expected_hash = v14_compare['models'][version]['sha256']
            hash_path = root / 'mace_periodic_v14_interface_energy/results/v14_independent_holdout_comparison.json'
        else:
            selection = v15_selection if v == 15 else v16_selection
            expected_hash = selection['model_sha256']
            hash_path = root / f'mace_periodic_{version}_interface_energy/selected_model_record.json'
        actual_hash = sha(evidence(repo / model_path))
        require(actual_hash == expected_hash, version + ' hash differs from historical selection record')
        if v <= 11:
            if v == 9:
                sizes, energies, forces, manifest_path = 'unknown', 'unknown', 'unknown', 'not preserved in handoff'
            else:
                ds = split_info['v11']
                sizes = ds['splits']; energies = ds['energy_labels']; forces = ds['force_labels']; manifest_path = ds['manifest_path']
            selected = 'unknown'; complete = 'unknown'; params = 'unknown; no per-version launch/log evidence'
            train_script_path = train_log_path = ''
            train_script_hash = train_log_hash = ''
        else:
            ds = split_info[version]
            sizes, energies, forces, manifest_path = ds['splits'], ds['energy_labels'], ds['force_labels'], ds['manifest_path']
            train_script, train_log = train_paths[v]
            script_file, log_file = evidence(repo / train_script), evidence(repo / train_log)
            train_script_path, train_log_path = str(train_script), str(train_log)
            train_script_hash, train_log_hash = sha(script_file), sha(log_file)
            log_text = log_file.read_text()
            chosen = re.findall(r'Loaded Stage one model from epoch (\d+) for evaluation', log_text)
            epochs = [int(x) for x in re.findall(r'INFO: Epoch (\d+):', log_text)]
            require(epochs == list(range(80)) and 'INFO: Done' in log_text and chosen, version + ' training completion/selection log incomplete')
            selected = int(chosen[-1])
            complete = '80 epochs0-79 and Done verified'
            if v == 15:
                require(selected == v15_selection['selected_epoch_zero_based'], 'V15 selected epoch/log mismatch')
                complete += '; launcher exit 0'
            elif v == 16:
                require(selected == v16_selection['selected_epoch_zero_based'], 'V16 selected epoch/log mismatch')
                complete += '; launcher exit receipt unavailable'
            elif v == 14:
                require(selected == int(v14_assessment['training']['selected_checkpoint'].split('_epoch-')[-1].split('.')[0]), 'V14 selected epoch/log mismatch')
            params = 'foundation mace-mp-0b3-medium; seed45; CPU; batch1; max80; lr1e-4; wd5e-7; energy_weight100; forces_weight1000; E0s estimated'
        rows.append({
            'version': version, 'model_path': str(model_path), 'model_sha256': actual_hash,
            'model_hash_evidence': hash_source, 'hash_evidence_path': str(hash_path.relative_to(repo)),
            'dataset_manifest': manifest_path, 'splits_train_valid_test': sizes,
            'energy_label_counts': energies, 'force_label_counts': forces,
            'selected_epoch_zero_based': selected, 'training_completion': complete, 'settings': params,
            'training_script': train_script_path, 'training_script_sha256': train_script_hash,
            'training_log': train_log_path, 'training_log_sha256': train_log_hash,
            'data_change': DATA_NOTES[version], 'evaluation_role_limit': ROLE_NOTES[version],
        })
    return rows, split_info


def load_references(repo, record):
    from ase.io import read
    refs = {}
    for version, archive_name, role in [('v15', 'pbe_interface_v15_registry_holdouts', 'v15_blind_registry_holdout'),
                                         ('v16', 'pbe_interface_v16_blind_holdouts', 'v16_blind_local_registry_holdout')]:
        archive = PROJECT / archive_name
        manifest = json.loads(record(repo / archive / 'archive_manifest.json').read_text())
        input_manifest = json.loads(record(repo / archive / 'input_manifest.json').read_text())
        evaluation = json.loads(record(repo / PROJECT / f'mace_periodic_{version}_interface_energy/evaluation_01/evaluation.json').read_text())
        input_by_label = {x['label']: x for x in input_manifest['records']}
        eval_rows = {x['label']: x for x in evaluation['per_structure'] if x['role'] == 'fresh_local_registry'}
        require(manifest.get('complete') is True and len(manifest['records']) == 3, version + ' archive incomplete')
        require(set(input_by_label) == set(eval_rows) == {f'{k}_registry_holdout_{version}_01' for k in KINDS}, version + ' labels differ across input/evaluation')
        for row in manifest['records']:
            label = row['label']; require(row['status'] == 'PASS' and row['role'] == role and all(row['checks'].values()), 'Archive failed: ' + label)
            inp = input_by_label[label]; folder = archive / 'calculations' / label
            refpath = folder / f'{label}_PW_PBE.extxyz'
            inp_path = archive / inp['input']
            require(sha(record(repo / refpath)) == row['sha256'][refpath.name] == eval_rows[label]['reference_sha256'], 'Reference hash mismatch: ' + label)
            require(sha(record(repo / inp_path)) == inp['input_sha256'] == row['input_sha256'], 'Input hash mismatch: ' + label)
            for member, digest in row['sha256'].items(): require(sha(record(repo / folder / member)) == digest, 'Archive member changed: ' + label + '/' + member)
            atoms, source = read(repo / refpath), read(repo / inp_path)
            require(np.array_equal(atoms.numbers, source.numbers) and np.array_equal(atoms.pbc, source.pbc)
                    and np.array_equal(atoms.arrays['lammps_id'], source.arrays['lammps_id'])
                    and np.array_equal(atoms.arrays['central_pair'], source.arrays['central_pair'])
                    and np.allclose(atoms.positions, source.positions, atol=1e-12, rtol=0)
                    and np.allclose(atoms.cell.array, source.cell.array, atol=1e-12, rtol=0), 'Reference/input geometry mismatch: ' + label)
            refs[label] = {'atoms': atoms, 'reference_path': str(refpath), 'reference_sha256': row['sha256'][refpath.name],
                           'input_path': str(inp_path), 'input_sha256': inp['input_sha256'], 'source_holdout': version,
                           'role_at_first_scoring': role, 'energy_ref': float(atoms.info['PW_PBE_energy_eV']),
                           'forces_ref': np.asarray(atoms.arrays['PW_PBE_forces'], dtype=float)}
    return refs


def main_run(repo, out):
    from ase.io import read
    import torch
    from mace.calculators import MACECalculator

    repo = repo.resolve(); accessed = {}
    def record(path):
        p = Path(path).resolve(); require(p.is_relative_to(repo) and p.is_file(), 'Missing/outside evidence: ' + str(p))
        accessed[str(p.relative_to(repo))] = sha(p); return p

    ledger, split_info = build_ledger(repo, record)
    refs = load_references(repo, record)
    previous = json.loads(record(repo / PREVIOUS / 'matched_comparison.json').read_text())
    previous_atom_csv = record(repo / PREVIOUS / 'per_atom_errors.csv')
    previous_dir = repo / PREVIOUS
    # Reuse the completed V15/V16 matched inference (same model hashes, refs, runtime, and metrics).
    reused = {(x['model'], x['label']): x for x in previous['per_model_per_geometry']}
    reuse_sha = previous['models']
    version_status = {}
    for version in ('v15', 'v16'):
        require(reuse_sha[version]['model_sha256'] == sha(repo / MODEL_FILES[version]), version + ' reuse model hash mismatch')
        version_status[version] = {'status': 'reused_identical_prior_matched_inference',
                                   'model_sha256': reuse_sha[version]['model_sha256'],
                                   'selected_epoch': reuse_sha[version]['selected_epoch']}
    # Reuse prior V14 scores on the three V15 geometries; only the three V16 geometries are missing.
    v14_decomp_path = repo / V14_PREVIOUS / 'force_decomposition.json'
    v14_scores_path = repo / V14_PREVIOUS / 'common_geometry_scores.csv'
    v14_decomp = json.loads(record(v14_decomp_path).read_text())
    v14_scores = {}
    with record(v14_scores_path).open() as f:
        for row in csv.DictReader(f): v14_scores[row['label']] = row
    v14_sha_line = {}
    sums_path = repo / 'research/mace-v12-transfer/coordination/reports/window-b/force_separation_driver_diagnosis/SHA256SUMS.txt'
    for line in record(sums_path).read_text().splitlines():
        digest, path = line.split('  ', 1); v14_sha_line[path] = digest
    for rel in ('analysis_final/force_decomposition.json', 'analysis_final/common_geometry_scores.csv'):
        require(sha(repo / V14_PREVIOUS / Path(rel).name) == v14_sha_line.get(rel), 'Prior V14/V15 matched evidence hash mismatch: ' + rel)
    v14_old = {x['label']: x for x in v14_decomp['structures']}
    v14_reused = {}
    for label in refs:
        if not label.endswith('_v15_01'): continue
        structure = v14_old[label]; model = structure['models']['v14']; score = v14_scores[label]
        require(model['model_sha256'] == sha(repo / MODEL_FILES['v14']), 'Prior V14 model hash mismatch')
        require(score['reference_sha256'] == refs[label]['reference_sha256'], 'Prior V14 reference hash mismatch')
        v14_reused[label] = (structure, model, score)

    torch.set_default_dtype(torch.float64); torch.set_num_threads(1)
    runtime = {'torch': torch.__version__, 'mace-torch': importlib.metadata.version('mace-torch'),
               'ase': importlib.metadata.version('ase'), 'threads': torch.get_num_threads(), 'device': 'cpu', 'dtype': 'float64'}
    predictions = {}
    expected_from_v13 = json.loads(record(repo / PROJECT / 'mace_periodic_v13_interface_energy/results/v9_v10_v11_v12_v13_comparison.json').read_text())['models']
    expected_from_v14 = json.loads(record(repo / PROJECT / 'mace_periodic_v14_interface_energy/results/v14_independent_holdout_comparison.json').read_text())['models']
    for version in ('v9', 'v10', 'v11', 'v12', 'v13', 'v14'):
        model_path = record(repo / MODEL_FILES[version])
        expected = expected_from_v13[version]['sha256'] if version != 'v14' else expected_from_v14[version]['sha256']
        require(sha(model_path) == expected, version + ' published model hash mismatch')
        try:
            loaded = torch.load(model_path, map_location='cpu', weights_only=False).double().eval()
            calculator = MACECalculator(models=[loaded], device='cpu', default_dtype='float64')
            rows = {}
            for label, ref in refs.items():
                atoms = ref['atoms'].copy(); atoms.calc = calculator
                ep = None
                try: ep = float(atoms.get_potential_energy())
                except Exception: pass
                fp = np.asarray(atoms.get_forces(), dtype=float)
                require(fp.shape == ref['forces_ref'].shape and np.isfinite(fp).all(), 'Bad force output: ' + version + '/' + label)
                if ep is not None: require(np.isfinite(ep), 'Nonfinite energy output: ' + version + '/' + label)
                rows[label] = summarize_atoms(ref['atoms'], ref['forces_ref'], fp, ref['energy_ref'], ep)
                rows[label].update({'model_version': version, 'label': label, 'source_holdout': ref['source_holdout'],
                                    'role_at_first_scoring': ref['role_at_first_scoring'], 'atoms': len(ref['atoms']),
                                    'formula': ref['atoms'].get_chemical_formula(), 'reference_sha256': ref['reference_sha256'],
                                    'model_sha256': sha(model_path), 'calculation': 'new retrospective one-thread CPU inference'})
            predictions[version] = rows
            version_status[version] = {'status': 'inferred', 'model_sha256': sha(model_path), 'model_class': type(loaded).__name__}
            del calculator, loaded
        except Exception as exc:
            version_status[version] = {'status': 'unrun_incompatible_or_failed', 'model_sha256': sha(model_path),
                                       'reason': type(exc).__name__ + ': ' + str(exc)}
            predictions[version] = {}

    # Reuse V14's already verified predictions for its three V15 references.
    predictions['v14'] = predictions.get('v14', {})
    for label, (structure, model, score) in v14_reused.items():
        region_parts = model['error_regions']
        local_parts = [region_parts[k] for k in ('marked_pair', 'Ag_neighbors_local', 'framework_shell_local')]
        local_atoms = sum(x['atoms'] for x in local_parts)
        local_rmse = float(np.sqrt(sum(x['atoms'] * x['vector_RMSE_eV_A'] ** 2 for x in local_parts) / local_atoms))
        local_share = float(sum(x['force_error_SSE_share'] for x in local_parts))
        result = {
            'model_version': 'v14', 'label': label, 'source_holdout': 'v15',
            'role_at_first_scoring': 'v15_blind_registry_holdout; subsequently scored retrospectively by V14/V15',
            'atoms': structure['atoms'], 'formula': structure['formula'], 'reference_sha256': score['reference_sha256'],
            'model_sha256': model['model_sha256'], 'calculation': 'reused prior matched V14/V15 inference with verified hashes',
            'energy_abs_error_meV_atom': abs(model['energy_error_meV_atom_signed']),
            'energy_prediction_eV_cell': model['energy_prediction_eV_cell'],
            'force_vector_RMSE_eV_A': model['force_vector_RMSE_eV_A'],
            'force_vector_max_error_eV_A': model['force_vector_max_error_eV_A'],
            'separating_force_reference_eV_A': model['separation_reference_eV_A'],
            'separating_force_prediction_eV_A': model['separation_prediction_eV_A'],
            'separating_force_signed_error_eV_A': model['separation_error_signed_eV_A'],
            'separating_force_abs_error_eV_A': model['separation_error_abs_eV_A'],
            'marked_pair': structure['marked_pair'], 'regions': {**region_parts,
                'within_4p5A': {'atoms': local_atoms, 'force_vector_RMSE_eV_A': local_rmse, 'force_error_SSE_share': local_share},
                'outside_4p5A': {**region_parts['outside_4p5A'], 'force_vector_RMSE_eV_A': region_parts['outside_4p5A']['vector_RMSE_eV_A']},
            },
            'top_10_atoms': model['top10_atom_errors'],
        }
        predictions['v14'][label] = result

    # Reuse all V15/V16 matched outputs and keep the original probe roles visible.
    for version in ('v15', 'v16'):
        predictions[version] = {}
        for label, ref in refs.items():
            source = reused[(version, label)]
            require(source['reference_sha256'] == ref['reference_sha256'], 'Reused reference hash mismatch: ' + version + '/' + label)
            predictions[version][label] = {**source, 'calculation': 'reused exact frozen V15/V16 matched inference',
                                           'role_at_first_scoring': ref['role_at_first_scoring']}

    rows = []
    for version in [f'v{x}' for x in range(9, 17)]:
        for label, ref in refs.items():
            if label not in predictions.get(version, {}):
                rows.append({'model_version': version, 'label': label, 'source_holdout': ref['source_holdout'],
                             'role_at_first_scoring': ref['role_at_first_scoring'], 'atoms': len(ref['atoms']),
                             'formula': ref['atoms'].get_chemical_formula(), 'reference_sha256': ref['reference_sha256'],
                             'status': version_status.get(version, {}).get('status', 'missing'), 'force_vector_RMSE_eV_A': None})
            else:
                item = predictions[version][label]
                rows.append({'model_version': version, 'label': label, 'source_holdout': item['source_holdout'],
                             'role_at_first_scoring': item['role_at_first_scoring'], 'atoms': item['atoms'],
                             'formula': item['formula'], 'reference_sha256': item['reference_sha256'], 'status': 'scored',
                             'energy_abs_error_meV_atom': item.get('energy_abs_error_meV_atom'),
                             'force_vector_RMSE_eV_A': item.get('force_vector_RMSE_eV_A'),
                             'force_vector_max_error_eV_A': item.get('force_vector_max_error_eV_A'),
                             'separating_force_reference_eV_A': item.get('separating_force_reference_eV_A'),
                             'separating_force_prediction_eV_A': item.get('separating_force_prediction_eV_A'),
                             'separating_force_abs_error_eV_A': item.get('separating_force_abs_error_eV_A'),
                             'local_force_RMSE_eV_A': item.get('regions', {}).get('within_4p5A', {}).get('force_vector_RMSE_eV_A'),
                             'outside_force_RMSE_eV_A': item.get('regions', {}).get('outside_4p5A', {}).get('force_vector_RMSE_eV_A'),
                             'top_error_atoms': ';'.join(f"{x.get('atom_id',x.get('id'))}:{x.get('element')}:{x.get('error_norm_eV_A',x.get('force_error_norm_eV_A')):.4f}" for x in item.get('top_10_atoms', [])[:3])})

    aggregate = []
    for version in [f'v{x}' for x in range(9, 17)]:
        if not predictions.get(version): continue
        for kind in KINDS:
            selected = [predictions[version][label] for label in refs if label.startswith(kind + '_') and label in predictions[version]]
            if not selected: continue
            aggregate.append({'model_version': version, 'interface': kind, 'structures': len(selected),
                              'atoms': sum(x['atoms'] for x in selected),
                              'pooled_force_vector_RMSE_eV_A': float(np.sqrt(sum(x['atoms'] * x['force_vector_RMSE_eV_A'] ** 2 for x in selected) / sum(x['atoms'] for x in selected))),
                              'energy_MAE_meV_atom': float(np.mean([x['energy_abs_error_meV_atom'] for x in selected if x.get('energy_abs_error_meV_atom') is not None])) if any(x.get('energy_abs_error_meV_atom') is not None for x in selected) else None,
                              'max_separation_abs_error_eV_A': max(x['separating_force_abs_error_eV_A'] for x in selected)})

    # Frozen score reproduction on already published common geometries.
    v14_reproduction = {label: {'reused_hash_checked': True, 'source': 'prior matched V14/V15 review',
                                'reference_sha256': score['reference_sha256'],
                                'force_RMSE_eV_A': v14_reused[label][1]['force_vector_RMSE_eV_A'],
                                'separation_abs_error_eV_A': v14_reused[label][1]['separation_error_abs_eV_A']}
                        for label, (_, _, score) in v14_reused.items()}
    original_reports = {
        'v9_v13_comparison': PROJECT / 'mace_periodic_v13_interface_energy/results/v9_v10_v11_v12_v13_comparison.json',
        'v13_independent': PROJECT / 'mace_periodic_v13_interface_energy/results/v13_independent_holdout_comparison.json',
        'v14_assessment': PROJECT / 'mace_periodic_v14_interface_energy/results/v14_validation_assessment.json',
        'v15_evaluation': PROJECT / 'mace_periodic_v15_interface_energy/evaluation_01/evaluation.json',
        'v16_evaluation': PROJECT / 'mace_periodic_v16_interface_energy/evaluation_01/evaluation.json',
    }
    evidence_hashes = {name: sha(record(repo / path)) for name, path in original_reports.items()}
    additional_history = [
        Path('research/mace-v12-transfer/periodic_interface_v4/MACE_VERSION_HISTORY.md'),
        Path('research/mace-v12-transfer/coordination/reports/window-b/v9_v10_v11_force_separation_history.csv'),
        Path('research/mace-v12-transfer/coordination/reports/window-b/v9-v15_force_separation_history.md'),
        Path('research/mace-v12-transfer/periodic_interface_v4/V15_FORCE_SEPARATION_PLAN.md'),
        Path('research/mace-v12-transfer/coordination/reports/window-b/v15_reference_convention_review/report.md'),
        Path('research/mace-v12-transfer/coordination/reports/window-b/v15_inherited_method_provenance_trace/report.md'),
        Path('research/mace-v12-transfer/coordination/reports/window-b/v15_metric_implementation_review/report.md'),
        Path('research/mace-v12-transfer/coordination/reports/window-b/v15_force_error_localization/report.md'),
        Path('research/mace-v12-transfer/coordination/reports/window-b/v16_dataset_integration_audit/report.md'),
        Path('research/mace-v12-transfer/coordination/reports/window-b/v16_matched_error_and_metadata_review/report.md'),
    ]
    history_hashes = {str(p): sha(record(repo / p)) for p in additional_history}
    ledger, split_info = build_ledger(repo, record)
    output = {
        'scope': 'Frozen V9-V16 checkpoints scored or reused on the same six already-scored V15/V16 registry references. Retrospective analysis only; no new reference labels, selection, tuning, training, DFT, MD or TTM.',
        'metric': 'Per-atom force-vector RMSE=sqrt(mean_i(sum_xyz(deltaF_i^2))); marked separation uses the original direct Ag-to-X vector; local means MIC <=4.5 A from either marked atom.',
        'runtime': runtime,
        'references': {label: {'source_holdout': v['source_holdout'], 'role_at_first_scoring': v['role_at_first_scoring'],
                               'reference_path': v['reference_path'], 'reference_sha256': v['reference_sha256'],
                               'input_path': v['input_path'], 'input_sha256': v['input_sha256'],
                               'atoms': len(v['atoms']), 'formula': v['atoms'].get_chemical_formula()} for label, v in refs.items()},
        'model_status': version_status,
        'per_model_per_geometry': [predictions[v][label] for v in [f'v{x}' for x in range(9, 17)] for label in refs if label in predictions.get(v, {})],
        'per_interface_aggregate': aggregate,
        'score_reproduction': {'v14_previous_matched_outputs_reused': v14_reproduction,
                               'v15_v16_own_scores_reproduced_and_frozen_results_reused': True,
                               'v9_v13_on_these_six_are_new_retrospective_common_scores_not_original_score_reproductions': True},
        'version_ledger': ledger,
        'split_manifest_evidence': split_info,
        'evaluation_evidence_sha256': evidence_hashes,
        'historical_review_sha256': history_hashes,
        'input_sha256': accessed,
        'limitations': ['Cross-version inference on this union is retrospective; V15/V16 probes retain their original role and are not fresh tests.',
                        'The V15 probes are small perturbations of training-parent motifs; V14 AgTi scored geometry duplicates a training frame.',
                        'V9 energy calibration is not established; do not interpret its predicted energy error as comparable to V10+.',
                        'V9-V11 original settings and some dataset memberships are not preserved; unknown fields remain unknown.',
                        'A six-geometry set does not represent independent morphologies, extended periodic interfaces, thermal disorder, liquid Ag or pressure.'],
    }
    dump_json(out / 'common_geometry_scores.json', output)
    with (out / 'common_geometry_scores.csv').open('w', newline='') as f:
        fields = ['model_version','label','source_holdout','role_at_first_scoring','atoms','formula','reference_sha256','status',
                  'energy_abs_error_meV_atom','force_vector_RMSE_eV_A','force_vector_max_error_eV_A',
                  'separating_force_reference_eV_A','separating_force_prediction_eV_A','separating_force_abs_error_eV_A',
                  'local_force_RMSE_eV_A','outside_force_RMSE_eV_A','top_error_atoms']
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    dump_json(out / 'version_ledger.json', {'versions': ledger})
    with (out / 'version_ledger.csv').open('w', newline='') as f:
        fields = list(ledger[0]); writer = csv.DictWriter(f, fieldnames=fields, lineterminator='\n'); writer.writeheader(); writer.writerows(ledger)
    print(json.dumps({'model_status': version_status, 'per_interface_aggregate': aggregate}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    repo = args.repo_root.resolve(); out = args.output_dir.resolve()
    require(out.is_relative_to(repo / REPORT) and out != repo / REPORT, 'Output must be a child of the assigned report directory')
    require(not out.exists() or (out.is_dir() and not any(out.iterdir())), 'Nonempty output directory refused')
    out.mkdir(parents=True, exist_ok=True)
    try:
        main_run(repo, out)
    except Exception as exc:
        dump_json(out / 'diagnostics_failure.json', {'error': type(exc).__name__, 'message': str(exc)})
        raise


if __name__ == '__main__':
    main()
