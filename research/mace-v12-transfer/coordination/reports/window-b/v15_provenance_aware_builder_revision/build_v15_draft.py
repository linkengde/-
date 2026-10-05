#!/usr/bin/env python3
"""Build an isolated V15 draft; never use V15 blind labels or write production data."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from ase.io import read, write

REPORT = Path('research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision')
PROJECT = Path('research/mace-v12-transfer/periodic_interface_v4')
EXCLUDED = {'AgC_primary_2.0A_LCAO', 'AgC_primary_2.4A_LCAO', 'AgTi_2p1A_LCAO', 'AgTi_2p7A_LCAO'}
RESIDUALS = {'AgC_residual_shell_v15_01', 'AgSi_residual_shell_v15_01', 'AgTi_residual_shell_v15_01'}
BLIND = {'AgC_registry_holdout_v15_01', 'AgSi_registry_holdout_v15_01', 'AgTi_registry_holdout_v15_01'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')


def ids(atoms):
    if 'lammps_id' not in atoms.arrays:
        return None
    values = [int(x) for x in atoms.arrays['lammps_id']]
    require(len(set(values)) == len(atoms), 'Duplicate persistent atom IDs')
    return {x: i for i, x in enumerate(values)}


def comparison(a, b):
    if len(a) != len(b) or Counter(a.numbers) != Counter(b.numbers):
        return None
    ma, mb = ids(a), ids(b)
    if ma is not None and mb is not None:
        if set(ma) != set(mb) or any(a[ma[x]].symbol != b[mb[x]].symbol for x in ma):
            return None
        order = sorted(ma)
        ia, ib = [ma[x] for x in order], [mb[x] for x in order]
        mapping = 'persistent_ID_and_element'
    elif np.array_equal(a.numbers, b.numbers):
        ia = ib = list(range(len(a)))
        mapping = 'ordered_species_fallback_no_persistent_ID_claim'
    else:
        return None
    delta = a.positions[ia] - b.positions[ib]
    centered = delta - delta.mean(axis=0)
    norms = np.linalg.norm(centered, axis=1)
    rms, maximum = float(np.sqrt(np.mean(norms**2))), float(norms.max())
    species = np.array([a[i].symbol for i in ia])
    ag, fw = species == 'Ag', species != 'Ag'
    registry_near = False
    ag_rms = framework_rms = None
    if ag.any() and fw.any():
        framework_aligned = delta - delta[fw].mean(axis=0)
        ag_rms = float(np.sqrt(np.mean(np.sum(framework_aligned[ag]**2, axis=1))))
        framework_rms = float(np.sqrt(np.mean(np.sum(framework_aligned[fw]**2, axis=1))))
        registry_near = ag_rms <= .15 + 1e-8 and framework_rms <= .15 + 1e-8
    same_cell = bool(np.allclose(a.cell, b.cell, atol=1e-8, rtol=0))
    same_pbc = bool(np.array_equal(a.pbc, b.pbc))
    exact = same_cell and maximum <= 1e-7
    near = same_cell and ((rms <= .05 + 1e-8 and maximum <= .15 + 1e-8) or registry_near)
    return {'mapping': mapping, 'cell_match': same_cell, 'pbc_match': same_pbc,
            'centered_RMS_A': rms, 'centered_max_A': maximum,
            'Ag_registry_RMS_A': ag_rms, 'framework_RMS_A': framework_rms,
            'exact_geometry': exact, 'near_geometry': near,
            'exact_geometry_and_boundary_input': exact and same_pbc,
            'independent_geometry': not (exact or near),
            'complete_calculator_input_identity': 'unknown; geometry/PBC identity does not establish method/settings identity'}


def identical_geometry(a, b, label):
    require(len(a) == len(b), label + ': atom count changed')
    require(np.array_equal(a.numbers, b.numbers), label + ': element/order changed')
    require(ids(a) is not None and ids(b) is not None, label + ': persistent IDs absent')
    require(np.array_equal(a.arrays['lammps_id'], b.arrays['lammps_id']), label + ': IDs/order changed')
    require(np.array_equal(a.pbc, b.pbc), label + ': PBC changed')
    require(np.allclose(a.cell, b.cell, atol=1e-8, rtol=0), label + ': cell changed')
    require(np.allclose(a.positions, b.positions, atol=1e-7, rtol=0), label + ': positions changed')
    for key in ('lammps_type', 'central_pair'):
        if key in a.arrays:
            require(key in b.arrays and np.array_equal(a.arrays[key], b.arrays[key]), label + ': ' + key + ' changed')


def finite_labels(a, energy_key='REF_energy', force_key='REF_forces'):
    require(energy_key in a.info and force_key in a.arrays, 'Missing energy/force label: ' + a.info.get('config_type', 'unknown'))
    energy = float(a.info[energy_key])
    forces = np.asarray(a.arrays[force_key], dtype=float)
    require(np.isfinite(energy) and forces.shape == (len(a), 3) and np.isfinite(forces).all(), 'Invalid energy/force labels')
    require(np.isfinite(a.positions).all() and np.isfinite(a.cell.array).all(), 'Nonfinite geometry')
    return energy, forces


def overlaps(left, right):
    hits = []
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            c = comparison(a, b)
            if c and (c['exact_geometry'] or c['near_geometry']):
                hits.append({'left_frame': i, 'left_config': a.info.get('config_type'),
                             'right_frame': j, 'right_config': b.info.get('config_type'), **c})
    return hits


def build(repo, out, sources, exclude_unresolved=False):
    project = repo / PROJECT
    blind_root = project / 'pbe_interface_v15_registry_holdouts'

    def record(path):
        path = path.resolve()
        require(path.is_relative_to(repo), 'Source escapes repository')
        if path.is_relative_to(blind_root):
            require(path == blind_root / 'input_manifest.json' or path.is_relative_to(blind_root / 'inputs'), 'V15 blind label/output access forbidden')
        require(path.is_file(), 'Source missing: ' + str(path))
        sources[str(path.relative_to(repo))] = digest(path)
        return path

    def load(path):
        return json.loads(record(path).read_text())

    def atoms(path):
        return read(record(path), index=':')

    # Pin all evidence paths explicitly. Never discover holdout calculation files.
    record(project / 'mace_periodic_v14_interface_energy/build_v14_dataset.py')
    record(project / 'V15_FORCE_SEPARATION_PLAN.md')
    record(repo / 'research/mace-v12-transfer/coordination/reports/window-b/v15_dataset_integrity_audit.md')
    trace_path = repo / 'research/mace-v12-transfer/coordination/reports/window-b/v15_inherited_method_provenance_trace/per_frame_ledger.json'
    trace = load(trace_path)
    trace_checksums = record(trace_path.parent / 'SHA256SUMS.txt')
    expected_trace = {name: h for h, name in (line.split('  ', 1) for line in trace_checksums.read_text().splitlines())}
    require(digest(trace_path) == expected_trace[str(trace_path.relative_to(repo))], 'Trace ledger checksum mismatch')
    require(len(trace['frames']) == 15, 'Expected exact 15-frame trace')
    traced = {(r['source_path'], r['source_frame'], r['source_sha256'], r['config_type']): r for r in trace['frames']}
    require(len(traced) == 15, 'Nonunique provenance trace identity')
    unresolved = [r for r in trace['frames'] if r['status'] == 'unresolved']
    require(len(unresolved) == 1 and unresolved[0]['config_type'] == 'AgTi_2p3871A_periodic_PW_PBE_force_only', 'Unexpected unresolved trace inventory')
    target = unresolved[0]
    require(target['source_frame'] == 13 and target['source_sha256'] == '3323cdbd3081fc8d91673662c5f83c8da93ba3b4d7a23428f6b70d0e7a4c9faf', 'Pinned unresolved source identity changed')
    matched_target = 0
    parent_root = project / 'mace_periodic_v14_interface_energy/data'
    old_manifest = load(parent_root / 'dataset_manifest.json')
    original = {}
    for split in ('train', 'valid', 'test'):
        path = record(parent_root / (split + '.extxyz'))
        require(digest(path) == old_manifest['output_sha256'][split], 'V14 split hash mismatch: ' + split)
        original[split] = atoms(path)
        require(len(original[split]) == old_manifest['split_sizes'][split], 'V14 split count mismatch: ' + split)
    require(len(original['train']) == 31, 'Expected 31 V14 training frames')
    ledger, exclusion, train = [], [], []
    for i, a in enumerate(original['train']):
        config = a.info.get('config_type')
        if config in EXCLUDED:
            require('REF_energy' not in a.info and 'REF_forces' in a.arrays, 'Excluded LCAO frame label inventory unexpectedly changed')
            exclusion.append({'frame': i, 'config_type': config, 'source_path': str((parent_root / 'train.extxyz').relative_to(repo)),
                              'source_sha256': sources[str((parent_root / 'train.extxyz').relative_to(repo))],
                              'reason': 'force-only LCAO; original method/basis/boundary/convergence provenance missing; not admissible as a unified PW-PBE reference'})
            continue
        finite_labels(a)
        identity = (str((parent_root / 'train.extxyz').relative_to(repo)), i,
                    sources[str((parent_root / 'train.extxyz').relative_to(repo))], config)
        target_identity = (target['source_path'], target['source_frame'], target['source_sha256'], target['config_type'])
        if identity == target_identity:
            matched_target += 1
            if exclude_unresolved:
                exclusion.append({'frame': i, 'config_type': config, 'source_path': identity[0],
                                  'source_sha256': identity[2], 'reason': 'explicit exclusion: sole unresolved original method/energy-force binding',
                                  'trace_evidence_path': str(trace_path.relative_to(repo)), 'trace_evidence_sha256': digest(trace_path),
                                  'matching_policy': 'config AND source frame AND source file SHA256; exactly one match'})
                continue
        status = 'inherited_explicit_method_metadata_not_reverified' if a.info.get('source_method') else 'unknown_inherited_method'
        ledger.append({'split': 'train', 'frame': len(train), 'source_frame': i,
                       'source_path': str((parent_root / 'train.extxyz').relative_to(repo)), 'provenance_status': status})
        train.append(a.copy())
    require(matched_target == 1, 'Unresolved exclusion selector must match exactly once')
    require({x['config_type'] for x in exclusion if x['config_type'] in EXCLUDED} == EXCLUDED and sum(x['config_type'] in EXCLUDED for x in exclusion) == 4, 'Four expected LCAO exclusions not found exactly once')
    acquire = project / 'pbe_interface_v15_targeted_acquisition'
    input_manifest = load(acquire / 'input_manifest.json')
    require(len(input_manifest['records']) == 3 and {r['label'] for r in input_manifest['records']} == RESIDUALS, 'Residual manifest must contain exactly the three approved acquisitions')
    residual_checks = []
    for entry in input_manifest['records']:
        label = entry['label']
        inp = record(acquire / entry['input'])
        require(inp.parent == acquire / 'inputs', 'Residual input path is outside approved inputs')
        require(digest(inp) == entry['input_sha256'], label + ': manifest input hash mismatch')
        parent_input = record(project / entry['source'])
        require(digest(parent_input) == entry['source_sha256'], label + ': source-parent input hash mismatch')
        folder = acquire / 'calculations' / label
        verify_path, summary_path, progress_path = folder / 'verification.json', folder / 'summary.json', folder / 'progress.json'
        verify, summary, progress = load(verify_path), load(summary_path), load(progress_path)
        require(verify.get('status') == 'PASS' and verify.get('checks') and all(x is True for x in verify['checks'].values()), label + ': archive verification not PASS')
        require(verify.get('input_sha256') == entry['input_sha256'] == summary.get('source_sha256') == progress.get('source_sha256'), label + ': source provenance hashes disagree')
        require(summary.get('label') == label and verify.get('label') == label, label + ': label identity mismatch')
        require(summary.get('scf_converged') is True and progress.get('status') == 'complete' and progress.get('scf_converged') is True, label + ': SCF not complete/converged')
        require(summary.get('mpi_ranks') == 4 and progress.get('mpi_ranks') == 4, label + ': expected four MPI ranks')
        method = summary.get('method', {})
        require(method.get('xc') == 'PBE' and method.get('basis') == 'plane wave' and method.get('cutoff_eV') == 500 and method.get('kpts') == [1, 1, 1], label + ': unexpected reference method')
        for name, expected in verify['sha256'].items():
            require(Path(name).name == name, label + ': unsafe archive member name')
            member = record(folder / name)
            require(digest(member) == expected, label + ': archive member hash mismatch: ' + name)
        output = folder / (label + '_PW_PBE.extxyz')
        require(output.name in verify['sha256'] and 'summary.json' in verify['sha256'] and 'progress.json' in verify['sha256'], label + ': compact archive hash coverage incomplete')
        input_frames, output_frames = atoms(inp), atoms(output)
        require(len(input_frames) == len(output_frames) == 1, label + ': expected one fixed-geometry frame')
        a = output_frames[0]
        identical_geometry(input_frames[0], a, label)
        energy, forces = finite_labels(a, 'PW_PBE_energy_eV', 'PW_PBE_forces')
        require(summary.get('atoms') == len(a) and summary.get('formula') == a.get_chemical_formula(), label + ': summary atom inventory differs')
        require(np.isclose(energy, summary['energy_eV_cell'], atol=1e-8, rtol=0), label + ': summary energy differs')
        require(not any((comparison(a, old) or {}).get('exact_geometry') for old in train), label + ': exact training geometry duplicate')
        a = a.copy(); a.calc = None
        a.info['REF_energy'] = energy; a.arrays['REF_forces'] = forces.copy()
        a.info.pop('PW_PBE_energy_eV', None); a.arrays.pop('PW_PBE_forces', None)
        a.info['config_type'] = label + '_periodic_PW_PBE_energy_force_v15_draft_train'
        a.info['source_method'] = 'GPAW ' + str(method['GPAW_version']) + ' PW-PBE 500 eV Gamma single point'
        a.info['dataset_role'] = 'v15_training_acquisition'
        ledger.append({'split': 'train', 'frame': len(train), 'source_frame': 0, 'source_path': str(output.relative_to(repo)),
                       'provenance_status': 'verified_residual_PW_PBE_archive', 'verification_path': str(verify_path.relative_to(repo)),
                       'summary_path': str(summary_path.relative_to(repo)), 'input_path': str(inp.relative_to(repo)), 'method': method})
        residual_checks.append({'label': label, 'status': 'PASS', 'input_output_IDs_elements_cell_PBC_positions_match': True,
                                'finite_labels': True, 'scf_converged': True, 'mpi_ranks': 4, 'all_declared_archive_hashes_match': True})
        train.append(a)
    require(len(train) == (29 if exclude_unresolved else 30), 'Expected 31 - 4 + 3 = 30 training frames; investigate mismatch')
    splits = {'train': train, 'valid': original['valid'], 'test': original['test']}
    held_manifest = load(blind_root / 'input_manifest.json')
    require(len(held_manifest['records']) == 3 and {r['label'] for r in held_manifest['records']} == BLIND, 'Expected all three frozen V15 blind inputs')
    held, held_records = [], []
    for entry in held_manifest['records']:
        p = record(blind_root / entry['input'])
        require(p.parent == blind_root / 'inputs', 'Holdout input outside frozen input directory')
        require(digest(p) == entry['input_sha256'], 'Frozen holdout input hash mismatch: ' + entry['label'])
        items = atoms(p); require(len(items) == 1, 'Expected one blind input geometry')
        a = items[0]
        require(a.calc is None and not any(k in a.info for k in ('energy', 'REF_energy', 'PW_PBE_energy_eV')) and not any(k in a.arrays for k in ('forces', 'REF_forces', 'PW_PBE_forces')), 'Frozen input unexpectedly contains labels')
        held.append(a); held_records.append({'label': entry['label'], 'path': str(p.relative_to(repo)), 'sha256': digest(p), 'label_content_read': False})
    blind_checks = {s: overlaps(frames, held) for s, frames in splits.items()}
    require(not blind_checks['train'] and not blind_checks['valid'], 'New V15 blind geometry overlaps train/valid, including PBC variants')
    require(not blind_checks['test'], 'New V15 blind geometry overlaps historical test')
    historical = {s: overlaps(original['train'], original[s]) for s in ('valid', 'test')}
    final_cross = {s: overlaps(train, splits[s]) for s in ('valid', 'test')}
    composition = np.asarray([[a.get_chemical_symbols().count(e) for e in ('Ag', 'C', 'Si', 'Ti')] for a in train], dtype=float)
    rank = int(np.linalg.matrix_rank(composition)); require(rank == 4, 'Energy composition matrix rank must be 4/4')
    for split in ('valid', 'test'):
        for i in range(len(splits[split])):
            ledger.append({'split': split, 'frame': i, 'source_frame': i, 'source_path': str((parent_root / (split + '.extxyz')).relative_to(repo)),
                           'provenance_status': 'inherited_explicit_method_metadata_not_reverified' if splits[split][i].info.get('source_method') else 'unknown_inherited_method'})
    for rec in ledger:
        a = splits[rec['split']][rec['frame']]
        finite_labels(a)
        rec.update({'config_type': a.info.get('config_type'), 'source_sha256': sources[rec['source_path']], 'source_method_as_stored': a.info.get('source_method'),
                    'atoms': len(a), 'formula': a.get_chemical_formula(), 'REF_energy_present': True, 'REF_forces_present': True,
                    'persistent_atom_ids': [int(x) for x in a.arrays['lammps_id']] if ids(a) is not None else None,
                    'cell_A': a.cell.tolist(), 'PBC': a.pbc.tolist(), 'source_original_settings_known': rec['provenance_status'] == 'verified_residual_PW_PBE_archive'})
    trace_matches = 0
    for rec in ledger:
        a = splits[rec['split']][rec['frame']]
        rec['source_method_as_stored'] = a.info.get('source_method')
        rec['dft_method_as_stored'] = a.info.get('dft_method')
        rec['method_declaration_evidence'] = {'path': rec['source_path'], 'frame': rec['source_frame'],
                                              'sha256': rec['source_sha256'], 'fields': ['source_method', 'dft_method']}
        key = (rec['source_path'], rec['source_frame'], rec['source_sha256'], rec['config_type'])
        recovered = traced.get(key)
        if recovered is not None:
            trace_matches += 1
            require(rec['source_method_as_stored'] == recovered['source_method_field'] and
                    rec['dft_method_as_stored'] == recovered['dft_method_declaration'], 'Raw declaration disagrees with trace')
            rec['provenance_status'] = recovered['status']
            rec['evidence_tier'] = recovered['status']
            rec['parameter_evidence'] = recovered['method_parameter_evidence']
            rec['original_run_evidence'] = {'status': 'not_verified', 'trace_path': str(trace_path.relative_to(repo)),
                  'trace_sha256': digest(trace_path), 'trace_frame': recovered['draft_train_frame'],
                  'original_output_identity': recovered['original_calculation_output_geometry_energy_force_identity'],
                  'declared_original_sources': recovered['original_calculation_claims'], 'missing_materials': recovered['missing_materials']}
            rec['energy_convention_evidence'] = {'selected_field': 'REF_energy', 'replacement_by_free_energy': False,
                  'serialized_alias_comparison': recovered['native_alias_identity'],
                  'independent_original_getter_verified': False}
        elif rec['provenance_status'] == 'verified_residual_PW_PBE_archive':
            rec['evidence_tier'] = 'verified_residual_archive'
            rec['parameter_evidence'] = {'stored_method': rec['method'], 'evidence_path': rec['summary_path'],
                                          'evidence_sha256': sources[rec['summary_path']]}
            rec['original_run_evidence'] = {'status': 'archive_verified', 'verification_path': rec['verification_path'],
                  'verification_sha256': sources[rec['verification_path']], 'summary_path': rec['summary_path'],
                  'summary_sha256': sources[rec['summary_path']]}
            rec['energy_convention_evidence'] = {'selected_field': 'PW_PBE_energy_eV -> REF_energy',
                  'replacement_by_free_energy': False, 'evidence_path': rec['source_path']}
        else:
            rec['evidence_tier'] = 'inherited_declaration_not_reverified'
            rec['parameter_evidence'] = {'source_method': rec['source_method_as_stored'], 'dft_method': rec['dft_method_as_stored'],
                                          'status': 'declaration_only; unparsed settings remain unknown'}
            rec['original_run_evidence'] = {'status': 'not_reverified', 'reason': 'No new original-run verification performed for explicitly declared historical frames'}
            rec['energy_convention_evidence'] = {'selected_field': 'REF_energy', 'replacement_by_free_energy': False,
                                                 'independent_original_getter_verified': False}
    require(trace_matches == (14 if exclude_unresolved else 15), 'Trace binding count mismatch')
    for a in train:
        a.calc = None
    write(out / 'train.extxyz', train, format='extxyz')
    for split in ('valid', 'test'):
        (out / (split + '.extxyz')).write_bytes((parent_root / (split + '.extxyz')).read_bytes())
    produced = {s: read(out / (s + '.extxyz'), index=':') for s in splits}
    for s in splits:
        require(len(produced[s]) == len(splits[s]), 'Output frame count changed')
        for before, after in zip(splits[s], produced[s]):
            finite_labels(after)
            require(np.array_equal(before.numbers, after.numbers) and np.allclose(before.positions, after.positions, atol=1e-8, rtol=0) and np.allclose(before.cell, after.cell, atol=1e-8, rtol=0) and np.array_equal(before.pbc, after.pbc), 'Output geometry changed')
            if ids(before) is not None:
                require(np.array_equal(before.arrays['lammps_id'], after.arrays['lammps_id']), 'Output atom IDs/order changed')
            require(np.isclose(before.info['REF_energy'], after.info['REF_energy'], atol=1e-10, rtol=0) and np.allclose(before.arrays['REF_forces'], after.arrays['REF_forces'], atol=1e-8, rtol=0), 'Output reference labels changed')
    for s in ('valid', 'test'):
        require(digest(out / (s + '.extxyz')) == sources[str((parent_root / (s + '.extxyz')).relative_to(repo))], 'Historical split bytes changed')
    diagnostics = {'engineering_status': 'PASS_WITH_PROVENANCE_AND_HISTORICAL_INDEPENDENCE_LIMITS', 'residual_archive_checks': residual_checks,
                   'frozen_V15_geometry_overlap_by_split': blind_checks, 'historical_V14_training_split_overlaps': historical,
                   'draft_V15_training_split_overlaps': final_cross, 'historical_AgTi_test_independent_geometry': False,
                   'blind_label_content_read': False, 'unknown_inherited_training_methods': sum(r['split'] == 'train' and r['provenance_status'] == 'unresolved' for r in ledger),
                   'physical_method_consistency': 'UNKNOWN; inherited records without original settings prevent a unified verified PW-PBE consistency claim',
                   'geometry_policy': 'Match IDs+elements; ordered-species fallback when IDs absent. Cell <=1e-8 A; remove common translation, no rotation. Exact geometry max <=1e-7 A independent of PBC; report PBC separately. Near all-atom RMS <=0.05/max <=0.15 A OR Ag/framework RMS <=0.15 A. All are review thresholds, not universal physics cutoffs.',
                   'production_dataset_modified': False, 'wide_span_proposals_included': False, 'DFT_or_training_run': False}
    manifest = {'version': 'V15 provenance-aware isolated draft; A integration pending', 'scenario': 'exclude_unresolved_29' if exclude_unresolved else 'baseline_30', 'explicit_unresolved_exclusion': exclude_unresolved, 'analyzed_repository_commit': subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip(),
                'split_sizes': {s: len(produced[s]) for s in produced}, 'energy_label_counts': {s: sum('REF_energy' in a.info for a in produced[s]) for s in produced},
                'force_label_counts': {s: sum('REF_forces' in a.arrays for a in produced[s]) for s in produced},
                'energy_composition_matrix': {'elements': ['Ag', 'C', 'Si', 'Ti'], 'rank': rank, 'matrix': composition.astype(int).tolist(), 'singular_values': np.linalg.svd(composition, compute_uv=False).tolist()},
                'excluded_frame_ledger': exclusion, 'per_frame_provenance': ledger, 'frozen_holdout_INPUTS_only': held_records,
                'sources_sha256': dict(sorted(sources.items())), 'output_sha256': {s: digest(out / (s + '.extxyz')) for s in produced},
                'valid_test_byte_identical_to_V14': True, 'provenance_status_counts_train': dict(Counter(r['provenance_status'] for r in ledger if r['split'] == 'train')),
                'evidence_tier_counts_train': dict(Counter(r['evidence_tier'] for r in ledger if r['split'] == 'train')),
                'physical_method_consistency': diagnostics['physical_method_consistency'], 'historical_split_independence': 'Historical AgTi test is not independent; retained byte-identical as regression evidence',
                'not_production_or_training_approval': True}
    save(out / 'dataset_manifest.json', manifest); save(out / 'verification_report.json', diagnostics)
    (out / 'SHA256SUMS.txt').write_text(''.join(f'{digest(p)}  {p.name}\n' for p in sorted(out.iterdir()) if p.is_file()))
    return {'split_sizes': manifest['split_sizes'], 'label_counts': manifest['energy_label_counts'], 'composition_rank': rank,
            'provenance_status_counts_train': manifest['provenance_status_counts_train'], 'engineering_status': diagnostics['engineering_status']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--exclude-unresolved-agti', action='store_true', help='Exclude exactly the pinned unresolved AgTi row; retain all 14 partial rows')
    args = parser.parse_args()
    repo, out = args.repo_root.resolve(), args.output_dir.resolve()
    require((repo / PROJECT).is_dir(), 'Invalid repository root')
    allowed = (repo / REPORT).resolve()
    require(out.is_relative_to(allowed) and out != allowed, 'Output must be an isolated child of B v15_provenance_aware_builder_revision; production paths refused')
    require(not out.exists() or out.is_dir(), 'Output exists and is not a directory')
    require(not out.exists() or not any(out.iterdir()), 'Refusing existing nonempty output directory')
    out.mkdir(parents=True, exist_ok=True)
    sources = {}
    try:
        result = build(repo, out, sources, args.exclude_unresolved_agti)
    except Exception as exc:
        save(out / 'diagnostics_failure.json', {'status': 'FAILED', 'exception': type(exc).__name__, 'blocker': str(exc), 'sources_sha256_read_before_failure': sources})
        raise
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print('Builder stopped; preserve diagnostics: ' + str(error), file=sys.stderr)
        sys.exit(1)
