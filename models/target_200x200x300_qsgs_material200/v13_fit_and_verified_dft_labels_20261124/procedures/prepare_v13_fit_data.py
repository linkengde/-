#!/workspace/venvs/gpaw/bin/python
"""Audit and stage only completed V13 DFT labels for POD fitting."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

WORK = Path('/workspace/potential-fit/round7_pod_snap/models/target_200x200x300_qsgs_material200/v13_source_and_fit_20261123')
BASE = WORK.parent / 'v13_candidate_sources_20261123'
TRAIN = BASE / 'training_candidate_manifest_v13.json'
VALID = BASE / 'validation_candidate_manifest_v13.json'
LABELS = BASE / 'dft_labels'
VALID_LABELS = BASE / 'validation_dft_labels'
FIT = WORK / 'fit_v13'
DESCRIPTOR_SOURCE = Path('/workspace/potential-fit/round7_pod_snap/models/final_v10/POD/AgTiSiC_POD.pod')


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def checked_rows(manifest_path: Path, label_root: Path, expected_count: int, role: str):
    manifest = json.loads(manifest_path.read_text())
    rows = manifest['candidates']
    if len(rows) != expected_count:
        raise RuntimeError(f'{role}: expected {expected_count} rows; found {len(rows)}')
    passed = []
    for row in rows:
        cid = row['candidate_id']
        label_json = label_root / cid / 'label.json'
        extxyz = label_root / cid / 'label.extxyz'
        if not label_json.is_file() or not extxyz.is_file():
            raise RuntimeError(f'{role}: missing converged DFT output for {cid}')
        result = json.loads(label_json.read_text())
        if result.get('run_status') != 'complete' or not result.get('scf_converged') or not result.get('finite_energy_force'):
            raise RuntimeError(f'{role}: DFT acceptance failed for {cid}: {result.get("run_status")}')
        if result.get('geometry_sha256') != row.get('geometry_file_sha256'):
            raise RuntimeError(f'{role}: geometry file hash mismatch for {cid}')
        atoms = read(extxyz)
        if len(atoms) != int(result['natoms']) or not np.isfinite(atoms.get_forces()).all():
            raise RuntimeError(f'{role}: invalid extxyz labels for {cid}')
        if not np.isfinite(float(atoms.get_potential_energy())):
            raise RuntimeError(f'{role}: missing/nonfinite fitting energy for {cid}')
        passed.append((row, result, extxyz))
    return manifest, passed


def main():
    if not TRAIN.is_file() or not VALID.is_file():
        raise SystemExit('Training or validation manifest is missing; no dataset staged.')
    train_manifest, train = checked_rows(TRAIN, LABELS, 27, 'training')
    valid_manifest, valid = checked_rows(VALID, VALID_LABELS, 6, 'validation')
    train_hashes = {r['canonical_geometry_sha256'] for r, _, _ in train}
    valid_hashes = {r['canonical_geometry_sha256'] for r, _, _ in valid}
    train_parents = {r['parent_source_id'] for r, _, _ in train}
    valid_parents = {r['parent_source_id'] for r, _, _ in valid}
    if train_hashes & valid_hashes:
        raise RuntimeError('Training/validation canonical geometry overlap detected')
    if train_parents & valid_parents:
        raise RuntimeError('Training/validation source parent overlap detected')

    for d in (FIT / 'training', FIT / 'validation'):
        d.mkdir(parents=True, exist_ok=True)
    staged = []
    for role, items in [('training', train), ('validation', valid)]:
        dst_dir = FIT / role
        for row, result, source in items:
            dst = dst_dir / f"{row['candidate_id']}.xyz"
            labelled = read(source)
            energy = float(labelled.get_potential_energy())
            expected_energy = float(result['energy_fit_target_eV'])
            if abs(energy - expected_energy) > 1.0e-8:
                raise RuntimeError(f"{row['candidate_id']}: extxyz energy differs from audited fit target")
            clean = Atoms(symbols=labelled.get_chemical_symbols(),
                          positions=labelled.get_positions(),
                          cell=labelled.cell, pbc=labelled.pbc)
            clean.calc = SinglePointCalculator(clean, energy=energy, forces=labelled.get_forces())
            write(dst, clean, format='extxyz')
            staged.append({'candidate_id': row['candidate_id'], 'role': role,
                           'parent_source_id': row['parent_source_id'],
                           'canonical_geometry_sha256': row['canonical_geometry_sha256'],
                           'source_label_sha256': sha(source), 'staged_file': str(dst),
                           'staged_file_sha256': sha(dst), 'energy_fit_target_eV': result['energy_fit_target_eV'],
                           'force_rmse_vector_norm_eV_A': result['force_rmse_vector_norm_eV_A'],
                           'max_force_norm_eV_A': result['max_force_norm_eV_A']})

    descriptor = FIT / 'AgTiSiC_POD_v13.pod'
    shutil.copy2(DESCRIPTOR_SOURCE, descriptor)
    fit_in = FIT / 'fitpod.in'
    fit_in.write_text('''clear
units metal
atom_style atomic
fitpod AgTiSiC_POD_v13.pod fitpod.data
''')
    fit_data = FIT / 'fitpod.data'
    fit_data.write_text(f'''file_format extxyz
file_extension xyz
path_to_training_data_set {FIT / 'training'}
path_to_test_data_set {FIT / 'validation'}
basename_for_output_files v13
fitting_weight_energy 100.0
fitting_weight_force 1.0
fitting_weight_stress 0.0
error_analysis_for_training_data_set 1
error_analysis_for_test_data_set 1
energy_force_calculation_for_training_data_set 0
energy_force_calculation_for_test_data_set 0
fraction_training_data_set 1.0
fraction_test_data_set 1.0
randomize_training_data_set 0
randomize_test_data_set 0
fitting_regularization_parameter 0.1
save_pod_descriptors 0
compute_pod_descriptors 0
''')
    record = {'status': 'staged_for_manual_fit_start_only_after_full_label_audit',
              'training_count': len(train), 'validation_count': len(valid),
              'training_manifest_sha256': sha(TRAIN), 'validation_manifest_sha256': sha(VALID),
              'descriptor_sha256': sha(descriptor), 'fit_input_sha256': sha(fit_in),
              'fit_data_sha256': sha(fit_data),
              'staged_labels': staged,
              'fit_protocol_note': 'validation set is used for model selection, not a final blind test; stress excluded; exploratory PW350 labels must be interpreted under their documented limits'}
    (FIT / 'fit_data_manifest.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'staged_training': len(train), 'staged_validation': len(valid), 'fit_input': str(fit_in)}, indent=2))


if __name__ == '__main__':
    main()
