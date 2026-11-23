#!/workspace/venvs/gpaw/bin/python
"""Audit available V13 training/validation DFT outputs and refresh manifests."""
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read

WORK = Path('/workspace/potential-fit/round7_pod_snap/models/target_200x200x300_qsgs_material200/v13_source_and_fit_20261123')
BASE = WORK.parent / 'v13_candidate_sources_20261123'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def canonical_hash(atoms):
    frac = atoms.get_scaled_positions(wrap=True)
    rows = sorted((s, *(round(float(x), 8) for x in q)) for s, q in zip(atoms.get_chemical_symbols(), frac))
    obj = {'cell': np.round(np.asarray(atoms.cell), 8).tolist(),
           'pbc': np.asarray(atoms.pbc, bool).tolist(), 'atoms': rows}
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def audit_manifest(manifest_path, label_root, role):
    manifest = json.loads(manifest_path.read_text())
    results = []
    for row in manifest['candidates']:
        cid = row['candidate_id']
        d = label_root / cid
        jp, xp = d / 'label.json', d / 'label.extxyz'
        item = {'candidate_id': cid, 'role': role, 'status': 'not_started'}
        if jp.exists():
            label = json.loads(jp.read_text())
            item.update({'status': label.get('run_status'), 'scf_converged': label.get('scf_converged'),
                         'finite_energy_force': label.get('finite_energy_force'),
                         'scf_iterations': label.get('scf_iterations'), 'max_force_norm_eV_A': label.get('max_force_norm_eV_A'),
                         'energy_fit_target_eV': label.get('energy_fit_target_eV'),
                         'geometry_file_sha_matches': label.get('geometry_sha256') == row.get('geometry_file_sha256'),
                         'canonical_hash_matches': label.get('canonical_geometry_sha256') == row.get('canonical_geometry_sha256')})
            if xp.exists():
                atoms = read(xp)
                e = float(atoms.get_potential_energy())
                f = np.asarray(atoms.get_forces(), float)
                item['extxyz_energy_finite'] = bool(np.isfinite(e))
                item['extxyz_forces_finite'] = bool(np.isfinite(f).all())
                item['extxyz_natoms_matches'] = len(atoms) == int(label.get('natoms', -1))
                item['label_extxyz_sha256'] = sha(xp)
                item['per_element_max_force_eV_A'] = {
                    s: float(np.linalg.norm(f[np.asarray(atoms.get_chemical_symbols()) == s], axis=1).max())
                    for s in sorted(set(atoms.get_chemical_symbols()))
                }
                item['acceptance_pass'] = bool(
                    label.get('run_status') == 'complete' and label.get('scf_converged') and label.get('finite_energy_force')
                    and item['geometry_file_sha_matches'] and item['canonical_hash_matches']
                    and item['extxyz_energy_finite'] and item['extxyz_forces_finite'] and item['extxyz_natoms_matches'])
            else:
                item['acceptance_pass'] = False
            if item['acceptance_pass']:
                row['DFT_label_status'] = 'complete; independent output audit passed'
                row['DFT_label_json'] = str(jp)
                row['DFT_label_extxyz'] = str(xp)
            else:
                row['DFT_label_status'] = f"review_required; runner status={label.get('run_status')}"
        elif (d / 'gpaw_scf.txt').exists():
            item['status'] = 'in_progress_or_incomplete'
            row['DFT_label_status'] = item['status']
        results.append(item)
    passed = sum(bool(x.get('acceptance_pass')) for x in results)
    active = sum(x['status'] == 'in_progress_or_incomplete' for x in results)
    not_started = sum(x['status'] == 'not_started' for x in results)
    manifest['completed_dft_labels'] = passed
    manifest['active_or_partial_dft_labels'] = active
    manifest['pending_dft_labels'] = not_started
    manifest['status'] = f'{passed} labels accepted; {active} in progress; {not_started} not started; no POD fit unless all split labels pass'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    return results, {'accepted': passed, 'active_or_partial': active, 'not_started': not_started}


def main():
    specs = [
        (BASE / 'training_candidate_manifest_v13.json', BASE / 'dft_labels', 'training'),
        (BASE / 'validation_candidate_manifest_v13.json', BASE / 'validation_dft_labels', 'validation'),
    ]
    all_results, counts = {}, {}
    for manifest, label_root, role in specs:
        if manifest.exists():
            all_results[role], counts[role] = audit_manifest(manifest, label_root, role)
    output = {'status': 'incremental label audit; fitting must wait until all training and validation labels pass',
              'counts': counts, 'labels': [x for group in all_results.values() for x in group]}
    out = BASE / 'dft_label_status_v13.json'
    out.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({'counts': counts, 'status_file': str(out)}, indent=2))


if __name__ == '__main__':
    main()
