#!/usr/bin/env python3
"""Four explicitly public B controls: frozen inference only; never DFT, MD or training."""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('MKL_NUM_THREADS', '2')
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'

import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

torch.set_num_threads(2)
torch.set_num_interop_threads(1)
torch.set_default_dtype(torch.float64)

ROOT = Path(os.environ.get('CLOUD_C_REPO_ROOT', Path(__file__).resolve().parents[6])).resolve()
BASE = ROOT / 'research/mace-v12-transfer'
PUBLIC = BASE / 'periodic_interface_v4/pbe_interface_v19_pure_phase_controls'
B_REPORT = BASE / 'coordination/reports/window-b/v19_Ti3SiC2_backbone_force_response_execution_window_b'
B_SHA = 'c34cfa6f6505fc6b385908b6fad79119354cf24d8647423a89863765f8955c4e'
B_COMMIT = '34b8be988015af96d3fd8791fc041f6cc811e2b7'
ROLE = 'numerical_pure_phase_control_not_training_or_independent_validation'
PARENT_LABEL = 'Ti3SiC2_baseline_k10x10x4_v19'
PARENT_INPUT = PUBLIC / 'inputs/Ti3SiC2_baseline_PROPOSAL.extxyz'
PARENT_SHA = '157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4'
METHOD_SHA = '702df71e586a61cd5785afc62653d7135a707032c91b72842b035b62198fd709'
MODELS = {
    'MACE-MP-0b3-medium': (
        BASE / 'periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model',
        '2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'),
    'V18': (
        BASE / 'periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/models/MACE_periodic_v18_clean_core.model',
        '757586671174db85acbf7971ccfe92de19cc58378ee2431a797efa6d6163cb62'),
}
SITES = {'Ti': 3, 'C': 9}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def csv_write(out, name, rows):
    assert rows
    with (out / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def rel(path):
    return str(path.relative_to(ROOT))


def identities(atoms):
    ids = np.asarray(atoms.arrays['local_new_id'], dtype=int)
    assert len(ids) == len(set(ids.tolist()))
    return ids


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    assert sha(B_REPORT / 'final_results.json') == B_SHA
    assert sha(PARENT_INPUT) == PARENT_SHA
    final = json.loads((B_REPORT / 'final_results.json').read_text())
    manifest = json.loads((PUBLIC / 'input_manifest.json').read_text())
    parent = read(PARENT_INPUT)
    assert len(parent) == 48 and parent.pbc.tolist() == [True, True, True]
    refs = {}
    verification = {}
    source_rows = []
    for label, published in final['points'].items():
        assert label.startswith('Ti3SiC2_') and label.endswith('_k10x10x4_v19')
        assert len(final['points']) == 4
        folder = PUBLIC / 'calculations' / label
        process = subprocess.run(
            [sys.executable, '-B', str(PUBLIC / 'verify_result.py'), label, str(folder)],
            cwd=ROOT, text=True, capture_output=True, check=True,
            env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
        current = json.loads(process.stdout)
        assert current['status'] == published['archive_verification'] == 'PASS'
        sums = {}
        for line in (folder / 'SHA256SUMS.txt').read_text().splitlines():
            expected, name = line.split('  ', 1)
            path = (folder / name).resolve()
            assert path.is_relative_to(folder.resolve()) and sha(path) == expected
            assert name not in sums
            sums[name] = expected
        for name, expected in published['archive_file_sha256'].items():
            assert sums[name] == expected and sha(folder / name) == expected
        record = next(item for item in manifest['records'] if item['label'] == label)
        source = PUBLIC / record['input']
        result = folder / (label + '_PW_PBE.extxyz')
        assert sha(source) == published['input_sha256']
        assert sha(ROOT / record['source_input']) == sha(source)
        src = read(source)
        atoms = read(result)
        summary = json.loads((folder / 'summary.json').read_text())
        assert summary['dataset_role'] == record['role'] == ROLE
        assert summary['scf_converged'] is True
        assert canonical_sha(summary['method']) == METHOD_SHA
        assert np.array_equal(src.positions, atoms.positions)
        assert np.array_equal(src.numbers, atoms.numbers)
        assert np.array_equal(src.cell.array, atoms.cell.array)
        assert np.array_equal(identities(src), identities(atoms))
        assert np.array_equal(identities(atoms), identities(parent))
        assert np.array_equal(atoms.numbers, parent.numbers)
        assert np.array_equal(atoms.cell.array, parent.cell.array)
        assert np.array_equal(atoms.pbc, parent.pbc)
        delta = atoms.positions - parent.positions
        changed = np.flatnonzero(np.linalg.norm(delta, axis=1) > 1e-12)
        assert len(changed) == 1
        index = int(changed[0])
        atom_id = int(identities(atoms)[index])
        element = str(atoms.symbols[index])
        assert atom_id == SITES[element] == published['moved_atom']['local_new_id']
        assert abs(abs(delta[index, 2]) - .02) < 1e-12
        assert abs(delta[index, 0]) < 1e-12 and abs(delta[index, 1]) < 1e-12
        forces = np.asarray(atoms.arrays['PW_PBE_forces'], dtype=float)
        assert forces.shape == (48, 3) and np.isfinite(forces).all()
        assert np.allclose(forces[index], published['moved_atom']['force_eV_A'], rtol=0, atol=1e-12)
        assert summary['energy_eV_cell'] == published['energy_eV_cell']
        assert summary['free_energy_eV_cell'] == published['free_energy_eV_cell']
        assert atoms.info['PW_PBE_energy_eV'] == summary['energy_eV_cell']
        sign = -1 if delta[index, 2] < 0 else 1
        refs[label] = {'atoms': atoms, 'forces': forces, 'summary': summary,
                       'index': index, 'atom_id': atom_id, 'element': element, 'sign': sign}
        identity_payload = {'symbols': atoms.get_chemical_symbols(), 'ids': identities(atoms).tolist()}
        geom_payload = {**identity_payload, 'positions_A': atoms.positions.tolist(),
                        'cell_A': atoms.cell.array.tolist(), 'pbc': atoms.pbc.tolist()}
        verification[label] = {
            'rerun_readonly_verifier': current, 'archive_sha256': sums,
            'archive_inventory_sha256': sha(folder / 'SHA256SUMS.txt'),
            'source_path': rel(source), 'source_sha256': sha(source),
            'result_path': rel(result), 'result_sha256': sha(result),
            'method_manifest_sha256': canonical_sha(summary['method']),
            'ordered_atom_identity_sha256': canonical_sha(identity_payload),
            'geometry_canonical_sha256': canonical_sha(geom_payload),
            'parent_input_sha256': PARENT_SHA, 'parent_family_id': 'Ti3SiC2_bulk_parent_v19',
            'ordered_ids': identities(atoms).tolist(), 'ordered_symbols': atoms.get_chemical_symbols(),
            'cell_A': atoms.cell.array.tolist(), 'pbc': atoms.pbc.tolist(),
            'target_local_new_id': atom_id, 'target_index_zero_based': index,
            'delta_z_A': float(delta[index, 2]), 'dataset_role': ROLE,
        }
        source_rows.append({
            'label': label, 'parent_family_id': 'Ti3SiC2_bulk_parent_v19', 'parent_label': PARENT_LABEL,
            'parent_input_sha256': PARENT_SHA, 'dataset_role': ROLE, 'source_path': rel(source),
            'input_sha256': sha(source), 'result_path': rel(result), 'result_sha256': sha(result),
            'geometry_canonical_sha256': canonical_sha(geom_payload),
            'ordered_atom_identity_sha256': canonical_sha(identity_payload), 'n_atoms': 48,
            'target_element': element, 'target_local_new_id': atom_id,
            'target_index_zero_based': index, 'delta_z_A': float(delta[index, 2]),
            'pbc': 'T T T', 'verification': 'PASS', 'scf_iterations': summary['scf_iterations'],
            'DFT_energy_eV_cell': summary['energy_eV_cell'], 'DFT_free_energy_eV_cell': summary['free_energy_eV_cell'],
            'energy_convention': summary['energy_convention'], 'method_manifest_sha256': METHOD_SHA,
        })
    predictions = {}
    atom_rows = []
    metric_rows = []
    element_rows = []
    model_records = {}
    for name, (path, expected) in MODELS.items():
        assert sha(path) == expected
        model_start = time.perf_counter()
        calc = MACECalculator(model_paths=str(path), device='cpu', default_dtype='float64')
        net = calc.models[0]
        net.eval()  # Explicit inference mode; does not modify the serialized weights.
        assert not net.training
        assert all(not p.requires_grad for p in net.parameters())
        supported = [int(value) for value in net.atomic_numbers.cpu().tolist()]
        assert set([6, 14, 22]).issubset(supported)
        assert set(str(p.dtype) for p in net.parameters()) == {'torch.float64'}
        model_records[name] = {'repo_relative_path': rel(path), 'sha256': expected,
                              'atomic_numbers_order': supported, 'r_max_A': float(net.r_max),
                              'num_interactions': int(net.num_interactions), 'dtype': 'float64',
                              'device': 'cpu', 'training': False, 'weight_gradients_enabled': False}
        for label, ref in refs.items():
            one_start = time.perf_counter()
            atoms = ref['atoms'].copy()
            atoms.calc = calc
            predicted = atoms.get_forces()
            assert predicted.shape == ref['forces'].shape and np.isfinite(predicted).all()
            predictions[name, label] = predicted.copy()
            delta = predicted - ref['forces']
            norms = np.linalg.norm(delta, axis=1)
            ids = identities(atoms)
            metric = {'label': label, 'model': name, 'model_sha256': expected,
                      'parent_family_id': 'Ti3SiC2_bulk_parent_v19', 'dataset_role': ROLE,
                      'n_atoms': len(atoms), 'force_vector_RMSE_eV_A': float(np.sqrt(np.mean(norms ** 2))),
                      'max_atom_force_vector_error_eV_A': float(norms.max()),
                      'x_RMSE_eV_A': float(np.sqrt(np.mean(delta[:, 0] ** 2))),
                      'y_RMSE_eV_A': float(np.sqrt(np.mean(delta[:, 1] ** 2))),
                      'z_RMSE_eV_A': float(np.sqrt(np.mean(delta[:, 2] ** 2))),
                      'max_error_local_new_id': int(ids[np.argmax(norms)]),
                      'target_local_new_id': ref['atom_id'], 'target_index_zero_based': ref['index'],
                      'target_DFT_Fz_eV_A': float(ref['forces'][ref['index'], 2]),
                      'target_model_Fz_eV_A': float(predicted[ref['index'], 2]),
                      'inference_elapsed_s': time.perf_counter() - one_start,
                      'prediction_origin': 'fresh_frozen_float64_CPU_2_threads'}
            metric_rows.append(metric)
            symbols = np.asarray(atoms.get_chemical_symbols())
            for element in sorted(set(symbols.tolist())):
                mask = symbols == element
                element_rows.append({'label': label, 'model': name, 'element': element,
                                     'n_atoms': int(mask.sum()),
                                     'force_vector_RMSE_eV_A': float(np.sqrt(np.mean(norms[mask] ** 2))),
                                     'max_atom_force_vector_error_eV_A': float(norms[mask].max()),
                                     'x_RMSE_eV_A': float(np.sqrt(np.mean(delta[mask, 0] ** 2))),
                                     'y_RMSE_eV_A': float(np.sqrt(np.mean(delta[mask, 1] ** 2))),
                                     'z_RMSE_eV_A': float(np.sqrt(np.mean(delta[mask, 2] ** 2)))})
            for i in range(len(atoms)):
                row = {'label': label, 'model': name, 'parent_family_id': 'Ti3SiC2_bulk_parent_v19',
                       'dataset_role': ROLE, 'index_zero_based': i, 'local_new_id': int(ids[i]),
                       'element': str(symbols[i]), 'x_A': float(atoms.positions[i, 0]),
                       'y_A': float(atoms.positions[i, 1]), 'z_A': float(atoms.positions[i, 2]),
                       'error_norm_eV_A': float(norms[i])}
                for k, axis in enumerate('xyz'):
                    row[f'DFT_F{axis}_eV_A'] = float(ref['forces'][i, k])
                    row[f'model_F{axis}_eV_A'] = float(predicted[i, k])
                    row[f'error_F{axis}_eV_A'] = float(delta[i, k])
                atom_rows.append(row)
            print(name, label, metric['force_vector_RMSE_eV_A'], metric['target_model_Fz_eV_A'], flush=True)
        assert sha(path) == expected
        model_records[name]['elapsed_s'] = time.perf_counter() - model_start
    response_rows = []
    for element, atom_id in SITES.items():
        minus = next(label for label, ref in refs.items() if ref['element'] == element and ref['sign'] == -1)
        plus = next(label for label, ref in refs.items() if ref['element'] == element and ref['sign'] == 1)
        idx = refs[minus]['index']
        assert idx == refs[plus]['index']
        span = float(refs[plus]['atoms'].positions[idx, 2] - refs[minus]['atoms'].positions[idx, 2])
        assert abs(span - .04) < 1e-12
        dft_m = float(refs[minus]['forces'][idx, 2])
        dft_p = float(refs[plus]['forces'][idx, 2])
        dft_k = -(dft_p - dft_m) / .04
        expected_k = final[f'{element}_4f_symmetric_z_force_response_eV_A2']
        assert abs(dft_k - expected_k) < 1e-12
        for name in MODELS:
            model_m = float(predictions[name, minus][idx, 2])
            model_p = float(predictions[name, plus][idx, 2])
            model_k = -(model_p - model_m) / .04
            response_rows.append({'site': element + ' 4f', 'local_new_id': atom_id,
                                  'index_zero_based': idx, 'parent_family_id': 'Ti3SiC2_bulk_parent_v19',
                                  'dataset_role': ROLE, 'minus_label': minus, 'plus_label': plus,
                                  'displacement_span_A': span, 'model': name,
                                  'DFT_Fz_minus_eV_A': dft_m, 'DFT_Fz_plus_eV_A': dft_p,
                                  'model_Fz_minus_eV_A': model_m, 'model_Fz_plus_eV_A': model_p,
                                  'DFT_Kz_eV_A2': dft_k, 'model_Kz_eV_A2': model_k,
                                  'Kz_signed_error_eV_A2': model_k - dft_k,
                                  'Kz_relative_signed_error_percent': 100 * (model_k - dft_k) / dft_k})
    csv_write(out, 'source_identity.csv', source_rows)
    csv_write(out, 'per_atom.csv', atom_rows)
    csv_write(out, 'per_configuration_metrics.csv', metric_rows)
    csv_write(out, 'per_element.csv', element_rows)
    csv_write(out, 'response_comparison.csv', response_rows)
    (out / 'source_verification.json').write_text(json.dumps(verification, indent=2) + '\n')
    metadata = {'source_B_commit': B_COMMIT, 'source_B_final_results_sha256': B_SHA,
                'checkout_HEAD_at_execution': subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
                'runtime': {'python': sys.version.split()[0], 'packages': {n: importlib.metadata.version(n) for n in ['torch', 'mace-torch', 'ase', 'numpy', 'scipy']},
                            'torch_threads': torch.get_num_threads(), 'torch_interop_threads': torch.get_num_interop_threads(),
                            'OMP_NUM_THREADS': os.environ['OMP_NUM_THREADS'], 'OPENBLAS_NUM_THREADS': os.environ['OPENBLAS_NUM_THREADS'],
                            'MKL_NUM_THREADS': os.environ['MKL_NUM_THREADS'], 'dtype': 'float64', 'device': 'cpu'},
                'models': model_records, 'total_elapsed_s': time.perf_counter() - start,
                'metric_definitions': {'force_vector_RMSE': 'sqrt(mean_atom(sum_xyz((F_model-F_DFT)^2)))',
                                       'max_atom_error': 'max_atom(norm_xyz(F_model-F_DFT))',
                                       'Kz': '-(Fz_plus-Fz_minus)/0.04 Angstrom',
                                       'relative_signed_error_percent': '100*(K_model-K_DFT)/K_DFT'},
                'scope': {'sealed_data_read': False, 'new_DFT': False, 'MD': False, 'training': False,
                          'model_weights_changed': False, 'data_roles_changed': False, 'archived_reports_changed': False},
                'prediction_policy': 'All eight frozen-model configurations recomputed uniformly; previous three-point diagnostics are cross-checks only.'}
    (out / 'inference_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(response_rows, indent=2), flush=True)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True, help='Independent C output directory; source archives are read-only.')
    args = ap.parse_args()
    main(args.output.resolve())
