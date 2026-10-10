#!/usr/bin/env python3
"""Validate saved public-control diagnostics without invoking any model or calculator."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(out):
    sources = read_csv(out / 'source_identity.csv')
    atoms = read_csv(out / 'per_atom.csv')
    metrics = read_csv(out / 'per_configuration_metrics.csv')
    responses = read_csv(out / 'response_comparison.csv')
    elements = read_csv(out / 'per_element.csv')
    verification = json.loads((out / 'source_verification.json').read_text())
    metadata = json.loads((out / 'inference_metadata.json').read_text())
    model_names = {'MACE-MP-0b3-medium', 'V18'}
    assert len(sources) == len({r['label'] for r in sources}) == 4
    labels = {r['label'] for r in sources}
    assert len(metrics) == 8 and {(r['label'], r['model']) for r in metrics} == {(l, m) for l in labels for m in model_names}
    assert len(atoms) == 384 and len(elements) == 24 and len(responses) == 4
    assert len({(r['label'], r['model'], r['local_new_id']) for r in atoms}) == 384
    assert {r['dataset_role'] for r in sources + atoms + metrics + responses} == {'numerical_pure_phase_control_not_training_or_independent_validation'}
    assert {r['parent_family_id'] for r in sources + atoms + metrics + responses} == {'Ti3SiC2_bulk_parent_v19'}
    assert all(value is False for value in metadata['scope'].values())
    assert metadata['runtime']['torch_threads'] == 2
    assert metadata['runtime']['dtype'] == 'float64' and metadata['runtime']['device'] == 'cpu'
    assert set(verification) == labels
    for source in sources:
        record = verification[source['label']]
        assert record['rerun_readonly_verifier']['status'] == 'PASS'
        assert all(record['rerun_readonly_verifier']['checks'].values())
        assert source['target_index_zero_based'] == str({'Ti': 2, 'C': 8}[source['target_element']])
        assert source['ordered_atom_identity_sha256'] == record['ordered_atom_identity_sha256']
    for metric in metrics:
        selected = [r for r in atoms if r['label'] == metric['label'] and r['model'] == metric['model']]
        assert len(selected) == 48
        norms = []
        for row in selected:
            square = 0
            for axis in 'xyz':
                reference = float(row[f'DFT_F{axis}_eV_A'])
                predicted = float(row[f'model_F{axis}_eV_A'])
                err = predicted - reference
                assert math.isfinite(err)
                assert abs(err - float(row[f'error_F{axis}_eV_A'])) < 1e-12
                square += err * err
            norms.append(math.sqrt(square))
            assert abs(norms[-1] - float(row['error_norm_eV_A'])) < 1e-12
        assert abs(math.sqrt(sum(x*x for x in norms) / 48) - float(metric['force_vector_RMSE_eV_A'])) < 1e-12
        assert abs(max(norms) - float(metric['max_atom_force_vector_error_eV_A'])) < 1e-12
    for response in responses:
        dft_k = -(float(response['DFT_Fz_plus_eV_A']) - float(response['DFT_Fz_minus_eV_A'])) / .04
        model_k = -(float(response['model_Fz_plus_eV_A']) - float(response['model_Fz_minus_eV_A'])) / .04
        assert abs(dft_k - float(response['DFT_Kz_eV_A2'])) < 1e-12
        assert abs(model_k - float(response['model_Kz_eV_A2'])) < 1e-12
        assert abs(100*(model_k - dft_k)/dft_k - float(response['Kz_relative_signed_error_percent'])) < 1e-12
    sums = out / 'SHA256SUMS.txt'
    hashes = {}
    if sums.exists():
        for line in sums.read_text().splitlines():
            digest, name = line.split('  ', 1)
            assert name not in hashes and name != 'SHA256SUMS.txt'
            assert (out/name).resolve().is_relative_to(out.resolve())
            assert sha(out/name) == digest
            hashes[name] = digest
        assert set(hashes) == {p.name for p in out.iterdir() if p.is_file() and p.name != 'SHA256SUMS.txt'}
    result = {'status': 'PASS', 'counts': {'public_structures': 4, 'model_configuration_pairs': 8, 'per_atom_rows': 384, 'element_rows': 24, 'response_rows': 4},
              'checks': {'unique_structure_labels': True, 'unique_atom_rows': True, 'roles_unchanged': True,
                         'same_parent_family': True, 'four_readonly_archive_verifiers_passed': True,
                         'target_ID_to_index_mapping': True, 'force_metrics_recomputed_from_saved_rows': True,
                         'symmetric_responses_recomputed_from_saved_rows': True,
                         'frozen_cpu_float64_two_threads': True, 'scope_flags_clear': True},
              'note': 'The SHA manifest is checked when present. This validator does not launch inference, DFT, MD or training.'}
    result_path = out / 'output_verification.json'
    # Once included in a frozen checksum inventory, verification reads all files without changing it.
    if result_path.exists() and hashes:
        assert json.loads(result_path.read_text()) == result
    else:
        result_path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
    main(parser.parse_args().output.resolve())
