#!/usr/bin/env python3
"""Integrate B's approved clean30 proposal; --check verifies without writing."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from ase.io import read

ENTRY = Path('research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core')
BASE = ENTRY.parent
REPORTS = Path('research/mace-v12-transfer/coordination/reports/window-b')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    proposal = repo / REPORTS / 'v17_clean_core_builder_draft'
    source = json.loads((proposal / 'manifest.json').read_text())
    require(source['status'] == 'PROPOSAL_ONLY_NOT_INTEGRATED', 'Unexpected proposal status')
    require(source['frames'] == 30 and len(source['rows']) == 30, 'Wrong proposal count')
    require(sha(proposal / 'train_candidate.extxyz') == source['output_sha256'], 'Proposal hash mismatch')
    pins = dict(source['sources_sha256'])
    for name in ['manifest.json', 'train_candidate.extxyz', 'build.py']:
        path = proposal / name
        pins[str(path.relative_to(repo))] = sha(path)
    dev = repo / BASE / 'mace_periodic_v17_interface_energy/data/valid.extxyz'
    v17manifest = json.loads((dev.parent / 'dataset_manifest.json').read_text())
    require(sha(dev) == v17manifest['development_validation_input_sha256'], 'V17 dev snapshot changed')
    pins[str(dev.relative_to(repo))] = sha(dev)
    review = repo / REPORTS / 'v17_iteration_review/report.md'
    pins[str(review.relative_to(repo))] = sha(review)
    for rel, expected in pins.items():
        path = (repo / rel).resolve()
        require(path.is_relative_to(repo) and path.is_file() and sha(path) == expected, f'Source mismatch: {rel}')
    train = read(proposal / 'train_candidate.extxyz', ':')
    valid = read(dev, ':')
    require(len(train) == 30 and len(valid) == 3, 'Serialized split count mismatch')
    for split, frames in [('train', train), ('development_validation', valid)]:
        for i, a in enumerate(frames):
            require('REF_energy' in a.info and 'REF_forces' in a.arrays, f'Missing labels {split}[{i}]')
            f = np.asarray(a.arrays['REF_forces'])
            require(f.shape == (len(a), 3) and np.isfinite(f).all() and np.isfinite(a.info['REF_energy']), f'Invalid labels {split}[{i}]')
            require(np.isfinite(a.positions).all() and np.isfinite(a.cell.array).all(), f'Invalid geometry {split}[{i}]')
    spec = importlib.util.spec_from_file_location('v17_preflight', repo / BASE / 'mace_periodic_v17_interface_energy/preflight_v17.py')
    rank_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rank_module)
    elements = ['Ag', 'C', 'Si', 'Ti']
    matrix = np.asarray([[a.get_chemical_symbols().count(e) for e in elements] for a in train])
    require(rank_module.exact_rank(matrix) == 4, 'Composition rank is not four')
    spec = importlib.util.spec_from_file_location('geometry', repo / REPORTS / 'v17_split_geometry_design/final_set/generate_split_candidates.py')
    geom = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(geom)
    for i, a in enumerate(train):
        for j, b in enumerate(valid):
            match = geom.compare_geometry(a, b)
            require(not match or not (match['exact_duplicate'] or match['near_duplicate']), f'Train/dev overlap {i}/{j}')
    foundation = repo / BASE / 'mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model'
    require(sha(foundation) == v17manifest['foundation_model_sha256'], 'Foundation hash mismatch')
    manifest = {
        'status': 'PROVISIONAL_CLEAN30_CONTROLLED_SCREENING',
        'roles': {'train': 30, 'development_validation': 3, 'test': 0},
        'train_input_sha256': sha(proposal / 'train_candidate.extxyz'),
        'development_validation_input_sha256': sha(dev),
        'foundation_model_sha256': sha(foundation),
        'source_files_sha256': pins,
        'composition_matrix_elements': elements,
        'composition_matrix': matrix.tolist(),
        'composition_matrix_rank_exact': 4,
        'rows': source['rows'],
        'provenance_counts': {'directly_archive_verified': 30},
        'energy_convention': source['energy_convention'],
        'frame0_kpoint_convergence': source['frame0_kpoint_convergence'],
        'withheld_labels_opened': False,
        'limitations': ['Removing 14 partial-source rows also removes coverage; adding frame0 changes membership and its declared mesh.', 'Archive support does not establish method uniformity or k-point convergence.', 'Development structures are reused; this comparison is not fresh blind validation.'],
    }
    data = repo / ENTRY / 'data'
    if args.check:
        require(json.loads((data / 'dataset_manifest.json').read_text()) == manifest, 'Integrated manifest changed')
        require(sha(data / 'train.extxyz') == manifest['train_input_sha256'], 'Integrated train hash mismatch')
        require(sha(data / 'valid.extxyz') == manifest['development_validation_input_sha256'], 'Integrated dev hash mismatch')
        require(not (data / 'test.extxyz').exists(), 'Unexpected test file')
    else:
        require(not data.exists(), 'Existing dataset refused; preserve prior outputs')
        data.mkdir()
        (data / 'train.extxyz').write_bytes((proposal / 'train_candidate.extxyz').read_bytes())
        (data / 'valid.extxyz').write_bytes(dev.read_bytes())
        (data / 'dataset_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        (data / 'SHA256SUMS.txt').write_text(''.join(f'{sha(data / n)}  {n}\n' for n in ['train.extxyz', 'valid.extxyz', 'dataset_manifest.json']))
    print(json.dumps({'status': 'PROVISIONAL_PREFLIGHT_PASS', 'roles': manifest['roles'], 'composition_rank_exact': 4, 'withheld_labels_opened': False}, indent=2))


if __name__ == '__main__':
    main()
