#!/usr/bin/env python3
"""Compare train43/train49 force residual patterns on frozen dev3 only."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

HERE = Path(__file__).resolve().parent
TASK = HERE.parent
DEV = HERE / 'dataset/valid.extxyz'
MODELS = {
    'train43': TASK / 'training_candidate_v19_train43_20261029_A/run/models/MACE_periodic_v19_residual_response.model',
    'train49': HERE / 'run/models/MACE_periodic_v19_residual_response.model',
}

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main() -> None:
    dev_hash = digest(DEV)
    manifest = json.loads((HERE / 'dataset/manifest.json').read_text())
    if dev_hash != manifest['dev_sha256']:
        raise SystemExit('Frozen dev3 hash mismatch')
    torch.set_default_dtype(torch.float64)
    calculators = {}
    model_hashes = {}
    for name, path in MODELS.items():
        model_hashes[name] = digest(path)
        model = torch.load(path, map_location='cpu', weights_only=False).double().eval()
        calculators[name] = MACECalculator(models=[model], device='cpu', default_dtype='float64')
    atoms_list = read(DEV, index=':')
    interface_rows = []
    all_errors = {name: [] for name in MODELS}
    for atoms in atoms_list:
        ref = np.asarray(atoms.arrays['REF_forces'], dtype=float)
        predicted = {}
        for name, calc in calculators.items():
            frame = atoms.copy()
            frame.calc = calc
            predicted[name] = np.asarray(frame.get_forces(), dtype=float)
        errors = {name: np.linalg.norm(force - ref, axis=1) for name, force in predicted.items()}
        for name in MODELS:
            all_errors[name].extend(errors[name].tolist())
        symbols = np.asarray(atoms.get_chemical_symbols())
        central = np.flatnonzero(np.asarray(atoms.arrays.get('central_pair', np.zeros(len(atoms))), dtype=bool))
        if len(central) != 2:
            raise SystemExit(f'Missing frozen central-pair marker: {atoms.info.get("config_type")}')
        distance = np.asarray([min(atoms.get_distance(i, int(j), mic=True) for j in central)
                               for i in range(len(atoms))])
        bins = {
            '<=3.5_A': distance <= 3.5,
            '3.5_to_5_A': (distance > 3.5) & (distance <= 5.0),
            '>5_A': distance > 5.0,
        }
        species_rows = {}
        for element in sorted(set(symbols.tolist())):
            mask = symbols == element
            vals = {name: float(np.sqrt(np.mean(errors[name][mask] ** 2))) for name in MODELS}
            species_rows[element] = {**vals, 'train49_minus_train43_eV_A': vals['train49'] - vals['train43']}
        distance_rows = {}
        for label, mask in bins.items():
            vals = {name: float(np.sqrt(np.mean(errors[name][mask] ** 2))) for name in MODELS}
            distance_rows[label] = {**vals, 'train49_minus_train43_eV_A': vals['train49'] - vals['train43']}
        ids = np.asarray(atoms.arrays.get('lammps_id', np.arange(len(atoms))))
        top = []
        for idx in np.argsort(errors['train49'])[-8:][::-1]:
            top.append({
                'atom_index_zero_based': int(idx), 'atom_id': int(ids[idx]), 'element': str(symbols[idx]),
                'z_angstrom': float(atoms.positions[idx, 2]),
                'nearest_central_pair_distance_angstrom': float(distance[idx]),
                'train43_error_norm_eV_A': float(errors['train43'][idx]),
                'train49_error_norm_eV_A': float(errors['train49'][idx]),
            })
        interface_rows.append({
            'config_type': atoms.info['config_type'],
            'aggregate_vector_rmse_eV_A': {name: float(np.sqrt(np.mean(errors[name] ** 2))) for name in MODELS},
            'species_vector_rmse_eV_A': species_rows,
            'nearest_contact_distance_bins_vector_rmse_eV_A': distance_rows,
            'largest_train49_atom_residuals': top,
        })
    result = {
        'status': 'TRAIN49_DEV3_RESIDUAL_DIAGNOSTIC_COMPLETE',
        'evaluation_scope': 'frozen dev3 only; no test or sealed data opened',
        'frozen_dev3_sha256': dev_hash,
        'model_sha256': model_hashes,
        'aggregate_force_vector_rmse_eV_A': {
            name: float(np.sqrt(np.mean(np.square(values)))) for name, values in all_errors.items()
        },
        'per_interface': interface_rows,
        'interpretation_limit': 'Residual correlations identify development error patterns; they do not prove a causal mechanism or authorize copying dev geometries into training.',
    }
    out = HERE / 'run/results/train49_residual_diagnostics.json'
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
