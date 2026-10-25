#!/usr/bin/env python3
"""Decompose run-0 forces into current versus diagnostic zero-Morse terms."""
from __future__ import annotations
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ELEMENTS = {1: 'Ti', 2: 'Si', 3: 'C', 4: 'Ag'}

def load_dump(name: str) -> np.ndarray:
    with gzip.open(HERE / name, 'rt') as stream:
        return np.loadtxt(stream, skiprows=9)

def stats(values: np.ndarray) -> dict[str, float]:
    return {key: float(np.percentile(values, percentile)) for key, percentile in
            [('p50', 50), ('p95', 95), ('p99', 99)]} | {'max': float(np.max(values))}

full = load_dump('initial_forces.dump.gz')
no_cross = load_dump('forces_no_cross_morse.dump.gz')
if not np.array_equal(full[:, :2], no_cross[:, :2]):
    raise RuntimeError('Atom IDs/types differ between force dumps')
if not np.allclose(full[:, 2:5], no_cross[:, 2:5], atol=1e-10, rtol=0):
    raise RuntimeError('Force dumps do not use identical coordinates')

symbols = full[:, 1].astype(int)
positions = full[:, 2:5].copy()
force_full = full[:, 5:8]
force_no_cross = no_cross[:, 5:8]
force_cross = force_full - force_no_cross
magnitude_full = np.linalg.norm(force_full, axis=1)
magnitude_no_cross = np.linalg.norm(force_no_cross, axis=1)
magnitude_cross = np.linalg.norm(force_cross, axis=1)
box = json.loads((HERE / 'geometry_report.json').read_text())['simulation_box_A']
positions[:, 0] %= box[0]
positions[:, 1] %= box[1]
trees = {
    'Ag': cKDTree(positions[symbols == 4], boxsize=[box[0], box[1], 1000.]),
    'skeleton': cKDTree(positions[symbols != 4], boxsize=[box[0], box[1], 1000.]),
}
report = {
    'diagnostic': 'Run0 vector-force difference (full current calibration-v0 minus same-coordinate diagnostic with all three Ag-skeleton Morse D0 values set to zero). Intra-phase potentials and coordinates are unchanged. Distances are nearest opposite-phase atom.',
    'atoms': int(len(full)),
    'by_element': {},
}
for atom_type, name in ELEMENTS.items():
    selected = symbols == atom_type
    phase_tree = trees['Ag' if atom_type != 4 else 'skeleton']
    distance = phase_tree.query(positions[selected], k=1)[0]
    atom_indices = np.flatnonzero(selected)
    maximum_full_index = atom_indices[np.argmax(magnitude_full[selected])]
    maximum_cross_index = atom_indices[np.argmax(magnitude_cross[selected])]
    within4 = distance <= 4.0
    beyond8 = distance > 8.0
    report['by_element'][name] = {
        'nearest_opposite_phase_distance_A': {
            'p50': float(np.percentile(distance, 50)),
            'p95': float(np.percentile(distance, 95)),
            'fraction_within_4A': float(within4.mean()),
            'fraction_within_8A': float(np.mean(distance <= 8.0)),
        },
        'full_force_eV_A': stats(magnitude_full[selected]),
        'without_cross_force_eV_A': stats(magnitude_no_cross[selected]),
        'cross_force_component_eV_A': stats(magnitude_cross[selected]),
        'fraction_cross_force_gt_1_eV_A': float(np.mean(magnitude_cross[selected] > 1.0)),
        'fraction_cross_force_gt_2_eV_A': float(np.mean(magnitude_cross[selected] > 2.0)),
        'no_cross_force_by_location': {
            'within_4A_of_interface': stats(magnitude_no_cross[selected][within4]),
            'beyond_8A_from_interface': stats(magnitude_no_cross[selected][beyond8]) if beyond8.any() else None,
        },
        'max_full_force_atom': {
            'id': int(full[maximum_full_index, 0]),
            'position_A': full[maximum_full_index, 2:5].tolist(),
            'force_eV_A': float(magnitude_full[maximum_full_index]),
            'cross_component_eV_A': float(magnitude_cross[maximum_full_index]),
            'opposite_phase_distance_A': float(distance[np.argmax(magnitude_full[selected])]),
        },
        'max_cross_force_atom': {
            'id': int(full[maximum_cross_index, 0]),
            'position_A': full[maximum_cross_index, 2:5].tolist(),
            'cross_force_eV_A': float(magnitude_cross[maximum_cross_index]),
            'opposite_phase_distance_A': float(distance[np.argmax(magnitude_cross[selected])]),
        },
    }
(HERE / 'cross_morse_force_decomposition.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
