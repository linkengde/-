#!/usr/bin/env python3
"""Prepare disjoint acquisition and reserved holdout geometries for window B."""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
INPUTS = ROOT / 'inputs'
MAIN_ROOT = PROJECT / 'pbe_interface_v14_main_holdouts'
MAIN_INPUTS = MAIN_ROOT / 'inputs'
if any(p.exists() and any(p.iterdir()) for p in (INPUTS, MAIN_INPUTS)):
    raise SystemExit('Inputs already exist; inspect their hashes instead of regenerating them.')
INPUTS.mkdir(parents=True, exist_ok=True)
MAIN_INPUTS.mkdir(parents=True, exist_ok=True)
SOURCES = {
    'AgC': PROJECT / 'pbe_cluster_transfer_pw_20261002_run3_mixer_retry2/AgC_cluster_PW_PBE.extxyz',
    'AgSi': PROJECT / 'pbe_external_v8_motif_holdouts/AgSi/retry1/AgSi_cluster_PW_PBE.extxyz',
    'AgTi': PROJECT / 'pbe_interface_energy_additions_v12/calculations/validation_distance_scans/AgTi_validation_d2p60/AgTi_validation_d2p60_PW_PBE.extxyz',
}
SPECS = [
    ('AgTi_registry_probe_v13_01', 'AgTi', 'acquisition', 0.22, 0.0, 0.015, 0.000, 14001),
    ('AgC_registry_strain_acq_v14_01', 'AgC', 'acquisition', 0.28, 0.7, 0.025, 0.005, 14002),
    ('AgSi_registry_strain_acq_v14_01', 'AgSi', 'acquisition', 0.28, -0.7, 0.025, -0.005, 14003),
    ('AgTi_registry_holdout_v14_01', 'AgTi', 'v14_holdout', 0.32, 1.6, 0.035, -0.010, 14101),
    ('AgC_registry_holdout_v14_01', 'AgC', 'v14_holdout', 0.32, -1.4, 0.035, -0.010, 14102),
    ('AgSi_registry_holdout_v14_01', 'AgSi', 'v14_holdout', 0.32, 1.5, 0.035, 0.010, 14103),
]
train = read(PROJECT / 'mace_periodic_v13_interface_energy/data/train.extxyz', index=':')
records = []
for label, kind, role, shift, angle, sigma, strain, seed in SPECS:
    source = SOURCES[kind]
    atoms = read(source)
    symbols = np.asarray(atoms.get_chemical_symbols())
    pair = np.flatnonzero(atoms.arrays['central_pair'])
    ia = int(pair[symbols[pair] == 'Ag'][0])
    ib = int(pair[symbols[pair] == kind[2:]][0])
    vector = atoms.get_distance(ia, ib, mic=True, vector=True)
    tangent1 = np.cross(vector, [0.0, 0.0, 1.0])
    if np.linalg.norm(tangent1) < 1e-8: tangent1 = np.cross(vector, [0.0, 1.0, 0.0])
    tangent1 /= np.linalg.norm(tangent1)
    tangent2 = np.cross(vector, tangent1); tangent2 /= np.linalg.norm(tangent2)
    translation = shift * (np.cos(angle) * tangent1 + np.sin(angle) * tangent2)
    atoms.set_cell(atoms.cell.array * (1.0 + strain), scale_atoms=True)
    ag = np.flatnonzero(symbols == 'Ag')
    atoms.positions[ag] += translation
    distances = atoms.get_all_distances(mic=True)
    local = np.flatnonzero(np.min(distances[:, pair], axis=1) < 4.2)
    non_ag = local[symbols[local] != 'Ag']
    rng = np.random.default_rng(seed)
    atoms.positions[non_ag] += rng.normal(0, sigma, (len(non_ag), 3))
    ag_noise = 0.020 if role == 'v14_holdout' else 0.008
    atoms.positions[ag] += rng.normal(0, ag_noise, (len(ag), 3))
    for key in ('REF_energy', 'PW_PBE_energy_eV', 'energy'): atoms.info.pop(key, None)
    for key in ('REF_forces', 'PW_PBE_forces', 'forces'): atoms.arrays.pop(key, None)
    atoms.calc = None
    atoms.info.update(config_type=label, dataset_role=role,
                      geometry_family=kind + '_existing_cluster_motif',
                      correlation_note='Derived from an existing motif; not an independent material morphology.')
    for prior in train:
        if (len(prior) == len(atoms) and prior.get_chemical_symbols() == atoms.get_chemical_symbols()
                and np.allclose(prior.cell.array, atoms.cell.array, atol=1e-10, rtol=0)
                and np.allclose(prior.positions, atoms.positions, atol=1e-10, rtol=0)):
            raise SystemExit('Generated geometry duplicates v13 training: ' + label)
    d = atoms.get_all_distances(mic=True); np.fill_diagonal(d, np.inf)
    i, j = map(int, np.unravel_index(d.argmin(), d.shape))
    min_ag_other = float(d[np.ix_(ag, np.flatnonzero(symbols != 'Ag'))].min())
    if min_ag_other < 1.8 or float(d.min()) < 1.5:
        raise SystemExit('Unexplained close contact in ' + label)
    target_inputs = MAIN_INPUTS if role == 'v14_holdout' and kind != 'AgTi' else INPUTS
    out = target_inputs / (label + '.extxyz')
    write(out, atoms, format='extxyz')
    records.append({
        'label': label, 'role': role, 'input': 'inputs/' + out.name,
        'input_sha256': sha256(out.read_bytes()).hexdigest(),
        'source': str(source.relative_to(PROJECT)), 'source_sha256': sha256(source.read_bytes()).hexdigest(),
        'atoms': len(atoms), 'formula': atoms.get_chemical_formula(),
        'central_pair': {'elements': 'Ag-' + kind[2:], 'Ag_id': int(atoms.arrays['lammps_id'][ia]),
                         'contact_id': int(atoms.arrays['lammps_id'][ib]),
                         'input_distance_A': float(atoms.get_distance(ia, ib, mic=True))},
        'minimum_pair': {'elements': str(symbols[i]) + '-' + str(symbols[j]), 'distance_A': float(d[i, j])},
        'minimum_Ag_nonAg_distance_A': min_ag_other,
        'perturbation': {'translation_A': translation.tolist(), 'strain': strain,
                         'nonAg_noise_sigma_A': sigma, 'Ag_noise_sigma_A': ag_noise, 'seed': seed},
        'scope': 'Periodic 28 A cluster-cell motif with new registry, small noise and optional strain; not a bulk interface or a temperature-labelled snapshot.',
    })
manifest = {'owner_task': 'window-b', 'purpose': 'Parallel labels without repeating window-a jobs.',
            'method': 'GPAW 26.7.0 PW-PBE 500 eV Gamma 0.1 eV; 4 MPI ranks with ScaLAPACK.',
            'reserved_holdout_policy': 'v14_holdout labels must stay outside v14 training.',
            'records': [r for r in records if r['role'] != 'v14_holdout' or r['label'].startswith('AgTi')]}
(ROOT / 'input_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
main_manifest = dict(manifest, owner_task='window-a', records=[r for r in records if r not in manifest['records']])
(MAIN_ROOT / 'input_manifest.json').write_text(json.dumps(main_manifest, indent=2) + '\n')
print(json.dumps({'prepared': len(records), 'records': [{k:r[k] for k in ('label','role','minimum_pair','minimum_Ag_nonAg_distance_A')} for r in records]}, indent=2))
