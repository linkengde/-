#!/workspace/venvs/dft-labels/bin/python
"""Compare DFT and POD net-force trends along targeted Ti-C pairs."""
from pathlib import Path
import json
import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
LABELS = [
    ('initial_pair1', 'tic_pair_initial_pw350_gamma'),
    ('cold_10fs_pair1', 'tic_pair_cold_10fs_pw350_gamma'),
    ('cold_10fs_pair2', 'tic_pair2_cold_10fs_pw350_gamma'),
]
MASS = {'Ti': 47.867, 'C': 12.011}
results = []
for name, stem in LABELS:
    path = HERE/'labels'/f'{stem}.extxyz'
    meta = HERE/'labels'/f'{stem}.json'
    if not path.is_file() or not meta.is_file():
        continue
    record = json.loads(meta.read_text())
    if record.get('status') != 'CONVERGED':
        continue
    atoms = read(path)
    ids = np.asarray(atoms.arrays['source_atom_id'], dtype=int)
    ti_id, c_id = map(int, record['center_atom_ids'])
    ti = int(np.flatnonzero(ids == ti_id)[0]); c = int(np.flatnonzero(ids == c_id)[0])
    if atoms[ti].symbol != 'Ti' or atoms[c].symbol != 'C':
        raise RuntimeError(f'Element/ID mismatch in {stem}')
    direction = atoms.positions[c] - atoms.positions[ti]
    distance = float(np.linalg.norm(direction))
    direction /= distance
    values = {}
    for model, key in [('DFT', 'REF_forces'), ('POD_v12g', 'POD_v12g_forces')]:
        f = np.asarray(atoms.arrays[key], dtype=float)
        projected = float(np.dot(f[c]/MASS['C'] - f[ti]/MASS['Ti'], direction))
        values[model] = {'relative_acceleration_projection_eV_A_amu': projected,
                         'interpretation': 'negative=net-force tendency to close; positive=net-force tendency to separate'}
    values['DFT_minus_POD_projection_eV_A_amu'] = values['DFT']['relative_acceleration_projection_eV_A_amu'] - values['POD_v12g']['relative_acceleration_projection_eV_A_amu']
    results.append({'name': name, 'center_atom_ids': [ti_id, c_id], 'pair_distance_A': distance,
                    'DFT_center_force_error_eV_A': record['center_POD_v12g_force_errors_eV_A'],
                    **values})
out = {'status': 'AUDITED', 'note': 'Uses full atomic net forces projected onto pair separation; this is not a pairwise bond-force decomposition.', 'pairs': results}
target = HERE/'labels/tic_pair_force_projection_audit.json'
target.write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
