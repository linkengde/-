#!/workspace/venvs/dft-labels/bin/python
"""Extract local DFT crops around four distinct POD v11 Si force hotspots."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import write

HERE = Path(__file__).resolve().parent
SOURCE = Path('/workspace/potential-fit/round7_pod_snap/models/v11_candidate/full_static_run0/POD/forces.dump')
CENTERS = [70588, 104359, 83874, 6223]
RADIUS = 5.5
CELL = 22.0
SYMBOL = {1: 'Ag', 2: 'C', 3: 'Si', 4: 'Ti'}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)
    lines = SOURCE.read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith('ITEM: ATOMS '))
    columns = lines[header].split()[2:]
    ix = {name: i for i, name in enumerate(columns)}
    data = np.loadtxt(lines[header + 1:], dtype=float, ndmin=2)
    data = data[np.argsort(data[:, ix['id']])]
    if len(data) != 347464:
        raise RuntimeError(f'Expected 347,464 rows, found {len(data)}')
    ids = data[:, ix['id']].astype(np.int64)
    types = data[:, ix['type']].astype(np.int32)
    xyz = data[:, [ix['x'], ix['y'], ix['z']]]
    forces = data[:, [ix['fx'], ix['fy'], ix['fz']]]
    id_to_i = {int(atom_id): i for i, atom_id in enumerate(ids)}
    frames, records = [], []
    for idx, center_id in enumerate(CENTERS):
        ci = id_to_i[center_id]
        if types[ci] != 3:
            raise RuntimeError(f'Center {center_id} is not mapped to Si')
        delta = xyz - xyz[ci]
        delta[:, 0] -= 200.0 * np.rint(delta[:, 0] / 200.0)
        delta[:, 1] -= 200.0 * np.rint(delta[:, 1] / 200.0)
        distance = np.linalg.norm(delta, axis=1)
        keep = np.flatnonzero(distance <= RADIUS)
        local = Atoms([SYMBOL[int(types[i])] for i in keep],
                      positions=delta[keep] + CELL / 2,
                      cell=np.eye(3) * CELL, pbc=True)
        local.new_array('source_atom_id', ids[keep].astype(np.int32))
        local.new_array('POD_v11_forces', forces[keep].copy())
        local.info.update({
            'config_type': f'POD_v11_initial_Si_hotspot_{idx:02d}',
            'geometry_source': str(SOURCE),
            'source_timestep': 55,
            'source_model': 'POD v11, full-model static run 0',
            'center_atom_id': int(center_id),
            'center_species': 'Si',
            'center_position_A': xyz[ci].tolist(),
            'center_POD_v11_force_eV_A': forces[ci].tolist(),
            'center_POD_v11_force_norm_eV_A': float(np.linalg.norm(forces[ci])),
            'crop_radius_A': RADIUS,
            'crop_cell_A': CELL,
            'source_pbc': 'x/y periodic, z open; local crop embedded in vacuum box',
            'label_role': 'targeted development DFT; not a blind or independent trajectory test',
            'training_split': 'holdout' if idx == 1 else 'training_candidate',
            'boundary_note': 'Use center force as DFT target; local crop energy and edge forces are excluded from fitting.',
        })
        frames.append(local)
        records.append({'index': idx, 'center_atom_id': int(center_id),
                        'center_species': 'Si', 'center_position_A': xyz[ci].tolist(),
                        'center_force_v11_eV_A': forces[ci].tolist(),
                        'center_force_norm_v11_eV_A': float(np.linalg.norm(forces[ci])),
                        'natoms': len(local), 'formula': local.get_chemical_formula(),
                        'role': 'holdout' if idx == 1 else 'training_candidate'})
    output = HERE / 'candidate_geometries.extxyz'
    write(output, frames, format='extxyz')
    manifest = {'status': 'extracted; DFT labels pending',
                'source_dump': str(SOURCE), 'source_sha256': sha256(SOURCE),
                'source_natoms': len(data), 'source_box_A': [200, 200, 300],
                'crop_radius_A': RADIUS, 'crop_cell_A': CELL,
                'center_ids': CENTERS, 'frames': records,
                'method': 'GPAW PBE/PW350 Gamma, center-force active learning',
                'role': 'three training candidates; ID 104359 held out for v12 development validation'}
    (HERE / 'extraction_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
