"""Extract a second independent cold-state short Ti-C contact crop."""
from pathlib import Path
import json
import numpy as np
from ase import Atoms
from ase.io import write

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
PAIR = (324028, 324034)  # independent Ti-C contact, 1.687 A after 10 fs
RADIUS = 5.5
CELL = 22.0
SOURCE_TYPE_TO_ELEMENT = {1: 'Ti', 2: 'Si', 3: 'C', 4: 'Ag'}
MODEL_TYPE_TO_ELEMENT = {1: 'Ag', 2: 'C', 3: 'Si', 4: 'Ti'}
path = MODEL/'cold_v12g_initial_458k_1fs/end_0p01ps.dump'
lines = path.read_text().splitlines()
h = next(i for i, s in enumerate(lines) if s.startswith('ITEM: ATOMS'))
cols = lines[h].split()[2:]
rows = np.loadtxt(lines[h+1:], dtype=float)
ix = {k: cols.index(k) for k in ('id', 'type', 'x', 'y', 'z', 'fx', 'fy', 'fz')}
data = {int(row[ix['id']]): row for row in rows}
if any(i not in data for i in PAIR):
    raise RuntimeError('Target Ti-C pair absent from cold dump')
pair_types = [int(data[i][ix['type']]) for i in PAIR]
if pair_types != [4, 2]:
    raise RuntimeError(f'Expected Ti-C (LAMMPS type 4,2), found {pair_types}')
p = np.asarray([[data[a][ix['x']], data[a][ix['y']], data[a][ix['z']]] for a in PAIR])
delta = p[1] - p[0]
delta[:2] -= 200*np.rint(delta[:2]/200)
center = p[0] + delta/2
ids = np.asarray(sorted(data), dtype=np.int64)
coords = np.asarray([[data[int(i)][ix['x']], data[int(i)][ix['y']], data[int(i)][ix['z']]] for i in ids])
d = coords - center
d[:, :2] -= 200*np.rint(d[:, :2]/200)
mask = np.linalg.norm(d, axis=1) <= RADIUS
kept = ids[mask]
if not set(PAIR).issubset(set(kept)):
    raise RuntimeError('Crop omitted one or both target atoms')
symbols = [MODEL_TYPE_TO_ELEMENT[int(data[int(i)][ix['type']])] for i in kept]
forces = np.asarray([[data[int(i)][ix['fx']], data[int(i)][ix['fy']], data[int(i)][ix['fz']]] for i in kept])
target = np.isin(kept, PAIR).astype(np.int8)
atoms = Atoms(symbols=symbols, positions=d[mask] + CELL/2, cell=[CELL]*3, pbc=True)
atoms.new_array('source_atom_id', kept.astype(np.int32))
atoms.new_array('POD_v12g_forces', forces)
atoms.new_array('is_DFT_force_target', target)
distance = float(np.linalg.norm(delta))
atoms.info.update({'config_type': 'full200_tic_pair2_cold_10fs', 'source_timestep': 10,
    'center_atom_id': PAIR[0], 'center_species': 'Ti', 'partner_atom_id': PAIR[1],
    'partner_species': 'C', 'center_partner_distance_A': distance,
    'crop_radius_A': RADIUS, 'crop_cell_A': CELL, 'geometry_source': str(path),
    'source_pbc': 'x/y periodic, z open; crop embedded in cubic vacuum box',
    'label_role': 'training candidate; SI01 center 104359 remains held out',
    'DFT_edge_force_note': 'Only target Ti-C pair forces are used for fitting; crop edge forces excluded.',
    'POD_v12g_force_pair_eV_A': json.dumps(forces[target == 1].tolist())})
out = HERE/'candidates/tic_pair2_cold_10fs.extxyz'
write(out, atoms, format='extxyz')
print(json.dumps({'pair_ids': PAIR, 'pair_types': ['Ti', 'C'], 'distance_A': distance,
                  'crop_atoms': len(atoms), 'target_rows': int(target.sum()), 'output': str(out)}, indent=2))
