#!/usr/bin/env python3
"""Read exported files back and check mass, atomic topology and periodic winding."""
import hashlib
import json
from pathlib import Path
import numpy as np
from ase.io import read
from scipy import sparse
from scipy.sparse.csgraph import breadth_first_order, connected_components
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
a = read(ROOT/'ag_ti3sic2_qsgs.extxyz')
b = read(ROOT/'ag_ti3sic2_qsgs.data', format='lammps-data', atom_style='atomic', Z_of_type={1:22, 2:14, 3:6, 4:47})
assert len(a) == len(b)
assert np.array_equal(a.numbers, b.numbers)
assert np.allclose(a.positions, b.positions, atol=6e-9, rtol=0)
assert np.allclose(a.cell, b.cell, atol=1e-10, rtol=0)
assert list(a.pbc) == [True, True, False]
size = np.array(a.cell.lengths()); size[2] = a.info['material_thickness_A']
symbols = np.array(a.get_chemical_symbols())
counts = {s: int(np.sum(symbols == s)) for s in ('Ti', 'Si', 'C', 'Ag')}
assert counts['Ti'] == 3*counts['Si'] and counts['C'] == 2*counts['Si']
mass = a.get_masses()
ratio = mass[symbols == 'Ag'].sum()/mass[symbols != 'Ag'].sum()
assert abs(ratio-1) < 5e-6
result = {'export_roundtrip': 'PASS', 'counts': counts, 'mass_ratio_Ag_over_Ti3SiC2': float(ratio), 'phases': {}}
for name, selected, cutoff in [('Ti3SiC2', symbols != 'Ag', 3.25), ('Ag', symbols == 'Ag', 3.15)]:
    pos = a.positions[selected]
    pairs = cKDTree(pos, boxsize=[size[0], size[1], 1000]).query_pairs(cutoff, output_type='ndarray')
    graph = sparse.coo_matrix((np.ones(len(pairs), np.uint8), (pairs[:, 0], pairs[:, 1])), shape=(len(pos), len(pos))).tocsr()
    number, _ = connected_components(graph, directed=False)
    assert number == 1
    order, parent = breadth_first_order(graph, 0, directed=False)
    # Lift a spanning tree to periodic images. Residual edge translations
    # detect noncontractible cycles, a stronger test than crossing a seam.
    images = np.zeros((len(pos), 2), np.int32)
    shifts = np.zeros_like(images)
    children = order[1:]
    shifts[children] = -np.rint((pos[children, :2]-pos[parent[children], :2])/size[:2]).astype(np.int32)
    for v in children:
        images[v] = images[parent[v]]+shifts[v]
    edge_shifts = -np.rint((pos[pairs[:, 1], :2]-pos[pairs[:, 0], :2])/size[:2]).astype(np.int32)
    cycles = images[pairs[:, 1]]-images[pairs[:, 0]]-edge_shifts
    cycles = np.unique(cycles[np.any(cycles != 0, axis=1)], axis=0)
    winding_rank = int(np.linalg.matrix_rank(cycles.astype(float))) if len(cycles) else 0
    assert winding_rank == 2, (name, cycles)
    zrange = [float(pos[:, 2].min()), float(pos[:, 2].max())]
    assert zrange[0] < 3.3 and zrange[1] > size[2]-3.3
    result['phases'][name] = {'atoms': len(pos), 'components': int(number), 'near_neighbour_cutoff_A': cutoff, 'periodic_xy_winding_rank': winding_rank, 'cycle_translations': cycles.tolist(), 'z_min_max_A': zrange, 'thickness_spanning': True}
    print(name, result['phases'][name], flush=True)
result['status'] = 'PASS_GEOMETRY_ONLY_NOT_THERMAL_STABILITY'
(ROOT/'independent_geometry_validation.json').write_text(json.dumps(result, indent=2)+'\n')
print('Export roundtrip, stoichiometry, equal mass and both-phase winding: PASS', flush=True)
