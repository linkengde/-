#!/usr/bin/env python3
"""Audit 81 force-ranked centers against the full MACE QSGS geometry.

This creates finite 6 A extraction clusters for inspection only. The clusters
have open boundaries and are explicitly not DFT-ready labels.
"""
from collections import Counter, defaultdict
from pathlib import Path
import csv, gzip, json
import numpy as np
from scipy.spatial import cKDTree
from ase import Atoms
from ase.io import write

ROOT = Path(__file__).resolve().parent
DATA = Path('/workspace/-/models/ag_ti3sic2_qsgs_1to1/ag_ti3sic2_qsgs.data.gz')
SHORTLIST = Path('/workspace/-/archive/mace-fit-ace-exploratory-20261119/followup_20261119/acquisition/label_design_center_shortlist.csv')
OUT_CSV = ROOT / 'local_environment_screen.csv'
OUT_XYZ = ROOT / 'screened_clusters_6A.extxyz'
OUT_JSON = ROOT / 'local_environment_screen.json'
BOX = np.array([299.635, 301.857484616085, 258.6115], dtype=float)
TYPE_TO_SYMBOL = {1: 'Ti', 2: 'Si', 3: 'C', 4: 'Ag'}
SYMBOLS = ['Ag', 'C', 'Si', 'Ti']
RADII = (2.5, 3.0, 4.0, 6.0, 8.0)
CUT = 6.0

# Allocate once from the atom count in the data header; row zero is unused.
atom_count = None
with gzip.open(DATA, 'rt') as f:
    for line in f:
        p = line.split()
        if len(p) == 2 and p[1] == 'atoms':
            try:
                atom_count = int(p[0])
                break
            except ValueError:
                pass
if atom_count is None:
    raise ValueError(f'could not read atom count from {DATA}')
coords = np.zeros((atom_count + 1, 3), dtype=np.float64)
types = np.zeros(atom_count + 1, dtype=np.int8)
with gzip.open(DATA, 'rt') as f:
    in_atoms = False
    maxid = 0
    for line in f:
        if line.startswith('Atoms #'):
            in_atoms = True
            continue
        if not in_atoms:
            continue
        p = line.split()
        if len(p) < 5:
            if maxid:
                continue
            continue
        try:
            atom_id, typ = int(p[0]), int(p[1])
            xyz = [float(p[2]), float(p[3]), float(p[4])]
        except ValueError:
            if maxid:
                break
            continue
        if atom_id > atom_count:
            raise ValueError(f'atom ID {atom_id} exceeds header count {atom_count}')
        coords[atom_id] = xyz
        types[atom_id] = typ
        maxid = max(maxid, atom_id)

valid = types > 0
coords = coords[:maxid + 1]
types = types[:maxid + 1]
ids = np.flatnonzero(valid)
xyz = coords[ids]
typ = types[ids]
# A 3Lz box with a shifted coordinate copy implements periodic x/y and open z.
query_xyz = xyz.copy()
query_xyz[:, 2] += BOX[2]
tree = cKDTree(query_xyz, boxsize=[BOX[0], BOX[1], 3 * BOX[2]], compact_nodes=True, balanced_tree=True)
id_to_row = np.full(maxid + 1, -1, dtype=np.int64)
id_to_row[ids] = np.arange(len(ids))

with SHORTLIST.open(encoding='utf-8-sig', newline='') as f:
    candidates = list(csv.DictReader(f))
results = []
cluster_frames = []
signatures = defaultdict(list)
coord_match_count = 0
for cand in candidates:
    aid = int(cand['atom_id'])
    row = int(id_to_row[aid])
    if row < 0:
        raise ValueError(f'candidate atom ID {aid} absent from source geometry')
    center = xyz[row]
    written_xyz = np.array([float(cand[k]) for k in ('x_A', 'y_A', 'z_A')])
    dcoord = np.abs(center - written_xyz)
    dcoord[0] = min(dcoord[0], BOX[0] - dcoord[0])
    dcoord[1] = min(dcoord[1], BOX[1] - dcoord[1])
    coord_match_count += int(np.max(dcoord) < 1e-3)
    qcenter = center.copy(); qcenter[2] += BOX[2]
    near_rows = np.asarray(tree.query_ball_point(qcenter, r=8.0), dtype=np.int64)
    # The tree rows are aligned with ids; reduce z copies to physical neighbors.
    neighbor_ids = ids[near_rows]
    neighbor_ids = np.unique(neighbor_ids)
    delta = coords[neighbor_ids] - center
    delta[:, 0] -= BOX[0] * np.rint(delta[:, 0] / BOX[0])
    delta[:, 1] -= BOX[1] * np.rint(delta[:, 1] / BOX[1])
    distance = np.linalg.norm(delta, axis=1)
    include = (neighbor_ids != aid) & (distance <= 8.0 + 1e-8)
    neighbor_ids, delta, distance = neighbor_ids[include], delta[include], distance[include]
    neighbor_symbols = np.array([TYPE_TO_SYMBOL[int(types[j])] for j in neighbor_ids])
    center_symbol = TYPE_TO_SYMBOL[int(types[aid])]
    shell_counts = {}
    for radius in RADII:
        mask = distance <= radius + 1e-8
        cnt = Counter(neighbor_symbols[mask])
        for symbol in SYMBOLS:
            shell_counts[f'n_{symbol}_le_{radius:g}A'] = cnt.get(symbol, 0)
    nearest_idx = int(np.argmin(distance)) if len(distance) else None
    nearest_d = float(distance[nearest_idx]) if nearest_idx is not None else None
    nearest_symbol = str(neighbor_symbols[nearest_idx]) if nearest_idx is not None else ''
    # Coarse radial-shell descriptor; identical descriptors are likely duplicates,
    # while distinct descriptors still require structural review.
    sig_items = [(str(sym), round(float(d), 1)) for sym, d in zip(neighbor_symbols, distance) if d <= CUT]
    sig = tuple(sorted(sig_items))
    signatures[(center_symbol, sig)].append(aid)
    entry = {
        'atom_id': aid, 'element': center_symbol,
        'model_force_eV_A': float(cand['model_force_eV_A']),
        'xyz_A': center.tolist(), 'z_bin': int(cand['z_bin']),
        'nearest_neighbor_A': nearest_d, 'nearest_neighbor_element': nearest_symbol,
        'neighbors_within_6A': int(np.sum(distance <= CUT)),
        'neighbors_in_5_to_6A_shell': int(np.sum((distance > 5.0) & (distance <= CUT))),
        'coarse_signature': len(sig), **shell_counts,
    }
    results.append(entry)
    # Inspection cluster: center translated to 10 A, 6 A radius, nonperiodic.
    local = distance <= CUT + 1e-8
    local_ids = np.concatenate(([aid], neighbor_ids[local]))
    local_delta = np.vstack((np.zeros(3), delta[local]))
    atoms = Atoms(symbols=[TYPE_TO_SYMBOL[int(types[i])] for i in local_ids],
                  positions=local_delta + 10.0, cell=np.eye(3) * 20.0, pbc=False)
    atoms.info.update({
        'source_atom_id': aid, 'center_element': center_symbol,
        'source_model_force_eV_A': float(cand['model_force_eV_A']),
        'extraction_radius_A': CUT, 'cluster_role': 'inspection only; open-boundary cutout, not DFT-ready',
    })
    atoms.arrays['source_atom_id'] = local_ids.astype(np.int64)
    atoms.arrays['is_center'] = np.array([1] + [0] * (len(local_ids) - 1), dtype=np.int8)
    cluster_frames.append(atoms)

fieldnames = list(results[0])
with OUT_CSV.open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader(); w.writerows(results)
write(OUT_XYZ, cluster_frames, format='extxyz')
unique_sigs = []
for (symbol, sig), atomids in signatures.items():
    unique_sigs.append({'center_element': symbol, 'atom_ids': atomids, 'n_centers': len(atomids),
                        'signature_shell_count': len(sig)})
unique_sigs.sort(key=lambda x: (-x['n_centers'], x['center_element'], x['atom_ids'][0]))
report = {
    'source_data': str(DATA), 'source_atoms': int(len(ids)),
    'candidate_count': len(candidates), 'candidate_coordinates_matching_dump': coord_match_count,
    'boundary_policy': 'periodic x/y, open z, matching the prior static force audit',
    'radii_A': list(RADII), 'local_cutout_radius_A': CUT,
    'local_cutout_note': 'Generated finite clusters are for neighbor inspection only. No DFT calculations should use these cutouts without an embedding/boundary design.',
    'unique_coarse_signatures': len(unique_sigs), 'signature_groups': unique_sigs,
    'center_element_counts': dict(Counter(x['element'] for x in results)),
    'outputs': {'screen_csv': str(OUT_CSV), 'inspection_extxyz': str(OUT_XYZ)},
}
OUT_JSON.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
