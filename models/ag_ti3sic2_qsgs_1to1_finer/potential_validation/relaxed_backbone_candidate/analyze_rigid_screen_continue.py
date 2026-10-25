#!/usr/bin/env python3
"""Summarize the 200-step Ag-only rigid-scaffold FIRE continuation."""
import gzip
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
BOX = np.asarray(json.loads((HERE / 'geometry_report.json').read_text())['simulation_box_A'], float)
BOX_FOR_TREE = [BOX[0], BOX[1], 1000.]
ELEMENT = {1: 'Ti', 2: 'Si', 3: 'C', 4: 'Ag'}

def load(name):
    with gzip.open(HERE / name, 'rt') as stream:
        return np.loadtxt(stream, skiprows=9)

def topology(position, ids, cutoff):
    pairs = cKDTree(position, boxsize=BOX_FOR_TREE).query_pairs(cutoff, output_type='ndarray')
    graph = sparse.coo_matrix((np.ones(len(pairs) * 2, dtype=np.uint8),
        (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
        shape=(len(position), len(position))).tocsr()
    count, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels)
    largest = int(np.argmax(sizes))
    detached = np.flatnonzero(labels != largest)
    return {'components': int(count), 'largest_component_atoms': int(sizes[largest]),
        'atoms_outside_largest': int(len(detached)), 'cutoff_A': cutoff,
        'detached_ids_sample': [int(x) for x in ids[detached[:20]]]}

initial = load('ag_only_rigid_skeleton_100.dump.gz')
final = load('ag_only_rigid_skeleton_continue_200.dump.gz')
if not np.array_equal(initial[:, :2], final[:, :2]):
    raise RuntimeError('Atom IDs/types changed across the relaxation')
if np.max(np.abs(initial[:, 4] - final[:, 4])) > 5.0:
    raise RuntimeError('Implausible z-coordinate change')
types = final[:, 1].astype(int)
pos = final[:, 2:5].copy()
pos[:, 0] %= BOX[0]
pos[:, 1] %= BOX[1]
forces = np.linalg.norm(final[:, 5:8], axis=1)
delta = final[:, 2:5] - initial[:, 2:5]
delta[:, 0] -= BOX[0] * np.rint(delta[:, 0] / BOX[0])
delta[:, 1] -= BOX[1] * np.rint(delta[:, 1] / BOX[1])
displacement = np.linalg.norm(delta, axis=1)
mobile = initial[:, 4] >= 20.0
sel_ag_mobile = (types == 4) & mobile
sel_skeleton_mobile = (types != 4) & mobile

def force_stats(mask):
    value = forces[mask]
    return {'n': int(len(value)), 'p50_eV_A': float(np.percentile(value, 50)),
        'p95_eV_A': float(np.percentile(value, 95)), 'p99_eV_A': float(np.percentile(value, 99)),
        'max_eV_A': float(np.max(value)), 'count_gt_2_eV_A': int(np.sum(value > 2.0))}

def displacement_stats(mask):
    value = displacement[mask]
    return {'p50_A': float(np.percentile(value, 50)), 'p95_A': float(np.percentile(value, 95)),
        'max_A': float(np.max(value))}

minimum_pairs = {}
for first, second in [('Ag', 'Ti'), ('Ag', 'Si'), ('Ag', 'C'), ('Ti', 'C')]:
    first_type = next(t for t, el in ELEMENT.items() if el == first)
    second_type = next(t for t, el in ELEMENT.items() if el == second)
    tree = cKDTree(pos[types == second_type], boxsize=BOX_FOR_TREE)
    minimum_pairs[first + '-' + second] = float(tree.query(pos[types == first_type], k=1)[0].min())

report = {
    'status': 'STATIC_200_STEP_AG_ONLY_RIGID_SCAFFOLD_CONTINUATION_NOT_A_COLD_BASELINE',
    'steps_in_this_segment': 200,
    'atoms': int(len(final)),
    'mobile_Ag_force_after_unfixing_skeleton': force_stats(sel_ag_mobile),
    'mobile_skeleton_force_after_unfixing_skeleton': force_stats(sel_skeleton_mobile),
    'mobile_Ag_displacement_from_100_step_endpoint': displacement_stats(sel_ag_mobile),
    'minimum_pair_distances_A': minimum_pairs,
    'connectivity': {
        'Ti3SiC2': topology(pos[types != 4], final[types != 4, 0].astype(int), 3.25),
        'Ag': topology(pos[types == 4], final[types == 4, 0].astype(int), 3.15),
    },
    'limitations': [
        'Ag alone was minimized with the full Ti3SiC2 scaffold fixed; no dynamics or heat source.',
        'FIRE reached the 200-iteration cap rather than the force tolerance.',
        'This is a coordinate-relaxation screen under the unvalidated calibration-v0 cross potential, not a production baseline.'
    ]
}
(HERE / 'rigid_scaffold_ag_minimize_continue_metrics.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
