#!/usr/bin/env python3
"""QSGS morphology and unrelaxed, equal-mass Ag/Ti3SiC2 atomic candidate.

All results are new geometry candidates, not validated finite-temperature states.
No dynamics, heating, potential fitting or changes to source research files.
"""
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from ase import Atoms
from ase.build import bulk, make_supercell
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write

OUT = Path(__file__).resolve().parent
REPO = next((p for p in OUT.parents if (p/'research/mace-v12-transfer').is_dir()), Path('/workspace/-'))
SOURCE = REPO/'research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif'
POT = REPO/'research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009'
BASIS_DATA = json.loads((OUT/'relaxed_parent_basis.json').read_text())
CELL_LENGTHS = np.array(BASIS_DATA['relaxed_unit_cell_lengths_A'])
REPEATS = tuple(BASIS_DATA['model_repeats'])
SIZE = CELL_LENGTHS * np.array(REPEATS)
GRID = (100, 100, 50)
SEED = 20261024
VACUUM = 100.
FIXED = 10.
CLEARANCE = 2.2
MASS = {s: float(atomic_masses[atomic_numbers[s]]) for s in ('Ag', 'Ti', 'Si', 'C')}
FORMULA_MASS = 3*MASS['Ti']+MASS['Si']+2*MASS['C']


def progress(message):
    print(message, flush=True)


def shifted(a, offset):
    b = np.roll(a, offset, axis=(0, 1, 2))
    dz = offset[2]
    if dz > 0:
        b[:, :, :dz] = False
    elif dz < 0:
        b[:, :, dz:] = False
    return b


def growth_ranks():
    """One reproducible QSGS growth history, not a Gaussian random field.

    Seed probability 0.0032; 26-neighbour attachment probabilities 0.18/r,
    where r is the voxel-neighbour distance. xy wrap; z does not wrap.
    The target fraction is subsequently selected from this growth history.
    Only one growing solid phase, so phase-interaction probability is inactive.
    """
    rng = np.random.default_rng(SEED)
    occupied = rng.random(GRID) < .0032
    rank = np.full(GRID, -1, np.int32)
    ids = np.flatnonzero(occupied)
    rng.shuffle(ids)
    rank.flat[ids] = np.arange(len(ids))
    count = len(ids)
    iteration = 0
    offsets = [v for v in itertools.product((-1, 0, 1), repeat=3) if any(v)]
    while count < occupied.size*.82:
        added = np.zeros(GRID, bool)
        for v in offsets:
            # Independent directional attachment from the previous generation.
            added |= shifted(occupied, v) & (rng.random(GRID) < .18/np.linalg.norm(v))
        added &= ~occupied
        ids = np.flatnonzero(added)
        rng.shuffle(ids)
        rank.flat[ids] = np.arange(count, count+len(ids))
        occupied |= added
        count += len(ids)
        iteration += 1
    rank[rank < 0] = rank.size
    progress(f'QSGS: {iteration} generations, {count} occupied voxels')
    return rank, iteration


def periodic_labels(mask):
    """Face-connected voxel labels with explicit xy seam union, open z."""
    labels, n = ndimage.label(mask)
    parent = np.arange(n+1)
    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for axis in (0, 1):
        a = np.take(labels, 0, axis=axis).ravel()
        b = np.take(labels, -1, axis=axis).ravel()
        pairs = np.unique(np.stack([a, b], axis=1), axis=0)
        for u, v in pairs:
            if u and v:
                ru, rv = root(u), root(v)
                if ru != rv:
                    parent[rv] = ru
    mapping = np.array([root(i) for i in range(n+1)])
    labels = mapping[labels]
    counts = np.bincount(labels.ravel(), minlength=n+1)
    counts[0] = 0
    return labels, counts


def clean_mask(mask):
    changes = []
    for _ in range(8):
        before = mask.copy()
        labels, counts = periodic_labels(mask)
        mask = labels == counts.argmax()
        labels, counts = periodic_labels(~mask)
        mask = ~(labels == counts.argmax())
        changes.append(int(np.count_nonzero(mask != before)))
        if changes[-1] == 0:
            break
    return mask, changes


def at_mask(mask, pos):
    indices = np.floor(pos / SIZE * np.array(GRID)).astype(int)
    indices = np.clip(indices, 0, np.array(GRID)-1)
    return mask[tuple(indices.T)]


def boundary_depth(mask, pos):
    # Tile xy so the depth ranking does not introduce false periodic surfaces.
    tiled = np.tile(mask, (3, 3, 1))
    depth = ndimage.distance_transform_edt(tiled, sampling=SIZE/GRID)
    depth = depth[GRID[0]:2*GRID[0], GRID[1]:2*GRID[1], :]
    indices = np.clip(np.floor(pos/SIZE*np.array(GRID)).astype(int), 0, np.array(GRID)-1)
    return depth[tuple(indices.T)]


def tree(pos):
    return cKDTree(pos, boxsize=[SIZE[0], SIZE[1], 1000.])


def component_filter(pos, cutoff=3.25):
    edges = tree(pos).query_pairs(cutoff, output_type='ndarray')
    graph = sparse.coo_matrix((np.ones(len(edges), np.uint8), (edges[:, 0], edges[:, 1])), shape=(len(pos), len(pos))).tocsr()
    n, labels = connected_components(graph, directed=False)
    counts = np.bincount(labels)
    largest = labels == counts.argmax()
    return largest, {'components': int(n), 'largest_atoms': int(counts.max()), 'discarded_atoms': int((~largest).sum()), 'cutoff_A': cutoff}


def stoichiometry(pos, symbols, mask):
    counts = {s: int(np.sum(symbols == s)) for s in ('Ti', 'Si', 'C')}
    units = min(counts['Ti']//3, counts['Si'], counts['C']//2)
    depth = boundary_depth(mask, pos)
    keep = np.ones(len(pos), bool)
    for s, number in [('Ti', 3), ('Si', 1), ('C', 2)]:
        ids = np.flatnonzero(symbols == s)
        excess = len(ids)-units*number
        if excess:
            # Remove excess atoms at phase boundaries; never replace species.
            order = np.lexsort((ids, depth[ids]))
            keep[ids[order[:excess]]] = False
    return pos[keep], symbols[keep], units, int((~keep).sum())


def main():
    OUT.mkdir(exist_ok=True)
    rank, generations = growth_ranks()
    np.save(OUT/'qsgs_growth_rank.npy', rank)
    cell = make_supercell(read(SOURCE), [[1, 0, 0], [1, 2, 0], [0, 0, 1]])
    cell.set_cell(np.diag(CELL_LENGTHS), scale_atoms=False)
    cell.set_scaled_positions(np.asarray(BASIS_DATA['relaxed_basis_fractional']))
    repeats = REPEATS
    skeleton = cell.repeat(repeats)
    unstrained_lengths = skeleton.cell.lengths().copy()
    target_cell = unstrained_lengths.copy()
    target_cell[:2] = SIZE[:2]
    skeleton.set_cell(target_cell, scale_atoms=True)
    valid = skeleton.positions[:, 2] < SIZE[2]
    ti_pos = skeleton.positions[valid].copy()
    ti_symbols = np.array(skeleton.get_chemical_symbols())[valid]
    silver = bulk('Ag', 'fcc', a=4.09, cubic=True).repeat((73, 74, 39))
    ag_lengths = silver.cell.lengths().copy()
    target_cell = ag_lengths.copy()
    target_cell[:2] = SIZE[:2]
    silver.set_cell(target_cell, scale_atoms=True)
    ag_pos = silver.positions[silver.positions[:, 2] < SIZE[2]].copy()
    for p in (ti_pos, ag_pos):
        p[:, :2] %= SIZE[:2]
    progress(f'Lattice pools: skeleton={len(ti_pos)}, silver={len(ag_pos)}')

    history = []
    low, high = .35, .70
    selected = None
    # Target a little excess available silver, which is trimmed only at boundaries.
    for iteration in range(12):
        fraction = (low+high)/2
        mask, changes = clean_mask(rank < int(rank.size*fraction))
        p = ti_pos[at_mask(mask, ti_pos)]
        sy = ti_symbols[at_mask(mask, ti_pos)]
        p, sy, units, stoich_trim = stoichiometry(p, sy, mask)
        candidates = ag_pos[~at_mask(mask, ag_pos)]
        distances, _ = tree(p).query(candidates, workers=2)
        candidates = candidates[distances >= CLEARANCE]
        ratio = len(candidates)*MASS['Ag']/(units*FORMULA_MASS)
        history.append({'requested_skeleton_voxel_fraction': fraction, 'cleaned_fraction': float(mask.mean()), 'available_mass_ratio_Ag_over_skeleton': ratio, 'stoichiometry_trimmed_atoms': stoich_trim, 'voxel_cleanup_changes': changes})
        progress(f'Mass search {iteration}: phi={mask.mean():.6f}, available mass ratio={ratio:.6f}')
        if ratio >= 1.08:
            low = fraction
            selected = (mask, p, sy, candidates, units)
        else:
            high = fraction
    if selected is None:
        raise RuntimeError('No mass-balanced morphology bracket')
    mask, p, sy, ag, units = selected
    cleanup_atoms = []
    for iteration in range(8):
        keep, info = component_filter(p)
        cleanup_atoms.append({'phase': 'Ti3SiC2', **info})
        p, sy = p[keep], sy[keep]
        p, sy, units, trimmed = stoichiometry(p, sy, mask)
        if info['discarded_atoms'] == 0 and trimmed == 0:
            break
    else:
        raise RuntimeError('Skeleton connectivity/stoichiometry did not stabilize')
    keep, info = component_filter(ag, cutoff=3.15)
    cleanup_atoms.append({'phase': 'Ag', **info})
    ag = ag[keep]
    desired_ag = round(units*FORMULA_MASS/MASS['Ag'])
    if len(ag) < desired_ag:
        raise RuntimeError(f'Connected silver insufficient: {len(ag)} < {desired_ag}')
    # Trim the closest interface Ag first; preserve deeper channel atoms.
    distances, _ = tree(p).query(ag, workers=2)
    original_ag = ag.copy()
    trim = len(ag)-desired_ag
    if trim:
        keep = np.ones(len(ag), bool)
        keep[np.argsort(distances, kind='stable')[:trim]] = False
        ag = ag[keep]
    keep, ag_info = component_filter(ag, cutoff=3.15)
    if not keep.all():
        # Start with the surviving connected interior and regrow only Ag sites
        # adjacent to it. This restores the exact count without isolated atoms.
        ag = ag[keep]
        current_tree = tree(ag)
        occupied = current_tree.query(original_ag, workers=2)[0] < 1.e-7
        original_depth = tree(p).query(original_ag, workers=2)[0]
        while occupied.sum() < desired_ag:
            candidates = np.flatnonzero(~occupied)
            neighbour_distance = tree(original_ag[occupied]).query(original_ag[candidates], workers=2)[0]
            candidates = candidates[neighbour_distance <= 3.15]
            if not len(candidates):
                raise RuntimeError('No connected Ag growth frontier')
            candidates = candidates[np.argsort(-original_depth[candidates], kind='stable')]
            occupied[candidates[:desired_ag-int(occupied.sum())]] = True
        ag = original_ag[occupied]
        keep, ag_info = component_filter(ag, cutoff=3.15)
        assert keep.all()
    keep, ti_info = component_filter(p)
    if not keep.all():
        raise RuntimeError('Skeleton final connectivity failure')

    positions = np.vstack([p, ag])
    symbols = sy.tolist()+['Ag']*len(ag)
    atoms = Atoms(symbols=symbols, positions=positions, cell=[SIZE[0], SIZE[1], SIZE[2]+VACUUM], pbc=[True, True, False])
    atoms.new_array('phase', np.r_[np.ones(len(p), np.int32), np.full(len(ag), 2, np.int32)])
    atoms.new_array('fixed_layer', (positions[:, 2] < FIXED).astype(np.int32))
    atoms.info.update({'status': 'TSCBOP_RELAXED_PARENT_QSGS_UNRELAXED_CANDIDATE', 'seed': SEED, 'material_thickness_A': float(SIZE[2]), 'vacuum_A': VACUUM, 'mass_ratio_target': 'Ag:Ti3SiC2=1:1'})
    write(OUT/'ag_ti3sic2_qsgs.extxyz', atoms)
    write(OUT/'ag_ti3sic2_qsgs.data', atoms, format='lammps-data', atom_style='atomic', specorder=['Ti', 'Si', 'C', 'Ag'], masses=True)
    np.savez_compressed(OUT/'qsgs_morphology.npz', skeleton=mask, size_A=SIZE, seed=SEED)
    pair_minima = {}
    ts = tree(positions)
    pairs = ts.query_pairs(3.3, output_type='ndarray')
    z = np.array(atoms.get_chemical_symbols())
    delta = positions[pairs[:, 0]]-positions[pairs[:, 1]]
    delta[:, :2] -= np.rint(delta[:, :2]/SIZE[:2])*SIZE[:2]
    lengths = np.linalg.norm(delta, axis=1)
    for i, s in enumerate(('Ti', 'Si', 'C', 'Ag')):
        for t in ('Ti', 'Si', 'C', 'Ag')[i:]:
            yes = ((z[pairs[:, 0]] == s)&(z[pairs[:, 1]] == t))|((z[pairs[:, 0]] == t)&(z[pairs[:, 1]] == s))
            if yes.any():
                pair_minima[s+'-'+t] = float(lengths[yes].min())
    connectivity = {}
    for name, pos, info in [('Ti3SiC2', p, ti_info), ('Ag', ag, ag_info)]:
        ed = tree(pos).query_pairs(info['cutoff_A'], output_type='ndarray')
        delta = np.abs(pos[ed[:, 0]]-pos[ed[:, 1]])
        # Single connected phase, each periodic seam crossed, and both z faces reached.
        connectivity[name] = {**info, 'xy_seam_edges': [int(np.sum(delta[:, k] > SIZE[k]/2)) for k in (0, 1)], 'z_min_max_A': [float(pos[:, 2].min()), float(pos[:, 2].max())], 'touches_bottom_top_within_A': 3.3}
        assert all(connectivity[name]['xy_seam_edges'])
        assert pos[:, 2].min() < 3.3 and pos[:, 2].max() > SIZE[2]-3.3
    achieved_ratio = len(ag)*MASS['Ag']/(units*FORMULA_MASS)
    report = {
        'status': 'TSCBOP_RELAXED_PARENT_UNRELAXED_GEOMETRY_CANDIDATE_NOT_A_300K_BASELINE',
        'material_dimensions_A': SIZE.tolist(), 'simulation_box_A': [float(SIZE[0]), float(SIZE[1]), float(SIZE[2]+VACUUM)],
        'pbc': [True, True, False], 'fixed_layer_z_A': [0, FIXED],
        'qsgs': {'seed': SEED, 'grid': GRID, 'voxel_size_A': (SIZE/GRID).tolist(), 'core_probability': .0032, 'growth_probability': '0.18 / neighbour-distance-in-voxels for 26 directions', 'generations': generations, 'interphase_growth_probability': 'inactive: one growing phase and complementary Ag fill', 'final_skeleton_voxel_fraction': float(mask.mean()), 'history': history},
        'atoms': len(atoms), 'counts': {s: symbols.count(s) for s in MASS}, 'Ti3SiC2_formula_units': units,
        'atomic_masses_u': MASS, 'mass_ratio_Ag_over_Ti3SiC2': achieved_ratio,
        'Ag_mass_fraction': achieved_ratio/(1+achieved_ratio),
        'skeleton_lattice': {'source': str(SOURCE), 'sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(), 'orthorhombic_repeats': repeats, 'xy_strain': (SIZE[:2]/unstrained_lengths[:2]-1).tolist(), 'orientation': 'Ti3SiC2 c axis along z; one crystal orientation', 'parent_basis_relaxed_with': 'repository TSCBOP at anisotropic zero normal stress', 'relaxed_unit_cell_lengths_A': CELL_LENGTHS.tolist(), 'parent_basis_max_displacement_A': BASIS_DATA['maximum_parent_basis_displacement_A']},
        'silver_lattice': {'fcc_reference_a_A': 4.09, 'xy_strain': (SIZE[:2]/ag_lengths[:2]-1).tolist(), 'orientation': '[100],[010],[001] along x,y,z'},
        'minimum_cross_phase_clearance_A': CLEARANCE, 'pair_minimum_distances_A': pair_minima,
        'atomic_connectivity': connectivity, 'atomic_cleanup': cleanup_atoms,
        'limitations': ['QSGS 26-neighbour growth is followed by face-connected voxel cleanup; final atomic mask differs from voxels near the interface.', 'Atomic boundary trimming enforces global Ti:Si:C=3:1:2 and equal phase mass; local surface coordination is not relaxed.', 'Finite vacuum 100 A and fixed layer 10 A are initial modelling defaults, not verified for later ejecta or thermal boundary convergence.', 'Global connectivity and periodic-seam crossing do not quantify channel transport or guarantee minimum neck thickness.', 'No minimization, 300 K equilibration, MD, heating, TTM or production potential validation has been performed.', 'Classical calibration-v0 is the repository candidate02 handoff, not authenticated Stage69.']}
    (OUT/'geometry_report.json').write_text(json.dumps(report, indent=2)+'\n')
    progress(f'DONE atoms={len(atoms)}, mass ratio={achieved_ratio:.10f}, minima={pair_minima}')


if __name__ == '__main__':
    main()
