#!/usr/bin/env python3
"""Build a 200 x 200 x 300 A Ag/Ti3SiC2 bicontinuous QSGS model.

Four QSGS controls: seed probability, growth probability, target phase mass
ratio, and RNG seed. Ag occupies the connected complementary volume.
"""
from __future__ import annotations
import gzip
import hashlib
import importlib.util
import json
import shutil
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write
from scipy import ndimage, sparse
from scipy.sparse.csgraph import connected_components

OUT = Path(__file__).resolve().parent
BASE_SCRIPT = Path('/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/build_model.py')
spec = importlib.util.spec_from_file_location('qsgs_helpers', BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

SIZE = np.array([200.0, 200.0, 200.0])
VACUUM = 100.0
BOX = np.array([200.0, 200.0, 300.0])
GRID = (100, 100, 100)
# Four stochastic QSGS controls.
SEED_PROBABILITY = 0.0064
GROWTH_PROBABILITY = 0.18
TARGET_AG_TO_SKELETON_MASS_RATIO = 1.0
RNG_SEED = 20261116
CLEARANCE_A = 2.2
base.SIZE = SIZE.copy()
base.GRID = GRID
base.SEED = RNG_SEED
base.OUT = OUT
FIXED_Z = (0.0, 20.0)
THERMOSTAT_Z = (20.0, 40.0)
NEWTON_Z = (40.0, 200.0)
SPUTTER_OBSERVATION_Z = (200.0, 250.0)
EVAPORATION_OBSERVATION_Z = (250.0, 300.0)
MASS = {s: float(atomic_masses[atomic_numbers[s]]) for s in ('Ag', 'Ti', 'Si', 'C')}


def qsgs_growth_ranks():
    """26-neighbor growth, periodic in x/y and open in z."""
    rng = np.random.default_rng(RNG_SEED)
    occupied = rng.random(GRID) < SEED_PROBABILITY
    rank = np.full(GRID, -1, np.int32)
    ids = np.flatnonzero(occupied)
    rng.shuffle(ids)
    rank.flat[ids] = np.arange(len(ids))
    count = len(ids)
    iteration = 0
    offsets = [v for v in __import__('itertools').product((-1, 0, 1), repeat=3) if any(v)]
    while count < occupied.size * 0.82:
        added = np.zeros(GRID, bool)
        for offset in offsets:
            distance = float(np.linalg.norm(offset))
            added |= base.shifted(occupied, offset) & (rng.random(GRID) < GROWTH_PROBABILITY / distance)
        added &= ~occupied
        ids = np.flatnonzero(added)
        rng.shuffle(ids)
        rank.flat[ids] = np.arange(count, count + len(ids))
        occupied |= added
        count += len(ids)
        iteration += 1
    rank[rank < 0] = rank.size
    return rank, iteration


def wrapped_graph(pos: np.ndarray, cutoff: float):
    tree = base.tree(pos)
    edges = tree.query_pairs(cutoff, output_type='ndarray')
    graph = sparse.coo_matrix((np.ones(len(edges), np.uint8),
                               (edges[:, 0], edges[:, 1])),
                              shape=(len(pos), len(pos))).tocsr()
    ncomp, labels = connected_components(graph, directed=False)
    counts = np.bincount(labels)
    return tree, edges, graph, int(ncomp), labels, counts


def spanning_metrics(pos: np.ndarray, cutoff: float):
    tree, edges, graph, ncomp, labels, counts = wrapped_graph(pos, cutoff)
    largest = int(counts.argmax())
    selected = labels == largest
    p = pos[selected]
    # Translate a spanning tree to periodic images and check both lateral winds.
    from scipy.sparse.csgraph import breadth_first_order
    order, parent = breadth_first_order(graph, int(np.flatnonzero(selected)[0]), directed=False)
    order = order[labels[order] == largest]
    images = np.zeros((len(pos), 2), np.int32)
    children = order[1:]
    shifts = -np.rint((pos[children, :2] - pos[parent[children], :2]) / SIZE[:2]).astype(np.int32)
    images[children] = images[parent[children]] + shifts
    edge_shifts = -np.rint((pos[edges[:, 1], :2] - pos[edges[:, 0], :2]) / SIZE[:2]).astype(np.int32)
    cycles = images[edges[:, 1]] - images[edges[:, 0]] - edge_shifts
    cycles = np.unique(cycles[np.any(cycles != 0, axis=1)], axis=0)
    winding_rank = int(np.linalg.matrix_rank(cycles.astype(float))) if len(cycles) else 0
    return p, {
        'components_before_largest_component_selection': ncomp,
        'atoms_in_largest_component': int(counts[largest]),
        'cutoff_A': cutoff,
        'xy_periodic_winding_rank': winding_rank,
        'z_extent_A': [float(p[:, 2].min()), float(p[:, 2].max())],
        'spans_material_thickness': bool(p[:, 2].min() < 3.3 and p[:, 2].max() > SIZE[2] - 3.3),
    }


def qsgs_morphology_metrics(mask: np.ndarray):
    """Voxel-scale line-intercept and interface-density measures."""
    spacing = SIZE / np.asarray(mask.shape)
    means = []
    for axis, periodic in ((0, True), (1, True), (2, False)):
        rows = np.moveaxis(mask, axis, -1).reshape(-1, mask.shape[axis])
        runs = []
        for row in rows:
            if periodic:
                if row.all():
                    runs.append(len(row))
                    continue
                start = (int(np.flatnonzero(~row)[0]) + 1) % len(row)
                row = np.roll(row, -start)
            padded = np.r_[False, row, False]
            edges = np.flatnonzero(padded[1:] != padded[:-1])
            runs.extend((edges[1::2] - edges[::2]).tolist())
        means.append(float(np.mean(runs) * spacing[axis]))
    transitions = sum(int(np.count_nonzero(mask != np.roll(mask, 1, axis=ax))) for ax in (0, 1))
    transitions += int(np.count_nonzero(mask[:, :, 1:] != mask[:, :, :-1]))
    face_area = float(np.prod(spacing[:2]))
    face_area_z = float(spacing[0] * spacing[1])
    interface_area = transitions * face_area
    domain_volume = float(np.prod(SIZE))
    return {
        'mean_skeleton_line_intercept_A_xyz': means,
        'voxel_interface_area_per_material_volume_per_A': float(interface_area / domain_volume),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rank, generations = qsgs_growth_ranks()
    source = base.SOURCE
    parent_cell = make_supercell(read(source), [[1, 0, 0], [1, 2, 0], [0, 0, 1]])
    # Near-commensurate repetitions; x/y are strained by less than 0.7%.
    skeleton_repeats = (65, 38, 12)
    skeleton = parent_cell.repeat(skeleton_repeats)
    raw_skel_cell = skeleton.cell.lengths().copy()
    target_cell = raw_skel_cell.copy()
    target_cell[:2] = SIZE[:2]
    skeleton.set_cell(target_cell, scale_atoms=True)
    skel_valid = skeleton.positions[:, 2] < SIZE[2]
    ti_positions = skeleton.positions[skel_valid].copy()
    ti_symbols = np.asarray(skeleton.get_chemical_symbols())[skel_valid]

    silver_repeats = (49, 49, 49)
    silver = bulk('Ag', 'fcc', a=4.09, cubic=True).repeat(silver_repeats)
    raw_ag_cell = silver.cell.lengths().copy()
    target_ag_cell = raw_ag_cell.copy()
    target_ag_cell[:2] = SIZE[:2]
    silver.set_cell(target_ag_cell, scale_atoms=True)
    ag_valid = silver.positions[:, 2] < SIZE[2]
    ag_positions = silver.positions[ag_valid].copy()
    for positions in (ti_positions, ag_positions):
        positions[:, :2] %= SIZE[:2]

    formula_mass = 3 * MASS['Ti'] + MASS['Si'] + 2 * MASS['C']
    history = []
    low, high = 0.30, 0.75
    best = None
    for _ in range(16):
        fraction = (low + high) / 2.0
        trial_mask, trial_cleanup = base.clean_mask(rank < int(rank.size * fraction))
        in_skeleton = base.at_mask(trial_mask, ti_positions)
        trial_pos, trial_symbols, trial_units, trial_trimmed = base.stoichiometry(
            ti_positions[in_skeleton], ti_symbols[in_skeleton], trial_mask)
        trial_ag = ag_positions[~base.at_mask(trial_mask, ag_positions)]
        distances, _ = base.tree(trial_pos).query(trial_ag, workers=2)
        trial_ag = trial_ag[distances >= CLEARANCE_A]
        mass_ratio = len(trial_ag) * MASS['Ag'] / (trial_units * formula_mass)
        history.append({
            'requested_skeleton_voxel_fraction': fraction,
            'cleaned_skeleton_voxel_fraction': float(trial_mask.mean()),
            'available_Ag_to_skeleton_mass_ratio': float(mass_ratio),
            'formula_units': int(trial_units),
            'Ag_candidates_after_clearance': int(len(trial_ag)),
            'voxel_cleanup_changes': trial_cleanup,
        })
        candidate = (trial_mask, trial_cleanup, trial_pos, trial_symbols,
                     trial_units, trial_trimmed, trial_ag, fraction, mass_ratio)
        if mass_ratio >= TARGET_AG_TO_SKELETON_MASS_RATIO * 1.10:
            if best is None or mass_ratio < best[-1]:
                best = candidate
            low = fraction
        else:
            high = fraction
    selected = best
    if selected is None:
        raise RuntimeError('QSGS mass-ratio search did not provide enough connected Ag candidates')
    (mask, cleanup, skel_positions, skel_symbols, formula_units,
     trimmed, ag_candidates, selected_fraction, available_mass_ratio) = selected
    np.savez_compressed(OUT / 'qsgs_bicontinuous_mask.npz',
                        skeleton=mask, size_A=SIZE, grid=np.array(GRID), seed=RNG_SEED)

    # Keep the largest connected atomic component in each phase. Both phases
    # must span the material thickness and wind around x and y periodically.
    for attempt in range(8):
        skel_keep, skel_conn = base.component_filter(skel_positions, cutoff=3.25)
        skel_positions, skel_symbols = skel_positions[skel_keep], skel_symbols[skel_keep]
        skel_positions, skel_symbols, formula_units, more_trimmed = base.stoichiometry(
            skel_positions, skel_symbols, mask)
        trimmed += more_trimmed
        if skel_keep.all() and more_trimmed == 0:
            break
    ag_connected, ag_available_conn = spanning_metrics(ag_candidates, 3.15)
    desired_ag = round(formula_units * formula_mass * TARGET_AG_TO_SKELETON_MASS_RATIO / MASS['Ag'])
    if len(ag_connected) < desired_ag:
        raise RuntimeError(
            f'Largest connected Ag phase has {len(ag_connected)} atoms, below 1:1 mass target {desired_ag}'
        )
    original_ag = ag_connected.copy()
    distances, _ = base.tree(skel_positions).query(original_ag, workers=2)
    trim_count = len(original_ag) - desired_ag
    if trim_count:
        keep = np.ones(len(original_ag), bool)
        keep[np.argsort(distances, kind='stable')[:trim_count]] = False
        ag_positions = original_ag[keep]
    else:
        ag_positions = original_ag
    ag_keep, _ = base.component_filter(ag_positions, cutoff=3.15)
    if not ag_keep.all():
        ag_positions = ag_positions[ag_keep]
        current_tree = base.tree(ag_positions)
        occupied = current_tree.query(original_ag, workers=2)[0] < 1.0e-7
        original_depth = base.tree(skel_positions).query(original_ag, workers=2)[0]
        while int(occupied.sum()) < desired_ag:
            candidates = np.flatnonzero(~occupied)
            nearest = base.tree(original_ag[occupied]).query(original_ag[candidates], workers=2)[0]
            candidates = candidates[nearest <= 3.15]
            if not len(candidates):
                raise RuntimeError('Cannot restore a connected Ag phase after enforcing 1:1 mass ratio')
            candidates = candidates[np.argsort(-original_depth[candidates], kind='stable')]
            occupied[candidates[:desired_ag - int(occupied.sum())]] = True
        ag_positions = original_ag[occupied]
    if len(ag_positions) != desired_ag:
        raise RuntimeError(f'Ag mass-ratio selection produced {len(ag_positions)} atoms, expected {desired_ag}')
    ag_positions, ag_conn = spanning_metrics(ag_positions, 3.15)
    if len(ag_positions) != desired_ag:
        raise RuntimeError('Ag trimming broke connectivity or periodic winding at the requested 1:1 mass ratio')
    skel_positions, skel_conn = spanning_metrics(skel_positions, 3.25)
    if skel_conn['xy_periodic_winding_rank'] != 2 or not skel_conn['spans_material_thickness']:
        raise RuntimeError(f'Ti3SiC2 phase is not bicontinuous: {skel_conn}')
    if ag_conn['xy_periodic_winding_rank'] != 2 or not ag_conn['spans_material_thickness']:
        raise RuntimeError(f'Ag phase is not bicontinuous: {ag_conn}')

    positions = np.vstack([skel_positions, ag_positions])
    symbols = skel_symbols.tolist() + ['Ag'] * len(ag_positions)
    atoms = Atoms(symbols=symbols, positions=positions,
                  cell=[BOX[0], BOX[1], BOX[2]], pbc=[True, True, False])
    atoms.new_array('phase', np.r_[np.ones(len(skel_positions), np.int32),
                                   np.full(len(ag_positions), 2, np.int32)])
    fixed_skeleton = np.r_[(skel_positions[:, 2] >= FIXED_Z[0]) & (skel_positions[:, 2] < FIXED_Z[1]),
                           np.zeros(len(ag_positions), bool)]
    thermostat_region = (positions[:, 2] >= THERMOSTAT_Z[0]) & (positions[:, 2] < THERMOSTAT_Z[1])
    atoms.new_array('fixed_skeleton', fixed_skeleton.astype(np.int32))
    atoms.new_array('thermostat_region', thermostat_region.astype(np.int32))
    atoms.info.update({
        'model': 'Ag-Ti3SiC2 bicontinuous random QSGS',
        'box_A': BOX.tolist(),
        'material_z_A': [0.0, 200.0],
        'vacuum_z_A': [200.0, 300.0],
        'qsgs_four_controls': {
            'seed_nucleation_probability': SEED_PROBABILITY,
            'growth_probability_prefactor': GROWTH_PROBABILITY,
            'target_Ag_to_Ti3SiC2_mass_ratio': TARGET_AG_TO_SKELETON_MASS_RATIO,
            'random_seed': RNG_SEED,
        },
        'interpretation': 'exploratory unrelaxed starting geometry; not a thermally equilibrated state',
    })
    data_path = OUT / 'Ag_Ti3SiC2_bicontinuous_200x200x300.data'
    xyz_path = OUT / 'Ag_Ti3SiC2_bicontinuous_200x200x300.extxyz'
    write(data_path, atoms, format='lammps-data', atom_style='atomic',
          specorder=['Ti', 'Si', 'C', 'Ag'], masses=True)
    write(xyz_path, atoms, format='extxyz')
    for path in (data_path, xyz_path):
        with path.open('rb') as inp, gzip.open(str(path) + '.gz', 'wb', compresslevel=6) as out:
            shutil.copyfileobj(inp, out)

    counts = {s: int(np.count_nonzero(np.asarray(symbols) == s)) for s in ('Ag', 'Ti', 'Si', 'C')}
    zcoord = positions[:, 2]
    phase_array = np.asarray(atoms.arrays['phase'], dtype=np.int32)
    def layer_counts(low, high):
        selected = (zcoord >= low) & (zcoord < high)
        return {
            'all_atoms': int(selected.sum()),
            'Ag': int(np.count_nonzero(selected & (phase_array == 2))),
            'Ti3SiC2': int(np.count_nonzero(selected & (phase_array == 1))),
        }
    morphology = qsgs_morphology_metrics(mask)
    previous_p0032_intercepts = [19.1346147381, 19.3553636367, 19.0785788975]
    morphology['same_target_p0_0p0032_reference'] = {
        'mean_skeleton_line_intercept_A_xyz': previous_p0032_intercepts,
        'line_intercept_change_fraction_xyz': [
            float(now / old - 1.0)
            for now, old in zip(morphology['mean_skeleton_line_intercept_A_xyz'], previous_p0032_intercepts)
        ],
        'note': 'Comparison against the prior 1:1-mass QSGS candidate with seed probability 0.0032. The prior voxel pitch was about 3 A and the new pitch is 2 A, so this is a morphology indicator rather than a grid-converged or atomistic ligament measurement.',
    }
    result = {
        'status': 'BUILT_GEOMETRY_ONLY_NOT_RELAXED',
        'box_A': BOX.tolist(), 'material_z_A': [0.0, 200.0], 'vacuum_z_A': [200.0, 300.0],
        'boundary': ['periodic', 'periodic', 'fixed_nonperiodic'],
        'layers_A': {
            'fixed_skeleton': list(FIXED_Z),
            '300K_ideal_heat_sink_thermostat': list(THERMOSTAT_Z),
            'newton_dynamics_heating_cooling': list(NEWTON_Z),
            'sputter_observation': list(SPUTTER_OBSERVATION_Z),
            'evaporation_observation': list(EVAPORATION_OBSERVATION_Z),
            'empty_part_of_newton_region_before_heating': [200.0, 200.0],
        },
        'layer_atom_counts_initial': {
            'fixed_slab_z0_20': layer_counts(0.0, 20.0),
            'fixed_skeleton_atoms_z0_20': int(fixed_skeleton.sum()),
            'thermostat_z20_40': layer_counts(20.0, 40.0),
            'Newton_z40_200': layer_counts(40.0, 200.0),
            'sputter_observation_z200_250': layer_counts(200.0, 250.0),
            'evaporation_observation_z250_300': layer_counts(250.0, 300.0),
        },
        'heat_sink_representation': 'Ideal 300 K thermostat slab on the existing Ag-Ti3SiC2 atoms; no explicit Al atoms because the current fitted POD/SNAP models contain Ag, Ti, Si, C only.',
        'qsgs_controls': {
            'seed_nucleation_probability': SEED_PROBABILITY,
            'growth_probability_prefactor': GROWTH_PROBABILITY,
            'directional_rule': '26 neighbors; probability prefactor divided by Euclidean neighbor distance; periodic x/y, open z',
            'target_Ag_to_Ti3SiC2_mass_ratio': TARGET_AG_TO_SKELETON_MASS_RATIO,
            'selected_skeleton_voxel_fraction_bracket': float(selected_fraction),
            'available_mass_ratio_before_Ag_trim': float(available_mass_ratio),
            'random_seed': RNG_SEED,
            'growth_generations': generations,
            'cleanup_voxel_changes': cleanup,
            'measured_cleaned_skeleton_voxel_fraction': float(mask.mean()),
            'morphology_metrics': morphology,
        },
        'atoms': len(atoms), 'counts': counts,
        'chemical_formula_Ti3SiC2_units': formula_units,
        'actual_Ag_to_Ti3SiC2_mass_ratio': float(
            len(ag_positions) * MASS['Ag'] / (formula_units * formula_mass)
        ),
        'target_Ag_to_Ti3SiC2_mass_ratio': TARGET_AG_TO_SKELETON_MASS_RATIO,
        'available_Ag_connectivity_before_mass_trim': ag_available_conn,
        'Ag_atoms_trimmed_to_mass_ratio': int(trim_count),
        'stoichiometry_trimmed_atoms': int(trimmed),
        'atomic_connectivity': {'Ti3SiC2': skel_conn, 'Ag': ag_conn},
        'lattice': {
            'Ti3SiC2_repeats': list(skeleton_repeats),
            'Ti3SiC2_xy_strain': (SIZE[:2] / raw_skel_cell[:2] - 1).tolist(),
            'Ag_repeats': list(silver_repeats),
            'Ag_xy_strain': (SIZE[:2] / raw_ag_cell[:2] - 1).tolist(),
            'Ag_fcc_a_A': 4.09,
            'minimum_Ag_skeleton_clearance_A': CLEARANCE_A,
        },
        'files': {
            'lammps_data_gz': data_path.name + '.gz',
            'extxyz_gz': xyz_path.name + '.gz',
            'qsgs_mask': 'qsgs_bicontinuous_mask.npz',
        },
        'sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (data_path, xyz_path)},
        'limitations': [
            'Geometry build only; no minimization, thermal equilibration, heating, or cooling performed.',
            'The ideal Al-Si heat sink is represented by the requested 300 K thermostat slab, not explicit Al-Si atoms.',
            'The random QSGS mask is bicontinuous at the voxel level; atomic phase connectivity is separately checked.',
            'The present 4-element POD/SNAP potential has not been validated on this new target geometry.',
        ],
    }
    (OUT / 'geometry_report.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('status', 'box_A', 'atoms', 'counts', 'atomic_connectivity')},
                     ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
