#!/usr/bin/env python3
"""Build a finer QSGS Ag/Ti3SiC2 geometry candidate; no MD is run."""
from __future__ import annotations

import hashlib
import heapq
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.csgraph import connected_components, breadth_first_order
from scipy.spatial import cKDTree
from ase import Atoms
from ase.build import bulk, make_supercell
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
REPO = next(p for p in ROOT.parents if (p / 'research/mace-v12-transfer').is_dir())
BASE = REPO / 'models/ag_ti3sic2_qsgs_1to1_finer/potential_validation/relaxed_backbone_candidate'
SOURCE_CIF = REPO / 'research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif'
REFERENCE_DATA = REPO / 'models/target_200x200x300_qsgs_material200/v14_exploratory_hybrid_ag_20261124/hybrid_458k_skeleton_only_cold/original.data'

SIZE = np.array([200.0, 200.0, 200.0])  # material extent
BOX_SIZE = np.array([200.0, 200.0, 500.0])  # 300 Å vacuum above the material
VACUUM_PADDING = 1000.0
GRID = np.array([200, 200, 200], dtype=int)  # 1 Å voxels for finer QSGS branches
SEED = 20261229
SEED_PROB = 0.0015
GROWTH_PROB = 0.050
BRANCH_PROB = 0.004
BACK_PROB = 0.0004
TARGET_FORMULAS = 38680  # adjusted for connectivity cleanup; about 300k atoms at 1:1 mass
MIN_CROSS_CLEARANCE_A = 2.2
FIXED_LAYER_A = 20.0
THERMOSTAT_LAYER_A = 40.0
AG_A0 = 4.09
TSC_CELL_A = np.array([3.0811793774974063, 5.336851145156566, 17.71017349395825])
TSC_REPEATS = tuple(np.ceil(SIZE / TSC_CELL_A).astype(int))
MASS = {s: float(atomic_masses[atomic_numbers[s]]) for s in ('Ag', 'Ti', 'Si', 'C')}
FORMULA_MASS = 3 * MASS['Ti'] + MASS['Si'] + 2 * MASS['C']
TARGET_AG = int(round(TARGET_FORMULAS * FORMULA_MASS / MASS['Ag']))


def p(message):
    print(message, flush=True)


def shifted(a: np.ndarray, offset: tuple[int, int, int]) -> np.ndarray:
    b = np.roll(a, offset, axis=(0, 1, 2))
    dz = offset[2]
    if dz > 0:
        b[:, :, :dz] = False
    elif dz < 0:
        b[:, :, dz:] = False
    return b


def growth_ranks(seed: int) -> tuple[np.ndarray, int]:
    """Four-parameter directional QSGS: seed, forward, branch, backward growth.

    Each newly added voxel inherits its incoming direction. Forward growth keeps
    a branch moving along its local axis; lateral growth creates daughter twigs.
    x/y are periodic and z is open.
    """
    rng = np.random.default_rng(seed)
    shape = tuple(GRID)
    occupied = rng.random(shape) < SEED_PROB
    offsets = [v for v in itertools.product((-1, 0, 1), repeat=3) if any(v)]
    unit = np.asarray(offsets, dtype=float)
    unit /= np.linalg.norm(unit, axis=1)[:, None]
    heading = np.full(shape, -1, dtype=np.int8)
    seed_ids = np.flatnonzero(occupied)
    heading.flat[seed_ids] = rng.integers(0, len(offsets), len(seed_ids), dtype=np.int8)
    rank = np.full(shape, -1, np.int32)
    rng.shuffle(seed_ids)
    count = len(seed_ids)
    rank.flat[seed_ids] = np.arange(count, dtype=np.int32)
    generation = 0
    target = int(occupied.size * 0.82)

    while count < target:
        added = np.zeros(shape, dtype=bool)
        new_heading = np.full(shape, -1, dtype=np.int8)
        for k, offset in enumerate(offsets):
            source = shifted(occupied, offset)
            source_heading = np.roll(heading, offset, axis=(0, 1, 2))
            dz = offset[2]
            if dz > 0:
                source_heading[:, :, :dz] = -1
            elif dz < 0:
                source_heading[:, :, dz:] = -1
            valid = source & (source_heading >= 0) & ~occupied
            if not np.any(valid):
                continue
            cosine = np.zeros(shape, dtype=np.float32)
            hv = source_heading[valid]
            cosine[valid] = unit[hv] @ unit[k]
            probability = np.where(cosine > 0.5, GROWTH_PROB,
                         np.where(cosine >= -0.5, BRANCH_PROB, BACK_PROB))
            success = valid & (rng.random(shape) < probability)
            take = success & ~added
            new_heading[take] = k
            added |= success
        ids = np.flatnonzero(added)
        if not len(ids):
            raise RuntimeError('Directional QSGS growth stopped before its rank pool target.')
        rng.shuffle(ids)
        rank.flat[ids] = np.arange(count, count + len(ids), dtype=np.int32)
        heading[added] = new_heading[added]
        occupied |= added
        count += len(ids)
        generation += 1
        if generation % 10 == 0:
            p(f'Directional QSGS generation {generation}: {count / occupied.size:.1%} occupied')
    rank[rank < 0] = rank.size
    p(f'Directional QSGS rank pool ready: {generation} generations, {count:,} voxels')
    return rank, generation


def periodic_labels(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    labels, n = ndimage.label(mask)
    parent = np.arange(n + 1)

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for axis in (0, 1):
        aa = np.take(labels, 0, axis=axis).ravel()
        bb = np.take(labels, -1, axis=axis).ravel()
        for u, v in np.unique(np.stack((aa, bb), axis=1), axis=0):
            if u and v:
                ru, rv = root(int(u)), root(int(v))
                if ru != rv:
                    parent[rv] = ru
    mapping = np.array([root(i) for i in range(n + 1)])
    labels = mapping[labels]
    counts = np.bincount(labels.ravel(), minlength=n + 1)
    counts[0] = 0
    return labels, counts


def clean_mask(mask: np.ndarray) -> np.ndarray:
    # Keep the largest connected QSGS skeleton and the largest connected pore.
    # The previous iterative complement operation filled enclosed channels and
    # made the skeleton look like bulky, solid clumps.
    labels, counts = periodic_labels(mask)
    counts[0] = 0
    largest_skeleton = int(np.argmax(counts))
    mask = labels == largest_skeleton
    void_labels, void_counts = periodic_labels(~mask)
    void_counts[0] = 0
    largest_void = int(np.argmax(void_counts))
    # Only isolated, non-percolating void pockets are assigned to the skeleton.
    # Preserve the connected pore network for the Ag phase.
    mask |= (void_labels != 0) & (void_labels != largest_void)
    return mask


def inside_mask(mask: np.ndarray, pos: np.ndarray) -> np.ndarray:
    ijk = np.floor(pos / SIZE * GRID).astype(np.int64)
    ijk = np.clip(ijk, 0, GRID - 1)
    return mask[ijk[:, 0], ijk[:, 1], ijk[:, 2]]


def tree(pos: np.ndarray) -> cKDTree:
    return cKDTree(pos, boxsize=[SIZE[0], SIZE[1], VACUUM_PADDING])


def connected_component(pos: np.ndarray, cutoff: float) -> tuple[np.ndarray, dict]:
    pairs = tree(pos).query_pairs(cutoff, output_type='ndarray')
    graph = sparse.coo_matrix((np.ones(len(pairs), np.uint8),
                               (pairs[:, 0], pairs[:, 1])),
                              shape=(len(pos), len(pos))).tocsr()
    graph = graph + graph.T
    n, labels = connected_components(graph, directed=False)
    counts = np.bincount(labels)
    keep = labels == int(np.argmax(counts))
    return keep, {'components': int(n), 'largest_atoms': int(counts.max()), 'cutoff_A': cutoff,
                  'edges': int(len(pairs))}


def winding_rank(pos: np.ndarray, cutoff: float) -> dict:
    pairs = tree(pos).query_pairs(cutoff, output_type='ndarray')
    graph = sparse.coo_matrix((np.ones(len(pairs), np.uint8),
                               (pairs[:, 0], pairs[:, 1])),
                              shape=(len(pos), len(pos))).tocsr()
    graph = graph + graph.T
    n, labels = connected_components(graph, directed=False)
    largest_label = int(np.argmax(np.bincount(labels)))
    keep = labels == largest_label
    sub = graph[keep][:, keep]
    points = pos[keep]
    order, parents = breadth_first_order(sub, 0, directed=False, return_predecessors=True)
    images = np.zeros((len(points), 2), dtype=np.int32)
    for v in order[1:]:
        par = parents[v]
        jump = -np.rint((points[v, :2] - points[par, :2]) / SIZE[:2]).astype(np.int32)
        images[v] = images[par] + jump
    # Rebuild minimal-image edge shifts for all edges inside the retained component.
    original_ids = np.flatnonzero(keep)
    local_map = np.full(len(pos), -1, dtype=np.int64)
    local_map[original_ids] = np.arange(len(original_ids))
    ep = pairs[keep[pairs[:, 0]] & keep[pairs[:, 1]]]
    ii, jj = local_map[ep[:, 0]], local_map[ep[:, 1]]
    shifts = -np.rint((pos[ep[:, 1], :2] - pos[ep[:, 0], :2]) / SIZE[:2]).astype(np.int32)
    cycles = images[jj] - images[ii] - shifts
    cycles = np.unique(cycles[np.any(cycles != 0, axis=1)], axis=0)
    rank = int(np.linalg.matrix_rank(cycles.astype(float))) if len(cycles) else 0
    return {'components': int(n), 'largest_atoms': int(keep.sum()), 'winding_rank_xy': rank,
            'z_range_A': [float(pos[keep, 2].min()), float(pos[keep, 2].max())],
            'spans_z': bool(pos[keep, 2].min() < 3.3 and pos[keep, 2].max() > SIZE[2] - 3.3)}


def build_tsc_lattice() -> tuple[np.ndarray, np.ndarray]:
    basis_doc = json.loads((BASE / 'relaxed_parent_basis.json').read_text())
    lengths = np.asarray(basis_doc['relaxed_unit_cell_lengths_A'], dtype=float)
    basis = make_supercell(read(SOURCE_CIF), [[1, 0, 0], [1, 2, 0], [0, 0, 1]])
    basis.set_cell(np.diag(lengths), scale_atoms=False)
    basis.set_scaled_positions(np.asarray(basis_doc['relaxed_basis_fractional'], dtype=float))
    repeated = basis.repeat(TSC_REPEATS)
    repeated.set_cell(np.diag(SIZE), scale_atoms=True)
    pos = repeated.positions.copy()
    syms = np.asarray(repeated.get_chemical_symbols())
    valid = np.all(pos < SIZE[None, :], axis=1)
    return pos[valid], syms[valid]


def boundary_depth(mask: np.ndarray, pos: np.ndarray) -> np.ndarray:
    tiled = np.tile(mask, (3, 3, 1))
    depth = ndimage.distance_transform_edt(tiled, sampling=SIZE / GRID)
    depth = depth[GRID[0]:2 * GRID[0], GRID[1]:2 * GRID[1], :]
    ijk = np.clip(np.floor(pos / SIZE * GRID).astype(np.int64), 0, GRID - 1)
    return depth[ijk[:, 0], ijk[:, 1], ijk[:, 2]]


def select_tsc(mask: np.ndarray, pool: np.ndarray, symbols: np.ndarray):
    inside = inside_mask(mask, pool)
    pos, syms = pool[inside], symbols[inside]
    dep = boundary_depth(mask, pos)
    kept = np.ones(len(pos), dtype=bool)
    target_by_species = {'Ti': 3 * TARGET_FORMULAS, 'Si': TARGET_FORMULAS, 'C': 2 * TARGET_FORMULAS}
    for symbol, target in target_by_species.items():
        ids = np.flatnonzero(syms == symbol)
        if len(ids) < target:
            return None, len(ids)
        excess = len(ids) - target
        if excess:
            # Preserve the fine outer network while trimming excess atoms from
            # thick interior regions (the earlier ascending order kept bulky cores).
            remove = ids[np.lexsort((ids, -dep[ids]))[:excess]]
            kept[remove] = False
    psel, ssel = pos[kept], syms[kept]
    return (psel, ssel), int(len(psel) // 6)


def select_ag(mask: np.ndarray, skeleton: np.ndarray, ag_pool: np.ndarray, target_ag: int):
    candidates = ag_pool[~inside_mask(mask, ag_pool)]
    dsk, _ = tree(skeleton).query(candidates, workers=4)
    candidates = candidates[dsk >= MIN_CROSS_CLEARANCE_A]
    if len(candidates) < target_ag:
        return None, {'candidate_count': int(len(candidates)), 'reason': 'not_enough_clear_sites'}

    # Start from the largest connected pore-space Ag candidate component.
    pairs = tree(candidates).query_pairs(3.15, output_type='ndarray')
    graph = sparse.coo_matrix((np.ones(len(pairs), np.uint8),
                               (pairs[:, 0], pairs[:, 1])),
                              shape=(len(candidates), len(candidates))).tocsr()
    graph = graph + graph.T
    ncomp, labels = connected_components(graph, directed=False)
    component_sizes = np.bincount(labels)
    largest_label = int(np.argmax(component_sizes))
    component = labels == largest_label
    original = candidates[component]
    if len(original) < target_ag:
        return None, {'candidate_count': int(len(candidates)), 'largest_component': int(len(original)),
                      'reason': 'largest_Ag_component_too_small'}
    dsk = tree(skeleton).query(original, workers=4)[0]

    # Build a small connected Ag backbone that spans z and preserves two
    # independent x/y periodic winding loops. Then grow a connected Ag fill
    # toward the Ti3SiC2 surfaces until the mass-ratio target is reached.
    ids = np.flatnonzero(component)
    local_map = np.full(len(candidates), -1, dtype=np.int64)
    local_map[ids] = np.arange(len(ids))
    edge_keep = component[pairs[:, 0]] & component[pairs[:, 1]]
    ep = local_map[pairs[edge_keep]]
    sub = sparse.coo_matrix((np.ones(len(ep), np.uint8), (ep[:, 0], ep[:, 1])),
                            shape=(len(original), len(original))).tocsr()
    sub = sub + sub.T
    order, parents = breadth_first_order(sub, 0, directed=False, return_predecessors=True)
    if len(order) != len(original):
        return None, {'candidate_count': int(len(candidates)), 'largest_component': int(len(original)),
                      'reason': 'largest_component_bfs_incomplete'}

    def tree_path(start: int, end: int) -> list[int]:
        ancestors = {}
        node = int(start)
        while node >= 0 and node not in ancestors:
            ancestors[node] = len(ancestors)
            par = int(parents[node])
            node = par if par != node else -1
        tail = []
        node = int(end)
        while node not in ancestors:
            if node < 0:
                return []
            tail.append(node)
            par = int(parents[node])
            node = par if par != node else -1
        lca = node
        head = []
        node = int(start)
        while node != lca:
            head.append(node)
            par = int(parents[node])
            node = par if par != node else -1
        head.append(lca)
        return head + list(reversed(tail))

    seed_nodes: set[int] = set()
    z_low = int(np.argmin(original[:, 2]))
    z_high = int(np.argmax(original[:, 2]))
    seed_nodes.update(tree_path(z_low, z_high))

    # Find independent periodic cycle vectors from the BFS tree and add their
    # fundamental paths, so the selected Ag component remains bicontinuous.
    images = np.zeros((len(original), 2), dtype=np.int32)
    for node in order[1:]:
        par = int(parents[node])
        jump = -np.rint((original[node, :2] - original[par, :2]) / SIZE[:2]).astype(np.int32)
        images[node] = images[par] + jump
    selected_cycles: list[np.ndarray] = []
    local_edges = ep
    shifts = -np.rint((original[local_edges[:, 1], :2] - original[local_edges[:, 0], :2]) /
                      SIZE[:2]).astype(np.int32)
    cycles = images[local_edges[:, 1]] - images[local_edges[:, 0]] - shifts
    for edge, cycle in zip(local_edges, cycles):
        if not np.any(cycle):
            continue
        trial = selected_cycles + [cycle]
        if np.linalg.matrix_rank(np.asarray(trial, dtype=float)) > len(selected_cycles):
            seed_nodes.update(tree_path(int(edge[0]), int(edge[1])))
            selected_cycles.append(cycle.copy())
            if len(selected_cycles) == 2:
                break
    if len(selected_cycles) < 2:
        return None, {'candidate_count': int(len(candidates)), 'largest_component': int(len(original)),
                      'reason': 'Ag_candidate_component_lacks_xy_winding'}
    if not seed_nodes or len(seed_nodes) > target_ag:
        return None, {'candidate_count': int(len(candidates)), 'largest_component': int(len(original)),
                      'seed_backbone_atoms': int(len(seed_nodes)), 'reason': 'Ag_seed_backbone_invalid'}

    occupied = np.zeros(len(original), dtype=bool)
    occupied[np.fromiter(seed_nodes, dtype=np.int64)] = True
    frontier: list[tuple[float, int]] = []
    queued = occupied.copy()

    def push_neighbors(nodes):
        for node in nodes:
            start, stop = sub.indptr[node], sub.indptr[node + 1]
            for neighbor in sub.indices[start:stop]:
                if not queued[neighbor]:
                    queued[neighbor] = True
                    heapq.heappush(frontier, (float(dsk[neighbor]), int(neighbor)))

    push_neighbors(np.flatnonzero(occupied))
    while int(occupied.sum()) < target_ag and frontier:
        _, node = heapq.heappop(frontier)
        if occupied[node]:
            continue
        occupied[node] = True
        push_neighbors([node])
    if int(occupied.sum()) != target_ag:
        return None, {'candidate_count': int(len(candidates)), 'largest_component': int(len(original)),
                      'seed_backbone_atoms': int(len(seed_nodes)), 'selected_atoms': int(occupied.sum()),
                      'reason': 'Ag_connected_growth_stalled'}
    selected = original[occupied]
    metrics = {'candidate_count': int(len(candidates)), 'largest_component': int(len(original)),
               'seed_backbone_atoms': int(len(seed_nodes)),
               'winding_seed_vectors': [v.tolist() for v in selected_cycles],
               'selected_atoms': int(len(selected)),
               'Ag_to_skeleton_distance_percentiles_A': {
                   str(q): float(np.percentile(dsk[occupied], q)) for q in (50, 75, 90, 95)}}
    return selected, metrics


def ref_nearest_phase_metrics(path: Path) -> dict:
    if not path.exists():
        return {}
    syms, chunks = [], []
    section = False
    with path.open() as f:
        for line in f:
            if line.startswith('Atoms'):
                section = True
                continue
            if not section or not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 5:
                if chunks:
                    break
                continue
            try:
                typ = int(fields[1])
                xyz = list(map(float, fields[2:5]))
            except ValueError:
                if chunks:
                    break
                continue
            syms.append(typ)
            chunks.append(xyz)
    types = np.asarray(syms)
    xyz = np.asarray(chunks)
    if not len(xyz):
        return {}
    ag = xyz[types == 4]
    tsc = xyz[types != 4]
    dist, _ = tree(ag).query(tsc, workers=4)
    return {'atoms': int(len(xyz)), 'TSC_atoms': int(len(tsc)), 'Ag_atoms': int(len(ag)),
            'TSC_to_nearest_Ag_distance_percentiles_A': {
                str(q): float(np.percentile(dist, q)) for q in (50, 75, 90, 95)}}


def clean_atomic_stoichiometry(pos: np.ndarray, symbols: np.ndarray):
    """Drop atom-isolated boundary debris, then restore exact Ti:Si:C stoichiometry."""
    history = []
    for iteration in range(8):
        largest, info = connected_component(pos, 3.25)
        if not largest.all():
            history.append({'iteration': iteration, 'components_before': info['components'],
                            'discarded_fragment_atoms': int((~largest).sum())})
            pos, symbols = pos[largest], symbols[largest]
        counts = {s: int(np.count_nonzero(symbols == s)) for s in ('Ti', 'Si', 'C')}
        formulas = min(counts['Ti'] // 3, counts['Si'], counts['C'] // 2)
        targets = {'Ti': 3 * formulas, 'Si': formulas, 'C': 2 * formulas}
        excess = {s: counts[s] - targets[s] for s in counts}
        excess = {s: n for s, n in excess.items() if n > 0}
        if not excess:
            final, final_info = connected_component(pos, 3.25)
            if final.all():
                return pos, symbols, formulas, history, final_info
            continue

        pairs = tree(pos).query_pairs(3.25, output_type='ndarray')
        graph = sparse.coo_matrix((np.ones(len(pairs), np.int8),
                                   (pairs[:, 0], pairs[:, 1])),
                                  shape=(len(pos), len(pos))).tocsr()
        graph = graph + graph.T
        degree = np.asarray(graph.sum(axis=1)).ravel()
        remove = np.zeros(len(pos), bool)
        for symbol, nremove in excess.items():
            candidates = np.flatnonzero(symbols == symbol)
            candidates = candidates[np.lexsort((candidates, degree[candidates]))]
            remove[candidates[:nremove]] = True
        history.append({'iteration': iteration, 'components_before': info['components'],
                        'removed_for_stoichiometry': {s: int(n) for s, n in excess.items()}})
        pos, symbols = pos[~remove], symbols[~remove]
    raise RuntimeError('Could not retain one connected, stoichiometric Ti3SiC2 atom network.')


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    p(f'Building TSC lattice pool: repeats={TSC_REPEATS}')
    tsc_pool, tsc_symbols = build_tsc_lattice()
    p(f'TSC lattice sites={len(tsc_pool):,}; target formulas={TARGET_FORMULAS:,}')
    ag_crystal = bulk('Ag', 'fcc', a=AG_A0, cubic=True).repeat((49, 49, 74))
    ag_crystal.set_cell(np.diag(SIZE), scale_atoms=True)
    ag_pool = ag_crystal.positions.copy()
    p(f'Ag FCC lattice sites={len(ag_pool):,}')

    rank, generations = growth_ranks(SEED)
    np.savez_compressed(ROOT / 'qsgs_rank_and_parameters.npz', rank=rank, grid=GRID,
                        size_A=SIZE, seed=SEED, seed_probability=SEED_PROB,
                        growth_probability=GROWTH_PROB, branch_probability=BRANCH_PROB, backward_probability=BACK_PROB)

    # Adjust the QSGS threshold to match the selected model's skeleton formula count.
    frac = TARGET_FORMULAS / (len(tsc_pool) / 6)
    best = None
    for attempt in range(12):
        mask = clean_mask(rank < int(rank.size * frac))
        chosen, units = select_tsc(mask, tsc_pool, tsc_symbols)
        if chosen is not None:
            skeleton, skeleton_syms = chosen
            best = (mask, skeleton, skeleton_syms, units, frac)
            p(f'QSGS threshold attempt {attempt + 1}: requested={frac:.5f}, cleaned_voxel_fraction={mask.mean():.5f}, formulas={units:,}')
            if units == TARGET_FORMULAS:
                break
            # The selected atom counts are exact when each species has enough sites;
            # this loop records topology and provides a smooth threshold fallback.
            break
        p(f'QSGS threshold attempt {attempt + 1}: fraction={frac:.5f} has fewer than target stoichiometric sites')
        frac = min(0.75, frac * 1.12 + 0.005)
    if best is None:
        raise RuntimeError('QSGS mask does not contain enough connected stoichiometric skeleton sites.')
    mask, skeleton, skeleton_syms, _, selected_frac = best

    # Preserve atomic stoichiometry and exact 1:1 phase mass. The QSGS mask is
    # cleaned to single connected components in both phases before atom mapping.
    skeleton, skeleton_syms, final_formulas, skeleton_cleanup, skeleton_info = clean_atomic_stoichiometry(skeleton, skeleton_syms)
    target_ag = int(round(final_formulas * FORMULA_MASS / MASS['Ag']))
    tsc_topology = winding_rank(skeleton, 3.25)
    p(f'Ti3SiC2 atomic topology: {tsc_topology}')
    if tsc_topology['components'] != 1 or tsc_topology['winding_rank_xy'] < 2 or not tsc_topology['spans_z']:
        raise RuntimeError('Fine QSGS mask did not produce a 3D-connected atomistic skeleton; candidate not exported.')

    ag, ag_selection = select_ag(mask, skeleton, ag_pool, target_ag)
    if ag is None:
        raise RuntimeError(f'Could not select connected Ag phase: {ag_selection}')
    ag_topology = winding_rank(ag, 3.15)
    p(f'Ag atomic topology: {ag_topology}')
    if ag_topology['components'] != 1 or ag_topology['winding_rank_xy'] < 2 or not ag_topology['spans_z']:
        raise RuntimeError('Selected Ag phase is not bicontinuous/spanning; candidate not exported.')

    positions = np.vstack((skeleton, ag))
    symbols = np.concatenate((skeleton_syms, np.full(len(ag), 'Ag')))
    atoms = Atoms(symbols=symbols, positions=positions,
                  cell=np.diag(BOX_SIZE), pbc=[True, True, False])
    phase = np.r_[np.ones(len(skeleton), np.int32), np.full(len(ag), 2, np.int32)]
    z = positions[:, 2]
    region_id = np.full(len(positions), 3, dtype=np.int32)
    region_id[z < FIXED_LAYER_A] = 0
    region_id[(z >= FIXED_LAYER_A) & (z < THERMOSTAT_LAYER_A)] = 1
    region_id[(z >= THERMOSTAT_LAYER_A) & (z < SIZE[2])] = 2
    atoms.new_array('phase', phase)
    atoms.new_array('z_region', region_id)
    # The geometry labels the 0-20 A slab for every species. The planned
    # mechanical fixed group includes only the ceramic skeleton; Ag stays mobile.
    atoms.new_array('fixed_skeleton', ((z < FIXED_LAYER_A) & (phase == 1)).astype(np.int32))
    atoms.new_array('thermostat_skeleton', ((z >= FIXED_LAYER_A) & (z < THERMOSTAT_LAYER_A) & (phase == 1)).astype(np.int32))
    atoms.new_array('ag_mobile', (phase == 2).astype(np.int32))
    atoms.new_array('newton_region', (region_id == 2).astype(np.int32))
    atoms.info.update({'status': 'GEOMETRY_CANDIDATE_NOT_MINIMIZED_OR_MD_VALIDATED',
                       'qsgs_seed': SEED, 'box_A': BOX_SIZE.tolist(),
                       'material_size_A': SIZE.tolist(),
                       'mass_ratio_Ag_to_Ti3SiC2': '1:1',
                       'z_regions_A': '0-20 fixed; 20-40 thermostat; 40-200 Newton; 200-500 vacuum',
                       'planned_group_policy': 'fixed and thermostat apply to Ti3SiC2 only; all Ag remains mobile'})
    write(ROOT / 'skeleton_microtwig_200x200x500.data', atoms, format='lammps-data',
          atom_style='atomic', specorder=['Ti', 'Si', 'C', 'Ag'], masses=True)
    write(ROOT / 'skeleton_microtwig_200x200x500.extxyz', atoms)
    np.savez_compressed(ROOT / 'qsgs_skeleton_microtwig.npz', skeleton=mask,
                        size_A=SIZE, seed=SEED, grid=GRID)

    masses = atoms.get_masses()
    phase = atoms.arrays['phase']
    ratio = float(masses[phase == 2].sum() / masses[phase == 1].sum())
    counts = {s: int(np.count_nonzero(symbols == s)) for s in ('Ti', 'Si', 'C', 'Ag')}
    z_regions = {
        'fixed_0_20_A': (z < FIXED_LAYER_A),
        'thermostat_20_40_A': ((z >= FIXED_LAYER_A) & (z < THERMOSTAT_LAYER_A)),
        'newton_40_200_A': ((z >= THERMOSTAT_LAYER_A) & (z < SIZE[2])),
        'vacuum_200_500_A': ((z >= SIZE[2]) & (z < BOX_SIZE[2])),
    }
    layer_counts = {}
    for name, selected in z_regions.items():
        layer_counts[name] = {
            'all_atoms': int(selected.sum()),
            'Ag_atoms': int(np.count_nonzero(selected & (phase == 2))),
            'Ti3SiC2_atoms': int(np.count_nonzero(selected & (phase == 1))),
        }
    nearest = tree(ag).query(skeleton, workers=4)[0]
    fine_metric = {str(q): float(np.percentile(nearest, q)) for q in (50, 75, 90, 95)}
    mask_dist = ndimage.distance_transform_edt(mask, sampling=SIZE / GRID)
    mask_radius = np.percentile(mask_dist[mask], [50, 75, 90, 95, 99])
    old_mask_doc = np.load(BASE / 'qsgs_morphology.npz')
    old_mask = old_mask_doc['skeleton'].astype(bool)
    old_spacing = old_mask_doc['size_A'] / np.asarray(old_mask.shape)
    old_dist = ndimage.distance_transform_edt(old_mask, sampling=old_spacing)
    old_radius = np.percentile(old_dist[old_mask], [50, 75, 90, 95, 99])
    report = {
        'status': 'GEOMETRY_CANDIDATE_ONLY_NOT_MINIMIZED_OR_MD_VALIDATED',
        'box_A': BOX_SIZE.tolist(), 'material_dimensions_A': SIZE.tolist(),
        'vacuum_above_material_A': float(BOX_SIZE[2] - SIZE[2]),
        'boundary_intended': 'periodic x/y; open z', 'atom_count': int(len(atoms)),
        'z_layers_A': {'fixed': [0.0, FIXED_LAYER_A],
                       'thermostat': [FIXED_LAYER_A, THERMOSTAT_LAYER_A],
                       'newton': [THERMOSTAT_LAYER_A, float(SIZE[2])],
                       'vacuum': [float(SIZE[2]), float(BOX_SIZE[2])]},
        'z_layer_atom_counts': layer_counts,
        'planned_group_policy': {'fixed': 'Ti3SiC2 atoms only in z=0-20 A; Ag remains mobile',
                                 'thermostat': 'Ti3SiC2 atoms only in z=20-40 A',
                                 'newton': 'all atoms in z=40-200 A; all Ag atoms remain mobile',
                                 'vacuum': 'z=200-500 A; no atoms initially'},
        'atom_counts': counts, 'TSC_formula_units': int(final_formulas),
        'mass_ratio_Ag_over_Ti3SiC2': ratio,
    'target_baseline_atom_count': 300000,
        'QSGS': {'seed': SEED, 'grid': GRID.tolist(), 'voxel_A': (SIZE / GRID).tolist(),
                 'seed_probability': SEED_PROB, 'growth_probability_rule': 'directional branch growth over 26 neighbours',
                 'forward_probability': GROWTH_PROB, 'lateral_branch_probability': BRANCH_PROB, 'backward_probability': BACK_PROB,
                 'growth_generations_to_82pct': generations, 'threshold_fraction': float(selected_frac),
                 'cleaned_skeleton_voxel_fraction': float(mask.mean()),
                 'both_phase_voxel_components_cleaned': True},
        'finer_skeleton_proxy': {'metric': 'distance from every Ti3SiC2 atom to nearest Ag atom',
                                 'percentiles_A': fine_metric,
                                 'actual_minimum_TSC_Ag_distance_A': float(nearest.min()),
                                 'comparison_reference': ref_nearest_phase_metrics(REFERENCE_DATA)},
        'QSGS_local_thickness_proxy': {
            'method': 'Euclidean distance transform inside the binary skeleton mask; reported as twice the distance-radius percentiles',
            'new_mask_radius_percentiles_A': {str(q): float(v) for q, v in zip((50, 75, 90, 95, 99), mask_radius)},
            'new_mask_diameter_P90_A': float(2 * mask_radius[2]),
            'repository_reference_radius_percentiles_A': {str(q): float(v) for q, v in zip((50, 75, 90, 95, 99), old_radius)},
            'repository_reference_diameter_P90_A': float(2 * old_radius[2]),
            'reference_note': 'Reference is the earlier 100x100x50 QSGS mask on a different box; this compares morphology scale, not an atom-relaxed neck width.'},
        'atomic_connectivity': {'Ti3SiC2': tsc_topology, 'Ag': ag_topology},
        'skeleton_cleanup': {'disconnected_fragment_and_stoichiometry_steps': skeleton_cleanup,
                             'last_component_report': skeleton_info},
        'Ag_selection': ag_selection,
        'minimum_requested_cross_clearance_A': MIN_CROSS_CLEARANCE_A,
        'files': {'lammps_data': 'skeleton_microtwig_200x200x500.data',
                  'extxyz': 'skeleton_microtwig_200x200x500.extxyz',
                  'qsgs_mask': 'qsgs_skeleton_microtwig.npz'},
        'limitations': ['This is a constructed geometry, not energy-minimized or thermally equilibrated.',
                        'No MD, heating, melting, cooling, or potential validation was run.',
                        'Finer morphology is exploratory; QSGS connectivity is checked geometrically only.',
                        'The fixed/thermostat/Newton layer labels describe z regions; the mechanical group policy keeps all Ag mobile.']}
    (ROOT / 'geometry_report.json').write_text(json.dumps(report, indent=2) + '\n')
    p(f'Finished: atoms={len(atoms):,}; Ag:skeleton mass={ratio:.8f}; counts={counts}')
    p(f'Skeleton-to-nearest-Ag distance percentiles: {fine_metric}')


if __name__ == '__main__':
    main()
