#!/usr/bin/env python3
"""Build a bicontinuous W-backbone/Ag-matrix QSGS geometry candidate."""

from __future__ import annotations

import argparse
import gzip
import itertools
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import write
from scipy import sparse
from scipy.sparse.csgraph import breadth_first_order, connected_components
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "models" / "ag_w_qsgs_200A"
MATERIAL_SIZE = np.array([200.0, 200.0, 200.0])
BOX_SIZE = np.array([200.0, 200.0, 300.0])
DEFAULT_GRID = (80, 80, 80)
SEED = 20261123
P_CORE = 0.0003
P_GROW = 0.145
MAX_GENERATIONS = 10
TARGET_W_FRACTION = 0.3545
LATTICE = {"W": ("bcc", 3.16), "Ag": ("fcc", 4.09)}
NEIGHBOR_CUTOFF = {"W": 2.95, "Ag": 3.15}


def shifted(mask: np.ndarray, offset: tuple[int, int, int]) -> np.ndarray:
    result = np.roll(mask, offset, axis=(0, 1, 2))
    if offset[2] > 0:
        result[:, :, : offset[2]] = False
    elif offset[2] < 0:
        result[:, :, offset[2] :] = False
    return result


def qsgs_mask(shape, seed, p_core, p_grow, max_generations, target_fraction):
    """Generate W through four controls: seed, growth, generations, target fraction."""
    rng = np.random.default_rng(seed)
    occupied = rng.random(shape) < p_core
    target = int(round(target_fraction * np.prod(shape)))
    offsets = [v for v in itertools.product((-1, 0, 1), repeat=3) if any(v)]
    history = [{"generation": 0, "occupied": int(occupied.sum()), "fraction": float(occupied.mean())}]
    for generation in range(1, max_generations + 1):
        added = np.zeros(shape, dtype=bool)
        for offset in offsets:
            added |= shifted(occupied, offset) & (rng.random(shape) < p_grow / np.linalg.norm(offset))
        added &= ~occupied
        count_needed = target - int(occupied.sum())
        ids = np.flatnonzero(added)
        if count_needed <= 0 or len(ids) >= count_needed:
            rng.shuffle(ids)
            occupied.flat[ids[:count_needed]] = True
            history.append({
                "generation": generation,
                "new_voxels_available": int(len(ids)),
                "new_voxels_added": int(count_needed),
                "occupied": int(occupied.sum()),
                "fraction": float(occupied.mean()),
                "target_reached": True,
            })
            return occupied, history
        occupied |= added
        history.append({
            "generation": generation,
            "new_voxels_available": int(len(ids)),
            "new_voxels_added": int(len(ids)),
            "occupied": int(occupied.sum()),
            "fraction": float(occupied.mean()),
            "target_reached": False,
        })
    raise ValueError("QSGS target fraction not reached within max_generations")


def make_lattice_sites(element: str):
    structure, a = LATTICE[element]
    repeat_xy = np.rint(MATERIAL_SIZE[:2] / a).astype(int)
    repeat = (*repeat_xy, int(np.ceil(MATERIAL_SIZE[2] / a)))
    crystal = bulk(element, structure, a=a, cubic=True).repeat(repeat)
    positions = crystal.positions[crystal.positions[:, 2] < MATERIAL_SIZE[2]].copy()
    strain = MATERIAL_SIZE[:2] / (repeat_xy * a) - 1.0
    positions[:, 0] = (positions[:, 0] * (1.0 + strain[0])) % BOX_SIZE[0]
    positions[:, 1] = (positions[:, 1] * (1.0 + strain[1])) % BOX_SIZE[1]
    positions[positions[:, 0] >= BOX_SIZE[0] - 1.0e-10, 0] = 0.0
    positions[positions[:, 1] >= BOX_SIZE[1] - 1.0e-10, 1] = 0.0
    return positions, strain.tolist()


def in_mask(mask, positions):
    index = np.floor(positions / MATERIAL_SIZE * np.asarray(mask.shape)).astype(int)
    index = np.clip(index, 0, np.asarray(mask.shape) - 1)
    return mask[tuple(index.T)]


def topology(positions, element):
    pairs = cKDTree(positions, boxsize=[BOX_SIZE[0], BOX_SIZE[1], 1000.0]).query_pairs(
        NEIGHBOR_CUTOFF[element], output_type="ndarray"
    )
    graph = sparse.coo_matrix(
        (np.ones(len(pairs), dtype=np.uint8), (pairs[:, 0], pairs[:, 1])),
        shape=(len(positions), len(positions)),
    ).tocsr()
    graph = graph + graph.T
    count, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels)
    keep = labels == sizes.argmax()
    pos = positions[keep]
    subgraph = graph[keep][:, keep]
    order, parents = breadth_first_order(subgraph, 0, directed=False)
    children = order[1:]
    delta = pos[children, :2] - pos[parents[children], :2]
    shifts = np.zeros((len(pos), 2), dtype=np.int32)
    shifts[children] = -np.rint(delta / BOX_SIZE[:2]).astype(np.int32)
    image = np.zeros_like(shifts)
    for child in children:
        image[child] = image[parents[child]] + shifts[child]
    left, right = subgraph.nonzero()
    edge = left < right
    left, right = left[edge], right[edge]
    jumps = -np.rint((pos[right, :2] - pos[left, :2]) / BOX_SIZE[:2]).astype(np.int32)
    winding = image[right] - image[left] - jumps
    cycles = np.unique(winding[np.any(winding, axis=1)], axis=0)
    rank = int(np.linalg.matrix_rank(cycles.astype(float))) if cycles.size else 0
    zmin, zmax = float(pos[:, 2].min()), float(pos[:, 2].max())
    return {
        "atoms": int(len(positions)),
        "components": int(count),
        "largest_component_atoms": int(keep.sum()),
        "largest_component_share": float(keep.mean()),
        "neighbor_cutoff_A": NEIGHBOR_CUTOFF[element],
        "xy_winding_rank": rank,
        "largest_component_z_extent_A": [zmin, zmax],
        "largest_component_spans_z": bool(zmin < 3.2 and zmax > 196.8),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--p-core", type=float, default=P_CORE)
    parser.add_argument("--p-grow", type=float, default=P_GROW)
    parser.add_argument("--generations", type=int, default=MAX_GENERATIONS)
    parser.add_argument("--target-w-fraction", type=float, default=TARGET_W_FRACTION)
    parser.add_argument("--grid", type=int, default=DEFAULT_GRID[0], help="cubic QSGS voxel count per side; default 80 preserves the original build")
    args = parser.parse_args()
    if args.grid < 8:
        parser.error("--grid must be at least 8 voxels per side")
    grid = (args.grid, args.grid, args.grid)
    args.output.mkdir(parents=True, exist_ok=True)
    for name in ("ag_w_qsgs.data.gz", "ag_w_qsgs.extxyz.gz", "qsgs_morphology.npz", "geometry_report.json"):
        if (args.output / name).exists():
            raise FileExistsError(f"Refusing to overwrite {args.output / name}")

    w_mask, history = qsgs_mask(
        grid, args.seed, args.p_core, args.p_grow, args.generations, args.target_w_fraction
    )
    w_pool, w_strain = make_lattice_sites("W")
    ag_pool, ag_strain = make_lattice_sites("Ag")
    w_positions = w_pool[in_mask(w_mask, w_pool)]
    ag_positions = ag_pool[~in_mask(w_mask, ag_pool)]
    symbols = ["W"] * len(w_positions) + ["Ag"] * len(ag_positions)
    positions = np.vstack((w_positions, ag_positions))
    atoms = Atoms(symbols=symbols, positions=positions, cell=BOX_SIZE, pbc=(True, True, False))

    z = atoms.positions[:, 2]
    layers = np.select(
        [z < 20.0, (z >= 20.0) & (z < 40.0), (z >= 40.0) & (z < 200.0)],
        [1, 2, 3], default=0,
    ).astype(np.int8)
    atoms.new_array("phase", np.where(np.asarray(symbols) == "W", 1, 2).astype(np.int8))
    atoms.new_array("layer", layers)
    atoms.new_array("fixed_layer", (layers == 1).astype(np.int8))
    atoms.new_array("thermostat_layer", (layers == 2).astype(np.int8))
    atoms.info.update({
        "status": "UNRELAXED_QSGS_GEOMETRY_ONLY",
        "seed": int(args.seed),
        "mass_ratio_target": "Ag:W approximately 1:1",
        "heating_method": "TTM requested; not run by the geometry builder",
    })

    n_w, n_ag = len(w_positions), len(ag_positions)
    masses = atoms.get_masses()
    ratio = float(masses[n_w:].sum() / masses[:n_w].sum())
    minimum_cross = float(cKDTree(w_positions, boxsize=[200.0, 200.0, 1000.0]).query(ag_positions, workers=-1)[0].min())
    phase_topology = {
        "W": topology(w_positions, "W"),
        "Ag": topology(ag_positions, "Ag"),
    }
    for element, result in phase_topology.items():
        if result["xy_winding_rank"] != 2 or not result["largest_component_spans_z"] or result["largest_component_share"] < 0.98:
            raise ValueError(f"Atomic phase {element} failed bicontinuity: {result}")

    layer_counts = {
        "fixed_0_20_A": int(np.count_nonzero(layers == 1)),
        "thermostat_20_40_A": int(np.count_nonzero(layers == 2)),
        "newton_40_200_A": int(np.count_nonzero(layers == 3)),
        "sputter_monitor_200_250_A": int(np.count_nonzero((z >= 200) & (z < 250))),
        "evaporation_monitor_250_300_A": int(np.count_nonzero((z >= 250) & (z < 300))),
    }
    with gzip.open(args.output / "ag_w_qsgs.data.gz", "wt") as stream:
        write(stream, atoms, format="lammps-data", atom_style="atomic", specorder=["W", "Ag"], masses=True)
    with gzip.open(args.output / "ag_w_qsgs.extxyz.gz", "wt") as stream:
        write(stream, atoms, format="extxyz")
    np.savez_compressed(
        args.output / "qsgs_morphology.npz",
        W_backbone=w_mask,
        Ag_matrix=~w_mask,
        material_size_A=MATERIAL_SIZE,
        grid=np.asarray(grid),
        seed=np.asarray(args.seed),
        p_core=np.asarray(args.p_core),
        p_grow=np.asarray(args.p_grow),
        generations=np.asarray(args.generations),
        target_w_fraction=np.asarray(args.target_w_fraction),
    )
    report = {
        "status": "PASS_UNRELAXED_BICONTINUOUS_GEOMETRY_ONLY",
        "material_size_A": MATERIAL_SIZE.tolist(),
        "box_size_A": BOX_SIZE.tolist(),
        "boundary": "p p f",
        "material_z_A": [0.0, 200.0],
        "vacuum_z_A": [200.0, 300.0],
        "layer_atom_counts": layer_counts,
        "total_atoms": int(len(atoms)),
        "atom_counts": {"W": int(n_w), "Ag": int(n_ag)},
        "atom_fractions": {"W": n_w / len(atoms), "Ag": n_ag / len(atoms)},
        "mass_ratio_Ag_over_W": ratio,
        "mass_ratio_relative_error": abs(ratio - 1.0),
        "minimum_initial_Ag_W_distance_A": minimum_cross,
        "initial_geometry_note": "Distinct perfect lattices have close cross-interface contacts; overlap relaxation is required before MD.",
        "qsgs": {
            "four_parameters": {
                "core_seed_probability": args.p_core,
                "growth_probability_base": args.p_grow,
                "maximum_generations": args.generations,
                "target_W_voxel_fraction": args.target_w_fraction,
            },
            "growth_rule": "26 directions; attachment probability p_grow/r; periodic x/y, open z",
            "seed": args.seed,
            "grid": list(grid),
            "voxel_size_A": (MATERIAL_SIZE / np.asarray(grid)).tolist(),
            "W_voxel_fraction_actual": float(w_mask.mean()),
            "Ag_voxel_fraction_actual": float((~w_mask).mean()),
            "generation_history": history,
        },
        "lattices": {
            "W": {"structure": "bcc", "lattice_parameter_A": LATTICE["W"][1], "xy_strain": w_strain},
            "Ag": {"structure": "fcc", "lattice_parameter_A": LATTICE["Ag"][1], "xy_strain": ag_strain},
        },
        "atomic_phase_topology": phase_topology,
        "TTM": {"requested": True, "status": "NOT_RUN", "note": "Stock fix ttm/mod does not support the free-z boundary; the workspace has an Ag-only finite-volume reference solver, but W electronic data and an Ag-W mixed-phase closure are not yet verified."},
    }
    (args.output / "geometry_report.json").write_text(json.dumps(report, indent=2) + "\n")
    (args.output / "in.geometry_read_check").write_text(
        "units metal\natom_style atomic\nboundary p p f\nread_data ag_w_qsgs.data.gz\n"
        "group tungsten type 1\ngroup silver type 2\n"
        "region fixed_zone block INF INF INF INF 0 20 units box\ngroup fixed region fixed_zone\n"
        "region thermo_zone block INF INF INF INF 20 40 units box\ngroup thermostat region thermo_zone\n"
        "region newton_zone block INF INF INF INF 40 200 units box\ngroup newton region newton_zone\n"
        "variable nW equal count(tungsten)\nvariable nAg equal count(silver)\n"
        "variable nFix equal count(fixed)\nvariable nThermo equal count(thermostat)\nvariable nNewton equal count(newton)\n"
        "print \"GEOMETRY_READ_OK W=${nW} Ag=${nAg} fixed=${nFix} thermostat=${nThermo} newton=${nNewton}\"\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
