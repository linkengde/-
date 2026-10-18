#!/usr/bin/env python3
"""Build an independent ideal-crystal Ag(111)/Ti3SiC2(0001) slab model."""
from pathlib import Path
import hashlib
import json
import numpy as np
from ase.build import make_supercell
from ase.io import read, write
from ase import Atoms
from ase.data import atomic_masses, atomic_numbers
from ase.neighborlist import neighbor_list
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parent
REPO = Path.cwd()
CIF = REPO / "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif"
MACE_MODEL = REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_03/models/MACE_periodic_v17_provisional_screening.model"
ORIGINAL = REPO / "research/mace-v12-transfer/coordination/reports/window-f/bicontinuous_ag_ti3sic2_ttm_20261015/exploratory_melt_cool_20261018/input/NEW_MODEL_ATOMS.data"

NTSC = 10
NAG = 11
N_AG_LAYERS = 6
AG_A0 = 4.09
INTERFACE_GAP = 2.70
VACUUM_TOP = 18.0


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def graph_report(atoms, phase, cutoff, phase_name):
    i, j = neighbor_list("ij", atoms, cutoff)
    keep = (phase[i] == phase_name) & (phase[j] == phase_name) & (i < j)
    pairs = np.column_stack((i[keep], j[keep]))
    nodes = np.flatnonzero(phase == phase_name)
    lookup = np.full(len(atoms), -1, dtype=int)
    lookup[nodes] = np.arange(len(nodes))
    if len(pairs):
        rr = lookup[pairs[:, 0]]
        cc = lookup[pairs[:, 1]]
        data = np.ones(2 * len(pairs), dtype=np.uint8)
        mat = coo_matrix((data, (np.r_[rr, cc], np.r_[cc, rr])), shape=(len(nodes), len(nodes))).tocsr()
    else:
        mat = coo_matrix((len(nodes), len(nodes))).tocsr()
    ncomp, labels = connected_components(mat, directed=False)
    sizes = np.bincount(labels) if len(labels) else np.array([], dtype=int)
    return {"phase": phase_name, "cutoff_A": cutoff, "components": int(ncomp),
            "largest_component_atoms": int(sizes.max()) if len(sizes) else 0,
            "largest_component_fraction": float(sizes.max() / len(nodes)) if len(nodes) else 0.0}


def build():
    parent = read(CIF)
    if len(parent) != 12 or sorted(parent.get_chemical_symbols()) != sorted(["Ti"] * 6 + ["Si"] * 2 + ["C"] * 4):
        raise RuntimeError("Unexpected COD 9009647 parent basis")
    tsc = make_supercell(parent, np.diag([NTSC, NTSC, 1]))
    tsc.pbc = [True, True, False]
    cell = np.array(tsc.cell.array, dtype=float)
    A, B = cell[0].copy(), cell[1].copy()
    z_top = float(tsc.positions[:, 2].max())
    ag_nn_ideal = AG_A0 / np.sqrt(2.0)
    ag_d111 = AG_A0 / np.sqrt(3.0)
    matched_primitive_a = np.linalg.norm(A) / NAG
    matched_primitive_b = np.linalg.norm(B) / NAG
    strain_a = matched_primitive_a / ag_nn_ideal - 1.0
    strain_b = matched_primitive_b / ag_nn_ideal - 1.0

    # Explicit FCC(111) triangular planes, aligned to the hexagonal basal plane.
    # Layer shifts follow ABC stacking. The 10:11 coincidence cell avoids a
    # nonperiodic lateral seam; all Ag sites are generated analytically.
    stacking = [(0.0, 0.0), (2.0 / 3.0, 1.0 / 3.0), (1.0 / 3.0, 2.0 / 3.0)]
    ag_positions = []
    ag_layer = []
    z0 = z_top + INTERFACE_GAP
    for k in range(N_AG_LAYERS):
        su, sv = stacking[k % 3]
        for i in range(NAG):
            for j in range(NAG):
                ag_positions.append(((i + su) / NAG) * A + ((j + sv) / NAG) * B + np.array([0.0, 0.0, z0 + k * ag_d111]))
                ag_layer.append(k)
    ag = Atoms("Ag" * len(ag_positions), positions=np.asarray(ag_positions), cell=cell, pbc=[True, True, False])
    combined = tsc + ag
    full_z = z0 + (N_AG_LAYERS - 1) * ag_d111 + VACUUM_TOP
    cell[2] = [0.0, 0.0, full_z]
    combined.set_cell(cell, scale_atoms=False)
    combined.pbc = [True, True, False]

    n_tsc = len(tsc)
    n_ag = len(ag)
    phase = np.array(["Ti3SiC2"] * n_tsc + ["Ag"] * n_ag, dtype="U8")
    symbols = np.array(combined.get_chemical_symbols())
    support = np.zeros(len(combined), dtype=bool)
    support[:n_tsc] = np.isclose(combined.positions[:n_tsc, 2], 0.0, atol=1e-6)
    ids = np.arange(1, len(combined) + 1, dtype=int)
    combined.arrays["id"] = ids
    combined.arrays["phase"] = phase
    combined.arrays["region"] = np.where(phase == "Ag", "Ag_film", "TSC_substrate")
    combined.arrays["support_fixed"] = support.astype(np.int32)
    combined.arrays["Ag111_layer"] = np.r_[np.full(n_tsc, -1, dtype=int), np.asarray(ag_layer, dtype=int)]
    combined.info["model_name"] = "independent_ideal_Ag111_Ti3SiC2_0001_slab"
    combined.info["source_CIF"] = str(CIF.relative_to(REPO))
    combined.info["Ag_lattice_parameter_A"] = AG_A0
    combined.info["Ag_in_plane_strain_x"] = float(strain_a)
    combined.info["Ag_in_plane_strain_y"] = float(strain_b)

    write(ROOT / "SMALL_MODEL_INITIAL.extxyz", combined, format="extxyz")
    # LAMMPS atomic-style data file for independent review / OVITO import.
    write(ROOT / "SMALL_MODEL_INITIAL.data", combined, format="lammps-data", atom_style="atomic", specorder=["Ti", "Si", "C", "Ag"])

    masses = {s: float(atomic_masses[atomic_numbers[s]]) for s in set(symbols)}
    mass_ag = n_ag * masses["Ag"]
    mass_tsc = sum(int(np.sum(symbols == s)) * masses[s] for s in ("Ti", "Si", "C"))
    topC = max(combined.positions[(symbols == "C") & (phase == "Ti3SiC2"), 2])
    # Cross-phase distances in the interfacial 5 A band.
    ii, jj, dd = neighbor_list("ijd", combined, 4.5)
    cross = (phase[ii] != phase[jj]) & (ii < jj)
    near_interface = cross & (combined.positions[ii, 2] >= topC - 5.0) & (combined.positions[jj, 2] >= topC - 5.0)
    cross_d = dd[near_interface]
    geom = {
        "atom_count": int(len(combined)), "TSC_atom_count": int(n_tsc), "Ag_atom_count": int(n_ag),
        "element_counts": {s: int(np.sum(symbols == s)) for s in ["Ti", "Si", "C", "Ag"]},
        "Ag_mass_fraction": float(mass_ag / (mass_ag + mass_tsc)),
        "TSC_mass_fraction": float(mass_tsc / (mass_ag + mass_tsc)),
        "cell_A": combined.cell.array.tolist(), "pbc": [True, True, False],
        "cell_volume_A3": float(combined.get_volume()), "material_z_min_A": float(combined.positions[:, 2].min()),
        "material_z_max_A": float(combined.positions[:, 2].max()), "top_vacuum_A": float(full_z - combined.positions[:, 2].max()),
        "interface_gap_to_topmost_TSC_C_A": INTERFACE_GAP,
        "interface_cross_pair_min_A": float(cross_d.min()), "interface_cross_pair_percentiles_A": np.percentile(cross_d, [5, 50, 95]).tolist(),
        "TSC_top_termination": "C-terminated Ti3SiC2(0001), one COD conventional c-period slab",
        "Ag_orientation": "FCC Ag(111), ABC stacking, 6 atomic planes",
        "TSC_lateral_repeats": [NTSC, NTSC], "Ag_lateral_repeats": [NAG, NAG],
        "ideal_Ag111_inplane_nn_A": float(ag_nn_ideal), "matched_Ag111_inplane_nn_A": [float(matched_primitive_a), float(matched_primitive_b)],
        "Ag_inplane_strain_percent": [float(100 * strain_a), float(100 * strain_b)],
        "Ag111_interplane_spacing_A": float(ag_d111), "support_fixed_atom_count": int(support.sum()),
        "support_definition": "lowest Ti plane at z=0 fixed; only 100 atoms, all other TSC and Ag atoms mobile",
        "hashes": {"CIF": sha256(CIF), "V17_weight": sha256(MACE_MODEL), "original_39486_data": sha256(ORIGINAL)},
        "original_model_audit": {"atom_count": 39486, "boundary": "f p p (X nonperiodic, Y/Z periodic)", "type_to_element": {"1": "Ti", "2": "Si", "3": "C", "4": "Ag"}},
    }
    geom["connectivity"] = [graph_report(combined, phase, 3.25, "Ti3SiC2"), graph_report(combined, phase, 3.10, "Ag")]
    (ROOT / "MODEL_GEOMETRY_AUDIT.json").write_text(json.dumps(geom, indent=2) + "\n")
    return combined, geom


if __name__ == "__main__":
    build()
