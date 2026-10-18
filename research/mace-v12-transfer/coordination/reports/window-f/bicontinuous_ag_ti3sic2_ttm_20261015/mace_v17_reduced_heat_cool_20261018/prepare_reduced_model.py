#!/usr/bin/env python3
"""Trim only remote ideal-crystal slab thickness from the 1,926-atom prototype."""
from pathlib import Path
import hashlib
import json
import numpy as np
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write
from ase.neighborlist import neighbor_list
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parent
SRC = ROOT.parent / "mace_v17_small_heating_cooling_20261018" / "SMALL_MODEL_INITIAL.extxyz"
SRC_AUDIT = ROOT.parent / "mace_v17_small_heating_cooling_20261018" / "MODEL_GEOMETRY_AUDIT.json"
MACE = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_03/models/MACE_periodic_v17_provisional_screening.model")
ORIGINAL = Path("research/mace-v12-transfer/coordination/reports/window-f/bicontinuous_ag_ti3sic2_ttm_20261015/exploratory_melt_cool_20261018/input/NEW_MODEL_ATOMS.data")
EXPECTED_MACE = "559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3"
EXPECTED_ORIGINAL = "e2b62a68c46e23037370a6ca20951aa0774edea2d728d3c3aadbcba4b3755674"
KEEP_TSC_Z_MIN = 8.812
KEEP_AG_LAYERS = {0, 1, 2, 3}
VACUUM_TOP = 18.0


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def components(atoms, mask, cutoff):
    i, j = neighbor_list("ij", atoms, cutoff)
    take = mask[i] & mask[j] & (i < j)
    ids = np.flatnonzero(mask)
    idx = np.full(len(atoms), -1, dtype=int)
    idx[ids] = np.arange(len(ids))
    if np.any(take):
        u, v = idx[i[take]], idx[j[take]]
        mat = coo_matrix((np.ones(2 * len(u)), (np.r_[u, v], np.r_[v, u])), shape=(len(ids), len(ids))).tocsr()
    else:
        mat = coo_matrix((len(ids), len(ids))).tocsr()
    n, labels = connected_components(mat, directed=False)
    sizes = np.bincount(labels) if len(labels) else np.array([], dtype=int)
    return {"components": int(n), "largest_component_atoms": int(sizes.max()) if len(sizes) else 0,
            "largest_component_fraction": float(sizes.max()/len(ids)) if len(ids) else 0.0}


def main():
    if sha(MACE) != EXPECTED_MACE or sha(ORIGINAL) != EXPECTED_ORIGINAL:
        raise SystemExit("V17 weight or original 39,486-atom model SHA mismatch")
    source = read(SRC)
    old_phase = source.arrays["phase"].astype(str)
    old_layer = np.asarray(source.arrays["Ag111_layer"], int)
    old_z = source.positions[:, 2]
    keep = ((old_phase == "Ti3SiC2") & (old_z >= KEEP_TSC_Z_MIN - 0.002)) | \
           ((old_phase == "Ag") & np.isin(old_layer, list(KEEP_AG_LAYERS)))
    old_ids = np.asarray(source.arrays["id"], int)
    atoms = source[keep]
    source_ids = old_ids[keep]
    phase = atoms.arrays["phase"].astype(str)
    symbols = np.asarray(atoms.get_chemical_symbols())
    old_layer = np.asarray(atoms.arrays["Ag111_layer"], int)
    atoms.positions[:, 2] -= float(atoms.positions[:, 2].min())
    z_top = float(atoms.positions[:, 2].max())
    cell = np.asarray(atoms.cell.array).copy()
    cell[2] = [0.0, 0.0, z_top + VACUUM_TOP]
    atoms.set_cell(cell, scale_atoms=False)
    atoms.set_pbc([True, True, False])
    tsc = phase == "Ti3SiC2"
    ag = phase == "Ag"
    # Preserve the full bottom Ti atomic plane as the only mechanical support.
    support = tsc & (symbols == "Ti") & np.isclose(atoms.positions[:, 2], 0.0, atol=1e-5)
    if int(support.sum()) != 100:
        raise RuntimeError(f"Expected one 100-atom Ti support plane, got {support.sum()}")
    counts = {e: int(np.sum(symbols == e)) for e in ("Ti", "Si", "C", "Ag")}
    if counts["Ti"] != 3 * counts["Si"] or counts["C"] != 2 * counts["Si"]:
        raise RuntimeError(f"Trimmed Ti3SiC2 repeat is not stoichiometric: {counts}")
    layer_ids = sorted(set(old_layer[ag].tolist()))
    if layer_ids != [0, 1, 2, 3]:
        raise RuntimeError(f"Ag layers are not the interface-side four consecutive planes: {layer_ids}")
    atoms.arrays["source_prototype_atom_id"] = source_ids.astype(int)
    atoms.arrays["id"] = np.arange(1, len(atoms) + 1, dtype=int)
    atoms.arrays["support_fixed"] = support.astype(np.int32)
    frac = atoms.positions @ np.linalg.inv(np.asarray(atoms.cell.array))
    duv = frac[:, :2] - np.array([0.5, 0.5])
    radial = np.full(len(atoms), np.inf)
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            xy = (duv + np.array([a, b])) @ np.asarray(atoms.cell.array)[:2, :2]
            radial = np.minimum(radial, np.linalg.norm(xy, axis=1))
    heat_region = ag & (old_layer >= 2) & (radial <= 10.0)
    cold_sink = tsc & ~support & (atoms.positions[:, 2] <= 4.5)
    if int(heat_region.sum()) < 50 or int(cold_sink.sum()) < 100:
        raise RuntimeError(f"Unexpected local reservoirs: heat={heat_region.sum()}, cold_sink={cold_sink.sum()}")
    atoms.arrays["heat_region_initial"] = heat_region.astype(np.int32)
    atoms.arrays["cold_sink"] = cold_sink.astype(np.int32)
    atoms.arrays["phase"] = phase.astype("U8")
    atoms.arrays["region"] = np.where(ag, "Ag_film", "TSC_substrate")
    atoms.info["model_name"] = "reduced_four_layer_Ag111_six_plane_stoichiometric_TSC_slab"
    atoms.info["source_model"] = str(SRC.relative_to(Path.cwd()))
    write(ROOT / "REDUCED_MODEL.extxyz", atoms, format="extxyz")
    write(ROOT / "REDUCED_MODEL.data", atoms, format="lammps-data", atom_style="atomic", specorder=["Ti", "Si", "C", "Ag"])

    masses = {e: float(atomic_masses[atomic_numbers[e]]) for e in counts}
    m_ag = counts["Ag"] * masses["Ag"]
    m_tsc = sum(counts[e] * masses[e] for e in ("Ti", "Si", "C"))
    i, j, d = neighbor_list("ijd", atoms, 4.5)
    cross = (phase[i] != phase[j]) & (i < j)
    old_audit = json.loads(SRC_AUDIT.read_text())
    tsc_z = np.unique(np.round(atoms.positions[tsc, 2], 3)).tolist()
    ag_z = np.unique(np.round(atoms.positions[ag, 2], 3)).tolist()
    report = {
        "status": "geometry_built_pending_static_force_and_MD",
        "atom_count": int(len(atoms)), "element_counts": counts,
        "TSC_atom_count": int(tsc.sum()), "Ag_atom_count": int(ag.sum()),
        "Ag_mass_fraction": float(m_ag / (m_ag + m_tsc)), "TSC_mass_fraction": float(m_tsc / (m_ag + m_tsc)),
        "TSC_stoichiometry": "Ti3SiC2 exactly; selected six ideal basal atomic planes form a complete 100-formula-unit repeat",
        "cell_A": atoms.cell.array.tolist(), "pbc": [True, True, False],
        "material_z_A": [float(atoms.positions[:, 2].min()), z_top], "top_vacuum_A": VACUUM_TOP,
        "TSC_kept_z_planes_A": tsc_z, "Ag_kept_layer_ids": layer_ids, "Ag_kept_z_planes_A": ag_z,
        "support_atom_count": int(support.sum()), "support_definition": "only bottom Ti plane fixed; all other Ti3SiC2 atoms and all Ag atoms mobile",
        "heat_source_atom_count": int(heat_region.sum()), "heat_source_definition": "fixed atom labels in upper two Ag(111) planes within 10 A lateral radius of XY cell center; Langevin 1900 K during heating only",
        "cold_sink_atom_count": int(cold_sink.sum()), "cold_sink_definition": "mobile TSC atoms with 0 < z <= 4.5 A; Langevin 300 K in heating and cooling",
        "TSC_network_cutoff_3p25_A": components(atoms, tsc, 3.25),
        "Ag_network_cutoff_3p10_A": components(atoms, ag, 3.10),
        "interface_cross_pair_min_A_within_4p5_A": float(np.min(d[cross])) if np.any(cross) else None,
        "interface_cross_pair_percentiles_A_within_4p5_A": np.percentile(d[cross], [5, 50, 95]).tolist() if np.any(cross) else [],
        "Ag111_in_plane_strain_percent": old_audit["Ag_inplane_strain_percent"],
        "Ag111_interplane_spacing_A": old_audit["Ag111_interplane_spacing_A"],
        "XY_cell_preserved_from_1926_prototype": True,
        "source_prototype_extxyz_sha256": sha(SRC), "MACE_V17_sha256": sha(MACE), "original_39486_data_sha256": sha(ORIGINAL),
        "source_atom_id_preserved_as_array": "source_prototype_atom_id",
        "reduction_rationale": "Keep all original lateral repeats and the Ag/TSC interface. Retain one complete stoichiometric Ti3SiC2 formula repeat (six basal planes; 7.54 A) and four consecutive FCC Ag(111) planes (9.44 A), preserving the Ag interface plane, subsurface Ag, and a free-surface region. Remove remote lower TSC repeat and two remote upper Ag planes only.",
        "difference_from_39486_bicontinuous_target": "This is a laterally periodic layered Ag(111)/Ti3SiC2(0001) interface slab, not the 39,486-atom 3D bicontinuous geometry; it represents a local interface patch and omits the original macroscopic surface/connected Ag and TSC networks.",
        "difference_from_1926_prototype": "Same XY coincidence cell, crystal orientation, ideal coordinates at the Ag/TSC interface, vacuum, and interface gap; reduced to six top TSC planes and four interface-side Ag planes. The thin slab and changed phase volume/mass fraction are finite-size limitations.",
    }
    (ROOT / "REDUCED_GEOMETRY_AUDIT.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
