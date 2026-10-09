"""Create and audit label-free rigid Ag-layer separation proposals; never runs DFT."""
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.data import covalent_radii
from ase.io import read, write

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[5]
PARENT = REPO / "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures"
SOURCE = json.loads((PARENT / "manifest.json").read_text())
hash_file = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
for line in (PARENT / "SHA256SUMS.txt").read_text().splitlines():
    digest, name = line.split("  ", 1)
    assert hash_file(PARENT / name) == digest, f"Parent source inventory mismatch: {name}"
normal_norm_tol = 1.0e-10
overlap_scale = 0.80
records = []
source_atoms = {}

for source in SOURCE["records"]:
    path = REPO / source["path"]
    assert path.is_file() and hash_file(path) == source["sha256"]
    atoms = read(path)
    assert len(atoms) == source["atoms"] == 52
    assert atoms.get_chemical_formula() == source["formula"]
    assert atoms.pbc.tolist() == [True, True, True]
    assert atoms.info.get("dataset_role") == "CIF_DERIVED_DEVELOPMENT_PROPOSAL"
    assert np.array_equal(atoms.arrays["lammps_id"],
                          np.arange(source["new_id_range"][0], source["new_id_range"][1] + 1))
    assert int(np.count_nonzero(atoms.arrays["central_pair"])) == 2
    assert set(atoms.get_chemical_symbols()[i] for i in np.flatnonzero(atoms.arrays["central_pair"]))
    # Derive the physical surface normal from the in-plane lattice vectors and orient it outward.
    normal = np.cross(atoms.cell[0], atoms.cell[1])
    normal /= np.linalg.norm(normal)
    if np.dot(normal, atoms.cell[2]) < 0:
        normal *= -1
    assert abs(np.dot(normal, atoms.cell[0])) < normal_norm_tol
    assert abs(np.dot(normal, atoms.cell[1])) < normal_norm_tol
    assert np.dot(normal, atoms.cell[2]) > 0
    symbols = np.asarray(atoms.get_chemical_symbols())
    ag = np.flatnonzero(symbols == "Ag")
    framework = np.flatnonzero(symbols != "Ag")
    assert len(ag) == 4 and len(framework) == 48
    proj = atoms.positions @ normal
    ag_plane = float(np.mean(proj[ag]))
    assert float(np.ptp(proj[ag])) < 1.0e-8
    top_plane = float(np.max(proj[framework]))
    gap = ag_plane - top_plane
    area = float(np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1])))
    cell_height = float(np.dot(atoms.cell[2], normal))
    parent_vacuum = cell_height - float(np.max(proj) - np.min(proj))
    assert abs(gap - source["gap_plane_A"]) < 1.0e-6
    source_atoms[source["interface"]] = atoms
    central_ids = source["marked_pair_ids"]
    id_to_index = {int(atom_id): i for i, atom_id in enumerate(atoms.arrays["lammps_id"])}
    pair_indices = [id_to_index[int(atom_id)] for atom_id in central_ids]
    assert all(atoms.arrays["central_pair"][pair_indices])
    pair_species = [symbols[i] for i in pair_indices]
    assert pair_species.count("Ag") == 1
    parent_distances = atoms.get_all_distances(mic=True)
    parent_marked_distance = float(parent_distances[pair_indices[0], pair_indices[1]])
    parent_species_min = {}
    for i in ag:
        for j in framework:
            key = f"Ag-{symbols[j]}"
            parent_species_min[key] = min(parent_species_min.get(key, float("inf")),
                                          float(parent_distances[i, j]))

    for delta in (-0.20, 0.20):
        variant = atoms.copy()
        translation = normal * delta
        variant.positions[ag] += translation
        # All non-Ag positions, cell, PBC, ordering and every atom array must stay byte-value identical.
        assert np.array_equal(variant.positions[framework], atoms.positions[framework])
        assert np.array_equal(variant.cell.array, atoms.cell.array)
        assert np.array_equal(variant.pbc, atoms.pbc)
        for key, value in atoms.arrays.items():
            if key == "positions":
                continue
            assert np.array_equal(variant.arrays[key], value), key
        moved = variant.positions - atoms.positions
        expected = np.zeros_like(moved)
        expected[ag] = translation
        assert np.allclose(moved, expected, atol=1.0e-10, rtol=0)
        distances = variant.get_all_distances(mic=True)
        minima = {}
        thresholds = {}
        unsafe = []
        for i in ag:
            for j in framework:
                key = f"Ag-{symbols[j]}"
                d = float(distances[i, j])
                minima[key] = min(minima.get(key, float("inf")), d)
                cutoff = overlap_scale * (covalent_radii[47] + covalent_radii[atoms[j].number])
                thresholds[key] = float(cutoff)
                if d < cutoff:
                    unsafe.append({"species_pair": key, "Ag_id": int(atoms.arrays["lammps_id"][i]),
                                   "framework_id": int(atoms.arrays["lammps_id"][j]),
                                   "distance_A": d, "minimum_allowed_A": float(cutoff)})
        assert not unsafe, f"Unsafe Ag/framework overlap for {source['interface']} delta={delta}: {unsafe}"
        marked_distance = float(distances[pair_indices[0], pair_indices[1]])
        termination = source["top_surface_species"]
        actual_top_plane = float(np.max(variant.positions[framework] @ normal))
        actual_ag_plane = float(np.min(variant.positions[ag] @ normal))
        actual_gap = actual_ag_plane - actual_top_plane
        variant_proj = variant.positions @ normal
        actual_vacuum = cell_height - float(np.max(variant_proj) - np.min(variant_proj))
        label = f"{source['interface']}_COD9009647_sep_{'p' if delta > 0 else 'm'}020A_PROPOSAL"
        variant.info = dict(atoms.info)
        variant.info.update({"config_type": label,
                             "dataset_role": "train_or_numerical_diagnostic_only",
                             "launch_enabled": False,
                             "owner": None,
                             "separation_delta_A": delta,
                             "separation_direction": "outward_along_geometry_derived_surface_normal" if delta > 0 else "inward_toward_framework_along_geometry_derived_surface_normal",
                             "separation_method": "rigid translation of all four Ag atoms; framework fixed"})
        name = label + ".extxyz"
        output = HERE / "inputs" / name
        output.parent.mkdir(parents=True, exist_ok=True)
        write(output, variant, format="extxyz")
        reread = read(output)
        assert hash_file(output)
        assert np.array_equal(reread.arrays["lammps_id"], atoms.arrays["lammps_id"])
        assert np.array_equal(reread.arrays["central_pair"], atoms.arrays["central_pair"])
        assert np.array_equal(reread.cell.array, atoms.cell.array) and np.array_equal(reread.pbc, atoms.pbc)
        assert np.allclose(reread.positions, variant.positions, atol=1.0e-12, rtol=0)
        assert not any(key in reread.info for key in ("energy", "free_energy", "REF_energy", "PW_PBE_energy_eV"))
        assert not any(key in reread.arrays for key in ("forces", "REF_forces", "PW_PBE_forces"))
        records.append({"interface": source["interface"], "label": label,
                        "path": str(output.relative_to(REPO)), "sha256": hash_file(output),
                        "role": "train_or_numerical_diagnostic_only", "launch_enabled": False,
                        "owner": None, "parent_label": source["label"], "parent_path": source["path"],
                        "parent_input_sha256": source["sha256"],
                        "source_manifest_path": str((PARENT / "manifest.json").relative_to(REPO)),
                        "source_manifest_sha256": hash_file(PARENT / "manifest.json"),
                        "source_cif_sha256": SOURCE["source_cif_sha256"],
                        "source_family": "COD9009647 Ti3SiC2 CIF; all three terminations share this parent",
                        "atom_count": len(atoms), "formula": atoms.get_chemical_formula(),
                        "ordered_symbols_sha256": hashlib.sha256("\n".join(symbols).encode()).hexdigest(),
                        "ordered_ids_sha256": hashlib.sha256(atoms.arrays["lammps_id"].tobytes()).hexdigest(),
                        "top_surface_species": termination, "surface_normal_outward": normal.tolist(),
                        "translation_vector_A": translation.tolist(), "translation_delta_A": delta,
                        "pbc": atoms.pbc.tolist(), "cell_A": atoms.cell.array.tolist(),
                        "interface_area_A2": area, "parent_gap_A": gap, "proposed_gap_A": actual_gap,
                        "parent_vacuum_A": parent_vacuum, "proposed_vacuum_A": actual_vacuum,
                        "parent_marked_pair_ids": [int(x) for x in central_ids],
                        "parent_marked_pair_species": pair_species,
                        "parent_marked_pair_distance_A": parent_marked_distance,
                        "proposed_marked_pair_distance_A": marked_distance,
                        "parent_Ag_framework_species_minima_A": parent_species_min,
                        "proposed_Ag_framework_species_minima_A": minima,
                        "overlap_rule": f"reject if periodic Ag-framework distance < {overlap_scale:.2f}*(ASE covalent_radius(Ag)+covalent_radius(X)); geometry screening only",
                        "overlap_thresholds_A": thresholds, "unsafe_overlaps": unsafe,
                        "changed_atoms": 4, "changed_species": "Ag only",
                        "energy_labels_present": False, "force_labels_present": False,
                        "method": "unassigned; common target pending A mesh/boundary/energy review"})

manifest = {"status": "GEOMETRY_ONLY_PROPOSALS_VERIFIED_PENDING_A_REVIEW",
            "source_manifest_path": str((PARENT / "manifest.json").relative_to(REPO)),
            "source_manifest_sha256": hash_file(PARENT / "manifest.json"),
            "source_cif_sha256": SOURCE["source_cif_sha256"],
            "source_family_count": 1,
            "independent_validation": False,
            "method": None,
            "records": records,
            "scope": "Six rigid whole-Ag-layer normal-translation proposals; no calculator, DFT, model inference, training, MD or heating executed.",
            "limitations": ["All three interfaces derive from one shared COD9009647 CIF parent; they are not three independent validation families.",
                            "Parent slabs are unrelaxed, one-sided thin slabs; proposals are not equilibrium gaps or interface strengths.",
                            "The 0.80 covalent-radius overlap threshold is a conservative geometry screen, not a chemical bond or stability criterion.",
                            "No common DFT target or launch owner is assigned; every record is launch-disabled and owner-null."]}
(HERE / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps({"proposals": len(records), "source_families": 1,
                  "minimum_proposed_Ag_framework_distance_A": min(min(x["proposed_Ag_framework_species_minima_A"].values()) for x in records),
                  "max_gap_A": max(x["proposed_gap_A"] for x in records)}, indent=2))
