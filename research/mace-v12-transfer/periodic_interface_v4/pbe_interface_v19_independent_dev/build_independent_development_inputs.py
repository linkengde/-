#!/usr/bin/env python3
"""Build and geometry-audit fixed, label-free V19 development inputs.

Only symbols, coordinates, atom IDs, PBC, and cells are read from existing
structures. No reference energy/force fields or sealed test files are opened.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cosine

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
P = REPO / "research/mace-v12-transfer/periodic_interface_v4"
SOURCE_ROOT = P / "v19_cif_parent_structures"
OUT = HERE / "inputs"
GEOMETRY_HELPER = REPO / "research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set/generate_split_candidates.py"
TRAIN = P / "mace_periodic_v19_residual_response/training_candidate/train.extxyz"
RESIDUAL = P / "pbe_interface_v19_residual_acquisition/inputs"
V17_DEV = REPO / "research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set/geometries"
AUDIT_SOURCE = REPO / "research/mace-v12-transfer/coordination/reports/window-b/v19_uploaded_CIF_rebuild_review/audit.json"
SOURCE_CIF = SOURCE_ROOT / "source/Ti3SiC2_COD_9009647.cif"
INTERFACES = {"AgC": "C", "AgSi": "Si", "AgTi": "Ti"}
AG_REGISTRY_SHIFT_FRAC = (-0.27, 0.19)
METHOD = {"xc": "PBE", "basis": "PW", "cutoff_eV": 500, "kpts": [6, 6, 1], "smearing_eV": 0.10}
PBC = [True, True, False]
POISSONSOLVER = {"dipolelayer": "xy"}
LAMMPS_TYPES = {"Ag": 1, "Ti": 2, "Si": 3, "C": 4}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_geometry_helper():
    spec = importlib.util.spec_from_file_location("v17_geometry_helpers", GEOMETRY_HELPER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def rdf_fingerprint(atoms: object) -> np.ndarray:
    # Descriptive normalized, species-pair-resolved MIC distance histogram.
    distances = atoms.get_all_distances(mic=True)
    symbols = atoms.get_chemical_symbols()
    bins = np.linspace(0.0, 8.0, 81)
    values: list[float] = []
    for e1 in ("Ag", "C", "Si", "Ti"):
        for e2 in ("Ag", "C", "Si", "Ti"):
            pair_distances = [
                float(distances[i, j])
                for i in range(len(atoms))
                for j in range(i + 1, len(atoms))
                if tuple(sorted((symbols[i], symbols[j]))) == tuple(sorted((e1, e2)))
            ]
            hist = np.histogram(pair_distances, bins=bins)[0].astype(float)
            if pair_distances:
                hist /= len(pair_distances)
            values.extend(hist.tolist())
    return np.asarray(values, dtype=float)


def geometry_sha(atoms: object) -> str:
    payload = {
        "symbols": atoms.get_chemical_symbols(),
        "positions_A": np.round(atoms.positions, 10).tolist(),
        "cell_A": np.round(atoms.cell.array, 10).tolist(),
        "pbc": atoms.pbc.astype(bool).tolist(),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def create_inputs(helper: object) -> list[dict]:
    from ase import Atoms

    OUT.mkdir(parents=True, exist_ok=True)
    audit_source = json.loads(AUDIT_SOURCE.read_text())
    train = helper.extxyz_geometry_frames(TRAIN)
    residual = []
    residual_records = []
    for interface, contact in INTERFACES.items():
        for sign, mode in (("m", {"AgC": "06", "AgSi": "03", "AgTi": "03"}[interface]),
                           ("p", {"AgC": "06", "AgSi": "03", "AgTi": "03"}[interface])):
            label = f"{interface}_residual_mode_{sign}_{mode}_v19_acq"
            path = RESIDUAL / f"{interface}_residual_mode_{sign}_{mode}_v19_proposal.extxyz"
            residual.extend(helper.extxyz_geometry_frames(path))
            residual_records.append({"label": label, "path": str(path.relative_to(REPO)), "sha256": sha(path)})
    v17_dev = []
    v17_records = []
    for interface in INTERFACES:
        label = f"{interface}_v17_dev_validation_01"
        path = V17_DEV / f"{label}.extxyz"
        v17_dev.extend(helper.extxyz_geometry_frames(path))
        v17_records.append({"label": label, "path": str(path.relative_to(REPO)), "sha256": sha(path)})

    source_hash = sha(SOURCE_CIF)
    entries = []
    atoms_by_interface = {}
    for index, (interface, contact) in enumerate(INTERFACES.items()):
        parent = SOURCE_ROOT / "inputs" / f"{interface}_cod9009647_slab_v19_dev_proposal.extxyz"
        parent_hash = sha(parent)
        frames = helper.extxyz_geometry_frames(parent)
        assert len(frames) == 1
        atoms = frames[0]
        ag_mask = np.asarray([symbol == "Ag" for symbol in atoms.get_chemical_symbols()])
        atoms.positions[ag_mask] += AG_REGISTRY_SHIFT_FRAC[0] * atoms.cell[0] + AG_REGISTRY_SHIFT_FRAC[1] * atoms.cell[1]
        atoms.set_pbc(PBC)
        atoms.wrap()
        symbol_list = atoms.get_chemical_symbols()
        if "lammps_type" in atoms.arrays:
            atoms.set_array("lammps_type", np.asarray([LAMMPS_TYPES[s] for s in symbol_list], dtype=int))
        else:
            atoms.new_array("lammps_type", np.asarray([LAMMPS_TYPES[s] for s in symbol_list], dtype=int))
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        target_mask = np.asarray([symbol == contact for symbol in atoms.get_chemical_symbols()])
        ag_indices = np.flatnonzero(ag_mask)
        target_indices = np.flatnonzero(target_mask)
        local = distances[np.ix_(ag_indices, target_indices)]
        ag_index, target_index = np.unravel_index(int(np.argmin(local)), local.shape)
        markers = np.zeros(len(atoms), dtype=int)
        markers[ag_indices[ag_index]] = 1
        markers[target_indices[target_index]] = 1
        if "central_pair" in atoms.arrays:
            atoms.set_array("central_pair", markers)
        else:
            atoms.new_array("central_pair", markers)
        # The inputs are written from this geometry object only; no REF columns exist.
        label = f"{interface}_v19_independent_dev_01"
        filename = f"{label}.extxyz"
        input_path = OUT / filename
        symbols = atoms.get_chemical_symbols()
        ids = atoms.arrays.get("lammps_id")
        types = atoms.arrays.get("lammps_type")
        assert ids is not None and types is not None
        assert "REF_forces" not in atoms.arrays and "REF_energy" not in atoms.info
        cell = " ".join(f"{x:.12f}" for x in atoms.cell.array.reshape(-1))
        pbc = " ".join("T" if value else "F" for value in PBC)
        header = (
            f'{len(atoms)}\nLattice="{cell}" '
            'Properties=species:S:1:pos:R:3:lammps_id:I:1:lammps_type:I:1:central_pair:I:1 '
            f'config_type="{label}" candidate_role="independent_development" '
            f'dataset_role="independent_development" source_path="{SOURCE_CIF.relative_to(REPO)}" '
            f'source_sha256={source_hash} parent_geometry_path="{parent.relative_to(REPO)}" '
            f'parent_geometry_sha256={parent_hash} role_frozen_before_labeling=T '
            f'lineage_review="APPROVED_NOT_REUSED_V17_DEVELOPMENT_FAMILY" '
            f'ag_registry_shift_primitive_frac="{AG_REGISTRY_SHIFT_FRAC[0]},{AG_REGISTRY_SHIFT_FRAC[1]}" '
            f'proposal_only_no_DFT_label=T pbc="{pbc}"\n'
        )
        lines = [header]
        for idx, atom in enumerate(atoms):
            pos = " ".join(f"{v:.12f}" for v in atom.position)
            lines.append(f"{atom.symbol} {pos} {int(ids[idx])} {int(types[idx])} {int(markers[idx])}\n")
        input_path.write_text("".join(lines), encoding="utf-8")
        reread = helper.extxyz_geometry_frames(input_path)[0]
        assert len(reread) == 52 and reread.pbc.tolist() == PBC
        atoms_by_interface[interface] = reread
        pair_minima: dict[str, float] = {}
        for i in range(len(reread)):
            for j in range(i + 1, len(reread)):
                pair = "-".join(sorted((symbols[i], symbols[j])))
                pair_minima[pair] = min(pair_minima.get(pair, float("inf")), float(distances[i, j]))
        entries.append({
            "interface": interface,
            "label": label,
            "input": str(input_path.relative_to(REPO)),
            "input_sha256": sha(input_path),
            "geometry_sha256": geometry_sha(reread),
            "parent_geometry": str(parent.relative_to(REPO)),
            "parent_geometry_sha256": parent_hash,
            "source_CIF": str(SOURCE_CIF.relative_to(REPO)),
            "source_CIF_sha256": source_hash,
            "atoms": len(reread),
            "formula": reread.get_chemical_formula(),
            "cell_A": reread.cell.array.tolist(),
            "pbc": reread.pbc.tolist(),
            "target_contact_pair": {
                "species": f"Ag-{contact}",
                "ag_lammps_id": int(ids[ag_indices[ag_index]]),
                "framework_lammps_id": int(ids[target_indices[target_index]]),
                "distance_A": float(distances[ag_indices[ag_index], target_indices[target_index]]),
            },
            "minimum_pair_distances_A": pair_minima,
            "role": "independent_development",
            "A_lineage_review": "APPROVED_NOT_REUSED_V17_DEVELOPMENT_FAMILY",
            "labels_read": False,
        })

    # Existing geometry comparator checks all candidates that share composition,
    # cell, and PBC. A separate normalized RDF descriptor is recorded for
    # dissimilar-sized families where the comparator correctly returns N/A.
    comparisons = {}
    refs = {"train36_geometry_only": train, "six_v19_residual_geometry_only": residual, "three_v17_dev_geometry_only": v17_dev}
    candidate_features = {name: rdf_fingerprint(atoms) for name, atoms in atoms_by_interface.items()}
    for name, atoms in atoms_by_interface.items():
        comparisons[name] = {}
        for family, structures in refs.items():
            pair_checks = [helper.compare_geometry(atoms, other) for other in structures]
            comparable = [result for result in pair_checks if result is not None]
            similarities = [1.0 - cosine(candidate_features[name], rdf_fingerprint(other)) for other in structures]
            comparisons[name][family] = {
                "n_geometries_checked": len(structures),
                "n_comparable_same_composition_cell_pbc": len(comparable),
                "exact_duplicates": sum(bool(result["exact_duplicate"]) for result in comparable),
                "near_duplicates": sum(bool(result["near_duplicate"]) for result in comparable),
                "max_species_pair_RDF_cosine_similarity_descriptive_only": float(max(similarities)),
                "labels_or_force_values_read": False,
            }
    candidate_pairs = {}
    names = list(atoms_by_interface)
    for i, left in enumerate(names):
        for right in names[i + 1:]:
            result = helper.compare_geometry(atoms_by_interface[left], atoms_by_interface[right])
            candidate_pairs[f"{left}__{right}"] = result

    report = {
        "status": "PASS_GEOMETRY_AND_ROLE_SCREEN_PENDING_EXTERNAL_B_AUDIT",
        "data_role_frozen_before_DFT": True,
        "test_data_opened": False,
        "training36_labels_or_force_values_read": False,
        "training36_geometry_only_manifest_sha256": sha(TRAIN),
        "source_CIF": str(SOURCE_CIF.relative_to(REPO)),
        "source_CIF_sha256": source_hash,
        "source_evidence_limit": "User-uploaded COD9009647 CIF is internally consistent and rebuild reproduces byte-identically; COD remote bytes/checksum remain unverified after the documented HTTP 403.",
        "shared_parent_limit": "All three development structures share one published Ti3SiC2 CIF substrate source but use distinct Ag-X terminations/contacts; they are independent from train36 and the six V19 residual-label geometry family, not three independent bulk-source families.",
        "fixed_DFT_method": {
            "xc": "PBE", "basis": "PW", "cutoff_eV": 500, "kpts": [6, 6, 1],
            "smearing_eV": 0.10, "pbc": PBC, "poissonsolver": POISSONSOLVER,
            "fixed_ions": True, "energy_convention": "native extrapolated total energy plus separately recorded free energy",
        },
        "reference_geometry_sets": {
            "train36": {"path": str(TRAIN.relative_to(REPO)), "sha256": sha(TRAIN), "frames": len(train), "reader": "geometry-only extxyz columns"},
            "six_v19_residual_inputs": residual_records,
            "three_public_v17_development_parents": v17_records,
        },
        "candidate_records": entries,
        "geometry_comparisons": comparisons,
        "pairwise_candidate_comparisons": candidate_pairs,
        "RDF_note": "Species-pair distance histogram cosine is descriptive only; exact/near classification follows the published V17 split comparator when formula, cell and PBC are comparable.",
    }
    (HERE / "geometry_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return entries


if __name__ == "__main__":
    helper = load_geometry_helper()
    rows = create_inputs(helper)
    print(json.dumps({"inputs": rows, "geometry_audit": str((HERE / "geometry_audit.json").relative_to(REPO))}, indent=2))
