#!/usr/bin/env python3
"""Build three unlabeled Ag/Ti3SiC2 interface candidates from the 1967 source.

The builder reads only the 1967 reconstructed CIF and geometry-only columns
from train/dev/acquisition references.  It does not read any DFT labels.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path

import numpy as np
from ase import Atom
from ase.io import read, write

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
P = REPO / "research/mace-v12-transfer/periodic_interface_v4"
SOURCE_CIF = REPO / "research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_jeitschko1967_source_audit_20261027/Ti3SiC2_Jeitschko1967_fulltext_Fourier_variant.cif"
SOURCE_AUDIT = REPO / "research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_jeitschko1967_source_audit_20261027/README.md"
SOURCE_PAPER_SHA256 = "6765f99b3cc83daf3e48e4674dbeba109e5fffef8fa66cd53e266733ac4a4a63"
TRAIN = P / "mace_periodic_v19_residual_response/data/train.extxyz"
DEV = P / "mace_periodic_v19_residual_response/data/valid.extxyz"
RESIDUAL = P / "pbe_interface_v19_residual_acquisition/inputs"
V17_DEV = REPO / "research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set/geometries"
GEOMETRY_HELPER = REPO / "research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set/generate_split_candidates.py"
OUT = HERE / "inputs"

CONTACTS = (("AgC", "C", 2.4, 20300000), ("AgSi", "Si", 2.6, 20301000), ("AgTi", "Ti", 2.6, 20302000))
REGISTRY_PRIMITIVE_FRAC = (0.37, 0.21)
PBC = [True, True, False]
TYPES = {"Ag": 1, "Ti": 2, "Si": 3, "C": 4}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_helper():
    spec = importlib.util.spec_from_file_location("geometry_helper", GEOMETRY_HELPER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def pair_minima(atoms):
    distances = atoms.get_all_distances(mic=True)
    symbols = atoms.get_chemical_symbols()
    minima = {}
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            pair = "-".join(sorted((symbols[i], symbols[j])))
            value = float(distances[i, j])
            if pair not in minima or value < minima[pair]["distance_A"]:
                minima[pair] = {"distance_A": value, "atom_indices": [i, j]}
    return minima


def geometry_sha(atoms):
    payload = {
        "symbols": atoms.get_chemical_symbols(),
        "positions_A": np.round(atoms.positions, 10).tolist(),
        "cell_A": np.round(atoms.cell.array, 10).tolist(),
        "pbc": atoms.pbc.astype(bool).tolist(),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def normalized_comparison(candidate, reference, helper):
    """Compare after scaling the small (<0.5%) lattice-parameter difference away."""
    normalized = candidate.copy()
    normalized.set_cell(reference.cell, scale_atoms=True)
    comparison = helper.compare_geometry(normalized, reference)
    return {
        "comparison_after_cell_normalization": comparison,
        "normalization_scale_factors": (reference.cell.lengths() / candidate.cell.lengths()).tolist(),
        "note": "A source-cell difference is normalized only for duplicate screening; the registered DFT input retains the 1967 source cell.",
    }


def main():
    helper = load_helper()
    bulk = read(SOURCE_CIF)
    assert bulk.get_chemical_formula() == "C4Si2Ti6" and len(bulk) == 12
    source_sha = sha(SOURCE_CIF)
    source_params = {
        "spacegroup": 194, "a_A": float(bulk.cell.lengths()[0]),
        "c_A": float(bulk.cell.lengths()[2]), "source_cif_sha256": source_sha,
    }
    OUT.mkdir(parents=True, exist_ok=True)

    # These readers select only geometry fields from extxyz. No energy, force,
    # validation metric, or sealed/test file is opened.
    train_frames = helper.extxyz_geometry_frames(TRAIN)
    dev_frames = helper.extxyz_geometry_frames(DEV)
    residual_frames = []
    residual_inputs = []
    for iface, contact, _, _ in CONTACTS:
        mode = "06" if iface == "AgC" else "03"
        for sign in ("m", "p"):
            path = RESIDUAL / f"{iface}_residual_mode_{sign}_{mode}_v19_proposal.extxyz"
            residual_frames.extend(helper.extxyz_geometry_frames(path))
            residual_inputs.append({"path": str(path.relative_to(REPO)), "sha256": sha(path)})
    v17_frames = []
    v17_inputs = []
    for iface in ("AgC", "AgSi", "AgTi"):
        path = V17_DEV / f"{iface}_v17_dev_validation_01.extxyz"
        v17_frames.extend(helper.extxyz_geometry_frames(path))
        v17_inputs.append({"path": str(path.relative_to(REPO)), "sha256": sha(path)})

    records = []
    atoms_by_label = {}
    for iface, contact, gap_A, id_base in CONTACTS:
        frac = bulk.get_scaled_positions(wrap=True)
        symbols = bulk.get_chemical_symbols()
        top = float(np.max(frac[np.asarray(symbols) == contact, 2]))
        cut = (top + 0.001) % 1.0
        unit = bulk.copy()
        frac[:, 2] = (frac[:, 2] - cut) % 1.0
        unit.set_scaled_positions(frac)

        slab = unit.repeat((2, 2, 1))
        slab.positions[:, 2] += 6.0
        top_z = float(slab.positions[:, 2].max())
        cell = slab.cell.array.copy()
        cell[2] = [0.0, 0.0, top_z + gap_A + 8.0]
        slab.set_cell(cell, scale_atoms=False)
        # Four Ag atoms form a fixed coherent monolayer in the registry used by
        # the original CIF-derived proposals. They are regenerated on the
        # Jeitschko cell; no old Ag coordinates or DFT labels are copied.
        for i in range(2):
            for j in range(2):
                xy = (i + REGISTRY_PRIMITIVE_FRAC[0]) * unit.cell[0] + (j + REGISTRY_PRIMITIVE_FRAC[1]) * unit.cell[1]
                slab.append(Atom("Ag", position=xy + np.array([0.0, 0.0, top_z + gap_A])))
        slab.set_pbc(PBC)
        slab.wrap()

        ids = np.arange(id_base, id_base + len(slab), dtype=int)
        slab.new_array("lammps_id", ids)
        slab.new_array("lammps_type", np.asarray([TYPES[s] for s in slab.get_chemical_symbols()], dtype=int))
        slab.new_array("central_pair", np.zeros(len(slab), dtype=int))
        symbols = slab.get_chemical_symbols()
        dist = slab.get_all_distances(mic=True)
        np.fill_diagonal(dist, np.inf)
        ag_indices = np.flatnonzero(np.asarray(symbols) == "Ag")
        contact_indices = np.flatnonzero(np.asarray(symbols) == contact)
        local = dist[np.ix_(ag_indices, contact_indices)]
        ai, ci = np.unravel_index(int(np.argmin(local)), local.shape)
        ag_i, contact_i = int(ag_indices[ai]), int(contact_indices[ci])
        markers = np.zeros(len(slab), dtype=int)
        markers[[ag_i, contact_i]] = 1
        slab.set_array("central_pair", markers)

        label = f"{iface}_JN1967_registryB_v21_acq01"
        slab.info = {
            "config_type": label,
            "candidate_role": "train_acquisition_candidate_pending_audit",
            "dataset_role": "train_acquisition_candidate_pending_audit",
            "source_doi": "10.1007/BF00899949",
            "source_family": "Jeitschko_Nowotny_1967_single_crystal",
            "source_cif_sha256": source_sha,
            "source_paper_sha256": SOURCE_PAPER_SHA256,
            "registry_primitive_fractional": f"{REGISTRY_PRIMITIVE_FRAC[0]},{REGISTRY_PRIMITIVE_FRAC[1]}",
            "gap_A": gap_A,
            "fixed_geometry_unlabeled": True,
            "full_occupancy_model_assumption": True,
            "pbc": "T T F",
        }
        path = OUT / f"{label}.extxyz"
        write(path, slab, format="extxyz")
        check = read(path)
        assert len(check) == 52 and check.get_chemical_formula() == "C16Ag4Si8Ti24"
        assert check.pbc.tolist() == PBC
        assert "PW_PBE_forces" not in check.arrays and "REF_forces" not in check.arrays
        minima = pair_minima(check)
        marked = np.flatnonzero(check.arrays["central_pair"])
        assert len(marked) == 2 and {check[i].symbol for i in marked} == {"Ag", contact}

        refs = {"train36_geometry_only": train_frames, "V19_dev3_geometry_only": dev_frames,
                "V19_residual_inputs_geometry_only": residual_frames,
                "V17_dev_geometry_only": v17_frames}
        screens = {}
        for group, structures in refs.items():
            values = []
            for ref in structures:
                if len(ref) != len(check) or Counter(ref.get_chemical_symbols()) != Counter(check.get_chemical_symbols()) or not np.array_equal(ref.pbc, check.pbc):
                    continue
                values.append(normalized_comparison(check, ref, helper))
            cmp_values = [v["comparison_after_cell_normalization"] for v in values if v["comparison_after_cell_normalization"] is not None]
            screens[group] = {
                "frames_geometry_only": len(structures),
                "same_composition_and_PBC_frames": sum(1 for ref in structures if len(ref) == len(check) and Counter(ref.get_chemical_symbols()) == Counter(check.get_chemical_symbols()) and np.array_equal(ref.pbc, check.pbc)),
                "cell_normalized_comparable_frames": len(cmp_values),
                "exact_hits": sum(bool(v["exact_duplicate"]) for v in cmp_values),
                "near_hits": sum(bool(v["near_duplicate"]) for v in cmp_values),
                "comparisons": values,
            }

        records.append({
            "label": label, "interface": iface, "contact_species": contact,
            "input": str(path.relative_to(REPO)), "input_sha256": sha(path),
            "geometry_sha256": geometry_sha(check), "atoms": len(check),
            "formula": check.get_chemical_formula(), "cell_A": check.cell.array.tolist(),
            "cell_lengths_A": check.cell.lengths().tolist(), "pbc": check.pbc.tolist(),
            "source_cif_sha256": source_sha, "source_paper_sha256": SOURCE_PAPER_SHA256,
            "surface_cut": {"fractional_cut": cut, "top_z_A": top_z, "gap_to_Ag_plane_A": gap_A},
            "Ag_registry_primitive_fractional": list(REGISTRY_PRIMITIVE_FRAC),
            "coherent_in_plane_strain_note": "Ag monolayer is fixed to the 1967 substrate in-plane cell; no relaxation or equilibrium-interface claim.",
            "marked_contact_pair": {"species": f"Ag-{contact}", "indices": [ag_i, contact_i], "distance_A": float(dist[ag_i, contact_i])},
            "minimum_pair_distances_A": {k: v["distance_A"] for k, v in minima.items()},
            "role": "train_acquisition_candidate_pending_audit",
            "launch_enabled": False,
            "geometry_screens": screens,
        })
        atoms_by_label[label] = check

    pairwise = {}
    labels = list(atoms_by_label)
    for i, left in enumerate(labels):
        for right in labels[i + 1:]:
            a = atoms_by_label[left].copy()
            b = atoms_by_label[right].copy()
            b.set_cell(a.cell, scale_atoms=True)
            pairwise[f"{left}__{right}"] = helper.compare_geometry(a, b)

    manifest = {
        "task_id": "user_directed_Jeitschko1967_AgTi3SiC2_interface_acquisition_v21",
        "status": "geometry_candidates_built_preflight_required",
        "authorization_basis": "User requested continued modeling and said to calculate needed data directly; DFT acquisition remains candidate-only pending slab memory/k-point preflight and data-role review.",
        "launch_enabled": False, "training_enabled": False, "dataset_integration_enabled": False,
        "source": {"doi": "10.1007/BF00899949", "source_family": "Jeitschko_Nowotny_1967_single_crystal", "source_cif": str(SOURCE_CIF.relative_to(REPO)), "source_cif_sha256": source_sha, "uploaded_paper_sha256": SOURCE_PAPER_SHA256, "geometry_variant": "full-text Fourier-estimate reconstruction, not an ICSD export"},
        "method_candidate": {"xc": "PBE", "basis": "plane wave", "cutoff_eV": 500, "smearing_eV": 0.10, "pbc": PBC, "poissonsolver": {"dipolelayer": "xy"}, "candidate_kpts": [6, 6, 1], "slab_control_kpts": [[6, 6, 1], [8, 8, 1]], "mpi_ranks": 4, "maximum_concurrent_jobs": 1, "GPAW_version": "26.7.0", "ASE_version": "3.29.0", "PAW_data_version": "1.2.1", "PAW_setup_sha256": {"Ag.PBE.gz": "bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97", "Ti.PBE.gz": "bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f", "Si.PBE.gz": "2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0", "C.PBE.gz": "fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983"}, "scf_thresholds": {"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5}, "mixer": {"type": "Pulay", "beta": 0.05, "nmaxold": 8, "weight": 100}},
        "preflight_jobs": [
            {"label": "AgSi_JN1967_registryB_k6x6x1_sigma0p10_preflight_v21", "geometry_label": "AgSi_JN1967_registryB_v21_acq01", "input": f"inputs/{Path(records[1]['input']).name}", "input_sha256": records[1]["input_sha256"], "kpts": [6, 6, 1], "smearing_eV": 0.10, "role": "train_acquisition_candidate_pending_audit", "launch_enabled": True},
            {"label": "AgSi_JN1967_registryB_k8x8x1_sigma0p10_control_v21", "geometry_label": "AgSi_JN1967_registryB_v21_acq01", "input": f"inputs/{Path(records[1]['input']).name}", "input_sha256": records[1]["input_sha256"], "kpts": [8, 8, 1], "smearing_eV": 0.10, "role": "numerical_kmesh_control_only_not_training", "launch_enabled": True},
        ],
        "energy_convention": "Record native extrapolated total energy and free energy separately; forces use Fermi-Dirac width 0.10 eV. No conversion or training target is assumed.",
        "runtime_guard": {"cgroup_memory_limit_bytes": 34359738368, "stop_threshold_bytes": 32212254720, "minimum_free_disk_bytes": 5000000000, "automatic_retries": 0, "checkpoints_archived": False},
        "role_policy": "All three 6x6x1 geometries are train-acquisition candidates pending B/source-family/geometry and energy-role audit. The 8x8x1 AgSi same-geometry endpoint is control-only and must not enter training.",
        "reference_geometry_hashes": {"train36": sha(TRAIN), "dev3": sha(DEV), "residual_inputs": residual_inputs, "v17_dev_inputs": v17_inputs},
        "records": records,
        "pairwise_candidate_comparisons_after_cell_normalization": pairwise,
        "limitations": ["Unrelaxed single Ag monolayer interface probes are not equilibrium morphology or stability claims.", "The source paper discusses possible nonstoichiometry; these inputs assume ideal full occupancy.", "The Ag registry and 2.4/2.6 A contact gaps are fixed acquisition choices, not source-measured values.", "The 1967 source family differs from V19 dev3, but all three candidates share the same 1967 parent family.", "Training remains disabled until method/energy convention and source-family/geometry audits are complete."],
    }
    (HERE / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    audit = {"status": "GEOMETRY_BUILD_COMPLETE_PRELAUNCH_AUDIT_PENDING", "source_parameters": source_params, "records": records, "pairwise_candidates": pairwise, "geometry_reader_scope": "train36, V19 dev3, V19 residual inputs and V17 dev geometries; geometry columns only", "labels_or_forces_read": False, "test_or_sealed_data_opened": False}
    (HERE / "geometry_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"task_id": manifest["task_id"], "records": [{"label": r["label"], "sha256": r["input_sha256"], "marked_pair_distance_A": r["marked_contact_pair"]["distance_A"], "near_hits": sum(v["near_hits"] for v in r["geometry_screens"].values())} for r in records], "pairwise": pairwise}, indent=2))


if __name__ == "__main__":
    main()
