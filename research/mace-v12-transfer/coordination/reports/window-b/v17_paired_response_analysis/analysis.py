#!/usr/bin/env python3
"""Matched V17 +/- diagnostic response analysis with frozen V16 CPU inference."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from ase.geometry import find_mic
from ase.io import read


PROJECT = Path("research/mace-v12-transfer/periodic_interface_v4")
V14_DATA = PROJECT / "mace_periodic_v14_interface_energy/data/train.extxyz"
A_DIR = PROJECT / "pbe_interface_v17_main_targeted_acquisition"
B_DIR = PROJECT / "pbe_interface_v17_parallel_targeted_acquisition"
PAIR_LABELS = {
    "AgC_framework_neighbor": {
        "interface": "AgC", "mode": "framework_neighbor", "step_A": 0.04,
        "positive": "AgC_framework_neighbor_p_acq_v17_01",
        "negative": "AgC_framework_neighbor_m_acq_v17_01",
    },
    "AgSi_transverse_Ag": {
        "interface": "AgSi", "mode": "transverse_Ag", "step_A": 0.15,
        "positive": "AgSi_transverse_Ag_p_acq_v17_01",
        "negative": "AgSi_transverse_Ag_m_acq_v17_01",
    },
    "AgSi_framework_neighbor": {
        "interface": "AgSi", "mode": "framework_neighbor", "step_A": 0.04,
        "positive": "AgSi_framework_neighbor_p_acq_v17_01",
        "negative": "AgSi_framework_neighbor_m_acq_v17_01",
    },
    "AgTi_framework_neighbor": {
        "interface": "AgTi", "mode": "framework_neighbor", "step_A": 0.04,
        "positive": "AgTi_framework_neighbor_p_acq_v17_01",
        "negative": "AgTi_framework_neighbor_m_acq_v17_01",
    },
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def require(ok: bool, message: str) -> None:
    if not ok:
        raise SystemExit(message)


def geometry_sha(atoms) -> str:
    # Local canonical digest; compare with the upstream parent digest and expose any mismatch.
    payload = {
        "symbols": atoms.get_chemical_symbols(),
        "positions_A": np.round(atoms.positions, 10).tolist(),
        "cell_A": np.round(atoms.cell.array, 10).tolist(),
        "pbc": atoms.pbc.astype(bool).tolist(),
        "lammps_id": atoms.arrays.get("lammps_id", np.array([], dtype=int)).tolist(),
    }
    return sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def get_force(atoms, candidates: tuple[str, ...]) -> np.ndarray:
    for key in candidates:
        if key in atoms.arrays:
            force = np.asarray(atoms.arrays[key], dtype=float)
            if force.shape == (len(atoms), 3):
                return force
    raise SystemExit(f"No force array among {candidates}")


def id_map(atoms) -> dict[int, int]:
    ids = np.asarray(atoms.arrays.get("lammps_id", []), dtype=int)
    require(len(ids) == len(atoms) and len(set(ids.tolist())) == len(ids), "Missing/duplicate persistent atom IDs")
    return {int(atom_id): index for index, atom_id in enumerate(ids)}


def direct_pair(atoms, pair_ids: dict[str, int], forces: np.ndarray) -> dict:
    mapping = id_map(atoms)
    ia, ix = mapping[int(pair_ids["Ag_id"])], mapping[int(pair_ids["contact_id"])]
    vector = np.asarray(atoms.positions[ix] - atoms.positions[ia], dtype=float)
    distance = float(np.linalg.norm(vector))
    require(distance > 0, "Zero direct marked-pair distance")
    unit = vector / distance
    force_difference = forces[ix] - forces[ia]
    longitudinal = float(np.dot(force_difference, unit))
    transverse = force_difference - longitudinal * unit
    return {
        "distance_A_direct": distance,
        "unit_Ag_to_X_direct": unit.tolist(),
        "Ag_projection_eV_A": float(np.dot(forces[ia], unit)),
        "X_projection_eV_A": float(np.dot(forces[ix], unit)),
        "separating_force_eV_A": longitudinal,
        "transverse_force_xyz_eV_A": transverse.tolist(),
        "transverse_force_magnitude_eV_A": float(np.linalg.norm(transverse)),
    }


def clean_geometry(atoms):
    result = atoms.copy()
    result.calc = None
    keep = {"numbers", "positions", "lammps_id", "lammps_type", "central_pair"}
    for name in list(result.arrays):
        if name not in keep:
            del result.arrays[name]
    result.info = {}
    return result


def force_stats(delta: np.ndarray) -> dict:
    norms = np.linalg.norm(delta, axis=1)
    return {
        "atoms": int(len(delta)),
        "mean_force_change_norm_eV_A": float(np.mean(norms)) if len(norms) else None,
        "rms_force_change_norm_eV_A": float(np.sqrt(np.mean(norms**2))) if len(norms) else None,
        "max_force_change_norm_eV_A": float(np.max(norms)) if len(norms) else None,
    }


def vector_rmse(delta: np.ndarray) -> float | None:
    return float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))) if len(delta) else None


def load_archives(repo: Path, manifests: dict[str, dict], archive_manifests: dict[str, dict], label: str) -> dict:
    owner = "A" if label.endswith("_p_acq_v17_01") else "B"
    base = repo / (A_DIR if owner == "A" else B_DIR)
    manifest = manifests[owner]
    record = next((r for r in manifest["records"] if r["label"] == label), None)
    require(record is not None, f"Missing manifest record: {label}")
    archive_manifest = archive_manifests[owner]
    require(archive_manifest.get("complete") is True, f"Archive manifest incomplete: {owner}")
    require(any(r.get("label") == label for r in archive_manifest.get("records", [])), f"Archive manifest missing {label}")
    folder = base / "calculations" / label
    verification = json.loads((folder / "verification.json").read_text())
    require(verification.get("status") == "PASS" and verification.get("label") == label, f"Verification failed: {label}")
    require(bool(verification.get("checks")) and all(verification["checks"].values()), f"False archive check: {label}")
    input_path = base / record["input"]
    output_path = folder / f"{label}_PW_PBE.extxyz"
    summary_path, progress_path, log_path = folder / "summary.json", folder / "progress.json", folder / "gpaw.log"
    summary = json.loads(summary_path.read_text())
    progress = json.loads(progress_path.read_text())
    require(digest(input_path) == record["input_sha256"] == summary.get("source_sha256"), f"Input hash mismatch: {label}")
    require(summary.get("scf_converged") is True and progress.get("status") == "complete", f"SCF incomplete: {label}")
    require(summary.get("mpi_ranks") == 4, f"Wrong MPI rank count: {label}")
    for filename, expected in verification.get("sha256", {}).items():
        candidate = folder / filename
        require(candidate.is_file() and digest(candidate) == expected, f"Archived file hash mismatch: {label}/{filename}")
    source_atoms = read(input_path, format="extxyz")
    result_atoms = read(output_path, format="extxyz")
    require(len(result_atoms) == len(source_atoms), f"Atom count mismatch: {label}")
    require(result_atoms.get_chemical_symbols() == source_atoms.get_chemical_symbols(), f"Element order mismatch: {label}")
    require(np.array_equal(result_atoms.arrays["lammps_id"], source_atoms.arrays["lammps_id"]), f"Atom ID mismatch: {label}")
    require(np.allclose(result_atoms.positions, source_atoms.positions, atol=2e-8, rtol=0), f"Position mismatch: {label}")
    require(np.allclose(result_atoms.cell.array, source_atoms.cell.array, atol=1e-10, rtol=0), f"Cell mismatch: {label}")
    require(np.array_equal(result_atoms.pbc, source_atoms.pbc), f"PBC mismatch: {label}")
    require(digest(input_path) == verification.get("input_sha256"), f"Verification input hash mismatch: {label}")
    native = float(summary["energy_eV_cell"])
    force = get_force(result_atoms, ("PW_PBE_forces", "forces"))
    require(np.isfinite(native) and np.isfinite(force).all(), f"Nonfinite DFT label: {label}")
    free = summary.get("free_energy_eV_cell")
    require(free is None or np.isfinite(float(free)), f"Nonfinite free energy: {label}")
    require(np.isclose(float(result_atoms.info["PW_PBE_energy_eV"]), native, atol=1e-10, rtol=0), f"Energy mismatch: {label}")
    return {
        "label": label,
        "owner": owner,
        "record": record,
        "summary": summary,
        "verification": verification,
        "archive_files_sha256": {
            "verification.json": digest(folder / "verification.json"),
            **verification.get("sha256", {}),
        },
        "input_path": input_path,
        "input_atoms": source_atoms,
        "output_atoms": result_atoms,
        "native_energy_eV": native,
        "free_energy_eV": None if free is None else float(free),
        "forces": force,
        "method": summary["method"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--selection-record", required=True, type=Path)
    args = parser.parse_args()
    repo, out = args.repo_root.resolve(), args.output_dir.resolve()
    require(out.is_relative_to(repo / "research/mace-v12-transfer/coordination/reports/window-b/v17_paired_response_analysis") and out != repo / "research/mace-v12-transfer/coordination/reports/window-b/v17_paired_response_analysis", "Output must be under the assigned B report directory")
    require(not out.exists() or (out.is_dir() and not any(out.iterdir())), "Nonempty output directory refused")

    a_base, b_base = repo / A_DIR, repo / B_DIR
    a_manifest = json.loads((a_base / "input_manifest.json").read_text())
    b_manifest = json.loads((b_base / "input_manifest.json").read_text())
    a_archives = json.loads((a_base / "archive_manifest.json").read_text())
    b_archives = json.loads((b_base / "archive_manifest.json").read_text())
    require(a_manifest.get("owner_task") == "window-a" and b_manifest.get("owner_task") == "window-b", "Unexpected diagnostic acquisition ownership")
    input_manifests = {"A": a_manifest, "B": b_manifest}
    archive_manifests = {"A": a_archives, "B": b_archives}

    # Verify frozen V16 selection provenance before inference.
    model = args.model.resolve()
    selection_path = args.selection_record.resolve()
    selection = json.loads(selection_path.read_text())
    model_hash, selection_hash = digest(model), digest(selection_path)
    require(selection.get("reviewed_by") == "window-a", "V16 model selection record not reviewed by A")
    require(selection.get("model_sha256") == model_hash, "Frozen V16 model hash mismatch")
    require(selection.get("selection_uses_only_training_validation") is True and selection.get("blind_labels_used_for_selection") is False, "V16 selection was not isolated from blind labels")
    require(selection.get("completed_epochs") == 80 and selection.get("selected_epoch") == 80, "V16 completion/selected epoch record mismatch")

    v14_path = repo / V14_DATA
    expected_parent_hash = a_manifest["records"][0]["source_parent"]["dataset_sha256"]
    require(digest(v14_path) == expected_parent_hash, "V14 source parent dataset SHA256 mismatch")
    parents: dict[tuple[str, int], object] = {}
    structures: dict[str, dict] = {}
    audit_rows = []
    for pair_name, pair_spec in PAIR_LABELS.items():
        positive = load_archives(repo, input_manifests, archive_manifests, pair_spec["positive"])
        negative = load_archives(repo, input_manifests, archive_manifests, pair_spec["negative"])
        pr, nr = positive["record"], negative["record"]
        pparent, nparent = pr["source_parent"], nr["source_parent"]
        for field in ("dataset_path", "dataset_sha256", "frame_index_zero_based", "config_type", "frame_geometry_sha256", "formula", "atoms", "control_atom_ids"):
            require(pparent[field] == nparent[field], f"+/- parent mismatch for {pair_name}: {field}")
        require(pr["design"]["mode"] == nr["design"]["mode"] == pair_spec["mode"], f"Unexpected perturbation axis: {pair_name}")
        require(pr["design"]["sign"] == 1 and nr["design"]["sign"] == -1, f"Sign mismatch: {pair_name}")
        pchanged = pr["design"]["changed_atoms"]
        nchanged = nr["design"]["changed_atoms"]
        require([int(x["id"]) for x in pchanged] == [int(x["id"]) for x in nchanged], f"Changed atom IDs/order mismatch: {pair_name}")
        for plus_delta, minus_delta in zip(pchanged, nchanged):
            require(np.allclose(np.asarray(plus_delta["delta_xyz_A"], float), -np.asarray(minus_delta["delta_xyz_A"], float), atol=1e-12, rtol=0), f"+/- axis vectors are not opposite: {pair_name}")
        parent_key = (str(pparent["dataset_path"]), int(pparent["frame_index_zero_based"]))
        if parent_key not in parents:
            parent = read(repo / parent_key[0], index=parent_key[1])
            require(parent.info.get("config_type") == pparent["config_type"], f"Parent config_type mismatch: {pair_name}")
            require(parent.get_chemical_formula() == pparent["formula"] and len(parent) == pparent["atoms"], f"Parent composition mismatch: {pair_name}")
            parents[parent_key] = parent
        parent = parents[parent_key]
        pid = id_map(parent)
        require(int(parent.arrays["lammps_id"][np.flatnonzero(parent.arrays["central_pair"])[0]]) in pid, "Parent marked pair IDs missing")
        expected_marked = {int(v) for v in pparent["control_atom_ids"].values()}
        marked_ids = {int(parent.arrays["lammps_id"][i]) for i in np.flatnonzero(np.asarray(parent.arrays["central_pair"], bool))}
        require(marked_ids == expected_marked, f"Parent marked-pair IDs mismatch: {pair_name}")
        require(int(pr["central_pair"]["Ag_id"]) == int(pparent["control_atom_ids"]["Ag_id"]), f"Marked Ag ID disagrees with parent controls: {pair_name}")
        require(int(pr["central_pair"]["contact_id"]) == int(pparent["control_atom_ids"]["X_id"]), f"Marked contact ID disagrees with parent controls: {pair_name}")

        # Verify each child is exactly the declared small perturbation from that parent.
        moved_ids = [int(x["id"]) for x in pchanged]
        moved_by_id = {int(x["id"]): np.asarray(x["delta_xyz_A"], float) for x in pchanged}
        parent_id_map = id_map(parent)
        child_frames = []
        for item, sign in ((positive, 1), (negative, -1)):
            child = item["input_atoms"]
            cm = id_map(child)
            require(set(cm) == set(parent_id_map), f"Parent/child atom-ID set mismatch: {item['label']}")
            require(child.get_chemical_symbols() == parent.get_chemical_symbols(), f"Parent/child symbol order mismatch: {item['label']}")
            require(np.allclose(child.cell.array, parent.cell.array, atol=1e-10, rtol=0) and np.array_equal(child.pbc, parent.pbc), f"Parent/child cell/PBC mismatch: {item['label']}")
            deltas = {}
            for atom_id, index in cm.items():
                pi = parent_id_map[atom_id]
                require(child[index].symbol == parent[pi].symbol, f"Element mismatch for persistent atom ID {atom_id}: {item['label']}")
                displacement, _ = find_mic(child.positions[index] - parent.positions[pi], parent.cell, pbc=parent.pbc)
                deltas[atom_id] = np.asarray(displacement, float)
                if atom_id in moved_by_id:
                    require(np.allclose(deltas[atom_id], sign * moved_by_id[atom_id], atol=3e-8, rtol=0), f"Declared displacement mismatch: {item['label']}/{atom_id}")
                else:
                    require(np.linalg.norm(deltas[atom_id]) <= 3e-8, f"Undeclared atom moved: {item['label']}/{atom_id}")
            item["displacements_from_parent"] = deltas
            child_frames.append(item)
            structures[item["label"]] = {"atoms": clean_geometry(item["output_atoms"]), "dft": item}
        parent_label = f"parent_{pair_name}"
        structures[parent_label] = {"atoms": clean_geometry(parent), "parent": parent, "source_parent": pparent}

        # Generalized coordinate uses one scalar q for the declared common displacement axis.
        axis_vectors = [moved_by_id[atom_id] for atom_id in moved_ids]
        q_plus_values = [float(np.linalg.norm(v)) for v in axis_vectors]
        require(max(q_plus_values) - min(q_plus_values) < 1e-9, f"Nonuniform + displacement magnitudes: {pair_name}")
        q = float(np.mean(q_plus_values))
        unit_by_id = {atom_id: moved_by_id[atom_id] / np.linalg.norm(moved_by_id[atom_id]) for atom_id in moved_ids}
        first_unit = next(iter(unit_by_id.values()))
        require(all(np.allclose(unit, first_unit, atol=1e-10, rtol=0) for unit in unit_by_id.values()), f"Perturbation axis not collective/uniform: {pair_name}")
        pvector = np.asarray(positive["input_atoms"].positions[id_map(positive["input_atoms"])[int(pr["central_pair"]["Ag_id"])]] - positive["input_atoms"].positions[id_map(positive["input_atoms"])[int(pr["central_pair"]["contact_id"])]])
        direct_distance = float(np.linalg.norm(pvector))
        require(np.isclose(direct_distance, float(pr["central_pair"]["input_distance_A"]), atol=2e-8, rtol=0), f"Direct marked-pair distance convention mismatch: {pair_name}")
        for field in ("elements", "Ag_id", "contact_id"):
            require(pr["central_pair"][field] == nr["central_pair"][field], f"Marked pair {field} mismatch: {pair_name}")
        audit_rows.append({
            "pair": pair_name,
            "interface": pair_spec["interface"],
            "mode": pair_spec["mode"],
            "parent_key": {"dataset_path": parent_key[0], "frame_index_zero_based": parent_key[1]},
            "parent_geometry_sha256": pparent["frame_geometry_sha256"],
            "recomputed_parent_geometry_sha256_v17_canonical": geometry_sha(parent),
            "upstream_geometry_hash_exact_match": geometry_sha(parent) == pparent["frame_geometry_sha256"],
            "source_dataset_sha256": pparent["dataset_sha256"],
            "changed_atom_ids": moved_ids,
            "step_magnitude_A": q,
            "axis_xyz": first_unit.tolist(),
            "plus_input_sha256": pr["input_sha256"],
            "minus_input_sha256": nr["input_sha256"],
            "plus_archive_sha256": positive["archive_files_sha256"],
            "minus_archive_sha256": negative["archive_files_sha256"],
            "plus_method": positive["method"],
            "minus_method": negative["method"],
            "archive_verification": "PASS",
            "parent_source_method": parent.info.get("source_method"),
            "parent_energy_convention_metadata": "REF_energy stored; source frame metadata does not distinguish native extrapolated versus free energy",
            "parent_smearing_metadata": "not declared in source_method field",
        })

    # Load the selected V16 model in the reviewed CPU/MACE environment and infer all 8 children + 3 unique parents.
    try:
        import torch
        from mace.calculators import MACECalculator
    except Exception as exc:
        raise SystemExit(f"Frozen V16 inference dependencies unavailable: {type(exc).__name__}: {exc}")
    torch.set_num_threads(1)
    torch.set_default_dtype(torch.float64)
    loaded = torch.load(model, map_location="cpu", weights_only=False).double().eval()
    calculator = MACECalculator(models=[loaded], device="cpu", default_dtype="float64")
    predictions: dict[str, dict] = {}
    for name, entry in structures.items():
        atoms = entry["atoms"].copy()
        atoms.calc = calculator
        energy_pred = float(atoms.get_potential_energy())
        force_pred = np.asarray(atoms.get_forces(), dtype=float)
        require(np.isfinite(energy_pred) and force_pred.shape == (len(atoms), 3) and np.isfinite(force_pred).all(), f"Nonfinite V16 prediction: {name}")
        predictions[name] = {"energy_eV": energy_pred, "forces": force_pred}

    # Load parent DFT labels from the verified V14 training source. No blind/test files are touched.
    for pair_name, pair_spec in PAIR_LABELS.items():
        parent_entry = structures[f"parent_{pair_name}"]
        parent = parent_entry["parent"]
        parent_entry["parent_forces"] = get_force(parent, ("REF_forces", "forces"))
        parent_entry["parent_energy"] = float(parent.info["REF_energy"])
        require(np.isfinite(parent_entry["parent_forces"]).all() and np.isfinite(parent_entry["parent_energy"]), f"Nonfinite parent reference: {pair_name}")

    pair_rows = []
    group_rows = []
    atom_rows = []
    inference_rows = []
    for pair_name, pair_spec in PAIR_LABELS.items():
        pos = predictions[pair_spec["positive"]]
        neg = predictions[pair_spec["negative"]]
        pos_dft = structures[pair_spec["positive"]]["dft"]
        neg_dft = structures[pair_spec["negative"]]["dft"]
        parent_entry = structures[f"parent_{pair_name}"]
        parent = parent_entry["parent"]
        parent_ref_f = parent_entry["parent_forces"]
        parent_pred = predictions[f"parent_{pair_name}"]
        plus_atoms = pos_dft["input_atoms"]
        minus_atoms = neg_dft["input_atoms"]
        mapping = id_map(parent)
        changed_ids = audit_rows[[x["pair"] for x in audit_rows].index(pair_name)]["changed_atom_ids"]
        plus_rec, minus_rec = pos_dft["record"], neg_dft["record"]
        qplus = float(np.mean([np.linalg.norm(np.asarray(x["delta_xyz_A"], float)) for x in plus_rec["design"]["changed_atoms"]]))
        qminus = float(np.mean([np.linalg.norm(np.asarray(x["delta_xyz_A"], float)) for x in minus_rec["design"]["changed_atoms"]]))
        axis = np.asarray(audit_rows[[x["pair"] for x in audit_rows].index(pair_name)]["axis_xyz"], float)
        unit_by_id = {int(x["id"]): axis for x in plus_rec["design"]["changed_atoms"]}

        def generalized(force):
            return float(sum(np.dot(force[mapping[atom_id]], unit) for atom_id, unit in unit_by_id.items()))

        g_parent = generalized(parent_ref_f)
        g_p_dft = generalized(pos_dft["forces"])
        g_m_dft = generalized(neg_dft["forces"])
        g_parent_model = generalized(parent_pred["forces"])
        g_p_model = generalized(pos["forces"])
        g_m_model = generalized(neg["forces"])
        denominator = qplus + qminus
        native_slope = (pos_dft["native_energy_eV"] - neg_dft["native_energy_eV"]) / denominator
        free_slope = None if pos_dft["free_energy_eV"] is None or neg_dft["free_energy_eV"] is None else (pos_dft["free_energy_eV"] - neg_dft["free_energy_eV"]) / denominator
        model_slope = (pos["energy_eV"] - neg["energy_eV"]) / denominator
        pair_id_info = plus_rec["central_pair"]
        pair_atoms = {}
        pair_force_rows = {}
        for tag, forces in (("parent", parent_ref_f), ("plus", pos_dft["forces"]), ("minus", neg_dft["forces"])):
            atoms = parent if tag == "parent" else pos_dft["output_atoms"] if tag == "plus" else neg_dft["output_atoms"]
            pair_force_rows[tag] = direct_pair(atoms, pair_id_info, forces)
        for tag, forces in (("parent", parent_pred["forces"]), ("plus", pos["forces"]), ("minus", neg["forces"])):
            atoms = parent if tag == "parent" else pos_dft["output_atoms"] if tag == "plus" else neg_dft["output_atoms"]
            pair_atoms[f"model_{tag}"] = direct_pair(atoms, pair_id_info, forces)

        # Local shell is defined on the parent: non-Ag atoms within 3.5 A MIC of any moved atom.
        distances = parent.get_all_distances(mic=True)
        moved_indices = [mapping[i] for i in changed_ids]
        shell_indices = set()
        for idx, atom in enumerate(parent):
            if atom.symbol == "Ag" or int(parent.arrays["lammps_id"][idx]) in changed_ids:
                continue
            if min(float(distances[idx, m]) for m in moved_indices) <= 3.5:
                shell_indices.add(idx)
        pair_indices = {mapping[int(pair_id_info["Ag_id"])], mapping[int(pair_id_info["contact_id"])]}
        groups = {}
        for i, atom in enumerate(parent):
            aid = int(parent.arrays["lammps_id"][i])
            if aid in changed_ids:
                group = "moved_atoms"
            elif i in pair_indices:
                group = "marked_pair_other_member"
            elif i in shell_indices:
                group = "framework_neighbor_shell_nonAg_within_3p5A"
            elif atom.symbol == "Ag":
                group = "other_Ag"
            else:
                group = "framework_outside_shell"
            groups[aid] = group

        row = {
            "pair": pair_name,
            "interface": pair_spec["interface"],
            "mode": pair_spec["mode"],
            "step_A": float((qplus + qminus) / 2),
            "changed_atom_count": len(changed_ids),
            "changed_atom_ids": ";".join(map(str, changed_ids)),
            "parent_dataset_sha256": plus_rec["source_parent"]["dataset_sha256"],
            "parent_frame_index": plus_rec["source_parent"]["frame_index_zero_based"],
            "parent_geometry_sha256": plus_rec["source_parent"]["frame_geometry_sha256"],
            "plus_input_sha256": plus_rec["input_sha256"],
            "minus_input_sha256": minus_rec["input_sha256"],
            "direct_pair_distance_parent_A": pair_force_rows["parent"]["distance_A_direct"],
            "direct_pair_distance_plus_A": pair_force_rows["plus"]["distance_A_direct"],
            "direct_pair_distance_minus_A": pair_force_rows["minus"]["distance_A_direct"],
            "DFT_native_E_parent_eV": parent_entry["parent_energy"],
            "DFT_native_E_plus_eV": pos_dft["native_energy_eV"],
            "DFT_native_E_minus_eV": neg_dft["native_energy_eV"],
            "DFT_native_E_plus_minus_eV": pos_dft["native_energy_eV"] - neg_dft["native_energy_eV"],
            "DFT_native_centered_slope_eV_A": native_slope,
            "DFT_free_E_plus_eV": pos_dft["free_energy_eV"],
            "DFT_free_E_minus_eV": neg_dft["free_energy_eV"],
            "DFT_free_E_plus_minus_eV": None if free_slope is None else pos_dft["free_energy_eV"] - neg_dft["free_energy_eV"],
            "DFT_free_centered_slope_eV_A": free_slope,
            "DFT_parent_generalized_force_eV_A": g_parent,
            "DFT_plus_generalized_force_eV_A": g_p_dft,
            "DFT_minus_generalized_force_eV_A": g_m_dft,
            "DFT_centered_native_slope_plus_mean_force_residual_eV_A": native_slope + 0.5 * (g_p_dft + g_m_dft),
            "DFT_centered_free_slope_plus_mean_force_residual_eV_A": None if free_slope is None else free_slope + 0.5 * (g_p_dft + g_m_dft),
            "V16_E_parent_eV": parent_pred["energy_eV"],
            "V16_E_plus_eV": pos["energy_eV"],
            "V16_E_minus_eV": neg["energy_eV"],
            "V16_centered_slope_eV_A": model_slope,
            "V16_parent_generalized_force_eV_A": g_parent_model,
            "V16_plus_generalized_force_eV_A": g_p_model,
            "V16_minus_generalized_force_eV_A": g_m_model,
            "V16_centered_slope_plus_mean_force_residual_eV_A": model_slope + 0.5 * (g_p_model + g_m_model),
            "DFT_separating_force_parent_eV_A": pair_force_rows["parent"]["separating_force_eV_A"],
            "DFT_separating_force_plus_eV_A": pair_force_rows["plus"]["separating_force_eV_A"],
            "DFT_separating_force_minus_eV_A": pair_force_rows["minus"]["separating_force_eV_A"],
            "V16_separating_force_parent_eV_A": pair_atoms["model_parent"]["separating_force_eV_A"],
            "V16_separating_force_plus_eV_A": pair_atoms["model_plus"]["separating_force_eV_A"],
            "V16_separating_force_minus_eV_A": pair_atoms["model_minus"]["separating_force_eV_A"],
            "DFT_pair_transverse_parent_eV_A": pair_force_rows["parent"]["transverse_force_magnitude_eV_A"],
            "DFT_pair_transverse_plus_eV_A": pair_force_rows["plus"]["transverse_force_magnitude_eV_A"],
            "DFT_pair_transverse_minus_eV_A": pair_force_rows["minus"]["transverse_force_magnitude_eV_A"],
            "V16_pair_transverse_plus_eV_A": pair_atoms["model_plus"]["transverse_force_magnitude_eV_A"],
            "V16_pair_transverse_minus_eV_A": pair_atoms["model_minus"]["transverse_force_magnitude_eV_A"],
            "framework_shell_atoms": len(shell_indices),
            "framework_outside_atoms": sum(1 for i, atom in enumerate(parent) if atom.symbol != "Ag" and i not in shell_indices and int(parent.arrays["lammps_id"][i]) not in changed_ids and i not in pair_indices),
            "parent_source_method": parent.info.get("source_method"),
        }
        pair_rows.append(row)

        # Per-structure V16 errors and per-atom vector changes; no validation gate is assigned here.
        for tag, item, pred in (("plus", pos_dft, pos), ("minus", neg_dft, neg)):
            force_error = pred["forces"] - item["forces"]
            element_rows = {}
            for element in sorted(set(item["output_atoms"].get_chemical_symbols())):
                inds = [i for i, atom in enumerate(item["output_atoms"]) if atom.symbol == element]
                element_rows[element] = vector_rmse(force_error[inds])
            inference_rows.append({
                "pair": pair_name,
                "label": item["label"],
                "interface": pair_spec["interface"],
                "mode": pair_spec["mode"],
                "sign": tag,
                "atoms": len(item["output_atoms"]),
                "DFT_native_energy_eV": item["native_energy_eV"],
                "DFT_free_energy_eV": item["free_energy_eV"],
                "V16_energy_eV": pred["energy_eV"],
                "V16_minus_DFT_native_eV": pred["energy_eV"] - item["native_energy_eV"],
                "V16_minus_DFT_free_eV": None if item["free_energy_eV"] is None else pred["energy_eV"] - item["free_energy_eV"],
                "force_vector_RMSE_eV_A": vector_rmse(force_error),
                "force_vector_max_eV_A": float(np.linalg.norm(force_error, axis=1).max()),
                "Ag_force_RMSE_eV_A": element_rows.get("Ag"),
                "nonAg_force_RMSE_eV_A": vector_rmse(force_error[[i for i, atom in enumerate(item["output_atoms"]) if atom.symbol != "Ag"]]),
                "separating_force_DFT_eV_A": pair_force_rows[tag]["separating_force_eV_A"],
                "separating_force_V16_eV_A": pair_atoms[f"model_{tag}"]["separating_force_eV_A"],
                "separating_force_error_V16_minus_DFT_eV_A": pair_atoms[f"model_{tag}"]["separating_force_eV_A"] - pair_force_rows[tag]["separating_force_eV_A"],
                "marked_pair_transverse_force_DFT_eV_A": pair_force_rows[tag]["transverse_force_magnitude_eV_A"],
                "marked_pair_transverse_force_V16_eV_A": pair_atoms[f"model_{tag}"]["transverse_force_magnitude_eV_A"],
            })

            # Each atom is aligned by persistent ID to the exact source parent.
            output_atoms = item["output_atoms"]
            output_map = id_map(output_atoms)
            parent_model_force = parent_pred["forces"]
            for atom_id, index in mapping.items():
                out_i = output_map[atom_id]
                dft_force = item["forces"][out_i]
                pred_force = pred["forces"][out_i]
                parent_force = parent_ref_f[index]
                parent_model = parent_model_force[index]
                atom_rows.append({
                    "pair": pair_name,
                    "label": item["label"],
                    "interface": pair_spec["interface"],
                    "mode": pair_spec["mode"],
                    "sign": tag,
                    "atom_id": atom_id,
                    "element": parent[index].symbol,
                    "response_group": groups[atom_id],
                    "DFT_fx_eV_A": dft_force[0], "DFT_fy_eV_A": dft_force[1], "DFT_fz_eV_A": dft_force[2],
                    "V16_fx_eV_A": pred_force[0], "V16_fy_eV_A": pred_force[1], "V16_fz_eV_A": pred_force[2],
                    "V16_minus_DFT_fx_eV_A": pred_force[0] - dft_force[0],
                    "V16_minus_DFT_fy_eV_A": pred_force[1] - dft_force[1],
                    "V16_minus_DFT_fz_eV_A": pred_force[2] - dft_force[2],
                    "DFT_delta_force_from_parent_norm_eV_A": float(np.linalg.norm(dft_force - parent_force)),
                    "V16_delta_force_from_parent_norm_eV_A": float(np.linalg.norm(pred_force - parent_model)),
                    "V16_DFT_force_error_norm_eV_A": float(np.linalg.norm(pred_force - dft_force)),
                })

            for group in sorted(set(groups.values())):
                ids = [atom_id for atom_id, group_name in groups.items() if group_name == group]
                indices_parent = [mapping[atom_id] for atom_id in ids]
                indices_out = [output_map[atom_id] for atom_id in ids]
                dft_delta = np.asarray([item["forces"][i] - parent_ref_f[j] for i, j in zip(indices_out, indices_parent)])
                v16_delta = np.asarray([pred["forces"][i] - parent_pred["forces"][j] for i, j in zip(indices_out, indices_parent)])
                force_error = np.asarray([pred["forces"][i] - item["forces"][i] for i in indices_out])
                group_rows.append({
                    "pair": pair_name, "label": item["label"], "interface": pair_spec["interface"],
                    "mode": pair_spec["mode"], "sign": tag, "response_group": group,
                    "atom_count": len(ids),
                    "DFT_delta_from_parent_mean_norm_eV_A": force_stats(dft_delta)["mean_force_change_norm_eV_A"],
                    "DFT_delta_from_parent_RMS_norm_eV_A": force_stats(dft_delta)["rms_force_change_norm_eV_A"],
                    "V16_delta_from_parent_mean_norm_eV_A": force_stats(v16_delta)["mean_force_change_norm_eV_A"],
                    "V16_delta_from_parent_RMS_norm_eV_A": force_stats(v16_delta)["rms_force_change_norm_eV_A"],
                    "V16_force_error_vector_RMSE_eV_A": vector_rmse(force_error),
                })

    # Paired +/- force-vector response and directional response error.
    directional_rows = []
    for pair_name, pair_spec in PAIR_LABELS.items():
        pos_label, neg_label = pair_spec["positive"], pair_spec["negative"]
        p = structures[pos_label]["dft"]; n = structures[neg_label]["dft"]
        p_pred, n_pred = predictions[pos_label], predictions[neg_label]
        pmap, nmap = id_map(p["output_atoms"]), id_map(n["output_atoms"])
        for atom_id in sorted(pmap):
            dft_response = p["forces"][pmap[atom_id]] - n["forces"][nmap[atom_id]]
            model_response = p_pred["forces"][pmap[atom_id]] - n_pred["forces"][nmap[atom_id]]
            directional_rows.append({
                "pair": pair_name,
                "interface": pair_spec["interface"],
                "mode": pair_spec["mode"],
                "atom_id": atom_id,
                "element": p["output_atoms"][pmap[atom_id]].symbol,
                "response_group": next(row["response_group"] for row in atom_rows if row["pair"] == pair_name and row["label"] == pos_label and row["atom_id"] == atom_id),
                "DFT_plus_minus_delta_fx_eV_A": dft_response[0],
                "DFT_plus_minus_delta_fy_eV_A": dft_response[1],
                "DFT_plus_minus_delta_fz_eV_A": dft_response[2],
                "V16_plus_minus_delta_fx_eV_A": model_response[0],
                "V16_plus_minus_delta_fy_eV_A": model_response[1],
                "V16_plus_minus_delta_fz_eV_A": model_response[2],
                "V16_minus_DFT_response_error_norm_eV_A": float(np.linalg.norm(model_response - dft_response)),
            })

    audit = {
        "task": "v17_paired_response_analysis",
        "scope": "eight matched V17 diagnostic training acquisitions, their three unique V14 training parents, and frozen V16 CPU inference; no holdout/test labels",
        "all_eight_archives_verified": True,
        "archive_verification_checks": {label: structures[label]["dft"]["verification"]["checks"] for spec in PAIR_LABELS.values() for label in (spec["positive"], spec["negative"])},
        "archive_records": {
            label: {
                "owner": structures[label]["dft"]["owner"],
                "input_sha256": structures[label]["dft"]["record"]["input_sha256"],
                "verification_sha256": structures[label]["dft"]["archive_files_sha256"]["verification.json"],
                "archive_files_sha256": structures[label]["dft"]["archive_files_sha256"],
                "method": structures[label]["dft"]["method"],
                "native_energy_eV": structures[label]["dft"]["native_energy_eV"],
                "free_energy_eV": structures[label]["dft"]["free_energy_eV"],
                "mpi_ranks": structures[label]["dft"]["summary"]["mpi_ranks"],
                "scf_converged": structures[label]["dft"]["summary"]["scf_converged"],
            }
            for spec in PAIR_LABELS.values() for label in (spec["positive"], spec["negative"])
        },
        "source_manifest_hashes": {
            "A_input_manifest_sha256": digest(a_base / "input_manifest.json"),
            "A_archive_manifest_sha256": digest(a_base / "archive_manifest.json"),
            "B_input_manifest_sha256": digest(b_base / "input_manifest.json"),
            "B_archive_manifest_sha256": digest(b_base / "archive_manifest.json"),
        },
        "metric_conventions": {
            "force_vector_RMSE": "sqrt(mean_i(sum_xyz((F_V16-F_DFT)^2))) in eV/A",
            "generalized_force": "sum over uniformly displaced atoms of F_i dot displacement_unit_i; plus/minus directions use the same positive axis",
            "central_secant": "(E_plus-E_minus)/(q_plus+q_minus); native extrapolated energy and free energy are evaluated separately",
            "separating_force": "(F_X-F_Ag) dot unit(r_Ag_to_X), using the direct Cartesian cluster vector, matching the V16 evaluator convention",
            "framework_neighbor_shell": "non-Ag, non-moved parent atoms within 3.5 A MIC of any moved atom",
        },
        "model": {
            "path": str(model.relative_to(repo)),
            "sha256": model_hash,
            "selection_record_path": str(selection_path.relative_to(repo)),
            "selection_record_sha256": selection_hash,
            "selected_epoch": selection["selected_epoch"],
            "selection_uses_only_training_validation": selection["selection_uses_only_training_validation"],
            "blind_labels_used_for_selection": selection["blind_labels_used_for_selection"],
            "inference_device": "cpu",
            "torch_version": torch.__version__,
            "mace_torch_version": __import__("importlib.metadata", fromlist=["version"]).version("mace-torch"),
        },
        "parent_dataset": {"path": str(v14_path.relative_to(repo)), "sha256": digest(v14_path), "read_frames_zero_based": sorted({int(r["source_parent"]["frame_index_zero_based"]) for m in (a_manifest, b_manifest) for r in m["records"]})},
        "energy_reference_caveat": "New +/- references record native extrapolated energy and free energy separately at FermiDirac 0.1 eV. Parent V14 REF_energy metadata says GPAW 26.7.0 PW-PBE 500 eV Gamma but does not declare smearing or native/free convention; parent-relative energy offsets and finite-difference-vs-parent force comparisons therefore remain method/convention-qualified.",
        "finite_difference_caveat": "0.04 A framework and 0.15 A collective transverse Ag central differences are finite-interval secants, not infinitesimal derivatives. Compare native and free slopes separately; GPAW finite-smearing forces are convention-sensitive.",
        "pair_audits": audit_rows,
        "predictions": {name: {"energy_eV": row["energy_eV"], "forces_hash": sha256(np.asarray(row["forces"], dtype="<f8").tobytes()).hexdigest()} for name, row in predictions.items()},
        "counts": {"pairs": len(PAIR_LABELS), "diagnostic_labels": 8, "unique_parent_controls": len(parents), "per_structure_force_rows": len(inference_rows), "per_atom_rows": len(atom_rows)},
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "analysis.json").write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")
    for filename, rows in (("paired_response.csv", pair_rows), ("per_structure_inference.csv", inference_rows), ("force_response_by_group.csv", group_rows), ("per_atom_force_vectors.csv", atom_rows), ("directional_force_response.csv", directional_rows)):
        if not rows:
            continue
        with (out / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
    print(json.dumps({"pairs": len(pair_rows), "inference_rows": len(inference_rows), "force_rows": len(atom_rows), "model_sha256": model_hash}, indent=2))


if __name__ == "__main__":
    main()
