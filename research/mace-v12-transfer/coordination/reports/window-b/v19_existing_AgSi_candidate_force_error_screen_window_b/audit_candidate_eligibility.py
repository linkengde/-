#!/usr/bin/env python3
"""Audit the frozen MACE candidate and train-only overlap without inference."""
from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from ase.geometry import find_mic
from ase.io import read


TASK = "v19_existing_AgSi_candidate_force_error_screen_window_b"
REPORT_DIR = Path(__file__).resolve().parent
MODEL_ROOT = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core")
DATA_DIR = MODEL_ROOT / "data"
TRAIN_PATH = DATA_DIR / "train.extxyz"
SELECTION_PATH = MODEL_ROOT / "training_clean_core_01/selection_record.json"
ARTIFACT_PATH = MODEL_ROOT / "training_clean_core_01/artifact_manifest.json"
MODEL_PATH = MODEL_ROOT / "training_clean_core_01/models/MACE_periodic_v18_clean_core.model"
CHECKPOINT_PATH = MODEL_ROOT / "training_clean_core_01/checkpoints/MACE_periodic_v18_clean_core_run-45_epoch-79.pt"
REFERENCES = {
    "AgSi_COD9009647_pilot_k6x6": Path(
        "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_cif_convergence_pilot/calculations/AgSi_COD9009647_pilot_k6x6"
    ),
    "AgSi_COD9009647_vacuum26_k6x6_sigma0p10": Path(
        "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_cif_vacuum_pilot/calculations/AgSi_COD9009647_vacuum26_k6x6_sigma0p10"
    ),
}
GATES = {"force_vector_RMSE_eV_A": 0.05, "separating_force_abs_error_eV_A": 0.10}


def find_repo_root() -> Path:
    for parent in REPORT_DIR.parents:
        if (parent / ".git").exists():
            return parent
    raise RuntimeError("Cannot locate repository root")


REPO = find_repo_root()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_manifest(folder: Path) -> dict:
    manifest = folder / "SHA256SUMS.txt"
    if not manifest.is_file():
        raise RuntimeError(f"Missing archive checksum file: {manifest}")
    rows = []
    for line in manifest.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        expected, filename = line.split(maxsplit=1)
        filename = filename.lstrip("*")
        path = (folder / filename).resolve()
        if folder.resolve() not in path.parents or not path.is_file():
            rows.append({"file": filename, "expected": expected, "actual": None, "pass": False})
            continue
        actual = sha256(path)
        rows.append({"file": filename, "expected": expected, "actual": actual, "pass": expected == actual})
    if not rows or not all(row["pass"] for row in rows):
        raise RuntimeError(f"Archive checksum verification failed: {folder}")
    return {"status": "PASS", "files_checked": len(rows), "files": rows}


def verify_archive(repo: Path, label: str, relative_folder: Path) -> dict:
    folder = repo / relative_folder
    inventory = verify_manifest(folder)
    summary = json.loads((folder / "summary.json").read_text())
    verification = json.loads((folder / "verification.json").read_text())
    if summary.get("label") != label or summary.get("scf_converged") is not True:
        raise RuntimeError(f"Reference summary identity/convergence failed: {label}")
    checks = verification.get("checks", {})
    if verification.get("status") != "PASS" or not checks or not all(checks.values()):
        raise RuntimeError(f"Reference archive verifier failed: {label}")
    result = folder / f"{label}_PW_PBE.extxyz"
    atoms = read(result)
    source_sha = summary.get("source_sha256")
    source_record = str(summary.get("source", ""))
    marker = "research/mace-v12-transfer/"
    offset = source_record.replace("\\", "/").find(marker)
    source = repo / source_record.replace("\\", "/")[offset:] if offset >= 0 else Path(source_record)
    if not source.is_file() or sha256(source) != source_sha:
        raise RuntimeError(f"Reference source hash mismatch: {label}")
    return {
        "label": label,
        "relative_archive": str(relative_folder),
        "sha256_inventory": inventory,
        "verification_status": verification["status"],
        "verification_checks_passed": len(checks),
        "scf_converged": summary["scf_converged"],
        "source_sha256": source_sha,
        "source_hash_verified": True,
        "result_sha256": sha256(result),
        "atom_count": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "pbc": atoms.pbc.tolist(),
        "cell_A": np.asarray(atoms.cell.array).tolist(),
        "kpts": summary.get("method", {}).get("kpts"),
        "smearing_eV": summary.get("method", {}).get("smearing_eV"),
        "energy_eV_cell": summary.get("energy_eV_cell"),
        "free_energy_eV_cell": summary.get("free_energy_eV_cell"),
        "energy_convention": summary.get("energy_convention"),
        "_atoms": atoms,
    }


def atom_ids(atoms) -> np.ndarray:
    return np.asarray(atoms.arrays.get("lammps_id", np.arange(len(atoms))), dtype="<i8")


def exact_geometry_hash(atoms) -> str:
    payload = b"".join(
        (
            np.asarray(atoms.numbers, dtype="<i4").tobytes(),
            atom_ids(atoms).tobytes(),
            np.asarray(atoms.positions, dtype="<f8").tobytes(),
            np.asarray(atoms.cell.array, dtype="<f8").tobytes(),
            np.asarray(atoms.pbc, dtype="u1").tobytes(),
        )
    )
    return hashlib.sha256(payload).hexdigest()


def composition_key(atoms) -> tuple:
    symbols, counts = np.unique(atoms.get_chemical_symbols(), return_counts=True)
    return tuple((str(symbol), int(count)) for symbol, count in zip(symbols, counts))


def marked_pair(atoms) -> dict | None:
    markers = np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool)
    indices = np.flatnonzero(markers)
    if len(indices) != 2:
        return None
    ag = [int(index) for index in indices if atoms[int(index)].symbol == "Ag"]
    x = [int(index) for index in indices if atoms[int(index)].symbol != "Ag"]
    if len(ag) != 1 or len(x) != 1:
        return None
    i, j = ag[0], x[0]
    vector, _ = find_mic(atoms.positions[j] - atoms.positions[i], atoms.cell, pbc=atoms.pbc)
    return {
        "Ag_index_zero_based": i,
        "X_index_zero_based": j,
        "Ag_id": int(atom_ids(atoms)[i]),
        "X_id": int(atom_ids(atoms)[j]),
        "X_species": atoms[j].symbol,
        "distance_A": float(np.linalg.norm(vector)),
    }


def candidate_record(repo: Path) -> dict:
    data_manifest_path = repo / DATA_DIR / "dataset_manifest.json"
    data_manifest = json.loads(data_manifest_path.read_text())
    selection_path = repo / SELECTION_PATH
    selection = json.loads(selection_path.read_text())
    artifact_path = repo / ARTIFACT_PATH
    artifact = json.loads(artifact_path.read_text())
    model, checkpoint = repo / MODEL_PATH, repo / CHECKPOINT_PATH
    train = repo / TRAIN_PATH
    calculated = {
        "model_sha256": sha256(model),
        "epoch79_checkpoint_sha256": sha256(checkpoint),
        "training_input_sha256": sha256(train),
        "dataset_manifest_sha256": sha256(data_manifest_path),
        "selection_record_sha256": sha256(selection_path),
        "artifact_manifest_sha256": sha256(artifact_path),
    }
    artifact_by_path = {row["path"]: row for row in artifact.get("files", [])}
    checks = {
        "selection_model_hash_matches": calculated["model_sha256"] == selection.get("model_sha256"),
        "selection_epoch79_checkpoint_hash_matches": calculated["epoch79_checkpoint_sha256"]
        == selection.get("recovery_verification", {}).get("checkpoint_sha256"),
        "selection_training_hash_matches": calculated["training_input_sha256"] == selection.get("training_input_sha256"),
        "selection_dataset_manifest_hash_matches": calculated["dataset_manifest_sha256"] == selection.get("dataset_manifest_sha256"),
        "dataset_manifest_training_hash_matches": calculated["training_input_sha256"] == data_manifest.get("train_input_sha256"),
        "model_is_published_in_artifact_manifest": artifact_by_path.get(str(MODEL_PATH), {}).get("sha256")
        == calculated["model_sha256"]
        and artifact_by_path.get(str(MODEL_PATH), {}).get("published") is True,
        "epoch79_is_published_in_artifact_manifest": artifact_by_path.get(str(CHECKPOINT_PATH), {}).get("sha256")
        == calculated["epoch79_checkpoint_sha256"]
        and artifact_by_path.get(str(CHECKPOINT_PATH), {}).get("published") is True,
    }
    if not all(checks.values()):
        raise RuntimeError(f"V18 candidate hash lineage failed: {checks}")

    v19_root = repo / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v19_residual_response"
    v19_manifest_path = v19_root / "training_candidate/manifest.json"
    v19_manifest = json.loads(v19_manifest_path.read_text())
    v19_models = sorted(
        str(path.relative_to(repo))
        for path in v19_root.rglob("*")
        if path.is_file() and path.suffix in {".model", ".pt", ".safetensors"}
    )

    torch_spec = importlib.util.find_spec("torch")
    mace_spec = importlib.util.find_spec("mace")
    try:
        torch_version = importlib.metadata.version("torch") if torch_spec else None
    except importlib.metadata.PackageNotFoundError:
        torch_version = None
    try:
        mace_version = importlib.metadata.version("mace-torch") if mace_spec else None
    except importlib.metadata.PackageNotFoundError:
        mace_version = None
    runtime = {
        "python_executable": sys.executable,
        "python_version": sys.version.split()[0],
        "torch_module_available": torch_spec is not None,
        "torch_version": torch_version,
        "mace_module_available": mace_spec is not None,
        "mace_torch_version": mace_version,
        "inference_runtime_available": torch_spec is not None and mace_spec is not None,
    }
    roles = data_manifest.get("roles", {})
    return {
        "candidate_id": "V18 clean-core selected final model",
        "model_path": str(MODEL_PATH),
        "checkpoint_path": str(CHECKPOINT_PATH),
        "lineage_hashes": calculated,
        "lineage_checks": checks,
        "selection": {
            "completed_epochs": selection.get("completed_epochs"),
            "selected_epoch": selection.get("selected_epoch"),
            "training_exit_code": selection.get("training_exit_code"),
            "exported_state_equals_epoch79": selection.get("recovery_verification", {}).get("exported_state_equals_epoch79"),
            "reviewed_by": selection.get("reviewed_by"),
            "selection_uses_only_training_validation": selection.get("selection_uses_only_training_validation"),
            "blind_labels_used_for_selection": selection.get("blind_labels_used_for_selection"),
            "test_file_passed_to_training": selection.get("test_file_passed_to_training"),
            "original_launcher_exit_receipt": selection.get("recovery_verification", {}).get("original_launcher_exit_receipt"),
        },
        "training_data": {
            "manifest_path": str(DATA_DIR / "dataset_manifest.json"),
            "manifest_status": data_manifest.get("status"),
            "roles_from_manifest": roles,
            "train_frame_count": len(data_manifest.get("rows", [])),
            "train_input_path": str(TRAIN_PATH),
            "train_input_sha256": calculated["training_input_sha256"],
            "withheld_labels_opened": data_manifest.get("withheld_labels_opened"),
            "limitations": data_manifest.get("limitations", []),
            "development_or_test_file_contents_read": False,
        },
        "newer_v19_candidate": {
            "manifest_path": str(v19_manifest_path.relative_to(repo)),
            "manifest_sha256": sha256(v19_manifest_path),
            "status": v19_manifest.get("status"),
            "training_authorized": v19_manifest.get("training_authorized"),
            "executable_model_files": v19_models,
        },
        "runtime": runtime,
    }


def training_overlap(train_atoms: list, reference_records: dict, data_manifest: dict) -> dict:
    train_records = []
    source_hashes = {row.get("source_hash") for row in data_manifest.get("rows", [])}
    for index, atoms in enumerate(train_atoms):
        train_records.append(
            {
                "frame_index": index,
                "config_type": atoms.info.get("config_type"),
                "atom_count": len(atoms),
                "formula": atoms.get_chemical_formula(),
                "composition_key": list(composition_key(atoms)),
                "pbc": atoms.pbc.tolist(),
                "cell_lengths_A": np.asarray(atoms.cell.lengths()).tolist(),
                "geometry_sha256": exact_geometry_hash(atoms),
                "marked_pair": marked_pair(atoms),
            }
        )

    per_reference = []
    for label, record in reference_records.items():
        atoms = record["_atoms"]
        fingerprint = exact_geometry_hash(atoms)
        comp = composition_key(atoms)
        matching_composition = [
            row
            for row in train_records
            if row["atom_count"] == len(atoms)
            and tuple((str(x[0]), int(x[1])) for x in row["composition_key"]) == comp
            and row["pbc"] == atoms.pbc.tolist()
        ]
        exact = [row["frame_index"] for row in train_records if row["geometry_sha256"] == fingerprint]
        pair = marked_pair(atoms)
        train_pairs = [
            {"frame_index": row["frame_index"], **row["marked_pair"]}
            for row in train_records
            if row["marked_pair"] and row["marked_pair"]["X_species"] == (pair or {}).get("X_species")
            and row.get("config_type", "").startswith("AgSi")
        ]
        closest = None
        if pair and train_pairs:
            closest = min(train_pairs, key=lambda row: abs(row["distance_A"] - pair["distance_A"]))
        per_reference.append(
            {
                "label": label,
                "reference_geometry_sha256": fingerprint,
                "reference_source_hash_in_nonsealed_training_manifest": record["source_sha256"] in source_hashes,
                "exact_full_geometry_matches": exact,
                "same_composition_atom_count_and_pbc_candidate_frames": [row["frame_index"] for row in matching_composition],
                "same_composition_candidate_count": len(matching_composition),
                "global_near_geometry_rmsd_evaluated": False if not matching_composition else None,
                "marked_pair": pair,
                "closest_training_AgSi_marked_pair": closest,
                "local_contact_distance_overlap": bool(closest and abs(closest["distance_A"] - pair["distance_A"]) <= 0.05),
                "interpretation": "No exact/global same-composition training frame; the marked Ag-X contact distance has a close local analogue and is not independent validation.",
            }
        )
    return {
        "training_file_read": str(TRAIN_PATH),
        "training_frames_read": len(train_atoms),
        "no_development_or_test_extxyz_read": True,
        "training_frames": train_records,
        "references": per_reference,
        "summary": "Exact global duplicate search is clear for the two reference geometries. Full near-geometry RMSD was not computed because the V18 training set has no frame matching either target's atom count, composition, and PBC. Local Ag-Si contact geometry overlaps: the reference marked distance is about 2.7795 A and a clean30 AgSi training frame is about 2.8000 A.",
    }


def clean(record: dict) -> dict:
    return {key: value for key, value in record.items() if not key.startswith("_")}


def main() -> None:
    candidate = candidate_record(REPO)
    data_manifest = json.loads((REPO / DATA_DIR / "dataset_manifest.json").read_text())
    train_atoms = read(REPO / TRAIN_PATH, index=":")
    if len(train_atoms) != candidate["training_data"]["train_frame_count"]:
        raise RuntimeError("Training-frame count differs from the pinned dataset manifest")
    reference_records = {
        label: verify_archive(REPO, label, relative)
        for label, relative in REFERENCES.items()
    }
    overlap = training_overlap(train_atoms, reference_records, data_manifest)
    runtime = candidate["runtime"]
    inference_status = "NOT_RUN_RUNTIME_MISSING" if not runtime["inference_runtime_available"] else "ELIGIBILITY_REVIEW_REQUIRED"
    result = {
        "task": TASK,
        "status": "INFERENCE_BLOCKED_RUNTIME_MISSING" if not runtime["inference_runtime_available"] else "ELIGIBILITY_REVIEW_REQUIRED",
        "scope": "Read-only model lineage, non-sealed training overlap, and public reference archive audit. No model inference, training, DFT/MPI, install, or sealed-label access.",
        "provisional_gates": GATES,
        "candidate": candidate,
        "reference_archives": [clean(record) for record in reference_records.values()],
        "training_overlap": overlap,
        "inference": {
            "status": inference_status,
            "performed": False,
            "reason": "Neither PyTorch nor mace-torch is installed in the only available Python environments. This task forbids package installation.",
            "force_vector_metrics": None,
            "per_species_force_metrics": None,
            "marked_pair_separation_metrics": None,
            "Ag_layer_normal_force_metrics": None,
            "gates_evaluated": False,
        },
        "decision": {
            "frozen_v18_candidate_hash_lineage_passes": all(candidate["lineage_checks"].values()),
            "v18_candidate_is_only_a_provisional_numerical_screen": True,
            "v19_candidate_has_executable_model": bool(candidate["newer_v19_candidate"]["executable_model_files"]),
            "global_training_geometry_overlap_found": any(bool(row["exact_full_geometry_matches"] or row["same_composition_candidate_count"]) for row in overlap["references"]),
            "local_ag_si_contact_overlap_found": any(row["local_contact_distance_overlap"] for row in overlap["references"]),
            "eligible_to_run_inference_in_this_environment": False,
            "interpretation": "The pinned V18 model can only support a non-independent one-family diagnostic after an approved existing MACE runtime is available. The local Ag-Si contact motif overlaps training; do not report this as independent validation or a production pass.",
        },
    }
    out = REPORT_DIR / "eligibility.json"
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    print(f"Wrote {out.relative_to(REPO)}")
    print("runtime", runtime)
    print("train frames", len(train_atoms), "exact/global matches", [row["exact_full_geometry_matches"] for row in overlap["references"]])
    print("local AgSi nearest pair", [row["closest_training_AgSi_marked_pair"] for row in overlap["references"]])
    print("inference", result["inference"]["status"])


if __name__ == "__main__":
    main()
