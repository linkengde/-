#!/usr/bin/env python3
"""Verify the pinned V18 candidate and score only the two public AgSi CIF references."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import time

# MACE imports matplotlib indirectly; keep its cache inside writable temporary storage.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mace-cpu-mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/mace-cpu-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import numpy as np
import torch
from ase.geometry import find_mic
from ase.io import read
from mace.calculators import MACECalculator


TASK = "v19_existing_AgSi_candidate_force_error_screen_install_and_infer_window_b"
MODEL_ROOT = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core")
MODEL_PATH = MODEL_ROOT / "training_clean_core_01/models/MACE_periodic_v18_clean_core.model"
CHECKPOINT_PATH = MODEL_ROOT / "training_clean_core_01/checkpoints/MACE_periodic_v18_clean_core_run-45_epoch-79.pt"
TRAIN_PATH = MODEL_ROOT / "data/train.extxyz"
DATASET_MANIFEST = MODEL_ROOT / "data/dataset_manifest.json"
SELECTION_RECORD = MODEL_ROOT / "training_clean_core_01/selection_record.json"
ARTIFACT_MANIFEST = MODEL_ROOT / "training_clean_core_01/artifact_manifest.json"
PREVIOUS_AUDIT = Path(
    "research/mace-v12-transfer/coordination/reports/window-b/"
    "v19_existing_AgSi_candidate_force_error_screen_window_b/eligibility.json"
)
REFERENCES = {
    "AgSi_COD9009647_pilot_k6x6": Path(
        "research/mace-v12-transfer/periodic_interface_v4/"
        "pbe_interface_v19_cif_convergence_pilot/calculations/AgSi_COD9009647_pilot_k6x6"
    ),
    "AgSi_COD9009647_vacuum26_k6x6_sigma0p10": Path(
        "research/mace-v12-transfer/periodic_interface_v4/"
        "pbe_interface_v19_cif_vacuum_pilot/calculations/"
        "AgSi_COD9009647_vacuum26_k6x6_sigma0p10"
    ),
}
GATES = {
    "force_vector_RMSE_eV_A": 0.05,
    "separating_force_abs_error_eV_A": 0.10,
}


def find_repo_root(start: Path) -> Path:
    for parent in (start, *start.parents):
        if (parent / ".git").exists():
            return parent
    raise RuntimeError("Cannot locate repository root")


REPO = find_repo_root(Path(__file__).resolve())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repo_file(relative: Path) -> Path:
    path = (REPO / relative).resolve()
    if REPO.resolve() not in path.parents or not path.is_file():
        raise RuntimeError(f"Missing or out-of-repository input: {relative}")
    return path


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


def composition_key(atoms) -> tuple[tuple[str, int], ...]:
    symbols, counts = np.unique(atoms.get_chemical_symbols(), return_counts=True)
    return tuple((str(symbol), int(count)) for symbol, count in zip(symbols, counts))


def marked_pair(atoms) -> dict:
    markers = np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool)
    indices = np.flatnonzero(markers)
    if len(indices) != 2:
        raise RuntimeError("Expected exactly two central_pair atoms in a public AgSi reference")
    ag = [int(index) for index in indices if atoms[int(index)].symbol == "Ag"]
    x = [int(index) for index in indices if atoms[int(index)].symbol != "Ag"]
    if len(ag) != 1 or len(x) != 1:
        raise RuntimeError("Expected one marked Ag and one marked contact atom")
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


def verify_model_lineage() -> tuple[dict, dict]:
    selection_path = repo_file(SELECTION_RECORD)
    artifact_path = repo_file(ARTIFACT_MANIFEST)
    dataset_path = repo_file(DATASET_MANIFEST)
    model_path = repo_file(MODEL_PATH)
    checkpoint_path = repo_file(CHECKPOINT_PATH)
    train_path = repo_file(TRAIN_PATH)
    selection = json.loads(selection_path.read_text())
    artifact = json.loads(artifact_path.read_text())
    dataset = json.loads(dataset_path.read_text())

    hashes = {
        "model_sha256": sha256(model_path),
        "epoch79_checkpoint_sha256": sha256(checkpoint_path),
        "training_input_sha256": sha256(train_path),
        "dataset_manifest_sha256": sha256(dataset_path),
        "selection_record_sha256": sha256(selection_path),
        "artifact_manifest_sha256": sha256(artifact_path),
        "previous_eligibility_audit_sha256": sha256(repo_file(PREVIOUS_AUDIT)),
    }
    artifact_by_path = {row["path"]: row for row in artifact.get("files", [])}
    checks = {
        "selected_model_hash_matches": hashes["model_sha256"] == selection.get("model_sha256"),
        "selected_epoch79_checkpoint_hash_matches": hashes["epoch79_checkpoint_sha256"]
        == selection.get("recovery_verification", {}).get("checkpoint_sha256"),
        "training_input_hash_matches_selection": hashes["training_input_sha256"]
        == selection.get("training_input_sha256"),
        "dataset_manifest_hash_matches_selection": hashes["dataset_manifest_sha256"]
        == selection.get("dataset_manifest_sha256"),
        "dataset_manifest_points_to_same_training_input": hashes["training_input_sha256"]
        == dataset.get("train_input_sha256"),
        "model_recorded_as_published": artifact_by_path.get(str(MODEL_PATH), {}).get("sha256")
        == hashes["model_sha256"]
        and artifact_by_path.get(str(MODEL_PATH), {}).get("published") is True,
        "checkpoint_recorded_as_published": artifact_by_path.get(str(CHECKPOINT_PATH), {}).get("sha256")
        == hashes["epoch79_checkpoint_sha256"]
        and artifact_by_path.get(str(CHECKPOINT_PATH), {}).get("published") is True,
        "selected_epoch_is_79": selection.get("selected_epoch") == 79,
        "export_matches_epoch79": selection.get("recovery_verification", {}).get("exported_state_equals_epoch79")
        is True,
        "blind_labels_not_used_for_selection": selection.get("blind_labels_used_for_selection") is False,
        "test_file_not_passed_to_training": selection.get("test_file_passed_to_training") is False,
        "nonsealed_training_manifest": dataset.get("withheld_labels_opened") is False,
    }
    if not all(checks.values()):
        raise RuntimeError(f"Pinned candidate lineage check failed: {checks}")
    return {"hashes": hashes, "checks": checks, "selection": selection, "dataset": dataset}, artifact


def verify_archive(label: str, folder_relative: Path) -> dict:
    folder = repo_file(folder_relative / "summary.json").parent
    checksum_path = repo_file(folder_relative / "SHA256SUMS.txt")
    rows = []
    for line in checksum_path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        expected, filename = line.split(maxsplit=1)
        filename = filename.lstrip("*")
        member = (folder / filename).resolve()
        if folder.resolve() not in member.parents or not member.is_file():
            raise RuntimeError(f"Invalid archive checksum member: {label}/{filename}")
        actual = sha256(member)
        rows.append({"file": filename, "expected_sha256": expected, "actual_sha256": actual, "pass": expected == actual})
    if not rows or not all(row["pass"] for row in rows):
        raise RuntimeError(f"Archive SHA256SUMS verification failed for {label}")

    summary = json.loads((folder / "summary.json").read_text())
    verification = json.loads((folder / "verification.json").read_text())
    if summary.get("label") != label or summary.get("scf_converged") is not True:
        raise RuntimeError(f"Archive label or convergence check failed for {label}")
    checks = verification.get("checks", {})
    if verification.get("status") != "PASS" or not checks or not all(checks.values()):
        raise RuntimeError(f"Recorded result verification failed for {label}")

    source_text = str(summary.get("source", "")).replace("\\", "/")
    marker = "research/mace-v12-transfer/"
    offset = source_text.find(marker)
    if offset < 0:
        raise RuntimeError(f"Cannot resolve public source structure for {label}")
    source_relative = Path(source_text[offset:])
    source_path = repo_file(source_relative)
    if sha256(source_path) != summary.get("source_sha256"):
        raise RuntimeError(f"Reference source checksum mismatch for {label}")

    result_relative = folder_relative / f"{label}_PW_PBE.extxyz"
    result_path = repo_file(result_relative)
    if verification.get("sha256", {}).get(result_path.name) != sha256(result_path):
        raise RuntimeError(f"Archived result hash differs from verification record for {label}")
    atoms = read(result_path, index=0)
    source_atoms = read(source_path, index=0)
    geometry_checks = {
        "atom_count": len(atoms) == len(source_atoms) == 52,
        "formula": atoms.get_chemical_formula() == source_atoms.get_chemical_formula() == "C16Ag4Si8Ti24",
        "atomic_numbers_and_order": np.array_equal(atoms.numbers, source_atoms.numbers),
        "atom_ids_and_order": np.array_equal(atom_ids(atoms), atom_ids(source_atoms)),
        "positions_exact": np.array_equal(atoms.positions, source_atoms.positions),
        "cell_exact": np.array_equal(atoms.cell.array, source_atoms.cell.array),
        "pbc_exact": np.array_equal(atoms.pbc, source_atoms.pbc),
        "fully_periodic_public_reference": atoms.pbc.tolist() == [True, True, True],
    }
    if not all(geometry_checks.values()):
        raise RuntimeError(f"Archived geometry / atom identity check failed for {label}: {geometry_checks}")
    if len(atoms) != int(summary.get("atoms", -1)) or atoms.get_chemical_formula() != summary.get("formula"):
        raise RuntimeError(f"Archive summary chemistry mismatch for {label}")

    return {
        "label": label,
        "archive_directory": str(folder_relative),
        "archive_checksum_file_sha256": sha256(checksum_path),
        "archive_files_checked": len(rows),
        "archive_checksum_checks": rows,
        "recorded_verification_status": verification["status"],
        "recorded_verification_checks": checks,
        "source_structure_path": str(source_relative),
        "source_structure_sha256": sha256(source_path),
        "archive_result_path": str(result_relative),
        "archive_result_sha256": sha256(result_path),
        "geometry_and_identity_checks": geometry_checks,
        "atom_count": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "pbc": atoms.pbc.tolist(),
        "cell_A": np.asarray(atoms.cell.array).tolist(),
        "energy_convention": summary.get("energy_convention"),
        "native_extrapolated_energy_eV_cell": float(summary["energy_eV_cell"]),
        "free_energy_eV_cell": float(summary["free_energy_eV_cell"]),
        "scf_converged": summary["scf_converged"],
        "method": summary.get("method"),
        "marked_pair": marked_pair(atoms),
        "_atoms": atoms,
    }


def verify_training_overlap(train_atoms: list, reference_records: dict, dataset: dict) -> dict:
    if len(train_atoms) != len(dataset.get("rows", [])):
        raise RuntimeError("Non-sealed training frame count differs from its dataset manifest")
    train_hashes = {row.get("source_hash") for row in dataset.get("rows", [])}
    training_records = []
    for index, atoms in enumerate(train_atoms):
        pair = None
        try:
            pair = marked_pair(atoms)
        except RuntimeError:
            pass
        training_records.append(
            {
                "frame_index_zero_based": index,
                "config_type": atoms.info.get("config_type"),
                "atom_count": len(atoms),
                "formula": atoms.get_chemical_formula(),
                "composition_key": composition_key(atoms),
                "pbc": atoms.pbc.tolist(),
                "geometry_sha256": exact_geometry_hash(atoms),
                "marked_pair": pair,
            }
        )

    references = []
    for label, record in reference_records.items():
        atoms = record["_atoms"]
        pair = record["marked_pair"]
        exact = [r["frame_index_zero_based"] for r in training_records if r["geometry_sha256"] == exact_geometry_hash(atoms)]
        matching = [
            r
            for r in training_records
            if r["atom_count"] == len(atoms)
            and r["composition_key"] == composition_key(atoms)
            and r["pbc"] == atoms.pbc.tolist()
        ]
        local_pairs = [
            {"frame_index_zero_based": r["frame_index_zero_based"], "config_type": r["config_type"], **r["marked_pair"]}
            for r in training_records
            if r["marked_pair"] is not None
            and r["marked_pair"]["X_species"] == pair["X_species"]
            and str(r.get("config_type") or "").startswith("AgSi")
        ]
        nearest = min(local_pairs, key=lambda row: abs(row["distance_A"] - pair["distance_A"])) if local_pairs else None
        references.append(
            {
                "label": label,
                "reference_geometry_sha256": exact_geometry_hash(atoms),
                "source_hash_appears_in_training_manifest": record["source_structure_sha256"] in train_hashes,
                "exact_full_geometry_training_matches": exact,
                "same_composition_count_pbc_training_frames": [r["frame_index_zero_based"] for r in matching],
                "near_geometry_rmsd_checked": False if not matching else None,
                "marked_pair": pair,
                "nearest_training_AgSi_pair": nearest,
                "local_contact_distance_overlap_within_0p05A": bool(
                    nearest is not None and abs(nearest["distance_A"] - pair["distance_A"]) <= 0.05
                ),
            }
        )
    if any(row["exact_full_geometry_training_matches"] or row["same_composition_count_pbc_training_frames"] for row in references):
        raise RuntimeError("Reference geometry overlaps the V18 training set globally; this task's eligibility premise changed")
    return {
        "training_file": str(TRAIN_PATH),
        "training_file_sha256": sha256(repo_file(TRAIN_PATH)),
        "training_frame_count": len(train_atoms),
        "development_or_test_structure_files_read": False,
        "withheld_labels_opened": False,
        "training_frames": training_records,
        "references": references,
        "interpretation": (
            "No exact or same-composition global training geometry match was found. The local Ag-Si contact distance "
            "overlaps the training coverage; these related CIF references support only a diagnostic screen."
        ),
    }


def signed_separation_force(forces: np.ndarray, atoms, pair: dict) -> float:
    i, j = pair["Ag_index_zero_based"], pair["X_index_zero_based"]
    vector, _ = find_mic(atoms.positions[j] - atoms.positions[i], atoms.cell, pbc=atoms.pbc)
    unit = vector / np.linalg.norm(vector)
    return float(np.dot(forces[j] - forces[i], unit))


def force_metrics(predicted: np.ndarray, reference: np.ndarray, atoms) -> dict:
    delta = predicted - reference
    vector_norms = np.linalg.norm(delta, axis=1)
    per_species = {}
    for species in sorted(set(atoms.get_chemical_symbols())):
        mask = np.asarray(atoms.get_chemical_symbols()) == species
        norms = vector_norms[mask]
        per_species[species] = {
            "atom_count": int(np.sum(mask)),
            "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(norms**2))),
            "maximum_atom_force_vector_error_eV_A": float(np.max(norms)),
            "mean_atom_force_vector_error_eV_A": float(np.mean(norms)),
        }

    pair = marked_pair(atoms)
    predicted_pair = signed_separation_force(predicted, atoms, pair)
    reference_pair = signed_separation_force(reference, atoms, pair)
    separation_error = predicted_pair - reference_pair
    ag_mask = np.asarray(atoms.get_chemical_symbols()) == "Ag"
    ag_fz_pred = float(np.sum(predicted[ag_mask, 2]))
    ag_fz_ref = float(np.sum(reference[ag_mask, 2]))
    ag_count = int(np.sum(ag_mask))
    return {
        "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
        "maximum_atom_force_vector_error_eV_A": float(np.max(vector_norms)),
        "mean_atom_force_vector_error_eV_A": float(np.mean(vector_norms)),
        "per_species": per_species,
        "marked_pair": pair,
        "marked_pair_separating_force_eV_A": {
            "definition": "dot(F_X - F_Ag, minimum_image_unit_vector_from_Ag_to_X)",
            "model_signed": predicted_pair,
            "dft_signed": reference_pair,
            "signed_error_model_minus_dft": separation_error,
            "absolute_error": abs(separation_error),
        },
        "all_Ag_layer_normal_force_eV_A": {
            "definition": "sum of F_z over all Ag atoms; positive is +z",
            "Ag_atom_count": ag_count,
            "model_net_signed": ag_fz_pred,
            "dft_net_signed": ag_fz_ref,
            "net_signed_error_model_minus_dft": ag_fz_pred - ag_fz_ref,
            "model_mean_per_Ag": ag_fz_pred / ag_count,
            "dft_mean_per_Ag": ag_fz_ref / ag_count,
            "mean_signed_error_model_minus_dft": (ag_fz_pred - ag_fz_ref) / ag_count,
        },
    }


def environment_record() -> dict:
    packages = {}
    for name in ("torch", "mace-torch", "ase", "numpy", "scipy", "e3nn", "matscipy"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "python_executable": sys.executable,
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": packages,
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "device": "cpu",
        "torch_num_threads": int(torch.get_num_threads()),
        "matplotlib_config_dir": os.environ["MPLCONFIGDIR"],
    }


def json_safe(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def main() -> None:
    global REPO
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO)
    parser.add_argument("--report-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    REPO = args.repo_root.resolve()
    report_dir = args.report_dir.resolve()
    allowed_root = (REPO / "research/mace-v12-transfer/coordination/reports/window-b").resolve()
    if allowed_root not in report_dir.parents:
        raise RuntimeError("Refusing to write outputs outside the registered Window B reports path")
    report_dir.mkdir(parents=True, exist_ok=True)

    if torch.cuda.is_available():
        raise RuntimeError("CPU-only inference task: unexpected CUDA availability")
    try:
        torch.set_num_threads(4)
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass

    lineage, _artifact = verify_model_lineage()
    dataset = lineage["dataset"]
    train_atoms = read(repo_file(TRAIN_PATH), index=":")
    reference_records = {label: verify_archive(label, relative) for label, relative in REFERENCES.items()}
    overlap = verify_training_overlap(train_atoms, reference_records, dataset)

    model_path = repo_file(MODEL_PATH)
    calculator = MACECalculator(
        model_paths=str(model_path),
        device="cpu",
        default_dtype="float64",
    )
    inference_rows = []
    start_all = time.perf_counter()
    smoke_pass = True
    for label, archive in reference_records.items():
        atoms = archive["_atoms"]
        # These archives store force labels under the explicit PW_PBE_forces array.
        dft_forces = np.asarray(atoms.arrays.get("PW_PBE_forces"), dtype=float).copy()
        atoms.calc = calculator
        started = time.perf_counter()
        predicted_energy = float(atoms.get_potential_energy())
        predicted_forces = np.asarray(atoms.get_forces(), dtype=float)
        elapsed = time.perf_counter() - started
        if predicted_forces.shape != (len(atoms), 3) or dft_forces.shape != (len(atoms), 3):
            raise RuntimeError(f"Unexpected force array shape for {label}")
        if not np.isfinite(predicted_energy) or not np.isfinite(predicted_forces).all() or not np.isfinite(dft_forces).all():
            raise RuntimeError(f"Nonfinite CPU inference or DFT label for {label}")
        smoke_pass = smoke_pass and predicted_forces.shape == (52, 3) and np.isfinite(predicted_forces).all()

        force_result = force_metrics(predicted_forces, dft_forces, atoms)
        native_dft = archive["native_extrapolated_energy_eV_cell"]
        free_dft = archive["free_energy_eV_cell"]
        inference_rows.append(
            {
                "label": label,
                "reference_archive_result_sha256": archive["archive_result_sha256"],
                "inference_device": "cpu",
                "elapsed_s": elapsed,
                "atom_count": len(atoms),
                "formula": atoms.get_chemical_formula(),
                "model_energy_eV_cell": predicted_energy,
                "energy_comparisons_separate_by_DFT_convention": {
                    "native_extrapolated": {
                        "dft_energy_eV_cell": native_dft,
                        "model_minus_dft_eV_cell": predicted_energy - native_dft,
                        "model_minus_dft_meV_atom": (predicted_energy - native_dft) * 1000.0 / len(atoms),
                    },
                    "free_energy": {
                        "dft_energy_eV_cell": free_dft,
                        "model_minus_dft_eV_cell": predicted_energy - free_dft,
                        "model_minus_dft_meV_atom": (predicted_energy - free_dft) * 1000.0 / len(atoms),
                    },
                    "energy_gate_applied": False,
                },
                "force_metrics": force_result,
                "provisional_gates": {
                    "force_vector_RMSE_eV_A": {
                        "limit": GATES["force_vector_RMSE_eV_A"],
                        "value": force_result["force_vector_RMSE_eV_A"],
                        "pass": force_result["force_vector_RMSE_eV_A"] <= GATES["force_vector_RMSE_eV_A"],
                    },
                    "separating_force_abs_error_eV_A": {
                        "limit": GATES["separating_force_abs_error_eV_A"],
                        "value": force_result["marked_pair_separating_force_eV_A"]["absolute_error"],
                        "pass": force_result["marked_pair_separating_force_eV_A"]["absolute_error"]
                        <= GATES["separating_force_abs_error_eV_A"],
                    },
                },
            }
        )
        atoms.calc = None

    gate_pass = all(all(gate["pass"] for gate in row["provisional_gates"].values()) for row in inference_rows)
    finished = time.perf_counter()
    env = environment_record()
    metrics = {
        "task": TASK,
        "status": "INFERENCE_COMPLETED",
        "scope": "CPU inference/evaluation only on the two specified public, archive-verified, non-sealed AgSi CIF references; no DFT/MPI, training, fine-tuning, model edit, or sealed-label access.",
        "candidate": {
            "name": "V18 clean-core selected final model, epoch 79",
            "model_path": str(MODEL_PATH),
            "checkpoint_path": str(CHECKPOINT_PATH),
            "model_and_data_lineage_hashes": lineage["hashes"],
            "lineage_checks": lineage["checks"],
            "selected_epoch": lineage["selection"].get("selected_epoch"),
            "training_frame_count": len(train_atoms),
            "development_or_test_structure_files_read": False,
        },
        "reference_archives": [
            {key: value for key, value in archive.items() if not key.startswith("_")}
            for archive in reference_records.values()
        ],
        "training_overlap": overlap,
        "runtime": env,
        "cpu_inference_smoke": {
            "status": "PASS" if smoke_pass else "FAIL",
            "device": "cpu",
            "structures_inferred": len(inference_rows),
            "all_predictions_finite": smoke_pass,
            "total_inference_elapsed_s": finished - start_all,
        },
        "provisional_gates": GATES,
        "per_reference_results": inference_rows,
        "screen_decision": {
            "provisional_force_and_separation_gates_both_pass_on_both_references": gate_pass,
            "status": "DIAGNOSTIC_SCREEN_PASS" if gate_pass else "DIAGNOSTIC_SCREEN_FAIL",
            "independent_validation": False,
            "production_pass": False,
            "interpretation": (
                "This is a one-family diagnostic only: both public structures derive from one CIF family and their local "
                "Ag-Si contact distances overlap V18 training coverage. It is not independent validation or a production pass."
            ),
        },
    }
    (report_dir / "metrics.json").write_text(json.dumps(json_safe(metrics), indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    (report_dir / "environment.json").write_text(json.dumps(env, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    print(json.dumps(json_safe({
        "status": metrics["status"],
        "cpu_inference_smoke": metrics["cpu_inference_smoke"],
        "screen_decision": metrics["screen_decision"],
        "per_reference": [
            {
                "label": row["label"],
                "force_vector_RMSE_eV_A": row["force_metrics"]["force_vector_RMSE_eV_A"],
                "separating_force_abs_error_eV_A": row["force_metrics"]["marked_pair_separating_force_eV_A"]["absolute_error"],
            }
            for row in inference_rows
        ],
        "outputs": [str(report_dir / "metrics.json"), str(report_dir / "environment.json")],
    }), indent=2))


if __name__ == "__main__":
    main()
