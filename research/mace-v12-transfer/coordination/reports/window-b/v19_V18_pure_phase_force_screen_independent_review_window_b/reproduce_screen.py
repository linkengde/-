#!/usr/bin/env python3
"""Independently reproduce the V18 pure-phase force screen on two public controls."""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mace-v18-pure-phase-review-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/mace-v18-pure-phase-review-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator


def find_repo_root(start: Path) -> Path:
    for parent in (start, *start.parents):
        if (parent / ".git").exists():
            return parent
    raise RuntimeError("Cannot locate repository root")


REPO = find_repo_root(Path(__file__).resolve())
MODEL_ROOT = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core")
MODEL_REL = MODEL_ROOT / "training_clean_core_01/models/MACE_periodic_v18_clean_core.model"
CHECKPOINT_REL = MODEL_ROOT / "training_clean_core_01/checkpoints/MACE_periodic_v18_clean_core_run-45_epoch-79.pt"
SELECTION_REL = MODEL_ROOT / "training_clean_core_01/selection_record.json"
ARTIFACT_MANIFEST_REL = MODEL_ROOT / "training_clean_core_01/artifact_manifest.json"
REFERENCE_METRICS_REL = Path(
    "research/mace-v12-transfer/coordination/reports/window-a/"
    "v19_V18_pure_phase_force_screen/metrics.json"
)
CONTROL_ROOT = Path(
    "research/mace-v12-transfer/periodic_interface_v4/"
    "pbe_interface_v19_pure_phase_controls/calculations"
)
REFERENCES = {
    "Ag_baseline_k8x8x8_v19": Path(
        "research/mace-v12-transfer/periodic_interface_v4/"
        "pbe_interface_v19_pure_phase_controls/inputs/Ag_baseline_PROPOSAL.extxyz"
    ),
    "Ti3SiC2_baseline_k10x10x2_v19": Path(
        "research/mace-v12-transfer/periodic_interface_v4/"
        "pbe_interface_v19_pure_phase_controls/inputs/Ti3SiC2_baseline_PROPOSAL.extxyz"
    ),
}
FORCE_GATE = 0.05
TOLERANCE = {"rtol": 1e-7, "atol": 1e-8}


def repo_file(relative: Path) -> Path:
    path = (REPO / relative).resolve()
    if REPO.resolve() not in path.parents or not path.is_file():
        raise RuntimeError(f"Missing or out-of-repository input: {relative}")
    return path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def geometry_hash(atoms) -> str:
    payload = b"".join(
        (
            np.asarray(atoms.numbers, dtype="<i4").tobytes(),
            np.asarray(atoms.positions, dtype="<f8").tobytes(),
            np.asarray(atoms.cell.array, dtype="<f8").tobytes(),
            np.asarray(atoms.pbc, dtype="u1").tobytes(),
        )
    )
    return hashlib.sha256(payload).hexdigest()


def close(actual: float, expected: float) -> bool:
    return bool(np.isclose(actual, expected, **TOLERANCE))


def verify_model_lineage(reference_metrics: dict) -> dict:
    model_path = repo_file(MODEL_REL)
    checkpoint_path = repo_file(CHECKPOINT_REL)
    selection_path = repo_file(SELECTION_REL)
    artifact_path = repo_file(ARTIFACT_MANIFEST_REL)
    selection = json.loads(selection_path.read_text())
    artifact = json.loads(artifact_path.read_text())
    selected = reference_metrics["candidate"]
    model_hash = sha256(model_path)
    checkpoint_hash = sha256(checkpoint_path)
    artifact_rows = {row["path"]: row for row in artifact.get("files", [])}
    checks = {
        "selected_epoch_is_79": selection.get("selected_epoch") == 79,
        "selected_model_hash_matches_selection": model_hash == selection.get("model_sha256"),
        "selected_checkpoint_hash_matches_selection": checkpoint_hash
        == selection.get("recovery_verification", {}).get("checkpoint_sha256"),
        "model_hash_matches_A_report": model_hash == selected.get("model_sha256"),
        "checkpoint_hash_matches_A_report": checkpoint_hash == selected.get("checkpoint_sha256"),
        "model_recorded_published": artifact_rows.get(str(MODEL_REL), {}).get("sha256") == model_hash
        and artifact_rows.get(str(MODEL_REL), {}).get("published") is True,
        "checkpoint_recorded_published": artifact_rows.get(str(CHECKPOINT_REL), {}).get("sha256")
        == checkpoint_hash
        and artifact_rows.get(str(CHECKPOINT_REL), {}).get("published") is True,
    }
    if not all(checks.values()):
        raise RuntimeError(f"V18 model/checkpoint lineage failed: {checks}")
    return {
        "model_path": str(MODEL_REL),
        "model_sha256": model_hash,
        "checkpoint_path": str(CHECKPOINT_REL),
        "checkpoint_sha256": checkpoint_hash,
        "selection_record_sha256": sha256(selection_path),
        "artifact_manifest_sha256": sha256(artifact_path),
        "selected_epoch": 79,
        "checks": checks,
    }


def verify_archive(label: str, expected_metrics: dict) -> dict:
    folder = repo_file(CONTROL_ROOT / label / "summary.json").parent
    checksums_path = repo_file(CONTROL_ROOT / label / "SHA256SUMS.txt")
    checksum_rows = []
    for line in checksums_path.read_text().splitlines():
        if not line.strip():
            continue
        expected_hash, filename = line.split(maxsplit=1)
        filename = filename.lstrip("*")
        member = (folder / filename).resolve()
        if folder.resolve() not in member.parents or not member.is_file():
            raise RuntimeError(f"Invalid SHA256 member for {label}: {filename}")
        actual_hash = sha256(member)
        checksum_rows.append(
            {"file": filename, "expected_sha256": expected_hash, "actual_sha256": actual_hash, "pass": expected_hash == actual_hash}
        )
    if not checksum_rows or not all(row["pass"] for row in checksum_rows):
        raise RuntimeError(f"Archive SHA256 verification failed for {label}")

    summary = json.loads((folder / "summary.json").read_text())
    verification = json.loads((folder / "verification.json").read_text())
    if summary.get("label") != label or summary.get("scf_converged") is not True:
        raise RuntimeError(f"Label or SCF convergence check failed for {label}")
    if verification.get("status") != "PASS" or not verification.get("checks") or not all(verification["checks"].values()):
        raise RuntimeError(f"Archive verification receipt failed for {label}")

    source = repo_file(REFERENCES[label])
    source_hash = sha256(source)
    source_text = str(summary.get("source", "")).replace("\\", "/")
    marker = "research/mace-v12-transfer/"
    offset = source_text.find(marker)
    if offset < 0 or Path(source_text[offset:]) != REFERENCES[label]:
        raise RuntimeError(f"Unexpected source geometry for {label}: {summary.get('source')}")
    if source_hash != summary.get("source_sha256") or source_hash != expected_metrics.get("DFT_source_sha256"):
        raise RuntimeError(f"Source geometry hash mismatch for {label}")

    result_name = f"{label}_PW_PBE.extxyz"
    result_path = folder / result_name
    result_hash = sha256(result_path)
    if result_hash != verification.get("sha256", {}).get(result_name):
        raise RuntimeError(f"Archived force result hash mismatch for {label}")
    if result_hash != expected_metrics.get("DFT_result_sha256"):
        raise RuntimeError(f"Result hash differs from A's reference metrics for {label}")

    atoms = read(result_path, index=0)
    source_atoms = read(source, index=0)
    id_key = "lammps_id"
    geometry_checks = {
        "atomic_numbers_and_order": np.array_equal(atoms.numbers, source_atoms.numbers),
        "positions_exact": np.array_equal(atoms.positions, source_atoms.positions),
        "cell_exact": np.array_equal(atoms.cell.array, source_atoms.cell.array),
        "pbc_exact": np.array_equal(atoms.pbc, source_atoms.pbc),
        "atom_ids_exact": id_key not in atoms.arrays
        and id_key not in source_atoms.arrays
        or id_key in atoms.arrays
        and id_key in source_atoms.arrays
        and np.array_equal(atoms.arrays[id_key], source_atoms.arrays[id_key]),
        "formula_matches_summary": atoms.get_chemical_formula() == summary.get("formula"),
        "atom_count_matches_summary": len(atoms) == summary.get("atoms"),
        "geometry_hash_matches_A_report": geometry_hash(atoms) == expected_metrics.get("result_geometry_sha256"),
    }
    if not all(geometry_checks.values()):
        raise RuntimeError(f"Geometry identity failed for {label}: {geometry_checks}")

    force_key = "PW_PBE_forces"
    if force_key not in atoms.arrays:
        raise RuntimeError(f"Missing archived force array {force_key} for {label}")
    reference_forces = np.asarray(atoms.arrays[force_key], dtype=float).copy()
    if not np.isfinite(reference_forces).all():
        raise RuntimeError(f"Non-finite archived force array for {label}")
    return {
        "atoms": atoms,
        "forces": reference_forces,
        "source_path": str(REFERENCES[label]),
        "source_sha256": source_hash,
        "result_path": str((CONTROL_ROOT / label / result_name)),
        "result_sha256": result_hash,
        "geometry_sha256": geometry_hash(atoms),
        "geometry_checks": geometry_checks,
        "archive_checksums": checksum_rows,
        "archive_verification_status": verification["status"],
        "scf_iterations": summary.get("scf_iterations"),
        "method": summary.get("method"),
    }


def evaluate(label: str, archived: dict, expected: dict, calculator) -> dict:
    atoms = archived["atoms"]
    reference = archived["forces"]
    atoms.calc = calculator
    prediction = np.asarray(atoms.get_forces(), dtype=float)
    atoms.calc = None
    if prediction.shape != reference.shape or not np.isfinite(prediction).all():
        raise RuntimeError(f"Invalid CPU MACE forces for {label}")

    delta = prediction - reference
    per_atom_norm = np.linalg.norm(delta, axis=1)
    component = {axis: float(np.sqrt(np.mean(delta[:, i] ** 2))) for i, axis in enumerate("xyz")}
    species_rows = {}
    symbols = atoms.get_chemical_symbols()
    for symbol in sorted(set(symbols)):
        mask = np.asarray([item == symbol for item in symbols])
        species_norm = per_atom_norm[mask]
        species_rows[symbol] = {
            "atom_count": int(mask.sum()),
            "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(species_norm**2))),
            "maximum_atom_vector_error_eV_A": float(np.max(species_norm)),
            "force_component_RMSE_eV_A": {
                axis: float(np.sqrt(np.mean(delta[mask, i] ** 2))) for i, axis in enumerate("xyz")
            },
        }

    rmse = float(np.sqrt(np.mean(per_atom_norm**2)))
    dft_rms = float(np.sqrt(np.mean(np.sum(reference**2, axis=1))))
    model_rms = float(np.sqrt(np.mean(np.sum(prediction**2, axis=1))))
    expected_per_species = expected["per_species"]
    comparisons = {
        "all_atom_force_vector_RMSE": close(rmse, expected["all_atom_force_vector_RMSE_eV_A"]),
        "maximum_atom_vector_error": close(per_atom_norm.max(), expected["maximum_atom_vector_error_eV_A"]),
        "DFT_force_vector_RMS": close(dft_rms, expected["DFT_force_vector_RMS_eV_A"]),
        "MACE_force_vector_RMS": close(model_rms, expected["MACE_force_vector_RMS_eV_A"]),
        "all_cartesian_component_RMSE": all(
            close(component[axis], expected["force_component_RMSE_eV_A"][axis]) for axis in "xyz"
        ),
        "per_species_force_vector_RMSE": all(
            close(species_rows[symbol]["force_vector_RMSE_eV_A"], expected_per_species[symbol]["force_vector_RMSE_eV_A"])
            for symbol in expected_per_species
        ),
    }
    if not all(comparisons.values()):
        raise RuntimeError(f"Independent reproduction differs from A's report for {label}: {comparisons}")

    return {
        "label": label,
        "formula": atoms.get_chemical_formula(),
        "atom_count": len(atoms),
        "source_path": archived["source_path"],
        "source_sha256": archived["source_sha256"],
        "result_path": archived["result_path"],
        "result_sha256": archived["result_sha256"],
        "result_geometry_sha256": archived["geometry_sha256"],
        "archive_verification_status": archived["archive_verification_status"],
        "archive_checksum_count": len(archived["archive_checksums"]),
        "archive_sha256_all_pass": all(row["pass"] for row in archived["archive_checksums"]),
        "geometry_identity_checks": archived["geometry_checks"],
        "method": archived["method"],
        "scf_iterations": archived["scf_iterations"],
        "DFT_force_vector_RMS_eV_A": dft_rms,
        "MACE_force_vector_RMS_eV_A": model_rms,
        "all_atom_force_vector_RMSE_eV_A": rmse,
        "maximum_atom_vector_error_eV_A": float(per_atom_norm.max()),
        "force_component_RMSE_eV_A": component,
        "per_species": species_rows,
        "provisional_force_gate": {"limit_eV_A": FORCE_GATE, "pass": rmse <= FORCE_GATE},
        "matches_A_report_within_tolerance": comparisons,
    }


def main() -> None:
    if torch.cuda.is_available():
        raise RuntimeError("This audit must use CPU inference")
    torch.set_num_threads(1)
    reference_path = repo_file(REFERENCE_METRICS_REL)
    reference_metrics = json.loads(reference_path.read_text())
    if reference_metrics.get("status") != "DIAGNOSTIC_COMPLETED":
        raise RuntimeError("A's public pure-phase screen is not marked complete")
    if reference_metrics.get("development_or_test_files_read") is not False or reference_metrics.get("withheld_labels_opened") is not False:
        raise RuntimeError("A's report scope flags do not pass")

    lineage = verify_model_lineage(reference_metrics)
    reference_by_label = {row["label"]: row for row in reference_metrics["references"]}
    if set(reference_by_label) != set(REFERENCES):
        raise RuntimeError("A's metrics do not contain exactly the two assigned public controls")
    archived = {
        label: verify_archive(label, reference_by_label[label]) for label in REFERENCES
    }

    model_path = repo_file(MODEL_REL)
    calculator = MACECalculator(model_paths=str(model_path), device="cpu", default_dtype="float64")
    results = [evaluate(label, archived[label], reference_by_label[label], calculator) for label in REFERENCES]
    versions = {
        package: importlib.metadata.version(package)
        for package in ("torch", "mace-torch", "ase", "numpy")
    }
    output = {
        "task": "v19_V18_pure_phase_force_screen_independent_review_window_b",
        "status": "INDEPENDENT_REPRODUCTION_PASS",
        "scope": "CPU-only force inference on the two named public pure-phase numerical controls. No train.extxyz, development/test/holdout structures or labels, DFT/MPI, or training was accessed or run.",
        "interpretation_limit": "The Ag and Ti3SiC2 structures are numerical controls, not independent validation. A reproduction match confirms the published diagnostic arithmetic and lineage only.",
        "model_lineage": lineage,
        "reference_metrics_sha256": sha256(reference_path),
        "runtime": {
            "python": sys.version.split()[0],
            "executable": sys.executable,
            "platform": platform.platform(),
            "device": "cpu",
            "torch_cuda_available": False,
            "torch_num_threads": 1,
            "packages": versions,
        },
        "force_gate_eV_A": FORCE_GATE,
        "comparison_tolerance": TOLERANCE,
        "references": results,
        "train_extxyz_opened": False,
        "development_or_test_files_read": False,
        "withheld_labels_opened": False,
    }

    out_dir = Path(__file__).resolve().parent
    (out_dir / "audit.json").write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n")
    with (out_dir / "reference_comparison.csv").open("w", newline="") as stream:
        columns = [
            "label", "atom_count", "DFT_force_vector_RMS_eV_A", "MACE_force_vector_RMS_eV_A",
            "all_atom_force_vector_RMSE_eV_A", "maximum_atom_vector_error_eV_A",
            "force_component_RMSE_x_eV_A", "force_component_RMSE_y_eV_A", "force_component_RMSE_z_eV_A",
            "force_gate_pass", "per_species_force_vector_RMSE_eV_A", "all_metrics_match_A_report",
        ]
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in results:
            writer.writerow(
                {
                    "label": row["label"],
                    "atom_count": row["atom_count"],
                    "DFT_force_vector_RMS_eV_A": row["DFT_force_vector_RMS_eV_A"],
                    "MACE_force_vector_RMS_eV_A": row["MACE_force_vector_RMS_eV_A"],
                    "all_atom_force_vector_RMSE_eV_A": row["all_atom_force_vector_RMSE_eV_A"],
                    "maximum_atom_vector_error_eV_A": row["maximum_atom_vector_error_eV_A"],
                    "force_component_RMSE_x_eV_A": row["force_component_RMSE_eV_A"]["x"],
                    "force_component_RMSE_y_eV_A": row["force_component_RMSE_eV_A"]["y"],
                    "force_component_RMSE_z_eV_A": row["force_component_RMSE_eV_A"]["z"],
                    "force_gate_pass": row["provisional_force_gate"]["pass"],
                    "per_species_force_vector_RMSE_eV_A": json.dumps(
                        {key: value["force_vector_RMSE_eV_A"] for key, value in row["per_species"].items()},
                        sort_keys=True,
                    ),
                    "all_metrics_match_A_report": all(row["matches_A_report_within_tolerance"].values()),
                }
            )
    print(
        json.dumps(
            {
                "status": output["status"],
                "model_sha256": lineage["model_sha256"],
                "references": [
                    {
                        "label": row["label"],
                        "force_RMSE_eV_A": row["all_atom_force_vector_RMSE_eV_A"],
                        "force_gate_pass": row["provisional_force_gate"]["pass"],
                    }
                    for row in results
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
