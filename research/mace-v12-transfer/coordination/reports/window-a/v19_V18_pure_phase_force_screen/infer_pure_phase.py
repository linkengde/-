#!/usr/bin/env python3
"""Screen the frozen V18 candidate on two verified pure-phase DFT controls."""
from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mace-pure-phase-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/mace-pure-phase-cache")
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
MODEL = MODEL_ROOT / "training_clean_core_01/models/MACE_periodic_v18_clean_core.model"
PREVIOUS_SCREEN = Path(
    "research/mace-v12-transfer/coordination/reports/window-b/"
    "v19_existing_AgSi_candidate_force_error_screen_install_and_infer_window_b/infer_and_evaluate.py"
)
CONTROL_ROOT = Path(
    "research/mace-v12-transfer/periodic_interface_v4/"
    "pbe_interface_v19_pure_phase_controls/calculations"
)
LABELS = ["Ag_baseline_k8x8x8_v19", "Ti3SiC2_baseline_k10x10x2_v19"]
FORCE_RMSE_GATE = 0.05


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def exact_geometry_hash(atoms) -> str:
    payload = b"".join(
        (
            np.asarray(atoms.numbers, dtype="<i4").tobytes(),
            np.asarray(atoms.positions, dtype="<f8").tobytes(),
            np.asarray(atoms.cell.array, dtype="<f8").tobytes(),
            np.asarray(atoms.pbc, dtype="u1").tobytes(),
        )
    )
    return hashlib.sha256(payload).hexdigest()


def load_lineage() -> dict:
    helper_path = REPO / PREVIOUS_SCREEN
    spec = importlib.util.spec_from_file_location("v18_screen_helpers", helper_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load the previously reviewed V18 lineage verifier")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    lineage, _ = helper.verify_model_lineage()
    return lineage


def verify_archive(label: str) -> dict:
    folder = REPO / CONTROL_ROOT / label
    if not folder.is_dir():
        raise RuntimeError(f"Missing pure-phase archive: {label}")
    checks = []
    for line in (folder / "SHA256SUMS.txt").read_text().splitlines():
        if not line.strip():
            continue
        expected, filename = line.split(maxsplit=1)
        filename = filename.lstrip("*")
        member = (folder / filename).resolve()
        if folder.resolve() not in member.parents or not member.is_file():
            raise RuntimeError(f"Invalid checksum member in {label}: {filename}")
        actual = sha256(member)
        checks.append({"file": filename, "expected_sha256": expected, "actual_sha256": actual, "pass": expected == actual})
    if not checks or not all(row["pass"] for row in checks):
        raise RuntimeError(f"Archive SHA256SUMS failed: {label}")

    summary = json.loads((folder / "summary.json").read_text())
    verification = json.loads((folder / "verification.json").read_text())
    if summary.get("label") != label or summary.get("scf_converged") is not True:
        raise RuntimeError(f"Archive label/convergence failed: {label}")
    if verification.get("status") != "PASS" or not verification.get("checks") or not all(verification["checks"].values()):
        raise RuntimeError(f"Recorded archive verification failed: {label}")

    source_text = str(summary.get("source", "")).replace("\\", "/")
    marker = "research/mace-v12-transfer/"
    offset = source_text.find(marker)
    if offset < 0:
        raise RuntimeError(f"Cannot resolve input source for {label}")
    source = (REPO / source_text[offset:]).resolve()
    if REPO.resolve() not in source.parents or not source.is_file() or sha256(source) != summary.get("source_sha256"):
        raise RuntimeError(f"Source geometry hash failed: {label}")

    result = folder / f"{label}_PW_PBE.extxyz"
    if not result.is_file() or verification.get("sha256", {}).get(result.name) != sha256(result):
        raise RuntimeError(f"Archived force result hash failed: {label}")
    atoms = read(result, index=0)
    source_atoms = read(source, index=0)
    geometry_checks = {
        "atomic_numbers_and_order_match_input": np.array_equal(atoms.numbers, source_atoms.numbers),
        "positions_match_input": np.array_equal(atoms.positions, source_atoms.positions),
        "cell_matches_input": np.array_equal(atoms.cell.array, source_atoms.cell.array),
        "pbc_matches_input": np.array_equal(atoms.pbc, source_atoms.pbc),
        "formula_matches_summary": atoms.get_chemical_formula() == summary.get("formula"),
        "atom_count_matches_summary": len(atoms) == summary.get("atoms"),
    }
    if not all(geometry_checks.values()):
        raise RuntimeError(f"Archive geometry/input identity failed: {label}: {geometry_checks}")
    force_key = "PW_PBE_forces"
    if force_key not in atoms.arrays:
        raise RuntimeError(f"Missing {force_key} in {label}")
    return {
        "folder": folder,
        "atoms": atoms,
        "reference_forces": np.asarray(atoms.arrays[force_key], dtype=float).copy(),
        "force_key": force_key,
        "checksum_checks": checks,
        "summary": summary,
        "verification": verification,
        "source_path": str(source.relative_to(REPO)),
        "source_sha256": summary["source_sha256"],
        "result_sha256": sha256(result),
        "geometry_sha256": exact_geometry_hash(atoms),
        "geometry_checks": geometry_checks,
    }


def main() -> None:
    if torch.cuda.is_available():
        raise RuntimeError("This inference screen is CPU-only")
    torch.set_num_threads(1)
    lineage = load_lineage()
    model_path = REPO / MODEL
    calculator = MACECalculator(model_paths=str(model_path), device="cpu", default_dtype="float64")
    refs = {label: verify_archive(label) for label in LABELS}
    rows = []
    for label, ref in refs.items():
        atoms = ref["atoms"]
        reference = ref["reference_forces"]
        atoms.calc = calculator
        prediction = np.asarray(atoms.get_forces(), dtype=float)
        atoms.calc = None
        if prediction.shape != reference.shape or not np.isfinite(prediction).all() or not np.isfinite(reference).all():
            raise RuntimeError(f"Invalid force arrays from inference: {label}")
        delta = prediction - reference
        norms = np.linalg.norm(delta, axis=1)
        symbols = atoms.get_chemical_symbols()
        species = {}
        for element in sorted(set(symbols)):
            indices = [i for i, symbol in enumerate(symbols) if symbol == element]
            species[element] = {
                "atom_count": len(indices),
                "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(norms[indices] ** 2))),
                "maximum_atom_error_eV_A": float(np.max(norms[indices])),
                "mean_error_vector_eV_A": np.mean(delta[indices], axis=0).tolist(),
            }
        rmse = float(np.sqrt(np.mean(norms**2)))
        rows.append(
            {
                "label": label,
                "formula": atoms.get_chemical_formula(),
                "atom_count": len(atoms),
                "cell_A": atoms.cell.array.tolist(),
                "pbc": atoms.pbc.tolist(),
                "method": ref["summary"]["method"],
                "dataset_role": ref["summary"].get("dataset_role"),
                "scf_iterations": ref["summary"].get("scf_iterations"),
                "DFT_source_path": ref["source_path"],
                "DFT_source_sha256": ref["source_sha256"],
                "DFT_result_sha256": ref["result_sha256"],
                "result_geometry_sha256": ref["geometry_sha256"],
                "archive_sha256_checks_all_pass": all(row["pass"] for row in ref["checksum_checks"]),
                "archive_verification_status": ref["verification"]["status"],
                "source_geometry_checks": ref["geometry_checks"],
                "DFT_force_vector_RMS_eV_A": float(np.sqrt(np.mean(np.sum(reference**2, axis=1)))),
                "MACE_force_vector_RMS_eV_A": float(np.sqrt(np.mean(np.sum(prediction**2, axis=1)))),
                "all_atom_force_vector_RMSE_eV_A": rmse,
                "maximum_atom_vector_error_eV_A": float(np.max(norms)),
                "force_component_RMSE_eV_A": {axis: float(np.sqrt(np.mean(delta[:, i] ** 2))) for i, axis in enumerate("xyz")},
                "per_species": species,
                "provisional_force_gate": {"limit_eV_A": FORCE_RMSE_GATE, "pass": rmse <= FORCE_RMSE_GATE},
            }
        )

    versions = {}
    for package in ("torch", "mace-torch", "ase", "numpy"):
        versions[package] = importlib.metadata.version(package)
    metrics = {
        "task": "v19_V18_pure_phase_force_screen_window_a",
        "status": "DIAGNOSTIC_COMPLETED",
        "scope": "CPU-only inference on two existing, verified pure-phase numerical-control archives. No new DFT/MPI, no training or model edits, no development/test structures, and no sealed labels.",
        "interpretation_limit": "These controls diagnose transfer to pure phases and are not independent validation or production readiness evidence. The Ti3SiC2 source itself has nonzero DFT forces; compare the model-minus-DFT residual, not only the model force norm.",
        "candidate": {
            "name": "V18 clean-core final model, epoch 79",
            "model_path": str(MODEL),
            "model_sha256": sha256(model_path),
            "checkpoint_sha256": lineage["hashes"]["epoch79_checkpoint_sha256"],
            "train_file_sha256_checked_without_loading_train_structures": lineage["hashes"]["training_input_sha256"],
            "lineage_checks": lineage["checks"],
        },
        "runtime": {"python": sys.version.split()[0], "executable": sys.executable, "platform": platform.platform(), "device": "cpu", "torch_cuda_available": False, "torch_num_threads": 1, "packages": versions},
        "provisional_force_RMSE_gate_eV_A": FORCE_RMSE_GATE,
        "references": rows,
        "development_or_test_files_read": False,
        "withheld_labels_opened": False,
    }
    output = Path(__file__).resolve().parent / "metrics.json"
    output.write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": metrics["status"], "model_sha256": metrics["candidate"]["model_sha256"], "references": [{"label": r["label"], "RMSE_eV_A": r["all_atom_force_vector_RMSE_eV_A"], "gate": r["provisional_force_gate"]} for r in rows]}, indent=2))


if __name__ == "__main__":
    main()
