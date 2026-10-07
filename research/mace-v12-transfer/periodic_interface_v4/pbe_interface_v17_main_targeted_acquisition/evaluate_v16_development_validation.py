#!/usr/bin/env python3
"""Score only the three authorized V17 development-validation DFT archives."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np


LABELS = (
    "AgC_v17_dev_validation_01",
    "AgSi_v17_dev_validation_01",
    "AgTi_v17_dev_validation_01",
)
GATES = {
    "energy_MAE_meV_atom": 10.0,
    "force_vector_RMSE_eV_A": 0.05,
    "separating_force_abs_error_eV_A": 0.10,
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def verify_archive(folder: Path, label: str) -> tuple[Path, dict, dict]:
    checksum_file = folder / "SHA256SUMS.txt"
    require(checksum_file.is_file(), f"Missing SHA256SUMS.txt: {label}")
    for line in checksum_file.read_text().splitlines():
        expected, name = line.split("  ", 1)
        candidate = folder / name
        require(candidate.is_file() and sha(candidate) == expected,
                f"Archive checksum mismatch: {label}/{name}")

    summary = json.loads((folder / "summary.json").read_text())
    verification = json.loads((folder / "verification.json").read_text())
    require(verification.get("status") == "PASS", f"Archive verification is not PASS: {label}")
    require(summary.get("label") == label and summary.get("role") == "development_validation",
            f"Wrong label or role: {label}")
    require(summary.get("scf_converged") is True and summary.get("finite") is True,
            f"DFT label is not converged/finite: {label}")
    require(summary.get("mpi_ranks") == 4, f"Unexpected MPI rank count: {label}")
    output = folder / f"{label}_PW_PBE.extxyz"
    require(output.is_file(), f"Missing compact DFT structure: {label}")
    require(verification.get("input_sha256") == summary.get("source_sha256"),
            f"Input hash disagreement: {label}")
    return output, summary, verification


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--archive-root", required=True, type=Path)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--selection-record", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    repo = args.repo_root.resolve()
    archive_root = args.archive_root.resolve()
    model_path = args.model.resolve()
    selection_path = args.selection_record.resolve()
    output_dir = args.output_dir.resolve()

    allowed_parent = (repo / "research/mace-v12-transfer/periodic_interface_v4/"
                      "pbe_interface_v17_main_targeted_acquisition").resolve()
    require(output_dir.is_relative_to(allowed_parent) and output_dir != allowed_parent,
            "Output must be a child of the V17 main acquisition directory")
    require(not output_dir.exists() or
            (output_dir.is_dir() and not any(output_dir.iterdir())),
            "Refusing to overwrite nonempty development evaluation output")
    require(model_path.is_file() and selection_path.is_file(), "Model/selection record missing")

    selection_bytes = selection_path.read_bytes()
    selection_hash = hashlib.sha256(selection_bytes).hexdigest()
    selection = json.loads(selection_bytes)
    model_hash = sha(model_path)
    require(selection.get("model_sha256") == model_hash, "Frozen V16 model hash mismatch")
    require(Path(selection.get("model_path", "")).resolve() == model_path,
            "Selection record does not identify this model")
    require(selection.get("reviewed_by") == "window-a" and
            selection.get("blind_labels_used_for_selection") is False and
            selection.get("selection_uses_only_training_validation") is True,
            "V16 selection record is not a frozen training/validation-only choice")
    require(selection.get("completed_epochs") == 80 and
            selection.get("selected_epoch_zero_based") == 79,
            "Unexpected V16 completion/selected checkpoint")
    for field in ("training_stdout", "epoch_completion_evidence"):
        evidence = selection[field]
        evidence_path = Path(evidence["path"]).resolve()
        require(evidence_path.is_file() and sha(evidence_path) == evidence["sha256"],
                f"V16 completion evidence hash mismatch: {field}")
    stdout = Path(selection["training_stdout"]["path"]).read_text()
    epochs = [int(x) for x in re.findall(r"INFO: Epoch (\d+):", stdout)]
    require(epochs == list(range(80)) and "INFO: Done" in stdout,
            "V16 80-epoch completion log is incomplete")
    chosen = re.findall(r"Loaded Stage one model from epoch (\d+) for evaluation", stdout)
    require(chosen and int(chosen[-1]) == selection["selected_epoch_zero_based"],
            "V16 selected checkpoint does not match completion log")

    verified_archives = []
    for label in LABELS:
        output, summary, verification = verify_archive(archive_root / label, label)
        verified_archives.append({
            "label": label,
            "archive_sha256_manifest": sha(archive_root / label / "SHA256SUMS.txt"),
            "output_sha256": sha(output),
            "source_sha256": summary["source_sha256"],
            "verification": verification["status"],
        })

    # Only explicitly named development archives are opened below. No blind/test paths
    # are traversed or passed to the calculator.
    from ase.io import read
    import torch
    from mace.calculators import MACECalculator

    torch.set_num_threads(1)
    torch.set_default_dtype(torch.float64)
    loaded = torch.load(model_path, map_location="cpu", weights_only=False).double().eval()
    calculator = MACECalculator(models=[loaded], device="cpu", default_dtype="float64")

    structure_rows = []
    atom_rows = []
    for label in LABELS:
        folder = archive_root / label
        output = folder / f"{label}_PW_PBE.extxyz"
        summary = json.loads((folder / "summary.json").read_text())
        atoms = read(output)
        require(len(atoms) == summary["atom_count"], f"Atom count mismatch: {label}")
        require(atoms.info.get("candidate_role") == "development_validation",
                f"Candidate role mismatch: {label}")
        require(atoms.info.get("PW_PBE_energy_eV") is not None,
                f"Native DFT energy missing: {label}")
        require("PW_PBE_forces" in atoms.arrays, f"DFT forces missing: {label}")
        reference_energy = float(summary["energy_eV_cell"])
        extxyz_energy = float(atoms.info["PW_PBE_energy_eV"])
        require(abs(reference_energy - extxyz_energy) <= 1e-10,
                f"Summary/extxyz energy disagreement: {label}")
        reference_forces = np.asarray(atoms.arrays["PW_PBE_forces"], dtype=float)
        require(reference_forces.shape == (len(atoms), 3) and np.isfinite(reference_forces).all(),
                f"Invalid DFT forces: {label}")

        predicted = atoms.copy()
        predicted.calc = calculator
        energy_prediction = float(predicted.get_potential_energy())
        force_prediction = np.asarray(predicted.get_forces(), dtype=float)
        require(np.isfinite(energy_prediction) and force_prediction.shape == reference_forces.shape
                and np.isfinite(force_prediction).all(), f"Invalid MACE prediction: {label}")

        delta = force_prediction - reference_forces
        force_norms = np.linalg.norm(delta, axis=1)
        central = np.flatnonzero(np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool))
        require(len(central) == 2, f"Marked Ag-X pair missing/duplicated: {label}")
        ag_idx = [int(i) for i in central if atoms[int(i)].symbol == "Ag"]
        x_idx = [int(i) for i in central if atoms[int(i)].symbol != "Ag"]
        require(len(ag_idx) == len(x_idx) == 1, f"Marked pair is not Ag-X: {label}")
        i, j = ag_idx[0], x_idx[0]
        vector = atoms.positions[j] - atoms.positions[i]
        distance = float(np.linalg.norm(vector))
        require(distance > 0, f"Invalid marked-pair distance: {label}")
        unit = vector / distance
        separation_ref = float(np.dot(reference_forces[j] - reference_forces[i], unit))
        separation_pred = float(np.dot(force_prediction[j] - force_prediction[i], unit))
        separation_error = abs(separation_pred - separation_ref)

        n = len(atoms)
        energy_error = abs(energy_prediction - reference_energy) * 1000.0 / n
        vector_rmse = float(np.sqrt(np.mean(np.sum(delta ** 2, axis=1))))
        max_index = int(np.argmax(force_norms))
        metrics = {
            "label": label,
            "interface": label.split("_", 1)[0],
            "role": "development_validation",
            "atoms": n,
            "formula": atoms.get_chemical_formula(),
            "reference_sha256": sha(output),
            "energy_DFT_native_eV_cell": reference_energy,
            "energy_DFT_free_eV_cell": summary.get("free_energy_eV_cell"),
            "energy_prediction_eV_cell": energy_prediction,
            "energy_abs_error_meV_atom": energy_error,
            "force_vector_RMSE_eV_A": vector_rmse,
            "force_vector_max_error_eV_A": float(force_norms[max_index]),
            "force_max_atom_index": max_index,
            "force_max_atom_id": int(atoms.arrays["lammps_id"][max_index]),
            "force_max_atom_symbol": atoms[max_index].symbol,
            "marked_Ag_id": int(atoms.arrays["lammps_id"][i]),
            "marked_X_symbol": atoms[j].symbol,
            "marked_X_id": int(atoms.arrays["lammps_id"][j]),
            "marked_pair_distance_A": distance,
            "separating_force_reference_eV_A": separation_ref,
            "separating_force_prediction_eV_A": separation_pred,
            "separating_force_abs_error_eV_A": separation_error,
        }
        structure_rows.append(metrics)
        ids = atoms.arrays.get("lammps_id", np.arange(len(atoms)))
        for idx, atom in enumerate(atoms):
            atom_rows.append({
                "label": label,
                "atom_index": idx,
                "atom_id": int(ids[idx]),
                "symbol": atom.symbol,
                "central_pair": bool(atoms.arrays.get("central_pair", np.zeros(len(atoms)))[idx]),
                "reference_fx_eV_A": float(reference_forces[idx, 0]),
                "reference_fy_eV_A": float(reference_forces[idx, 1]),
                "reference_fz_eV_A": float(reference_forces[idx, 2]),
                "prediction_fx_eV_A": float(force_prediction[idx, 0]),
                "prediction_fy_eV_A": float(force_prediction[idx, 1]),
                "prediction_fz_eV_A": float(force_prediction[idx, 2]),
                "error_fx_eV_A": float(delta[idx, 0]),
                "error_fy_eV_A": float(delta[idx, 1]),
                "error_fz_eV_A": float(delta[idx, 2]),
                "error_vector_norm_eV_A": float(force_norms[idx]),
            })

    by_interface = {}
    for kind in ("AgC", "AgSi", "AgTi"):
        row = next(r for r in structure_rows if r["interface"] == kind)
        gates = {
            "energy": row["energy_abs_error_meV_atom"] <= GATES["energy_MAE_meV_atom"],
            "force_vector": row["force_vector_RMSE_eV_A"] <= GATES["force_vector_RMSE_eV_A"],
            "separation": row["separating_force_abs_error_eV_A"] <= GATES["separating_force_abs_error_eV_A"],
        }
        by_interface[kind] = {"metrics": row, "provisional_gates": gates,
                              "all_provisional_gates_met": all(gates.values())}
    all_gates = all(x["all_provisional_gates_met"] for x in by_interface.values())
    status = "DEV_SCREEN_PASS_LIMITED" if all_gates else "DEV_SCREEN_FAIL"

    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "overall_status": status,
        "scope": "Only the three hash-verified V17 development-validation archives were scored.",
        "model_path": str(model_path.relative_to(repo)),
        "model_sha256": model_hash,
        "selection_record_sha256": selection_hash,
        "selected_epoch": selection["selected_epoch"],
        "archive_inputs": verified_archives,
        "thresholds": GATES,
        "per_structure": structure_rows,
        "by_interface": by_interface,
        "limitations": "Development screen only; these candidates share small-cluster parent motifs. Passing narrow gates would not establish general transferability or authorize MD/TTM. V16 is frozen; no held-out labels or outputs were opened.",
    }
    (output_dir / "evaluation.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    with (output_dir / "per_structure.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(structure_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(structure_rows)
    with (output_dir / "per_atom.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(atom_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(atom_rows)

    files = sorted(p for p in output_dir.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt")
    (output_dir / "SHA256SUMS.txt").write_text(
        "".join(f"{sha(p)}  {p.name}\n" for p in files))
    print(json.dumps({"overall_status": status,
                      "by_interface": {k: v["metrics"] for k, v in by_interface.items()}},
                     indent=2))


if __name__ == "__main__":
    main()
