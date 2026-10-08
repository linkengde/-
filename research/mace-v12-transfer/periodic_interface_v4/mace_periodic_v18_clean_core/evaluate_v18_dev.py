#!/usr/bin/env python3
"""Score only the three V17 development structures with a frozen V17 model."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.io import read

ENTRY = Path("research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core")
GATES = {"energy_abs_error_meV_atom": 10.0, "force_vector_RMSE_eV_A": 0.05, "separating_force_abs_error_eV_A": 0.10}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--model", required=True, type=Path)
    ap.add_argument("--selection-record", required=True, type=Path)
    args = ap.parse_args()
    repo = args.repo_root.resolve()
    entry = repo / ENTRY
    out = args.output_dir.resolve()
    model = args.model.resolve()
    selection_path = args.selection_record.resolve()
    require(out.is_relative_to(entry) and out != entry, "Evaluation output must be a V18 directory child")
    require(not out.exists() or (out.is_dir() and not any(out.iterdir())), "Nonempty evaluation output refused")
    require(model.is_file() and selection_path.is_file(), "Model or selection record missing")
    manifest = json.loads((entry / "data/dataset_manifest.json").read_text())
    require(manifest["status"].startswith("PROVISIONAL"), "Refuse evaluation if data status is not explicit")
    require(sha(entry / "data/train.extxyz") == manifest["train_input_sha256"] and sha(entry / "data/valid.extxyz") == manifest["development_validation_input_sha256"], "Dataset changed since build")
    selection_bytes = selection_path.read_bytes()
    selection = json.loads(selection_bytes)
    model_hash = sha(model)
    require(selection.get("dataset_manifest_sha256") == sha(entry / "data/dataset_manifest.json"), "Selection/dataset manifest snapshot mismatch")
    require(selection.get("training_input_sha256") == sha(entry / "data/train.extxyz") and selection.get("development_validation_input_sha256") == sha(entry / "data/valid.extxyz"), "Selection/data snapshot mismatch")
    require(selection.get("model_sha256") == model_hash and Path(selection.get("model_path", "")).resolve() == model, "Selection/model hash mismatch")
    require(selection.get("training_exit_code") == 0 and selection.get("completed_epochs") == 80 and selection.get("selected_epoch") == 79, "A must verify complete fixed 80-epoch run and final epoch selection")
    require(selection.get("selection_uses_only_training_validation") is True and selection.get("blind_labels_used_for_selection") is False and selection.get("test_file_passed_to_training") is False, "Invalid split/selection policy")
    require(selection.get("reviewed_by") == "window-a" and selection.get("selection_reason"), "A selection record required")
    train_log = Path(selection["training_stdout"]["path"]).resolve()
    require(sha(train_log) == selection["training_stdout"]["sha256"], "Training stdout hash mismatch")
    log = train_log.read_text()
    require([int(x) for x in re.findall(r"INFO: Epoch (\d+):", log)] == list(range(80)), "Training log does not show all 80 epochs")
    require("INFO: Training complete" in log and "INFO: Done" in log, "Training completion evidence absent")
    require(sha(selection_path) == hashlib.sha256(selection_bytes).hexdigest(), "Selection record changed during evaluation")

    # No test archive or test label path is opened here. The only references are the
    # three already authorized V17 development-validation frames.
    frames = read(entry / "data/valid.extxyz", index=":")
    require(len(frames) == 3 and {str(a.info.get("config_type", "")).split("_v17_dev_validation")[0] for a in frames} == {"AgC", "AgSi", "AgTi"}, "Unexpected V17 development set")
    for a in frames:
        require("REF_energy" in a.info and "REF_forces" in a.arrays, "Dev reference labels missing")
    import torch
    from mace.calculators import MACECalculator

    torch.set_default_dtype(torch.float64)
    loaded = torch.load(model, map_location="cpu", weights_only=False).double().eval()
    calculator = MACECalculator(models=[loaded], device="cpu", default_dtype="float64")
    rows = []
    atom_rows = []
    for atoms in frames:
        name = str(atoms.info["config_type"])
        predicted = atoms.copy()
        predicted.calc = calculator
        ep = float(predicted.get_potential_energy())
        fp = np.asarray(predicted.get_forces(), dtype=float)
        eref = float(atoms.info["REF_energy"])
        fref = np.asarray(atoms.arrays["REF_forces"], dtype=float)
        require(np.isfinite(ep) and fp.shape == fref.shape and np.isfinite(fp).all(), f"Invalid prediction: {name}")
        delta = fp - fref
        central = np.flatnonzero(np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool))
        require(len(central) == 2, f"Marked pair missing: {name}")
        ag = [int(i) for i in central if atoms[i].symbol == "Ag"]
        contact = [int(i) for i in central if atoms[i].symbol != "Ag"]
        require(len(ag) == len(contact) == 1, f"Marked pair is not Ag-X: {name}")
        i, j = ag[0], contact[0]
        vector = atoms.positions[j] - atoms.positions[i]
        distance = float(np.linalg.norm(vector))
        require(distance > 0, f"Zero marked-pair distance: {name}")
        direction = vector / distance
        radial_ref = float(np.dot(fref[j] - fref[i], direction))
        radial_pred = float(np.dot(fp[j] - fp[i], direction))
        energy_error = abs(ep - eref) * 1000.0 / len(atoms)
        force_rmse = float(np.sqrt(np.mean(np.sum(delta**2, axis=1))))
        max_force_error = float(np.linalg.norm(delta, axis=1).max())
        separating_error = abs(radial_pred - radial_ref)
        require(all(np.isfinite(x) for x in (energy_error, force_rmse, max_force_error, separating_error)), f"Nonfinite metric: {name}")
        interface = name.split("_v17_dev_validation")[0]
        for k, atom in enumerate(atoms):
            longitudinal = float(np.dot(delta[k], direction))
            transverse = float(np.linalg.norm(delta[k] - longitudinal * direction))
            atom_rows.append({"label": name, "atom_index": k, "lammps_id": int(atoms.arrays["lammps_id"][k]), "species": atom.symbol, "central_pair": bool(atoms.arrays["central_pair"][k]), "x_A": float(atom.position[0]), "y_A": float(atom.position[1]), "z_A": float(atom.position[2]), "ref_fx": float(fref[k,0]), "ref_fy": float(fref[k,1]), "ref_fz": float(fref[k,2]), "pred_fx": float(fp[k,0]), "pred_fy": float(fp[k,1]), "pred_fz": float(fp[k,2]), "error_norm_eV_A": float(np.linalg.norm(delta[k])), "error_projected_pair_eV_A": longitudinal, "error_transverse_pair_eV_A": transverse})
        rows.append({"label": name, "interface": interface, "role": "development_validation_only", "atoms": len(atoms), "formula": atoms.get_chemical_formula(), "model_sha256": model_hash, "reference_sha256": atoms.info["v17_reference_sha256"], "energy_abs_error_meV_atom": energy_error, "force_vector_RMSE_eV_A": force_rmse, "force_vector_max_error_eV_A": max_force_error, "marked_pair_distance_A": distance, "separating_force_abs_error_eV_A": separating_error})
    per_interface = {}
    for row in rows:
        passed = all(row[key] <= GATES[key] for key in GATES)
        per_interface[row["interface"]] = {**{key: row[key] for key in GATES}, "provisional_gates_met": passed}
    overall = "FAIL" if any(not row["provisional_gates_met"] for row in per_interface.values()) else "UNDETERMINED"
    result = {"status": overall, "cycle": "V18 clean30 controlled numerical screening", "model_path": str(model), "model_sha256": model_hash, "selection_record_sha256": sha(selection_path), "development_validation_only": True, "withheld_test_labels_opened": False, "gates": GATES, "per_structure": rows, "by_interface": per_interface, "limitations": ["All 30 rows have archive support; frame0 Gamma convergence and uniform method consistency remain unresolved.", "The three dev structures were already used to diagnose V16 and are not fresh final blind tests.", "All eight added labels are paired local training diagnostics, not independent tests.", "A gate pass would not establish thermal, pressure, morphology, or production transferability."]}
    out.mkdir(parents=True, exist_ok=True)
    (out / "evaluation.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    with (out / "per_structure.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with (out / "per_atom.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(atom_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(atom_rows)
    (out / "SHA256SUMS.txt").write_text("".join(f"{sha(p)}  {p.name}\n" for p in sorted(out.iterdir()) if p.is_file()))
    print(json.dumps({"status": overall, "by_interface": per_interface}, indent=2))


if __name__ == "__main__":
    main()
