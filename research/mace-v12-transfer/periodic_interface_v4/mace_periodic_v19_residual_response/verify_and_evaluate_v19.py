#!/usr/bin/env python3
"""Verify all fixed V19 epochs and score only the approved dev3 split."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
DEFAULT_RUN = ROOT / "runs/v19_train36_dev3_20261023_A"
TAG = "MACE_periodic_v19_residual_response_run-45"
AUDIT_JSON = ROOT.parents[1] / "coordination/reports/window-b/v19_dataset_entry_independent_audit/audit.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_path(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN)
    args = parser.parse_args()
    run = args.run_dir.resolve()
    log = run / f"logs/{TAG}.log"
    stdout = run / "train.stdout"
    checkpoint_dir = run / "checkpoints"
    checkpoints = {int(m.group(1)): p for p in checkpoint_dir.glob(f"{TAG}_epoch-*.pt") if (m := re.fullmatch(rf"{re.escape(TAG)}_epoch-(\d+)\.pt", p.name))}
    log_text = log.read_text()
    epoch_lines = [int(x) for x in re.findall(r"INFO: Epoch (\d+):", log_text)]
    metric_pattern = re.compile(
        r"INFO: Epoch (\d+): head: Default, loss=\s*([0-9.eE+-]+), "
        r"RMSE_E_per_atom=\s*([0-9.eE+-]+) meV, RMSE_F=\s*([0-9.eE+-]+) meV / A"
    )
    epoch_metrics = [
        {
            "epoch": int(m.group(1)),
            "validation_loss": float(m.group(2)),
            "validation_energy_RMSE_meV_per_atom": float(m.group(3)),
            "validation_force_component_RMSE_meV_A": float(m.group(4)),
        }
        for m in metric_pattern.finditer(log_text)
    ]
    marker = "Loaded Stage one model from epoch 79 for evaluation"
    receipt = json.loads((run / "launcher_completion.json").read_text())
    assert epoch_lines == list(range(80)), f"Expected exactly epochs 0-79; got {epoch_lines}"
    assert [row["epoch"] for row in epoch_metrics] == list(range(80)), "Expected validation metrics for every epoch 0-79"
    assert sorted(checkpoints) == list(range(80)), f"Expected checkpoints 0-79; got {sorted(checkpoints)}"
    assert marker in stdout.read_text()
    assert receipt["training_command_exit_code"] == 0 and receipt["export_epoch79_marker_verified"] is True

    epoch79 = checkpoints[79]
    exported = run / "models/MACE_periodic_v19_residual_response.model"
    checkpoint_export = checkpoint_dir / f"{TAG}.model"
    assert exported.is_file() and checkpoint_export.is_file()
    torch.set_default_dtype(torch.float64)
    ckpt = torch.load(epoch79, map_location="cpu", weights_only=False)
    model = torch.load(exported, map_location="cpu", weights_only=False).double().eval()
    state = model.state_dict()
    assert set(ckpt["model"]) == set(state)
    differing = [key for key in state if not torch.equal(ckpt["model"][key].cpu(), state[key].cpu())]
    assert not differing, f"Exported model differs from epoch79 checkpoint tensors: {differing[:10]}"

    dev_path = ROOT / "data/valid.extxyz"
    dev = read(dev_path, index=":")
    assert len(dev) == 3
    audit = json.loads(AUDIT_JSON.read_text())
    calculator = MACECalculator(models=[model], device="cpu", default_dtype="float64")
    rows = []
    species_errors: dict[str, list[float]] = {}
    aggregate_sq = 0.0
    aggregate_atoms = 0
    separation_errors = []
    for atoms in dev:
        ref_e = float(atoms.info["REF_energy"])
        ref_f = np.asarray(atoms.arrays["REF_forces"], dtype=float)
        atoms.calc = calculator
        pred_e = float(atoms.get_potential_energy())
        pred_f = np.asarray(atoms.get_forces(), dtype=float)
        delta = pred_f - ref_f
        vector_norm = np.linalg.norm(delta, axis=1)
        vector_rmse = float(np.sqrt(np.mean(vector_norm**2)))
        component_rmse = float(np.sqrt(np.mean(delta**2)))
        aggregate_sq += float(np.sum(delta**2))
        aggregate_atoms += len(atoms)
        for element in sorted(set(atoms.get_chemical_symbols())):
            mask = np.asarray(atoms.get_chemical_symbols()) == element
            species_errors.setdefault(element, []).extend(vector_norm[mask].tolist())
        central = np.flatnonzero(np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool))
        assert len(central) == 2, f"Missing two-atom central-pair marker: {atoms.info.get('config_type')}"
        ia, ib = map(int, central)
        # The slab is periodic in x/y; use the minimum-image contact vector.
        vec = atoms.get_distance(ia, ib, mic=True, vector=True)
        label = atoms.info["config_type"]
        expected_contact = audit["candidate_entries"][label]["marked_contact"]["distance_A"]
        assert np.isclose(np.linalg.norm(vec), expected_contact, rtol=0.0, atol=1e-5)
        unit = vec / np.linalg.norm(vec)
        ref_sep = float(np.dot(ref_f[ib] - ref_f[ia], unit))
        pred_sep = float(np.dot(pred_f[ib] - pred_f[ia], unit))
        sep_error = abs(pred_sep - ref_sep)
        separation_errors.append(sep_error)
        rows.append({
            "config_type": atoms.info["config_type"],
            "interface": atoms.info["config_type"].split("_v19_", 1)[0],
            "atoms": len(atoms),
            "energy_error_eV_per_atom": (pred_e - ref_e) / len(atoms),
            "energy_abs_error_meV_per_atom": abs(pred_e - ref_e) * 1000.0 / len(atoms),
            "force_vector_RMSE_eV_A": vector_rmse,
            "force_component_RMSE_eV_A": component_rmse,
            "force_component_MAE_eV_A": float(np.mean(np.abs(delta))),
            "force_vector_max_error_eV_A": float(vector_norm.max()),
            "central_pair_elements": f"{atoms[ia].symbol}-{atoms[ib].symbol}",
            "central_pair_distance_A": float(np.linalg.norm(vec)),
            "central_pair_distance_matches_audit": True,
            "separation_force_reference_eV_A": ref_sep,
            "separation_force_prediction_eV_A": pred_sep,
            "separation_force_abs_error_eV_A": sep_error,
        })

    for element, errors in species_errors.items():
        species_errors[element] = float(np.sqrt(np.mean(np.square(errors))))
    total_force_rmse = float(np.sqrt(aggregate_sq / aggregate_atoms))
    force_threshold = 0.05
    separation_threshold = 0.10
    force_pass = all(row["force_vector_RMSE_eV_A"] <= force_threshold for row in rows)
    separation_pass = all(err <= separation_threshold for err in separation_errors)
    result = {
        "status": "TRAINING_COMPLETE_DEVELOPMENT_SCREEN_FAIL" if not (force_pass and separation_pass) else "TRAINING_COMPLETE_DEVELOPMENT_SCREEN_PASS",
        "verification": {
            "epoch_log_coverage": epoch_lines,
            "epoch_log_coverage_pass": True,
            "epoch_checkpoint_count": len(checkpoints),
            "epoch_checkpoint_coverage_pass": True,
            "epoch79_loaded_for_export_marker": True,
            "epoch79_checkpoint": repo_path(epoch79),
            "epoch79_checkpoint_sha256": digest(epoch79),
            "exported_model": repo_path(exported),
            "exported_model_sha256": digest(exported),
            "checkpoint_exported_model_sha256": digest(checkpoint_export),
            "epoch79_tensor_key_count": len(state),
            "epoch79_tensor_identity_pass": not differing,
            "launcher_completion": receipt,
        },
        "epoch_validation_metrics": epoch_metrics,
        "epoch_validation_metric_convention": {
            "RMSE_F": "MACE native per-force-component RMSE in meV/A",
            "project_force_vector_RMSE": "sqrt(3) times the MACE component RMSE for isotropic component aggregation; per-structure vector metrics are also reported below",
        },
        "evaluation_scope": {
            "dataset": "approved development validation only",
            "frames": len(dev),
            "test_or_sealed_data_opened": False,
            "model": repo_path(exported),
            "model_sha256": digest(exported),
            "marked_contact_distance_source": repo_path(AUDIT_JSON),
            "minimum_image_used_for_periodic_contact_axis": True,
        },
        "thresholds": {
            "per_interface_force_vector_RMSE_eV_A": {"max": force_threshold, "operator": "<="},
            "per_interface_separation_force_abs_error_eV_A": {"max": separation_threshold, "operator": "<="},
        },
        "aggregate_force_vector_RMSE_eV_A": total_force_rmse,
        "per_species_force_vector_RMSE_eV_A": species_errors,
        "per_interface": rows,
        "gates": {
            "force_vector_pass": force_pass,
            "separation_force_pass": separation_pass,
            "overall_pass": force_pass and separation_pass,
        },
        "limitations": [
            "The three development structures share one COD 9009647 Ti3SiC2 source family.",
            "This development screen is not blind test evidence or production validation.",
            "Slab k-point, smearing, and vacuum convergence remain unresolved.",
        ],
    }
    out = run / "results/v19_training_verification.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    curve = run / "results/v19_training_epoch_metrics.csv"
    curve.write_text(
        "epoch,validation_loss,validation_energy_RMSE_meV_per_atom,validation_force_component_RMSE_meV_A\n"
        + "".join(
            f"{row['epoch']},{row['validation_loss']},{row['validation_energy_RMSE_meV_per_atom']},{row['validation_force_component_RMSE_meV_A']}\n"
            for row in epoch_metrics
        )
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
