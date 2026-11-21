#!/usr/bin/env python3
"""Refit the B1616 ACE candidate with audited v24 and center-only v11 labels."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read
from pyace import ACEBBasisSet, PyACECalculator, create_multispecies_basis_config
from pyace.linearacefit import LinearACEFit, LinearACEDataset
from sklearn.linear_model import Ridge

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "fit_candidate_acquisition"
BASIS_CONFIG = Path("/workspace/potential-comparison/mace_fit_data_20261118/ace_mace_cutoff6/basis_config.json")
BASELINE_MODEL = Path("/workspace/-/archive/mace-fit-ace-exploratory-20261119/followup_20261119/ACE/AgCSiTi_B1616_E10000_F1_a0.001_Ti_gap_weight1.ace.yaml")
OUT = ROOT / "ace_acquisition_refits"
RCUT = 6.0
ENERGY_WEIGHT = 10000.0
ALPHA = 0.001
HOTSPOT_FORCE_SCALES = [1.0, 10.0, 100.0]
ENERGY_GATE_MEV_ATOM = 10.0
FORCE_COMPONENT_GATE = 0.2
HOTSPOT_CENTER_GATE = 0.5


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_dataset(path: Path) -> pd.DataFrame:
    rows = []
    for source_atoms in read(path, ":"):
        energy = float(source_atoms.info["REF_energy"])
        forces = np.asarray(source_atoms.arrays["REF_forces"], dtype=float).copy()
        if forces.shape != (len(source_atoms), 3) or not np.isfinite(forces).all():
            raise ValueError(f"invalid reference forces: {source_atoms.info.get('config_type')}")
        atoms = source_atoms.copy()
        if atoms.pbc.any() and not atoms.pbc.all():
            old_distances = atoms.get_all_distances(mic=True)
            cell = np.asarray(atoms.cell, dtype=float).copy()
            for axis in range(3):
                if not atoms.pbc[axis]:
                    height = float(np.linalg.norm(cell[axis]))
                    span = float(np.ptp(atoms.get_scaled_positions(wrap=False)[:, axis]) * height)
                    if height - span < RCUT + 2.0:
                        cell[axis] *= (span + RCUT + 2.0) / height
            atoms.set_cell(cell, scale_atoms=False)
            atoms.pbc = True
            new_distances = atoms.get_all_distances(mic=True)
            mask = (old_distances < RCUT + 1e-8) | (new_distances < RCUT + 1e-8)
            if not np.allclose(old_distances[mask], new_distances[mask], atol=1e-8):
                raise ValueError(f"vacuum padding changed cutoff neighbors: {atoms.info.get('config_type')}")
        rows.append({
            "ase_atoms": atoms,
            "energy": energy,
            "energy_corrected": energy,
            "forces": forces,
            "NUMBER_OF_ATOMS": len(atoms),
            "config_type": str(atoms.info.get("config_type", "unknown")),
            "dataset_role": str(atoms.info.get("dataset_role", "existing_candidate")),
            "training_label_policy": str(atoms.info.get("training_label_policy", "all_atom_energy_force")),
        })
    return pd.DataFrame(rows)


def evaluate_full_label_model(bbasis, dataframe: pd.DataFrame) -> dict:
    calculator = PyACECalculator(bbasis)
    energy_errors = []
    force_errors = []
    per_frame = []
    for _, row in dataframe.iterrows():
        atoms = row["ase_atoms"].copy()
        atoms.calc = calculator
        delta_energy = (float(atoms.get_potential_energy()) - float(row["energy"])) / len(atoms)
        delta_forces = np.asarray(atoms.get_forces(), dtype=float) - np.asarray(row["forces"], dtype=float)
        energy_errors.append(delta_energy)
        force_errors.extend(delta_forces.ravel())
        per_frame.append({
            "config_type": row["config_type"],
            "natoms": len(atoms),
            "energy_error_meV_atom": float(1000 * delta_energy),
            "force_component_rmse_eV_A": float(np.sqrt(np.mean(delta_forces ** 2))),
        })
    return {
        "energy_rmse_meV_atom": float(1000 * np.sqrt(np.mean(np.square(energy_errors)))),
        "force_component_rmse_eV_A": float(np.sqrt(np.mean(np.square(force_errors)))),
        "worst_frame_force_component_rmse_eV_A": float(max(x["force_component_rmse_eV_A"] for x in per_frame)),
        "per_configuration": per_frame,
    }


def evaluate_hotspot_center(bbasis, path: Path) -> dict:
    atoms = read(path)
    center_id = int(atoms.info["center_atom_id"])
    source_ids = np.asarray(atoms.arrays["source_atom_id"], dtype=int)
    center = np.flatnonzero(source_ids == center_id)
    if len(center) != 1:
        raise ValueError(f"center atom is ambiguous in {path.name}")
    center_idx = int(center[0])
    if atoms.get_chemical_symbols()[center_idx] != "Si":
        raise ValueError(f"registered center is not Si in {path.name}")
    reference = np.asarray(atoms.arrays["REF_forces"][center_idx], dtype=float)
    atoms.calc = PyACECalculator(bbasis)
    predicted = np.asarray(atoms.get_forces()[center_idx], dtype=float)
    error = float(np.linalg.norm(predicted - reference))
    return {
        "config_type": atoms.info.get("config_type"),
        "center_atom_id": center_id,
        "reference_force_eV_A": reference.tolist(),
        "predicted_force_eV_A": predicted.tolist(),
        "center_vector_error_eV_A": error,
        "gate_eV_A": HOTSPOT_CENTER_GATE,
        "pass": error <= HOTSPOT_CENTER_GATE,
    }


def fit_weights(train: LinearACEDataset, hotspot_scale: float) -> np.ndarray:
    energy_rows = len(train.df)
    weights = np.ones(len(train.get_target_vector()), dtype=float)
    weights[:energy_rows] = ENERGY_WEIGHT
    force_start = energy_rows
    for _, row in train.df.iterrows():
        atoms = row["ase_atoms"]
        config_type = str(atoms.info.get("config_type", ""))
        if config_type.startswith("POD_v11_initial_Si_hotspot_"):
            weights[force_start:force_start + 3 * len(atoms)] = 0.0
            center_id = int(atoms.info["force_target_center_atom_id"])
            source_ids = np.asarray(atoms.arrays["source_atom_id"], dtype=int)
            center = np.flatnonzero(source_ids == center_id)
            if len(center) != 1:
                raise ValueError(f"registered center is ambiguous: {config_type}")
            center_rows = force_start + 3 * int(center[0])
            weights[center_rows:center_rows + 3] = hotspot_scale
            weights[energy_rows + train.df.index.get_loc(row.name)] = 0.0
        force_start += 3 * len(atoms)
    if force_start != len(weights):
        raise RuntimeError(f"weight-vector length mismatch: expected {force_start}, got {len(weights)}")
    return weights


def main() -> None:
    manifest = json.loads((DATA / "candidate_manifest.json").read_text())
    for name, digest in manifest["files"].items():
        path = DATA / name
        if sha(path) != digest:
            raise SystemExit(f"prepared acquisition dataset hash mismatch: {name}")
    if not BASIS_CONFIG.is_file() or not BASELINE_MODEL.is_file():
        raise SystemExit("registered B1616 basis or baseline model is missing")

    training = load_dataset(DATA / "train_acquisition.extxyz")
    periodic_development = load_dataset(DATA / "dev_periodic.extxyz")
    hotspot_development = DATA / "dev_hotspot01_center_force.extxyz"
    basis_config = json.loads(BASIS_CONFIG.read_text())
    basis = create_multispecies_basis_config(potential_config=basis_config)
    train = LinearACEDataset(basis, training)
    development = LinearACEDataset(basis, periodic_development)
    train.get_design_matrix(max_workers=2, verbose=True)
    development.get_design_matrix(max_workers=2, verbose=True)

    OUT.mkdir(parents=True, exist_ok=True)
    baseline = ACEBBasisSet(str(BASELINE_MODEL))
    baseline_result = {
        "model_sha256": sha(BASELINE_MODEL),
        "full_label_training": evaluate_full_label_model(baseline, training[~training["config_type"].str.startswith("POD_v11_initial_Si_hotspot_")]),
        "periodic_development": evaluate_full_label_model(baseline, periodic_development),
        "hotspot_center_development": evaluate_hotspot_center(baseline, hotspot_development),
        "status": "baseline reference; not retrained with v11 center labels",
    }
    (OUT / "baseline_metrics.json").write_text(json.dumps(baseline_result, indent=2) + "\n")

    results = []
    for scale in HOTSPOT_FORCE_SCALES:
        weights = fit_weights(train, scale)
        fit = LinearACEFit(train_dataset=train, model=Ridge(alpha=ALPHA))
        fit.fit(sample_weight=weights)
        bbasis = fit.get_bbasis()
        model_name = f"AgCSiTi_B1616_E10000_F1_a{ALPHA:g}_hotspot_center_weight{scale:g}.ace.yaml"
        model_path = OUT / model_name
        bbasis.save(str(model_path))
        full_train = training[~training["config_type"].str.startswith("POD_v11_initial_Si_hotspot_")]
        metrics = {
            "model": model_name,
            "model_sha256": sha(model_path),
            "basis_config_sha256": sha(BASIS_CONFIG),
            "ridge_alpha": ALPHA,
            "energy_weight_full_label_frames": ENERGY_WEIGHT,
            "full_force_weight": 1.0,
            "center_Si_force_weight": scale,
            "crop_energy_weight": 0.0,
            "crop_noncenter_force_weight": 0.0,
            "full_label_train_frames": len(full_train),
            "center_only_train_frames": 3,
            "full_label_training": evaluate_full_label_model(bbasis, full_train),
            "periodic_development": evaluate_full_label_model(bbasis, periodic_development),
            "hotspot_center_training": [
                evaluate_hotspot_center(bbasis, p)
                for p in sorted((ROOT / "v11_hotspot_reference/archive/round7-v11-si-hotspot-labels-20261118/round7_pod_snap/models/v11_candidate/active_si_hotspots/labels").glob("si_hotspot_0[023]_pw350_gamma.extxyz"))
            ],
            "hotspot_center_development": evaluate_hotspot_center(bbasis, hotspot_development),
            "thresholds": {
                "energy_rmse_meV_atom": ENERGY_GATE_MEV_ATOM,
                "force_component_rmse_eV_A": FORCE_COMPONENT_GATE,
                "hotspot_center_vector_error_eV_A": HOTSPOT_CENTER_GATE,
            },
            "gate_pass": False,
            "status": "exploratory; development label is not a blind test",
        }
        train_center_pass = all(
            item["center_vector_error_eV_A"] <= HOTSPOT_CENTER_GATE
            for item in metrics["hotspot_center_training"]
        )
        metrics["gate_components"] = {
            "train_energy_pass": metrics["full_label_training"]["energy_rmse_meV_atom"] <= ENERGY_GATE_MEV_ATOM,
            "train_force_component_pass": metrics["full_label_training"]["force_component_rmse_eV_A"] <= FORCE_COMPONENT_GATE,
            "development_energy_pass": metrics["periodic_development"]["energy_rmse_meV_atom"] <= ENERGY_GATE_MEV_ATOM,
            "development_force_component_pass": metrics["periodic_development"]["force_component_rmse_eV_A"] <= FORCE_COMPONENT_GATE,
            "hotspot_train_centers_pass": train_center_pass,
            "hotspot_development_center_pass": metrics["hotspot_center_development"]["center_vector_error_eV_A"] <= HOTSPOT_CENTER_GATE,
        }
        metrics["gate_pass"] = all(metrics["gate_components"].values())
        metric_path = OUT / f"metrics_hotspot_center_weight{scale:g}.json"
        metric_path.write_text(json.dumps(metrics, indent=2) + "\n")
        results.append({"model": model_name, "metrics": metric_path.name,
                        "gate_pass": metrics["gate_pass"],
                        "energy_rmse_meV_atom": metrics["periodic_development"]["energy_rmse_meV_atom"],
                        "force_component_rmse_eV_A": metrics["periodic_development"]["force_component_rmse_eV_A"],
                        "hotspot_center_vector_error_eV_A": metrics["hotspot_center_development"]["center_vector_error_eV_A"]})
        print(json.dumps(results[-1], indent=2), flush=True)

    summary = {
        "status": "candidate fits completed; exploratory only",
        "train_manifest_sha256": sha(DATA / "candidate_manifest.json"),
        "baseline": baseline_result,
        "models": results,
        "passing_models": [x["model"] for x in results if x["gate_pass"]],
        "blind_test_opened": False,
    }
    (OUT / "fit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"status": summary["status"], "passing_models": summary["passing_models"]}, indent=2))


if __name__ == "__main__":
    main()
