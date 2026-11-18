#!/usr/bin/env python3
"""Reproduce the ACE force-weight sensitivity fit for the V17 Ti-gap row."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read
from pyace import PyACECalculator, create_multispecies_basis_config
from pyace.linearacefit import LinearACEFit, LinearACEDataset
from sklearn.linear_model import Ridge

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "sensitivity_Ti_gap_weight"
FRAME_TYPE = "Ti_gap_2p3_rank_recovery_v17_01"
RCUT = 6.0


def load_dataset(path):
    rows = []
    for atoms in read(path, ":"):
        energy = float(atoms.info["REF_energy"])
        forces = np.asarray(atoms.arrays["REF_forces"], dtype=float).copy()
        atoms = atoms.copy()
        if atoms.pbc.any() and not atoms.pbc.all():
            old = atoms.get_all_distances(mic=True)
            cell = np.asarray(atoms.cell).copy()
            for axis in range(3):
                if not atoms.pbc[axis]:
                    h = float(np.linalg.norm(cell[axis]))
                    span = float(np.ptp(atoms.get_scaled_positions(wrap=False)[:, axis]) * h)
                    if h - span < RCUT + 2:
                        cell[axis] *= (span + RCUT + 2) / h
            atoms.set_cell(cell, scale_atoms=False)
            atoms.pbc = True
            new = atoms.get_all_distances(mic=True)
            mask = (old < RCUT + 1e-8) | (new < RCUT + 1e-8)
            if not np.allclose(old[mask], new[mask], atol=1e-8):
                raise ValueError(f"Vacuum padding changed cutoff neighbors: {path}")
        rows.append({
            "ase_atoms": atoms,
            "energy": energy,
            "energy_corrected": energy,
            "forces": forces,
            "NUMBER_OF_ATOMS": len(atoms),
        })
    return pd.DataFrame(rows)


def evaluate(bbasis, dataset):
    calculator = PyACECalculator(bbasis)
    energy_errors = []
    force_errors = []
    si_vector_errors = []
    per_configuration = []
    for _, row in dataset.iterrows():
        atoms = row["ase_atoms"].copy()
        atoms.calc = calculator
        predicted_energy = float(atoms.get_potential_energy())
        predicted_forces = atoms.get_forces()
        delta = predicted_forces - np.asarray(row["forces"], dtype=float)
        energy_error = (predicted_energy - float(row["energy"])) / len(atoms)
        energy_errors.append(energy_error)
        force_errors.extend(delta.ravel())
        symbols = np.asarray(atoms.get_chemical_symbols())
        si_vector_errors.extend(np.linalg.norm(delta[symbols == "Si"], axis=1))
        per_configuration.append({
            "config_type": atoms.info.get("config_type"),
            "energy_error_eV_atom": float(energy_error),
            "force_component_rmse_eV_A": float(np.sqrt(np.mean(delta**2))),
        })
    return {
        "energy_rmse_meV_atom": float(1000 * np.sqrt(np.mean(np.square(energy_errors)))),
        "force_component_rmse_eV_A": float(np.sqrt(np.mean(np.square(force_errors)))),
        "si_site_vector_error_max_eV_A": float(max(si_vector_errors)),
        "per_configuration": per_configuration,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--force-weight-scale", type=float, default=0.05)
    args = parser.parse_args()
    if args.force_weight_scale < 0:
        raise ValueError("force-weight scale must be nonnegative")

    basis_config = json.loads((ROOT / "ace_mace_cutoff6/basis_config.json").read_text())
    basis = create_multispecies_basis_config(potential_config=basis_config)
    train = LinearACEDataset(basis, load_dataset(ROOT / "train.extxyz"))
    development = LinearACEDataset(basis, load_dataset(ROOT / "dev.extxyz"))
    train.get_design_matrix(max_workers=2, verbose=True)
    development.get_design_matrix(max_workers=2, verbose=True)

    targets = train.get_target_vector()
    energy_rows = len(train.df)
    frame_indices = [
        i for i, row in train.df.iterrows()
        if row["ase_atoms"].info.get("config_type") == FRAME_TYPE
    ]
    if len(frame_indices) != 1:
        raise ValueError(f"Expected exactly one {FRAME_TYPE}, found {len(frame_indices)}")
    frame_index = frame_indices[0]
    force_start = energy_rows + sum(
        3 * len(train.df.iloc[i]["ase_atoms"]) for i in range(frame_index)
    )
    force_rows = 3 * len(train.df.iloc[frame_index]["ase_atoms"])
    weights = np.ones(len(targets))
    weights[:energy_rows] = 10000.0
    weights[force_start:force_start + force_rows] *= args.force_weight_scale

    fit = LinearACEFit(train_dataset=train, model=Ridge(alpha=0.001))
    fit.fit(sample_weight=weights)
    OUT.mkdir(exist_ok=True)
    suffix = f"{args.force_weight_scale:g}"
    model_name = f"AgCSiTi_B1616_E10000_F1_a0.001_Ti_gap_weight{suffix}.ace.yaml"
    model_path = OUT / model_name
    fit.get_bbasis().save(str(model_path))
    results = {
        "frame_type": FRAME_TYPE,
        "force_weight_scale": args.force_weight_scale,
        "energy_weight": 10000.0,
        "force_weight_other_frames": 1.0,
        "ridge_alpha": 0.001,
        "model": model_name,
        "train": evaluate(fit.get_bbasis(), train.df),
        "development": evaluate(fit.get_bbasis(), development.df),
        "status": "exploratory sensitivity candidate; k-point convergence of selected frame unresolved",
    }
    (OUT / f"metrics_weight{suffix}.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
