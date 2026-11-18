#!/usr/bin/env python3
"""Refit ACE against the 1x6x3 Ti-gap label at two force-weight scales."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read
from pyace import PyACECalculator, create_multispecies_basis_config
from pyace.linearacefit import LinearACEFit, LinearACEDataset
from sklearn.linear_model import Ridge


HERE = Path(__file__).resolve().parent
DATA = HERE / "fit_candidate"
FRAME = "Ti_gap_2p3_rank_recovery_v17_01"
RCUT = 6.0


def load_dataset(path):
    rows = []
    for atoms in read(path, ":"):
        energy = float(atoms.info["REF_energy"])
        forces = np.asarray(atoms.arrays["REF_forces"], dtype=float).copy()
        atoms = atoms.copy()
        if atoms.pbc.any() and not atoms.pbc.all():
            old_distances = atoms.get_all_distances(mic=True)
            cell = np.asarray(atoms.cell).copy()
            for axis in range(3):
                if not atoms.pbc[axis]:
                    height = float(np.linalg.norm(cell[axis]))
                    span = float(np.ptp(atoms.get_scaled_positions(wrap=False)[:, axis]) * height)
                    if height - span < RCUT + 2:
                        cell[axis] *= (span + RCUT + 2) / height
            atoms.set_cell(cell, scale_atoms=False)
            atoms.pbc = True
            new_distances = atoms.get_all_distances(mic=True)
            mask = (old_distances < RCUT + 1e-8) | (new_distances < RCUT + 1e-8)
            if not np.allclose(old_distances[mask], new_distances[mask], atol=1e-8):
                raise ValueError(f"Vacuum padding changed cutoff neighbors: {path}")
        rows.append({"ase_atoms": atoms, "energy": energy, "energy_corrected": energy,
                     "forces": forces, "NUMBER_OF_ATOMS": len(atoms)})
    return pd.DataFrame(rows)


def evaluate(bbasis, dataframe):
    calculator = PyACECalculator(bbasis)
    energy_errors, force_errors, si_site_errors, per_frame = [], [], [], []
    for _, row in dataframe.iterrows():
        atoms = row["ase_atoms"].copy()
        atoms.calc = calculator
        delta_energy = (float(atoms.get_potential_energy()) - float(row["energy"])) / len(atoms)
        delta_forces = atoms.get_forces() - np.asarray(row["forces"], dtype=float)
        energy_errors.append(delta_energy)
        force_errors.extend(delta_forces.ravel())
        si_mask = np.asarray(atoms.get_chemical_symbols()) == "Si"
        if si_mask.any():
            si_site_errors.extend(np.linalg.norm(delta_forces[si_mask], axis=1))
        per_frame.append({
            "config_type": atoms.info.get("config_type"),
            "natoms": len(atoms),
            "energy_error_meV_atom": 1000 * float(delta_energy),
            "force_component_rmse_eV_A": float(np.sqrt(np.mean(delta_forces ** 2))),
            "si_site_vector_error_max_eV_A": float(np.linalg.norm(delta_forces[si_mask], axis=1).max()) if si_mask.any() else None,
        })
    return {
        "energy_rmse_meV_atom": float(1000 * np.sqrt(np.mean(np.square(energy_errors)))),
        "force_component_rmse_eV_A": float(np.sqrt(np.mean(np.square(force_errors)))),
        "si_site_vector_error_max_eV_A_proxy": float(max(si_site_errors)) if si_site_errors else None,
        "si_center_metric_status": "proxy over all Si sites; registered hotspot-center IDs are absent",
        "worst_frame_force_component_rmse_eV_A": float(max(x["force_component_rmse_eV_A"] for x in per_frame)),
        "per_configuration": per_frame,
    }


def main():
    basis_config = json.loads((HERE.parent / "ace_mace_cutoff6" / "basis_config.json").read_text())
    basis = create_multispecies_basis_config(potential_config=basis_config)
    train = LinearACEDataset(basis, load_dataset(DATA / "train_kconv.extxyz"))
    development = LinearACEDataset(basis, load_dataset(DATA / "dev.extxyz"))
    train.get_design_matrix(max_workers=2, verbose=True)
    development.get_design_matrix(max_workers=2, verbose=True)
    energy_rows = len(train.df)
    frame_indices = [i for i, row in train.df.iterrows()
                     if row["ase_atoms"].info.get("config_type") == FRAME]
    if len(frame_indices) != 1:
        raise ValueError(f"Expected one {FRAME}, found {len(frame_indices)}")
    frame_index = frame_indices[0]
    force_start = energy_rows + sum(3 * len(train.df.iloc[i]["ase_atoms"]) for i in range(frame_index))
    force_rows = 3 * len(train.df.iloc[frame_index]["ase_atoms"])
    output = HERE / "ace_kconv_refits"
    output.mkdir(exist_ok=True)
    for scale in (1.0, 0.05):
        weights = np.ones(len(train.get_target_vector()))
        weights[:energy_rows] = 10000.0
        weights[force_start:force_start + force_rows] *= scale
        fit = LinearACEFit(train_dataset=train, model=Ridge(alpha=0.001))
        fit.fit(sample_weight=weights)
        model_name = f"AgCSiTi_B1616_E10000_F1_a0.001_Ti_gap_weight{scale:g}.ace.yaml"
        fit.get_bbasis().save(str(output / model_name))
        result = {
            "label_mesh": "1x6x3 GPAW PBE/PW500; see DFT comparison.json",
            "force_weight_scale_for_Ti_gap_frame": scale,
            "energy_weight": 10000.0,
            "force_weight_other_frames": 1.0,
            "ridge_alpha": 0.001,
            "model": model_name,
            "train": evaluate(fit.get_bbasis(), train.df),
            "development": evaluate(fit.get_bbasis(), development.df),
            "status": "exploratory candidate; evaluate only after DFT convergence and label review",
        }
        (output / f"metrics_weight{scale:g}.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
