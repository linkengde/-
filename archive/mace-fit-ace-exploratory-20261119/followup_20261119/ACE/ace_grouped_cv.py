#!/usr/bin/env python3
"""Exploratory grouped ACE validation; folds keep related source families together."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read
from pyace import PyACECalculator, create_multispecies_basis_config
from pyace.linearacefit import LinearACEFit, LinearACEDataset
from sklearn.linear_model import Ridge


HERE = Path(__file__).resolve().parent
DATA = HERE / "fit_candidate" / "train_kconv.extxyz"
FRAME = "Ti_gap_2p3_rank_recovery_v17_01"
RCUT = 6.0


def group_name(atoms):
    name = str(atoms.info.get("config_type", "unknown"))
    lower = name.lower()
    if atoms.info.get("source_family"):
        return "JN1967_source_family"
    if "v17" in lower or "ti_gap" in lower:
        return "v17_framework_and_ti_gap"
    if "v21" in lower or "v22" in lower or "v23" in lower or "v24" in lower:
        return "JN1967_source_family"
    if "v25" in lower or "ti3sic2" in lower or "bulk" in lower:
        return "Ti3SiC2_bulk_v25"
    for pair in ("AgC", "AgSi", "AgTi"):
        if pair.lower() in lower:
            return f"MACE_{pair}_interface_family"
    raise ValueError(f"Cannot assign a source family: {name}")


def load_dataset(atoms_list):
    rows = []
    for source in atoms_list:
        energy = float(source.info["REF_energy"])
        forces = np.asarray(source.arrays["REF_forces"], dtype=float).copy()
        atoms = source.copy()
        if atoms.pbc.any() and not atoms.pbc.all():
            old = atoms.get_all_distances(mic=True)
            cell = np.asarray(atoms.cell).copy()
            for axis in range(3):
                if not atoms.pbc[axis]:
                    height = float(np.linalg.norm(cell[axis]))
                    span = float(np.ptp(atoms.get_scaled_positions(wrap=False)[:, axis]) * height)
                    if height - span < RCUT + 2:
                        cell[axis] *= (span + RCUT + 2) / height
            atoms.set_cell(cell, scale_atoms=False)
            atoms.pbc = True
            new = atoms.get_all_distances(mic=True)
            mask = (old < RCUT + 1e-8) | (new < RCUT + 1e-8)
            if not np.allclose(old[mask], new[mask], atol=1e-8):
                raise ValueError(f"Vacuum padding changed cutoff neighbors: {atoms.info.get('config_type')}")
        rows.append({"ase_atoms": atoms, "energy": energy, "energy_corrected": energy,
                     "forces": forces, "NUMBER_OF_ATOMS": len(atoms)})
    return pd.DataFrame(rows)


def evaluate(bbasis, dataframe):
    calculator = PyACECalculator(bbasis)
    energy_errors, force_errors, si_errors, per = [], [], [], []
    for _, row in dataframe.iterrows():
        atoms = row["ase_atoms"].copy()
        atoms.calc = calculator
        de = (float(atoms.get_potential_energy()) - float(row["energy"])) / len(atoms)
        df = atoms.get_forces() - np.asarray(row["forces"], dtype=float)
        energy_errors.append(de)
        force_errors.extend(df.ravel())
        si = np.asarray(atoms.get_chemical_symbols()) == "Si"
        if si.any():
            si_errors.extend(np.linalg.norm(df[si], axis=1))
        per.append({"config_type": atoms.info.get("config_type"),
                    "natoms": len(atoms),
                    "energy_error_meV_atom": 1000 * float(de),
                    "force_component_rmse_eV_A": float(np.sqrt(np.mean(df ** 2))),
                    "si_site_vector_error_max_eV_A": float(np.linalg.norm(df[si], axis=1).max()) if si.any() else None})
    return {"energy_rmse_meV_atom": float(1000 * np.sqrt(np.mean(np.square(energy_errors)))),
            "force_component_rmse_eV_A": float(np.sqrt(np.mean(np.square(force_errors)))),
            "si_site_vector_error_max_eV_A_proxy": float(max(si_errors)) if si_errors else None,
            "per_configuration": per}


def main():
    all_atoms = read(DATA, ":")
    groups = [group_name(a) for a in all_atoms]
    unique = sorted(set(groups))
    if len(unique) < 3:
        raise ValueError(f"Need at least three groups; got {unique}")
    config_path = HERE.parent / "ace_mace_cutoff6" / "basis_config.json"
    basis = create_multispecies_basis_config(potential_config=json.loads(config_path.read_text()))
    output = HERE / "ace_grouped_cv"
    output.mkdir(exist_ok=True)
    results = {"status": "exploratory grouped validation; folds are not blind test",
               "group_definition": "JN1967 source family; v17 perturbation family; Ti3SiC2 bulk; otherwise pair-family AgC/AgSi/AgTi",
               "groups": {g: sum(x == g for x in groups) for g in unique}, "folds": {}}
    for heldout in unique:
        keep = [i for i, g in enumerate(groups) if g != heldout]
        test = [i for i, g in enumerate(groups) if g == heldout]
        train_ds = LinearACEDataset(basis, load_dataset([all_atoms[i] for i in keep]))
        valid_ds = LinearACEDataset(basis, load_dataset([all_atoms[i] for i in test]))
        train_ds.get_design_matrix(max_workers=2, verbose=True)
        valid_ds.get_design_matrix(max_workers=2, verbose=True)
        targets = train_ds.get_target_vector()
        energy_rows = len(train_ds.df)
        weights = np.ones(len(targets))
        weights[:energy_rows] = 10000.0
        fit = LinearACEFit(train_dataset=train_ds, model=Ridge(alpha=0.001))
        fit.fit(sample_weight=weights)
        results["folds"][heldout] = {
            "training_frames": len(keep), "heldout_frames": len(test),
            "train": evaluate(fit.get_bbasis(), train_ds.df),
            "heldout": evaluate(fit.get_bbasis(), valid_ds.df),
        }
        (output / f"ACE_without_{heldout}.ace.yaml").parent.mkdir(exist_ok=True)
        fit.get_bbasis().save(str(output / f"ACE_without_{heldout}.ace.yaml"))
        (output / "metrics.json").write_text(json.dumps(results, indent=2) + "\n")
        print(json.dumps({heldout: results["folds"][heldout]}, indent=2), flush=True)
    (output / "metrics.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps({"groups": results["groups"], "folds": list(results["folds"])}, indent=2))


if __name__ == "__main__":
    main()
