#!/usr/bin/env python3
"""Compare v9-v13 on frozen and external frames excluded from v13 training."""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
MODELS = {
    "v9": PROJECT / "mace_periodic_v9_forceonly/models/Ag_Ti_Si_C_MACE_periodic_v9_forceonly.model",
    "v10": PROJECT / "mace_periodic_v10_energycal/checkpoints/MACE_periodic_v10_energycal_run-43.model",
    "v11": PROJECT / "mace_periodic_v11_energyfocus/checkpoints/MACE_periodic_v11_energyfocus_run-44.model",
    "v12": PROJECT / "mace_periodic_v12_interface_energy/checkpoints/MACE_periodic_v12_interface_energy_run-45.model",
    "v13": ROOT / "checkpoints/MACE_periodic_v13_interface_energy_run-45.model",
}
TEST = ROOT / "data/test.extxyz"
HOLDOUTS = {
    "AgSi_midpoint": PROJECT / "pbe_external_v8_motif_holdouts/AgSi/retry1/AgSi_cluster_PW_PBE.extxyz",
    "AgC_midpoint": PROJECT / "pbe_cluster_transfer_pw_20261002_run3_mixer_retry2/AgC_cluster_PW_PBE.extxyz",
}
if any(not path.is_file() for path in [*MODELS.values(), TEST, *HOLDOUTS.values()]):
    missing = [str(path) for path in [*MODELS.values(), TEST, *HOLDOUTS.values()] if not path.is_file()]
    raise SystemExit("Required model/data file missing:\n" + "\n".join(missing))

torch.set_default_dtype(torch.float64)
calculators = {}
for name, path in MODELS.items():
    model = torch.load(path, map_location="cpu", weights_only=False).double().eval()
    calculators[name] = MACECalculator(models=[model], device="cpu", default_dtype="float64")


def reference(atoms):
    energy_key = next((key for key in ("REF_energy", "PW_PBE_energy_eV", "energy") if key in atoms.info), None)
    force_key = next((key for key in ("REF_forces", "PW_PBE_forces", "forces") if key in atoms.arrays), None)
    energy = atoms.info.get(energy_key) if energy_key else None
    forces = atoms.arrays.get(force_key) if force_key else None
    results = getattr(atoms.calc, "results", {}) if atoms.calc is not None else {}
    if energy is None and "energy" in results:
        energy, energy_key = results["energy"], "SinglePointCalculator.energy"
    if forces is None and "forces" in results:
        forces, force_key = results["forces"], "SinglePointCalculator.forces"
    if energy is None or forces is None:
        raise SystemExit(f"Missing reference labels for {atoms.info.get('config_type')}")
    forces = np.asarray(forces, dtype=float)
    if forces.shape != (len(atoms), 3) or not np.isfinite(forces).all() or not np.isfinite(float(energy)):
        raise SystemExit(f"Invalid reference labels for {atoms.info.get('config_type')}")
    return float(energy), forces, energy_key, force_key


def pair_force(atoms, forces):
    central = np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool)
    indices = np.flatnonzero(central)
    if len(indices) != 2:
        return {}
    ia, ib = map(int, indices)
    symbols = atoms.get_chemical_symbols()
    if symbols[ia] == symbols[ib]:
        return {}
    vector = atoms.positions[ib] - atoms.positions[ia]
    distance = float(np.linalg.norm(vector))
    return {
        "central_pair_elements": f"{symbols[ia]}-{symbols[ib]}",
        "central_pair_distance_A": distance,
        "separating_force_eV_A": float(np.dot(forces[ib] - forces[ia], vector / distance)),
    }


def evaluate(source, label, split):
    eref, fref, energy_key, force_key = reference(source)
    atoms = source.copy()
    row = {
        "split": split,
        "label": label,
        "config_type": atoms.info.get("config_type", ""),
        "formula": atoms.get_chemical_formula(),
        "atoms": len(atoms),
        "DFT_energy_eV_cell": eref,
        "DFT_energy_key": energy_key,
        "DFT_force_key": force_key,
    }
    row.update({f"DFT_{key}": value for key, value in pair_force(atoms, fref).items()})
    for name, calc in calculators.items():
        atoms.calc = calc
        epred = float(atoms.get_potential_energy())
        fpred = np.asarray(atoms.get_forces(), dtype=float)
        delta = fpred - fref
        row[f"{name}_energy_error_DFT_minus_MACE_meV_atom"] = (eref - epred) * 1000.0 / len(atoms)
        row[f"{name}_force_vector_RMSE_eV_A"] = float(np.sqrt(np.mean(np.sum(delta**2, axis=1))))
        row[f"{name}_force_component_MAE_eV_A"] = float(np.mean(np.abs(delta)))
        row[f"{name}_force_vector_max_error_eV_A"] = float(np.linalg.norm(delta, axis=1).max())
        row[f"{name}_Fmax_eV_A"] = float(np.linalg.norm(fpred, axis=1).max())
        row.update({f"{name}_{key}": value for key, value in pair_force(atoms, fpred).items() if key not in {"central_pair_elements", "central_pair_distance_A"}})
    return row


rows = []
for atoms in read(TEST, index=":"):
    rows.append(evaluate(atoms, atoms.info.get("config_type", ""), "frozen_test"))
for label, path in HOLDOUTS.items():
    rows.append(evaluate(read(path), label, "external_holdout"))

summary = {
    "evaluation_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "models": {name: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for name, path in MODELS.items()},
    "test_file_sha256": hashlib.sha256(TEST.read_bytes()).hexdigest(),
    "holdout_file_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in HOLDOUTS.items()},
    "per_structure": rows,
    "scope_note": "The three distance points added to v13 training are excluded from this score. The frozen test and two midpoint structures remain outside the v13 training file; the midpoint geometries are related to the new training motifs and are not sufficient by themselves for production validation.",
}
for split in ("frozen_test", "external_holdout"):
    split_rows = [row for row in rows if row["split"] == split]
    summary[split] = {"structures": len(split_rows), "atoms": sum(row["atoms"] for row in split_rows), "models": {}}
    for name in MODELS:
        energies = [abs(row[f"{name}_energy_error_DFT_minus_MACE_meV_atom"]) for row in split_rows]
        force_sq = sum(row["atoms"] * row[f"{name}_force_vector_RMSE_eV_A"] ** 2 for row in split_rows)
        summary[split]["models"][name] = {
            "energy_MAE_meV_atom": float(np.mean(energies)),
            "energy_max_abs_error_meV_atom": float(max(energies)),
            "force_vector_RMSE_eV_A": float(np.sqrt(force_sq / sum(row["atoms"] for row in split_rows))),
        }

out = ROOT / "results"
out.mkdir(exist_ok=True)
(out / "v9_v10_v11_v12_v13_comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
with (out / "v9_v10_v11_v12_v13_comparison.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=sorted(set().union(*(row.keys() for row in rows))), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(summary, indent=2))
