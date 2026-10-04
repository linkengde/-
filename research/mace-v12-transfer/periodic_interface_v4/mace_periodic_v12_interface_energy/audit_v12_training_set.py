#!/usr/bin/env python3
"""Re-evaluate the selected v12 model on each labeled training frame."""
from pathlib import Path
import csv
import hashlib
import json
import re

import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parent
MODEL = ROOT / "checkpoints/MACE_periodic_v12_interface_energy_run-45.model"
TRAIN = ROOT / "data/train.extxyz"
TRAIN_LOG = ROOT / "train.stdout"
OUT = ROOT / "results"

torch.set_default_dtype(torch.float64)
model = torch.load(MODEL, map_location="cpu", weights_only=False).double().eval()
calc = MACECalculator(models=[model], device="cpu", default_dtype="float64")
rows = []
for index, source in enumerate(read(TRAIN, index=":")):
    atoms = source.copy()
    atoms.calc = calc
    epred = float(atoms.get_potential_energy())
    fpred = np.asarray(atoms.get_forces(), dtype=float)
    row = {
        "frame": index,
        "config_type": source.info.get("config_type", ""),
        "atoms": len(source),
        "has_reference_energy": "REF_energy" in source.info,
        "has_reference_forces": "REF_forces" in source.arrays,
        "predicted_energy_eV_cell": epred,
    }
    if "REF_energy" in source.info:
        delta = float(source.info["REF_energy"]) - epred
        row["energy_error_DFT_minus_MACE_meV_atom"] = delta * 1000.0 / len(source)
    if "REF_forces" in source.arrays:
        delta_f = fpred - np.asarray(source.arrays["REF_forces"], dtype=float)
        row["force_vector_RMSE_eV_A"] = float(np.sqrt(np.mean(np.sum(delta_f**2, axis=1))))
    rows.append(row)

energy_rows = [r for r in rows if "energy_error_DFT_minus_MACE_meV_atom" in r]
force_rows = [r for r in rows if "force_vector_RMSE_eV_A" in r]
energy_errors = np.asarray([r["energy_error_DFT_minus_MACE_meV_atom"] for r in energy_rows])
atom_counts = np.asarray([r["atoms"] for r in energy_rows], dtype=float)
force_squared = sum(r["atoms"] * r["force_vector_RMSE_eV_A"] ** 2 for r in force_rows)
reported = None
match = re.search(r"\|\s*train_Default\s*\|\s*([0-9.]+)", TRAIN_LOG.read_text(errors="replace"))
if match:
    reported = float(match.group(1))

OUT.mkdir(exist_ok=True)
csv_path = OUT / "v12_training_frame_audit.csv"
with csv_path.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=sorted(set().union(*(r.keys() for r in rows))))
    writer.writeheader()
    writer.writerows(rows)

summary = {
    "model_sha256": hashlib.sha256(MODEL.read_bytes()).hexdigest(),
    "train_file_sha256": hashlib.sha256(TRAIN.read_bytes()).hexdigest(),
    "train_frames": len(rows),
    "energy_labeled_frames": len(energy_rows),
    "force_labeled_frames": len(force_rows),
    "standalone_energy_MAE_meV_atom": float(np.mean(np.abs(energy_errors))),
    "standalone_energy_RMSE_meV_atom_unweighted_by_atoms": float(np.sqrt(np.mean(energy_errors**2))),
    "standalone_energy_RMSE_meV_atom_atom_weighted": float(np.sqrt(np.sum(atom_counts * energy_errors**2) / np.sum(atom_counts))),
    "standalone_force_vector_RMSE_eV_A_atom_weighted": float(np.sqrt(force_squared / sum(r["atoms"] for r in force_rows))),
    "mace_training_stdout_reported_RMSE_E_meV_atom": reported,
    "assessment": "The standalone selected-checkpoint training-set RMSE does not reproduce the MACE training stdout summary; investigate before treating the model as validated.",
    "per_frame_csv": csv_path.name,
}
(OUT / "v12_training_frame_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
