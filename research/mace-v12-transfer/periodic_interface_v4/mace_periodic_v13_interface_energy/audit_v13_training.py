#!/usr/bin/env python3
"""Audit v13 checkpoint errors on the training frames and mask missing energies."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parent
MODEL = ROOT / "checkpoints/MACE_periodic_v13_interface_energy_run-45.model"
TRAIN = ROOT / "data/train.extxyz"
torch.set_default_dtype(torch.float64)
model = torch.load(MODEL, map_location="cpu", weights_only=False).double().eval()
calculator = MACECalculator(models=[model], device="cpu", default_dtype="float64")
frames = read(TRAIN, index=":")
rows = []
energy_errors = []
energy_errors_atom_weighted = []
zero_filled_errors = []
force_squared_atom_sum = 0.0
force_atom_count = 0

for index, source in enumerate(frames):
    atoms = source.copy()
    atoms.calc = calculator
    predicted_energy = float(atoms.get_potential_energy())
    predicted_forces = np.asarray(atoms.get_forces(), dtype=float)
    has_energy = "REF_energy" in source.info
    has_forces = "REF_forces" in source.arrays
    energy_error_meV_atom = None
    force_vector_rmse = None
    if has_energy:
        energy_error_meV_atom = (float(source.info["REF_energy"]) - predicted_energy) * 1000.0 / len(source)
        energy_errors.append(energy_error_meV_atom)
        energy_errors_atom_weighted.append((len(source), energy_error_meV_atom))
        zero_filled_errors.append(energy_error_meV_atom)
    else:
        # This intentionally reproduces MACE's unmasked summary if missing labels are zero-filled.
        zero_filled_errors.append(-predicted_energy * 1000.0 / len(source))
    if has_forces:
        force_delta = predicted_forces - np.asarray(source.arrays["REF_forces"], dtype=float)
        force_vector_rmse = float(np.sqrt(np.mean(np.sum(force_delta**2, axis=1))))
        force_squared_atom_sum += float(np.sum(force_delta**2))
        force_atom_count += len(source)
    rows.append({
        "frame": index,
        "config_type": source.info.get("config_type", ""),
        "atoms": len(source),
        "has_reference_energy": has_energy,
        "has_reference_forces": has_forces,
        "energy_error_DFT_minus_MACE_meV_atom": energy_error_meV_atom,
        "force_vector_RMSE_eV_A": force_vector_rmse,
        "predicted_energy_eV_cell": predicted_energy,
    })

unweighted_energy_rmse = float(np.sqrt(np.mean(np.square(energy_errors))))
atom_weighted_energy_rmse = float(np.sqrt(
    sum(natoms * error**2 for natoms, error in energy_errors_atom_weighted)
    / sum(natoms for natoms, _ in energy_errors_atom_weighted)
))
zero_filled_energy_rmse = float(np.sqrt(np.mean(np.square(zero_filled_errors))))
summary = {
    "model_sha256": hashlib.sha256(MODEL.read_bytes()).hexdigest(),
    "train_file_sha256": hashlib.sha256(TRAIN.read_bytes()).hexdigest(),
    "train_frames": len(frames),
    "energy_labelled_frames": len(energy_errors),
    "force_labelled_frames": sum(row["has_reference_forces"] for row in rows),
    "force_labelled_atoms": force_atom_count,
    "labelled_energy_RMSE_meV_atom_unweighted_by_atoms": unweighted_energy_rmse,
    "labelled_energy_RMSE_meV_atom_atom_weighted": atom_weighted_energy_rmse,
    "force_vector_RMSE_eV_A_atom_weighted": float(np.sqrt(force_squared_atom_sum / force_atom_count)),
    "zero_filled_all_frame_energy_RMSE_meV_atom": zero_filled_energy_rmse,
    "interpretation": "Only frames with REF_energy are valid energy targets. The zero-filled all-frame metric is calculated solely to diagnose the unmasked MACE summary metric.",
    "per_frame_csv": "v13_training_frame_audit.csv",
}
out = ROOT / "results"
out.mkdir(exist_ok=True)
with (out / "v13_training_frame_audit.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
(out / "v13_training_frame_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
