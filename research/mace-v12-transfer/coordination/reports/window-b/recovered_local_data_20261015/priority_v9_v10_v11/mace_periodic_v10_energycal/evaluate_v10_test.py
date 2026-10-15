#!/usr/bin/env python3
"""Evaluate an already trained v10 model on the frozen test split."""
from pathlib import Path
import csv, json
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parent
MODEL = ROOT / "models/MACE_periodic_v10_energycal.model"
TEST = ROOT / "data/test.extxyz"
if not MODEL.exists():
    raise SystemExit(f"Trained model not found: {MODEL}")
torch.set_default_dtype(torch.float64)
model = torch.load(MODEL, map_location="cpu", weights_only=False).double().eval()
calc = MACECalculator(models=[model], device="cpu", default_dtype="float64")
frames = read(TEST, index=":")
rows = []
for atoms in frames:
    atoms.calc = calc
    eref = float(atoms.info["REF_energy"])
    fref = np.asarray(atoms.arrays["REF_forces"], dtype=float)
    epred = float(atoms.get_potential_energy())
    fpred = np.asarray(atoms.get_forces(), dtype=float)
    diff = fpred - fref
    row = {
        "config_type": atoms.info.get("config_type"),
        "formula": atoms.get_chemical_formula(),
        "atoms": len(atoms),
        "energy_DFT_eV_cell": eref,
        "energy_MACE_eV_cell": epred,
        "energy_error_DFT_minus_MACE_meV_atom": (eref - epred) * 1000.0 / len(atoms),
        "force_RMSE_eV_A": float(np.sqrt(np.mean(np.sum(diff**2, axis=1)))),
        "force_component_MAE_eV_A": float(np.mean(np.abs(diff))),
        "force_vector_max_error_eV_A": float(np.linalg.norm(diff, axis=1).max()),
        "DFT_Fmax_eV_A": float(np.linalg.norm(fref, axis=1).max()),
        "MACE_Fmax_eV_A": float(np.linalg.norm(fpred, axis=1).max()),
    }
    central = atoms.arrays.get("central_pair")
    if central is not None and np.asarray(central, dtype=bool).sum() == 2:
        syms = atoms.get_chemical_symbols()
        idx = np.flatnonzero(np.asarray(central, dtype=bool))
        pair = {syms[i]: int(i) for i in idx}
        if "Ag" in pair and "Ti" in pair:
            i, j = pair["Ag"], pair["Ti"]
            vec = atoms.positions[j] - atoms.positions[i]
            dist = float(np.linalg.norm(vec))
            unit = vec / dist
            row["central_AgTi_distance_A"] = dist
            row["central_pair_separating_force_DFT_eV_A"] = float(np.dot(fref[j] - fref[i], unit))
            row["central_pair_separating_force_MACE_eV_A"] = float(np.dot(fpred[j] - fpred[i], unit))
    rows.append(row)

all_diff = []
all_force_sq = []
for atoms in read(TEST, index=":"):
    atoms.calc = calc
    all_diff.extend([float(atoms.info["REF_energy"]) - float(atoms.get_potential_energy())])
    delta = np.asarray(atoms.get_forces()) - np.asarray(atoms.arrays["REF_forces"])
    all_force_sq.extend(np.sum(delta**2, axis=1).tolist())
summary = {
    "model": str(MODEL),
    "model_sha256": __import__("hashlib").sha256(MODEL.read_bytes()).hexdigest(),
    "test_file": str(TEST),
    "test_frames": len(rows),
    "test_atoms": sum(row["atoms"] for row in rows),
    "energy_MAE_meV_atom": float(np.mean([abs(row["energy_error_DFT_minus_MACE_meV_atom"]) for row in rows])),
    "energy_max_abs_error_meV_atom": float(max(abs(row["energy_error_DFT_minus_MACE_meV_atom"]) for row in rows)),
    "force_RMSE_eV_A": float(np.sqrt(np.mean(all_force_sq))),
    "per_structure": rows,
    "interpretation": "This is a frozen test split. AgTi 2.70 A is a new periodic distance in a known local motif; it does not validate liquid, TTM, or production-scale behavior.",
}
(ROOT / "v10_test_metrics.json").write_text(json.dumps(summary, indent=2) + "\n")
with (ROOT / "v10_test_metrics.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print(json.dumps(summary, indent=2))
