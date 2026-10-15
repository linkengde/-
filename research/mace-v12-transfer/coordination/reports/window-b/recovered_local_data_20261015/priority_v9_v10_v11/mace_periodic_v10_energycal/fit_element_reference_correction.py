#!/usr/bin/env python3
"""Fit an energy-only per-element correction on v10 train labels, then score holdouts."""
from pathlib import Path
import csv
from copy import deepcopy
import hashlib
import json
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
MODEL_IN = ROOT / "checkpoints/MACE_periodic_v10_energycal_run-43.model"
OUTDIR = ROOT / "energy_reference_correction"
MODEL_OUT = OUTDIR / "MACE_periodic_v10_refcorr.model"
if MODEL_OUT.exists():
    raise SystemExit(f"Refusing to overwrite {MODEL_OUT}")
OUTDIR.mkdir(parents=True, exist_ok=True)

Z_TO_SYMBOL = {6: "C", 14: "Si", 22: "Ti", 47: "Ag"}
torch.set_default_dtype(torch.float64)
raw_model = torch.load(MODEL_IN, map_location="cpu", weights_only=False).double().eval()
z_order = [int(z) for z in raw_model.atomic_numbers.tolist()]
if set(z_order) != set(Z_TO_SYMBOL):
    raise SystemExit(f"Unexpected model element set: {z_order}")
calc_raw = MACECalculator(models=[raw_model], device="cpu", default_dtype="float64")

def energy_label(atoms):
    if "REF_energy" in atoms.info:
        return float(atoms.info["REF_energy"])
    return float(atoms.get_potential_energy())

def counts(atoms):
    numbers = np.asarray(atoms.numbers)
    return np.array([np.count_nonzero(numbers == z) for z in z_order], dtype=float)

train_path = ROOT / "data/train.extxyz"
train = [a for a in read(train_path, index=":") if "REF_energy" in a.info]
if len(train) != 15:
    raise SystemExit(f"Expected 15 labeled train frames; found {len(train)}")
A = np.stack([counts(a) for a in train])
y = []
for atoms in train:
    atoms.calc = calc_raw
    y.append(energy_label(atoms) - float(atoms.get_potential_energy()))
y = np.asarray(y, dtype=float)
delta, residuals, rank, singular_values = np.linalg.lstsq(A, y, rcond=None)
if rank != 4:
    raise SystemExit(f"Element-reference fit is rank deficient: {rank}/4")
train_residual = y - A @ delta

corrected_model = deepcopy(raw_model)
z_index = {z: i for i, z in enumerate(z_order)}
with torch.no_grad():
    buf = corrected_model.atomic_energies_fn.atomic_energies
    for z, correction in zip(z_order, delta):
        buf[0, z_index[z]] += float(correction)
calc_corr = MACECalculator(models=[corrected_model], device="cpu", default_dtype="float64")
torch.save(corrected_model, MODEL_OUT)

sources = {
    "train": [train_path],
    "valid": [ROOT / "data/valid.extxyz"],
    "test": [ROOT / "data/test.extxyz"],
    "AgSi_external": [PROJECT / "pbe_external_v8_motif_holdouts/AgSi/retry1/AgSi_cluster_PW_PBE.extxyz"],
    "AgC_method_transfer": [PROJECT / "pbe_cluster_transfer_pw_20261002_run3_mixer_retry2/AgC_cluster_PW_PBE.extxyz"],
}
rows = []
for split, paths in sources.items():
    frames = [a for path in paths for a in read(path, index=":")]
    for atoms in frames:
        if split == "train" and "REF_energy" not in atoms.info:
            continue
        eref = energy_label(atoms)
        n = counts(atoms)
        atoms.calc = calc_raw
        eraw = float(atoms.get_potential_energy())
        atoms.calc = calc_corr
        ecorr = float(atoms.get_potential_energy())
        rows.append({
            "split": split,
            "config_type": atoms.info.get("config_type", split),
            "formula": atoms.get_chemical_formula(),
            "atoms": len(atoms),
            "DFT_energy_eV_cell": eref,
            "raw_model_minus_DFT_meV_atom": (eraw - eref) * 1000.0 / len(atoms),
            "refcorr_model_minus_DFT_meV_atom": (ecorr - eref) * 1000.0 / len(atoms),
            "element_counts": json.dumps({Z_TO_SYMBOL[z]: int(c) for z, c in zip(z_order, n)}),
        })

metrics = {}
for split in sources:
    rr = [r for r in rows if r["split"] == split]
    if not rr:
        continue
    for key, suffix in [("raw_model_minus_DFT_meV_atom", "raw"), ("refcorr_model_minus_DFT_meV_atom", "refcorr")]:
        err = np.asarray([r[key] for r in rr], dtype=float)
        metrics.setdefault(split, {})[f"{suffix}_energy_MAE_meV_atom"] = float(np.mean(np.abs(err)))
        metrics[split][f"{suffix}_energy_RMSE_meV_atom"] = float(np.sqrt(np.mean(err**2)))

summary = {
    "fit_split": "train only; validation, test, and external frames are not used to fit corrections",
    "fit_energy_frames": len(train),
    "composition_matrix_rank": int(rank),
    "element_order": [Z_TO_SYMBOL[z] for z in z_order],
    "correction_eV_per_atom": {Z_TO_SYMBOL[z]: float(v) for z, v in zip(z_order, delta)},
    "train_residual_RMSE_eV_cell": float(np.sqrt(np.mean(train_residual**2))),
    "train_residual_RMSE_meV_atom": float(np.sqrt(np.mean((train_residual / np.array([len(a) for a in train]))**2)) * 1000.0),
    "raw_model_sha256": hashlib.sha256(MODEL_IN.read_bytes()).hexdigest(),
    "corrected_model": str(MODEL_OUT),
    "corrected_model_sha256": hashlib.sha256(MODEL_OUT.read_bytes()).hexdigest(),
    "metrics": metrics,
    "per_structure": rows,
    "interpretation": "This is a four-element, energy-only calibration fit on labeled training compositions. It leaves forces unchanged and must be judged by untouched validation, frozen test, and external structures. It is not a physical validation or a production potential endorsement.",
}
(OUTDIR / "energy_reference_correction_metrics.json").write_text(json.dumps(summary, indent=2) + "\n")
with (OUTDIR / "energy_reference_correction_metrics.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(summary, indent=2))
