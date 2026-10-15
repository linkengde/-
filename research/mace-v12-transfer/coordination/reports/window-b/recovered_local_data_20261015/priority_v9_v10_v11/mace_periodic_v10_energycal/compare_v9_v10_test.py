#!/usr/bin/env python3
"""Compare v9 and v10 on the same frozen v10 test structures and DFT labels."""
from pathlib import Path
import csv, hashlib, json
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
TEST = ROOT / "data/test.extxyz"
MODEL_PATHS = {
    "v9": PROJECT / "mace_periodic_v9_forceonly/models/Ag_Ti_Si_C_MACE_periodic_v9_forceonly.model",
    "v10": ROOT / "checkpoints/MACE_periodic_v10_energycal_run-43.model",
    "v11": PROJECT / "mace_periodic_v11_energyfocus/checkpoints/MACE_periodic_v11_energyfocus_run-44.model",
}
if any(not p.exists() for p in MODEL_PATHS.values()):
    raise SystemExit("Need both v9 and final v10 model artifacts before scoring")
torch.set_default_dtype(torch.float64)
calcs = {}
for name, path in MODEL_PATHS.items():
    model = torch.load(path, map_location="cpu", weights_only=False).double().eval()
    calcs[name] = MACECalculator(models=[model], device="cpu", default_dtype="float64")
frames = read(TEST, index=":")
rows = []
for source in frames:
    atoms = source.copy()
    eref = float(atoms.info["REF_energy"])
    fref = np.asarray(atoms.arrays["REF_forces"], dtype=float)
    row = {"config_type": atoms.info.get("config_type"), "formula": atoms.get_chemical_formula(), "atoms": len(atoms), "DFT_energy_eV_cell": eref}
    for name, calc in calcs.items():
        atoms.calc = calc
        epred = float(atoms.get_potential_energy())
        fpred = np.asarray(atoms.get_forces(), dtype=float)
        diff = fpred - fref
        row.update({
            f"{name}_energy_error_DFT_minus_MACE_meV_atom": (eref - epred) * 1000.0 / len(atoms),
            f"{name}_force_RMSE_eV_A": float(np.sqrt(np.mean(np.sum(diff**2, axis=1)))),
            f"{name}_force_component_MAE_eV_A": float(np.mean(np.abs(diff))),
            f"{name}_force_vector_max_error_eV_A": float(np.linalg.norm(diff, axis=1).max()),
            f"{name}_Fmax_eV_A": float(np.linalg.norm(fpred, axis=1).max()),
        })
        central = atoms.arrays.get("central_pair")
        if central is not None and np.asarray(central, dtype=bool).sum() == 2:
            syms = atoms.get_chemical_symbols()
            idx = np.flatnonzero(np.asarray(central, dtype=bool))
            pair = {syms[i]: int(i) for i in idx}
            if "Ag" in pair and "Ti" in pair:
                i, j = pair["Ag"], pair["Ti"]
                vec = atoms.positions[j] - atoms.positions[i]
                dist = float(np.linalg.norm(vec)); unit = vec / dist
                row[f"AgTi_distance_A"] = dist
                row[f"{name}_AgTi_separating_force_eV_A"] = float(np.dot(fpred[j] - fpred[i], unit))
                row["DFT_AgTi_separating_force_eV_A"] = float(np.dot(fref[j] - fref[i], unit))
    rows.append(row)

summary = {"test_file": str(TEST), "test_frames": len(rows), "test_atoms": sum(r["atoms"] for r in rows), "models": {}, "per_structure": rows,
           "scope_note": "Same frozen DFT test labels for all versions. Includes a same-motif AgTi distance holdout; this is not a liquid, long-MD, or TTM validation."}
for name in MODEL_PATHS:
    energy_errors = [abs(r[f"{name}_energy_error_DFT_minus_MACE_meV_atom"]) for r in rows]
    force_sq = []
    force_comp_abs = []
    for atoms, row in zip(frames, rows):
        atoms.calc = calcs[name]
        d = np.asarray(atoms.get_forces()) - np.asarray(atoms.arrays["REF_forces"])
        force_sq.extend(np.sum(d**2, axis=1).tolist())
        force_comp_abs.extend(np.abs(d).ravel().tolist())
    summary["models"][name] = {
        "model": str(MODEL_PATHS[name]),
        "model_sha256": hashlib.sha256(MODEL_PATHS[name].read_bytes()).hexdigest(),
        "energy_MAE_meV_atom": float(np.mean(energy_errors)),
        "energy_max_abs_error_meV_atom": float(max(energy_errors)),
        "force_RMSE_eV_A": float(np.sqrt(np.mean(force_sq))),
        "force_component_MAE_eV_A": float(np.mean(force_comp_abs)),
    }
(ROOT / "v9_v10_v11_frozen_test_comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
with (ROOT / "v9_v10_v11_frozen_test_comparison.csv").open("w", newline="") as f:
    fields = sorted(set().union(*(row.keys() for row in rows)))
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
print(json.dumps(summary, indent=2))
