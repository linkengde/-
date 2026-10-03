#!/usr/bin/env python3
"""Score v9-v12 on the frozen test split and independent interface holdouts."""
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
    "v12": ROOT / "checkpoints/MACE_periodic_v12_interface_energy_run-45.model",
}
TEST = ROOT / "data/test.extxyz"
HOLDOUTS = {
    "AgSi_midpoint": PROJECT / "pbe_external_v8_motif_holdouts/AgSi/retry1/AgSi_cluster_PW_PBE.extxyz",
    "AgC_midpoint": PROJECT / "pbe_cluster_transfer_pw_20261002_run3_mixer_retry2/AgC_cluster_PW_PBE.extxyz",
}
if any(not p.is_file() for p in [*MODELS.values(), TEST, *HOLDOUTS.values()]):
    missing = [str(p) for p in [*MODELS.values(), TEST, *HOLDOUTS.values()] if not p.is_file()]
    raise SystemExit("Required model/data file missing:\n" + "\n".join(missing))

torch.set_default_dtype(torch.float64)
calculators = {}
for name, path in MODELS.items():
    model = torch.load(path, map_location="cpu", weights_only=False).double().eval()
    calculators[name] = MACECalculator(models=[model], device="cpu", default_dtype="float64")


def reference(atoms):
    energy = atoms.info.get("REF_energy", atoms.info.get("PW_PBE_energy_eV"))
    forces = atoms.arrays.get("REF_forces", atoms.arrays.get("PW_PBE_forces"))
    if energy is None or forces is None:
        raise SystemExit(f"Missing reference energy/forces for {atoms.info.get('config_type')}")
    return float(energy), np.asarray(forces, dtype=float)


def pair_force(atoms, forces):
    central = np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool)
    idx = np.flatnonzero(central)
    if len(idx) != 2:
        return {}
    symbols = atoms.get_chemical_symbols()
    ia, ib = int(idx[0]), int(idx[1])
    if symbols[ia] == symbols[ib]:
        return {}
    vec = atoms.positions[ib] - atoms.positions[ia]
    distance = float(np.linalg.norm(vec))
    unit = vec / distance
    return {
        "central_pair_elements": f"{symbols[ia]}-{symbols[ib]}",
        "central_pair_distance_A": distance,
        "separating_force_eV_A": float(np.dot(forces[ib] - forces[ia], unit)),
    }


def evaluate(source, label, split):
    atoms = source.copy()
    eref, fref = reference(atoms)
    row = {
        "split": split,
        "label": label,
        "config_type": atoms.info.get("config_type", ""),
        "formula": atoms.get_chemical_formula(),
        "atoms": len(atoms),
        "DFT_energy_eV_cell": eref,
    }
    row.update({f"DFT_{k}": v for k, v in pair_force(atoms, fref).items()})
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
        row.update({f"{name}_{k}": v for k, v in pair_force(atoms, fpred).items() if k != "central_pair_elements" and k != "central_pair_distance_A"})
    return row


rows = []
for atoms in read(TEST, index=":"):
    rows.append(evaluate(atoms, atoms.info.get("config_type", ""), "frozen_test"))
for label, path in HOLDOUTS.items():
    rows.append(evaluate(read(path), label, "external_holdout"))

summary = {
    "models": {
        name: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for name, path in MODELS.items()
    },
    "test_file_sha256": hashlib.sha256(TEST.read_bytes()).hexdigest(),
    "holdout_file_sha256": {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in HOLDOUTS.items()},
    "per_structure": rows,
    "scope_note": "Small fixed-geometry checks only. AgC is a method-transfer check on a previously studied geometry; these metrics do not establish liquid, TTM, long-MD, or production validation.",
}
for split in ["frozen_test", "external_holdout"]:
    split_rows = [r for r in rows if r["split"] == split]
    summary[split] = {"structures": len(split_rows), "atoms": sum(r["atoms"] for r in split_rows), "models": {}}
    for name in MODELS:
        e = [abs(r[f"{name}_energy_error_DFT_minus_MACE_meV_atom"]) for r in split_rows]
        # Recompute force aggregation from per-structure RMSE, weighted by atom count.
        f2 = sum(r["atoms"] * r[f"{name}_force_vector_RMSE_eV_A"] ** 2 for r in split_rows)
        summary[split]["models"][name] = {
            "energy_MAE_meV_atom": float(np.mean(e)),
            "energy_max_abs_error_meV_atom": float(max(e)),
            "force_vector_RMSE_eV_A": float(np.sqrt(f2 / sum(r["atoms"] for r in split_rows))),
        }

out = ROOT / "results"
out.mkdir(exist_ok=True)
(out / "v9_v10_v11_v12_candidate_comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
fields = sorted(set().union(*(r.keys() for r in rows)))
with (out / "v9_v10_v11_v12_candidate_comparison.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(summary, indent=2))
