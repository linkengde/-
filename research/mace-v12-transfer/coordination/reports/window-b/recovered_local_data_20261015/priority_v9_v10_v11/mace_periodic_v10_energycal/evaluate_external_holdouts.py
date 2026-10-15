#!/usr/bin/env python3
"""Compare frozen v9 and v10 models on independent Ag-Si and Ag-C PBE frames."""
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
HOLDOUTS = {
    "AgSi": PROJECT / "pbe_external_v8_motif_holdouts/AgSi/retry1/AgSi_cluster_PW_PBE.extxyz",
    "AgC": PROJECT / "pbe_cluster_transfer_pw_20261002_run3_mixer_retry2/AgC_cluster_PW_PBE.extxyz",
}
MODEL_PATHS = {
    "v9": PROJECT / "mace_periodic_v9_forceonly/models/Ag_Ti_Si_C_MACE_periodic_v9_forceonly.model",
    "v10": ROOT / "checkpoints/MACE_periodic_v10_energycal_run-43.model",
    "v11": PROJECT / "mace_periodic_v11_energyfocus/checkpoints/MACE_periodic_v11_energyfocus_run-44.model",
}
if any(not p.exists() for p in MODEL_PATHS.values()):
    raise SystemExit("Need both v9 and final v10 model artifacts before scoring")
if any(not p.exists() for p in HOLDOUTS.values()):
    raise SystemExit("One or more independent holdout structures are missing")

torch.set_default_dtype(torch.float64)
calcs = {}
for name, path in MODEL_PATHS.items():
    model = torch.load(path, map_location="cpu", weights_only=False).double().eval()
    calcs[name] = MACECalculator(models=[model], device="cpu", default_dtype="float64")

rows = []
for contact, path in HOLDOUTS.items():
    atoms = read(path)
    # ASE maps extxyz energy/forces to atoms.info['energy'] and atoms.arrays['forces'].
    eref = float(atoms.info["REF_energy"]) if "REF_energy" in atoms.info else float(atoms.get_potential_energy())
    fref = np.asarray(atoms.arrays["REF_forces"], dtype=float) if "REF_forces" in atoms.arrays else np.asarray(atoms.get_forces(), dtype=float)
    base = {
        "contact": contact,
        "source_file": str(path),
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "formula": atoms.get_chemical_formula(),
        "atoms": len(atoms),
        "DFT_energy_eV_cell": eref,
    }
    central = np.asarray(atoms.arrays.get("central_pair", np.zeros(len(atoms))), dtype=bool)
    pair_indices = np.flatnonzero(central)
    if len(pair_indices) == 2:
        syms = atoms.get_chemical_symbols()
        pair = {syms[i]: int(i) for i in pair_indices}
    else:
        pair = {}
    for name, calc in calcs.items():
        atoms.calc = calc
        epred = float(atoms.get_potential_energy())
        fpred = np.asarray(atoms.get_forces(), dtype=float)
        delta = fpred - fref
        base[f"{name}_energy_error_meV_atom"] = (epred - eref) * 1000.0 / len(atoms)
        base[f"{name}_force_vector_RMSE_eV_A"] = float(np.sqrt(np.mean(np.sum(delta**2, axis=1))))
        base[f"{name}_force_component_RMSE_eV_A"] = float(np.sqrt(np.mean(delta**2)))
        base[f"{name}_force_vector_max_error_eV_A"] = float(np.linalg.norm(delta, axis=1).max())
        if contact == "AgSi" and "Ag" in pair and "Si" in pair:
            i, j = pair["Ag"], pair["Si"]
        elif contact == "AgC" and "Ag" in pair and "C" in pair:
            i, j = pair["Ag"], pair["C"]
        else:
            i = j = None
        if i is not None:
            vec = atoms.positions[j] - atoms.positions[i]
            dist = float(np.linalg.norm(vec))
            unit = vec / dist
            base["central_pair_distance_A"] = dist
            base["DFT_pair_separating_force_eV_A"] = float(np.dot(fref[j] - fref[i], unit))
            base[f"{name}_pair_separating_force_eV_A"] = float(np.dot(fpred[j] - fpred[i], unit))
    rows.append(base)

summary = {
    "models": {
        name: {
            "path": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        for name, path in MODEL_PATHS.items()
    },
    "holdout_count": len(rows),
    "per_structure": rows,
    "scope_note": (
        "These independent structures were excluded from v10 training. Ag-Si is an external periodic PBE motif; "
        "Ag-C is a method-transfer check on a previously studied geometry, not an unseen geometry. "
        "Passing these small frozen tests does not validate liquid, TTM, or production-scale behavior."
    ),
}
(ROOT / "v9_v10_v11_external_holdouts.json").write_text(json.dumps(summary, indent=2) + "\n")
with (ROOT / "v9_v10_v11_external_holdouts.csv").open("w", newline="") as f:
    fields = sorted(set().union(*(row.keys() for row in rows)))
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(summary, indent=2))
