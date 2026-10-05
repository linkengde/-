#!/usr/bin/env python3
"""Localize per-atom force residuals on the three scored v14 blind holdouts."""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

HERE = Path(__file__).resolve().parent
TRANSFER = HERE.parents[2]
PROJECT = TRANSFER / "periodic_interface_v4"
EVAL = PROJECT / "mace_periodic_v14_interface_energy/results/v14_independent_holdout_comparison.json"
comparison = json.loads(EVAL.read_text())
MODEL = PROJECT / comparison["models"]["v14"]["path"]
model_hash = sha256(MODEL.read_bytes()).hexdigest()
if comparison["models"]["v14"]["sha256"] != model_hash:
    raise SystemExit("v14 checkpoint differs from the checkpoint used for the published evaluation.")

SOURCES = [
    (PROJECT / "pbe_interface_v14_main_holdouts", "AgC_registry_holdout_v14_01"),
    (PROJECT / "pbe_interface_v14_main_holdouts", "AgSi_registry_holdout_v14_01"),
    (PROJECT / "pbe_interface_v14_parallel_acquisition", "AgTi_registry_holdout_v14_01"),
]
torch.set_default_dtype(torch.float64)
model = torch.load(MODEL, map_location="cpu", weights_only=False).double().eval()
calc = MACECalculator(models=[model], device="cpu", default_dtype="float64")
records = []

for folder, label in SOURCES:
    expected = comparison["independent_geometry_holdouts"][label]
    result = folder / "calculations" / label
    verification = json.loads((result / "verification.json").read_text())
    if verification.get("status") != "PASS" or verification.get("role") != "v14_holdout":
        raise SystemExit("Unverified or incorrectly assigned reference: " + label)
    output = result / (label + "_PW_PBE.extxyz")
    output_hash = sha256(output.read_bytes()).hexdigest()
    if output_hash != verification["sha256"][output.name] or output_hash != expected["sha256"]:
        raise SystemExit("Reference hash mismatch: " + label)
    atoms = read(output)
    ref = np.asarray(atoms.arrays["PW_PBE_forces"], dtype=float)
    atoms.calc = calc
    pred = np.asarray(atoms.get_forces(), dtype=float)
    delta = pred - ref
    error = np.linalg.norm(delta, axis=1)
    dist = atoms.get_all_distances(mic=True)
    central = np.flatnonzero(np.asarray(atoms.arrays["central_pair"], dtype=bool))
    if len(central) != 2:
        raise SystemExit("Expected two marked central-pair atoms: " + label)
    local_distance = np.min(dist[:, central], axis=1)
    local = local_distance <= 4.5
    symbols = atoms.get_chemical_symbols()
    ids = np.asarray(atoms.arrays["lammps_id"], dtype=int)
    i, j = map(int, central)
    vec = atoms.positions[j] - atoms.positions[i]
    separation = vec / np.linalg.norm(vec)
    ref_sep = float(np.dot(ref[j] - ref[i], separation))
    pred_sep = float(np.dot(pred[j] - pred[i], separation))
    top = []
    for atom_i in np.argsort(error)[::-1][:10]:
        neighbors = [int(k) for k in np.argsort(dist[atom_i]) if k != atom_i][:4]
        denom = np.linalg.norm(ref[atom_i]) * np.linalg.norm(pred[atom_i])
        cosine = None if denom == 0 else float(np.dot(ref[atom_i], pred[atom_i]) / denom)
        top.append({
            "index": int(atom_i), "id": int(ids[atom_i]), "element": symbols[atom_i],
            "position_A": atoms.positions[atom_i].tolist(),
            "force_error_vector_eV_A": delta[atom_i].tolist(),
            "force_error_norm_eV_A": float(error[atom_i]),
            "dft_force_norm_eV_A": float(np.linalg.norm(ref[atom_i])),
            "mace_force_norm_eV_A": float(np.linalg.norm(pred[atom_i])),
            "force_direction_cosine": cosine,
            "distance_to_marked_contact_A": float(local_distance[atom_i]),
            "nearest_neighbors": [
                {"id": int(ids[k]), "element": symbols[k], "distance_A": float(dist[atom_i, k])}
                for k in neighbors
            ],
        })
    records.append({
        "label": label, "role": verification["role"], "atoms": len(atoms),
        "formula": atoms.get_chemical_formula(), "reference_sha256": output_hash,
        "force_vector_RMSE_all_eV_A": float(np.sqrt(np.mean(error ** 2))),
        "force_vector_RMSE_within_4p5A_of_marked_pair_eV_A": float(np.sqrt(np.mean(error[local] ** 2))),
        "force_vector_RMSE_outside_4p5A_eV_A": None if not np.any(~local) else float(np.sqrt(np.mean(error[~local] ** 2))),
        "max_force_error_eV_A": float(error.max()),
        "marked_pair": {"ids": [int(ids[i]), int(ids[j])], "elements": [symbols[i], symbols[j]],
                        "DFT_separating_force_eV_A": ref_sep, "MACE_separating_force_eV_A": pred_sep,
                        "absolute_error_eV_A": abs(pred_sep - ref_sep)},
        "top_10_atom_errors": top,
    })

report = {
    "purpose": "Post-hoc diagnosis of the published v14 holdout force failures; identify atoms and neighbors for future DFT sampling.",
    "v14_model_sha256": model_hash,
    "evaluation_sha256": sha256(EVAL.read_bytes()).hexdigest(),
    "local_region_definition": "Atoms within 4.5 Angstrom of either marked central-pair atom, using minimum-image distances.",
    "sources": records,
    "limitations": "Only one related small-cell motif per interface is analyzed. This diagnostic is not another validation pass and does not authorize production MD/TTM.",
}
out = HERE / "v14_force_localization.json"
out.write_text(json.dumps(report, indent=2) + "\n")
lines = [
    "# v14 force-error localization", "", f"Model SHA256: `{model_hash}`", "",
    "Post-hoc diagnosis of the published blind-holdout screen; local means within 4.5 Å of either marked contact atom.", "",
    "| Label | All-atom vector RMSE (eV/Å) | Local RMSE (eV/Å) | Outside RMSE (eV/Å) | Max atom error (eV/Å) | Contact separation error (eV/Å) |",
    "|---|---:|---:|---:|---:|---:|",
]
for row in records:
    outside = "n/a" if row["force_vector_RMSE_outside_4p5A_eV_A"] is None else f"{row['force_vector_RMSE_outside_4p5A_eV_A']:.4f}"
    lines.append(f"| {row['label']} | {row['force_vector_RMSE_all_eV_A']:.4f} | {row['force_vector_RMSE_within_4p5A_of_marked_pair_eV_A']:.4f} | {outside} | {row['max_force_error_eV_A']:.4f} | {row['marked_pair']['absolute_error_eV_A']:.4f} |")
lines += ["", "Top-error atom IDs, force directions and neighbors are in the JSON. Use this only to guide follow-up DFT; the v14 screen remains failed.", ""]
(HERE / "v14_force_localization.md").write_text("\n".join(lines))
print(json.dumps({"report": str(out), "records": len(records), "model_sha256": model_hash}, indent=2))
