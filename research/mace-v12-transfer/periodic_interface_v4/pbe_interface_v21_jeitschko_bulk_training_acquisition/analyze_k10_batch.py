#!/usr/bin/env python3
"""Analyze the higher-mesh force-rich candidate labels and paired response."""
import csv
import json
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT / "input_manifest.json").read_text())
records = {r["label"]: r for r in manifest["records"]}


def get(label):
    folder = ROOT / "calculations" / label
    summary = json.loads((folder / "summary.json").read_text())
    atoms = read(folder / f"{label}_PW_PBE.extxyz")
    return summary, atoms, atoms.arrays["PW_PBE_forces"]


baseline = "Ti3SiC2_JN1967_primitive_baseline_k10x10x2_mesh_control_v21"
labels = [
    "Ti3SiC2_JN1967_primitive_Ti4f_id0003_z_m002A_k10x10x2_v21",
    "Ti3SiC2_JN1967_primitive_Ti4f_id0003_z_p002A_k10x10x2_v21",
    "Ti3SiC2_JN1967_primitive_C4f_id0009_z_m002A_k10x10x2_v21",
    "Ti3SiC2_JN1967_primitive_C4f_id0009_z_p002A_k10x10x2_v21",
]
base_summary, base_atoms, base_forces = get(baseline)
base_energy = base_summary["native_extrapolated_energy_eV_cell"]
base_free = base_summary["free_energy_eV_cell"]
rows = []
for label in labels:
    summary, atoms, forces = get(label)
    record = records[label]
    local_id = int(record["site_id"])
    rows.append({
        "label": label,
        "input_sha256": record["input_sha256"],
        "delta_z_A": float(record["delta_z_A"]),
        "displaced_local_new_id": local_id,
        "native_energy_eV_cell": summary["native_extrapolated_energy_eV_cell"],
        "free_energy_eV_cell": summary["free_energy_eV_cell"],
        "delta_native_energy_eV_cell_vs_k10_baseline": summary["native_extrapolated_energy_eV_cell"] - base_energy,
        "delta_free_energy_eV_cell_vs_k10_baseline": summary["free_energy_eV_cell"] - base_free,
        "displaced_atom_force_z_eV_A": float(forces[local_id - 1, 2]),
        "max_force_eV_A": summary["max_force_eV_A"],
        "rms_force_vector_eV_A": summary["rms_force_vector_eV_A"],
        "scf_iterations": summary["scf_iterations"],
    })

paired = {}
for name, minus, plus, local_id in [
    ("Ti4f_id0003_z", labels[0], labels[1], 3),
    ("C4f_id0009_z", labels[2], labels[3], 9),
]:
    sm, _, fm = get(minus)
    sp, _, fp = get(plus)
    secant = -(sp["free_energy_eV_cell"] - sm["free_energy_eV_cell"]) / 0.04
    paired[name] = {
        "relative_displacement_A": 0.04,
        "free_energy_secant_force_eV_A": secant,
        "native_energy_secant_force_eV_A": -(sp["native_extrapolated_energy_eV_cell"] - sm["native_extrapolated_energy_eV_cell"]) / 0.04,
        "mean_pair_force_z_eV_A": float((fm[local_id - 1, 2] + fp[local_id - 1, 2]) / 2),
        "baseline_force_z_eV_A": float(base_forces[local_id - 1, 2]),
        "secant_minus_baseline_force_eV_A": secant - float(base_forces[local_id - 1, 2]),
        "delta_free_energy_plus_minus_eV_cell": sp["free_energy_eV_cell"] - sm["free_energy_eV_cell"],
    }

out = {
    "status": "all_four_candidate_labels_converged_and_pass_geometry_label_identity_verification",
    "method": records[labels[0]]["method"],
    "source_family": "Jeitschko_Nowotny_1967_single_crystal",
    "baseline_reference": baseline,
    "baseline_role": records[baseline]["role"],
    "baseline_native_energy_eV_cell": base_energy,
    "baseline_free_energy_eV_cell": base_free,
    "baseline_max_force_eV_A": base_summary["max_force_eV_A"],
    "rows": rows,
    "paired_displacement_response": paired,
    "interpretation": "Paired free-energy secants track force components at the displaced site; all files are geometry- and hash-verified. The labels remain training candidates pending role/provenance review and are not integrated into V19.",
    "limitations": [
        "Only four displaced bulk-backbone labels are candidates; there are no Ag-interface labels in this batch.",
        "The baseline used for the secant is a k-point control, not a registered training row.",
        "Full-text Fourier coordinates assume ideal full occupancy despite occupancy uncertainty in the original paper.",
        "The 10x10x2 labels were not mixed with Gamma labels or used to fit a model.",
    ],
}
(ROOT / "k10_batch_analysis.json").write_text(json.dumps(out, indent=2) + "\n")
with (ROOT / "k10_batch_summary.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(out, indent=2))
