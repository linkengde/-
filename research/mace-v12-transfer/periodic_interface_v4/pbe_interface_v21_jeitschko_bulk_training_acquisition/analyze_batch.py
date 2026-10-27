#!/usr/bin/env python3
"""Summarize the completed small-cell labels without integrating them into training."""
import csv
import json
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT / "input_manifest.json").read_text())
records = {r["label"]: r for r in manifest["records"]}


def result(label):
    folder = ROOT / "calculations" / label
    summary = json.loads((folder / "summary.json").read_text())
    atoms = read(folder / f"{label}_PW_PBE.extxyz")
    return summary, atoms, atoms.arrays["PW_PBE_forces"]


labels = [
    "Ti3SiC2_JN1967_primitive_baseline_k6x6x2_v21",
    "Ti3SiC2_JN1967_primitive_Ti4f_id0003_z_m002A_k6x6x2_v21",
    "Ti3SiC2_JN1967_primitive_Ti4f_id0003_z_p002A_k6x6x2_v21",
    "Ti3SiC2_JN1967_primitive_C4f_id0009_z_m002A_k6x6x2_v21",
    "Ti3SiC2_JN1967_primitive_C4f_id0009_z_p002A_k6x6x2_v21",
]
data = {label: result(label) for label in labels}
base_label = labels[0]
base_summary, base_atoms, base_forces = data[base_label]
base_energy = base_summary["native_extrapolated_energy_eV_cell"]
base_free = base_summary["free_energy_eV_cell"]

rows = []
for label in labels:
    summary, atoms, forces = data[label]
    record = records[label]
    delta = float(record["delta_z_A"])
    displaced_raw = record.get("displaced_local_new_id", record.get("site_id", -1))
    displaced_id = -1 if displaced_raw is None else int(displaced_raw)
    target_force = None if displaced_id < 1 else float(forces[displaced_id - 1, 2])
    rows.append({
        "label": label,
        "input_sha256": record["input_sha256"],
        "native_energy_eV_cell": summary["native_extrapolated_energy_eV_cell"],
        "free_energy_eV_cell": summary["free_energy_eV_cell"],
        "delta_native_energy_eV_cell_vs_baseline": summary["native_extrapolated_energy_eV_cell"] - base_energy,
        "delta_free_energy_eV_cell_vs_baseline": summary["free_energy_eV_cell"] - base_free,
        "delta_z_A": delta,
        "displaced_local_new_id": displaced_id,
        "displaced_atom_force_z_eV_A": target_force,
        "max_force_eV_A": summary["max_force_eV_A"],
        "rms_force_vector_eV_A": summary["rms_force_vector_eV_A"],
        "scf_iterations": summary["scf_iterations"],
    })

paired = {}
for name, minus, plus, local_id in [
    ("Ti4f_id0003_z", labels[1], labels[2], 3),
    ("C4f_id0009_z", labels[3], labels[4], 9),
]:
    sm, am, fm = data[minus]
    sp, ap, fp = data[plus]
    h = 0.02
    secant_free_force = -(sp["free_energy_eV_cell"] - sm["free_energy_eV_cell"]) / (2 * h)
    secant_native_force = -(sp["native_extrapolated_energy_eV_cell"] - sm["native_extrapolated_energy_eV_cell"]) / (2 * h)
    paired[name] = {
        "minus_label": minus,
        "plus_label": plus,
        "relative_displacement_A": 2 * h,
        "free_energy_secant_force_eV_A": secant_free_force,
        "native_energy_secant_force_eV_A": secant_native_force,
        "mean_displaced_atom_force_z_eV_A": float((fm[local_id - 1, 2] + fp[local_id - 1, 2]) / 2),
        "baseline_displaced_atom_force_z_eV_A": float(base_forces[local_id - 1, 2]),
        "secant_minus_baseline_force_eV_A": secant_free_force - float(base_forces[local_id - 1, 2]),
        "delta_free_energy_plus_minus_eV_cell": sp["free_energy_eV_cell"] - sm["free_energy_eV_cell"],
    }

species_metrics = {}
for element in sorted(set(base_atoms.get_chemical_symbols())):
    idx = [i for i, symbol in enumerate(base_atoms.get_chemical_symbols()) if symbol == element]
    magnitudes = np.linalg.norm(base_forces, axis=1)
    species_metrics[element] = {
        "max_force_eV_A": float(magnitudes[idx].max()),
        "rms_force_vector_eV_A": float(np.sqrt(np.mean(magnitudes[idx] ** 2))),
    }

mesh6 = data[base_label]
mesh10 = result("Ti3SiC2_JN1967_primitive_baseline_k10x10x2_mesh_control_v21")
abstract10 = result("Ti3SiC2_JN1967_abstract_primitive_baseline_k10x10x2_mesh_control_v21")
mesh_compare = {
    "same_geometry_k6x6x2_vs_k10x10x2": {
        "native_delta_eV_cell": mesh10[0]["native_extrapolated_energy_eV_cell"] - mesh6[0]["native_extrapolated_energy_eV_cell"],
        "free_delta_eV_cell": mesh10[0]["free_energy_eV_cell"] - mesh6[0]["free_energy_eV_cell"],
        "max_force_k6_eV_A": mesh6[0]["max_force_eV_A"],
        "max_force_k10_eV_A": mesh10[0]["max_force_eV_A"],
        "force_vector_rms_difference_eV_A": float(np.sqrt(np.mean(np.sum((mesh10[2] - mesh6[2]) ** 2, axis=1)))),
    },
    "same_k10x10x2_abstract_vs_fulltext_fourier_coordinates": {
        "native_delta_eV_cell_abstract_minus_fourier": abstract10[0]["native_extrapolated_energy_eV_cell"] - mesh10[0]["native_extrapolated_energy_eV_cell"],
        "free_delta_eV_cell_abstract_minus_fourier": abstract10[0]["free_energy_eV_cell"] - mesh10[0]["free_energy_eV_cell"],
        "max_force_fourier_eV_A": mesh10[0]["max_force_eV_A"],
        "max_force_abstract_eV_A": abstract10[0]["max_force_eV_A"],
        "force_vector_rms_difference_eV_A": float(np.sqrt(np.mean(np.sum((abstract10[2] - mesh10[2]) ** 2, axis=1)))),
    },
}

out = {
    "status": "geometry_and_label_checks_pass; training_integration_not_performed",
    "source_doi": "10.1007/BF00899949",
    "source_family": "Jeitschko_Nowotny_1967_single_crystal",
    "underconverged_mesh_control_method": records[base_label]["method"],
    "formula_units_per_cell": 2,
    "baseline_species_force_metrics": species_metrics,
    "rows": rows,
    "paired_displacement_response": paired,
    "method_and_coordinate_controls": mesh_compare,
    "limitations": [
        "These are five fixed-cell DFT labels from one parent source family, not independent validation data.",
        "All five 6x6x2 results are underconverged-mesh controls and must not be integrated into training.",
        "Native extrapolated energies and Fermi-width free energies are kept separate; ASE forces at 0.10 eV width are free-energy derivatives.",
        "No Ag-Ti, Ag-Si, or Ag-C interface labels were calculated in this batch.",
        "No model was fit or evaluated and no official dataset was changed.",
    ],
}
(ROOT / "dft_batch_analysis.json").write_text(json.dumps(out, indent=2) + "\n")
with (ROOT / "dft_batch_summary.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(out, indent=2))
