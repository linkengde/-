#!/usr/bin/env python3
"""Create auditable not-run status records; never launches GPAW."""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[5]
INPUTS = REPO / "research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_ag_interface_fit_dft_inputs_20261014"
INPUT_COMMIT = "1fec365ecedb8555ae2c857509d1018953639b86"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, fields, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


manifest = json.loads((INPUTS / "INPUT_MANIFEST.json").read_text())
geometry = json.loads((INPUTS / "GEOMETRY_VALIDATION.json").read_text())
geometry_by_id = {x["task_id"]: x for x in geometry["items"]}
tasks = manifest["backbone_tasks"] + manifest["interface_scan_tasks"]
structure_rows, label_rows = [], []
for task in tasks:
    p = INPUTS / task["structure_path"]
    atoms = read(p)
    symbols = atoms.get_chemical_symbols()
    labels = {k: {"count": len(atoms.arrays[k]), "unique": len(set(atoms.arrays[k].tolist())) == len(atoms)} for k in atoms.arrays if k not in ("numbers", "positions")}
    atom_label_status = all(v["unique"] for k, v in labels.items() if k in ("atom_id", "local_new_id"))
    family = task["parent_family"]
    interface = family != "TSC_KISI_COD9009647_BASELINE_DERIVED_FAMILY"
    method = task["method"]
    final_status = "HOLD_RESOURCE_GATE"
    if interface:
        final_status += ";KPOINT_CONVERGENCE_PENDING"
    structure_rows.append({
        "task_id": task["task_id"],
        "structure_id": task["structure_id"],
        "group": "B_INTERFACE" if interface else "A_TSC_BACKBONE",
        "parent_family": family,
        "input_package_commit": INPUT_COMMIT,
        "structure_path": str((INPUTS / task["structure_path"]).relative_to(REPO)),
        "geometry_sha256": task["structure_sha256"],
        "config_path": str((INPUTS / task["config_path"]).relative_to(REPO)),
        "config_sha256": task["config_sha256"],
        "natoms": len(atoms),
        "composition": atoms.get_chemical_formula(),
        "ordered_symbols_sha256": hashlib.sha256(" ".join(symbols).encode()).hexdigest(),
        "atom_labels": json.dumps(labels, sort_keys=True),
        "atom_label_unique": atom_label_status,
        "cell_A": json.dumps(atoms.cell.array.tolist(), separators=(",", ":")),
        "pbc": "".join("T" if bool(x) else "F" for x in atoms.pbc),
        "minimum_periodic_distance_A": geometry_by_id[task["task_id"]].get("minimum_periodic_distance_A", ""),
        "kpoints": "x".join(map(str, method["kpoints"])),
        "kpoint_convention": method["kpoint_convention"],
        "cutoff_eV": method["cutoff_eV"],
        "smearing_eV": method["occupations"]["width_eV"],
        "scf_thresholds": json.dumps(method["convergence"], sort_keys=True),
        "PAW_setup_sha256": json.dumps(method["paw_setup_sha256"], sort_keys=True),
        "input_audit": "PASS" if not interface else "PASS_WITH_KPOINT_CONVERGENCE_HOLD",
        "SCF_status": "NOT_STARTED",
        "result_status": final_status,
    })
    label_rows.append({
        "task_id": task["task_id"],
        "structure_sha256": task["structure_sha256"],
        "native_total_energy_eV": "",
        "force_consistent_free_energy_eV": "",
        "forces_artifact": "",
        "stress_eV_per_A3_artifact": "",
        "SCF_converged": "",
        "SCF_iterations": "",
        "status": "NOT_RUN_RESOURCE_GATE",
        "note": "No DFT calculation was started; blank labels are missing, not zero.",
    })

write_csv(ROOT / "STRUCTURE_RESULTS.csv", list(structure_rows[0]), structure_rows)
write_csv(ROOT / "ENERGY_FORCE_STRESS_SUMMARY.csv", list(label_rows[0]), label_rows)
print(f"Prepared pending status rows for {len(tasks)} inputs; no DFT started.")
