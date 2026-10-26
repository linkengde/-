#!/usr/bin/env python3
"""Independent compact-archive checks for one B1 mesh endpoint."""
import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
TASK = json.loads((HERE / "task_manifest.json").read_text())
parser = argparse.ArgumentParser()
parser.add_argument("label")
args = parser.parse_args()
record = next((r for r in TASK["records"] + TASK["conditional_records"] if r["label"] == args.label), None)
assert record is not None, f"Unregistered label: {args.label}"
RUN = HERE / TASK["outputs_directory"] / args.label
source = ROOT / record["input"]
assert hashlib.sha256(source.read_bytes()).hexdigest() == record["input_sha256"]
input_atoms, output_atoms = read(source), read(RUN / "result.extxyz")
summary = json.loads((RUN / "summary.json").read_text())
progress = json.loads((RUN / "progress.json").read_text())
log = (RUN / "gpaw.log").read_text(errors="replace")
forces = np.asarray(output_atoms.arrays["PW_PBE_forces"], dtype=float)
checks = {
    "label_and_source": summary["label"] == progress["label"] == args.label and summary["source_sha256"] == progress["source_sha256"] == record["input_sha256"] and output_atoms.info["source_sha256"] == record["input_sha256"],
    "same_geometry": np.array_equal(input_atoms.numbers, output_atoms.numbers) and np.array_equal(input_atoms.positions, output_atoms.positions) and np.array_equal(input_atoms.cell.array, output_atoms.cell.array) and np.array_equal(input_atoms.pbc, output_atoms.pbc),
    "scf": summary["scf_converged"] is True and progress["scf_converged"] is True and progress["status"] == "complete" and f"Converged in {summary['scf_iterations']} steps" in log,
    "mpi": summary["mpi_ranks"] == progress["mpi_ranks"] == 4,
    "method": summary["method"] == record["method"] and summary["method"]["kpts"] == record["kpts"],
    "finite_values": np.isfinite(summary["native_energy_eV_cell"]) and np.isfinite(summary["free_energy_eV_cell"]) and forces.shape == (len(input_atoms),3) and np.isfinite(forces).all(),
    "energy_consistency": np.isclose(float(output_atoms.info["PW_PBE_energy_eV"]), summary["native_energy_eV_cell"], rtol=0, atol=1e-12) and np.isclose(float(output_atoms.info["free_energy_eV"]), summary["free_energy_eV_cell"], rtol=0, atol=1e-12),
    "control_only": summary["dataset_role"] == output_atoms.info["dataset_role"] == "NUMERICAL_CONTROL_ONLY",
}
assert all(checks.values()), checks
checksum_file = RUN / "SHA256SUMS.txt"
if checksum_file.exists():
    subprocess.run(["sha256sum", "-c", "SHA256SUMS.txt"], cwd=RUN, check=True, stdout=subprocess.DEVNULL)
result = {"status":"PASS" if all(checks.values()) else "FAIL", "label":args.label, "checks":checks}
print(json.dumps(result, indent=2))
