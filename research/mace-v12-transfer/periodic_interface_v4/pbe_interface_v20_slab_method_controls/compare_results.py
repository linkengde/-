"""Compare two verified fixed-geometry slab numerical controls."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from ase.geometry import find_mic
from ase.io import read

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
if len(sys.argv) != 3 or sys.argv[1] == sys.argv[2]:
    raise SystemExit("usage: compare_results.py LABEL_A LABEL_B")

new = json.loads((ROOT / "input_manifest.json").read_text())["records"]
labels = {r["label"]: r for r in new}
old_label = "AgSi_rigid_Ag_z_parent_k6x6x1_sigma0p10_dipoleXY_v19"
old_root = REPO / "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_gap_scan"
dev_label = "AgSi_v19_independent_dev_01"
dev_root = REPO / "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_independent_dev"


def get(label):
    if label in labels:
        folder = ROOT / "calculations" / label
        subprocess.run([sys.executable, str(ROOT / "verify_result.py"), label, str(folder)], check=True, stdout=subprocess.DEVNULL)
        summary = json.loads((folder / "summary.json").read_text())
        atoms = read(folder / (label + "_PW_PBE.extxyz"))
        method = summary["method"]
        return summary, atoms, method, folder, labels[label]["input_sha256"]
    if label == old_label:
        folder = old_root / "calculations" / label
        subprocess.run([sys.executable, str(old_root / "verify_result.py"), label, str(folder)], check=True, stdout=subprocess.DEVNULL)
        summary = json.loads((folder / "summary.json").read_text())
        atoms = read(folder / (label + "_PW_PBE.extxyz"))
        source_manifest = json.loads((old_root / "input_manifest.json").read_text())
        record = next(r for r in source_manifest["records"] if r["label"] == label)
        return summary, atoms, summary["method"], folder, record.get("input_sha256", record.get("source_sha256"))
    if label == dev_label:
        folder = dev_root / "calculations" / label
        subprocess.run([sys.executable, str(dev_root / "verify_result.py"), label, str(folder)], check=True, stdout=subprocess.DEVNULL)
        summary = json.loads((folder / "summary.json").read_text())
        atoms = read(folder / (label + "_PW_PBE.extxyz"))
        source_manifest = json.loads((dev_root / "input_manifest.json").read_text())
        record = next(r for r in source_manifest["records"] if r["label"] == label)
        return summary, atoms, summary["method"], folder, record["input_sha256"]
    raise SystemExit("Label is not in this scan or its registered baseline: " + label)


(sa, aa, ma, da, ha), (sb, ab, mb, db, hb) = get(sys.argv[1]), get(sys.argv[2])
assert np.array_equal(aa.numbers, ab.numbers)
assert np.array_equal(aa.positions, ab.positions) and np.array_equal(aa.cell.array, ab.cell.array)
assert np.array_equal(aa.pbc, ab.pbc) and np.array_equal(aa.arrays["lammps_id"], ab.arrays["lammps_id"])
symbols = np.asarray(aa.get_chemical_symbols())
marked = aa.arrays["central_pair"].astype(bool)
ia, ix = np.flatnonzero(marked & (symbols == "Ag")), np.flatnonzero(marked & (symbols == "Si"))
assert len(ia) == len(ix) == 1
ia, ix = int(ia[0]), int(ix[0])
unit = find_mic(aa.positions[ix] - aa.positions[ia], aa.cell, aa.pbc)[0]
unit /= np.linalg.norm(unit)
delta = ab.arrays["PW_PBE_forces"] - aa.arrays["PW_PBE_forces"]
metrics = {
    "native_energy_delta_meV_atom": (sb["energy_eV_cell"] - sa["energy_eV_cell"]) * 1000 / len(aa),
    "free_energy_delta_meV_atom": (sb["free_energy_eV_cell"] - sa["free_energy_eV_cell"]) * 1000 / len(aa),
    "force_vector_difference_RMS_eV_A": float(np.sqrt(np.mean(np.sum(delta * delta, axis=1)))),
    "maximum_atom_force_vector_difference_eV_A": float(np.linalg.norm(delta, axis=1).max()),
    "marked_Ag_Si_pair_projection_delta_eV_A": float(np.dot(delta[ix] - delta[ia], unit)),
}
budgets = json.loads((ROOT / "input_manifest.json").read_text())["numerical_budgets"]
checks = {
    "native_energy": abs(metrics["native_energy_delta_meV_atom"]) <= budgets["native_energy_meV_per_atom"],
    "free_energy": abs(metrics["free_energy_delta_meV_atom"]) <= budgets["free_energy_meV_per_atom"],
    "force_vector_RMS": metrics["force_vector_difference_RMS_eV_A"] <= budgets["force_vector_rms_eV_A"],
    "marked_pair_projection": abs(metrics["marked_Ag_Si_pair_projection_delta_eV_A"]) <= budgets["marked_pair_projection_eV_A"],
}
result = {
    "status": "PASS_WITHIN_PROVISIONAL_BUDGET" if all(checks.values()) else "SENSITIVITY_EXCEEDS_BUDGET",
    "labels": [sys.argv[1], sys.argv[2]],
    "methods": [ma, mb],
    "input_geometry_sha256": ha if ha == hb else None,
    "input_geometry_file_sha256": {sys.argv[1]: ha, sys.argv[2]: hb},
    "metrics": metrics,
    "budgets": budgets,
    "checks": checks,
    "archives_sha256": {
        label: {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.iterdir() if p.is_file() and p.name not in ("verification.json", "SHA256SUMS.txt") and p.suffix != ".gpw"}
        for label, folder in ((sys.argv[1], da), (sys.argv[2], db))
    },
    "role": "numerical_method_control_only_not_training_or_independent_validation",
}
out = ROOT / "comparisons"
out.mkdir(exist_ok=True)
path = out / ("comparison_" + sys.argv[1] + "__" + sys.argv[2] + ".json")
path.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
