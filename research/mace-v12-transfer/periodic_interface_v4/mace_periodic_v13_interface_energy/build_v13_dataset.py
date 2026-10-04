#!/usr/bin/env python3
"""Add the three independently evaluated PW-PBE distance labels for v13."""
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
V12 = ROOT.parent / "mace_periodic_v12_interface_energy"
ADD = ROOT.parent / "pbe_interface_energy_additions_v12"
ARCHIVE = ADD / "calculations/validation_distance_scans"
DATA = ROOT / "data"
LABELS = [
    "AgC_validation_d2p30",
    "AgSi_validation_d2p60",
    "AgTi_validation_d2p60",
]
ELEMENTS = ["Ag", "C", "Si", "Ti"]

if DATA.exists() and any(DATA.iterdir()):
    raise SystemExit(f"Refusing to overwrite existing v13 dataset files in {DATA}")
DATA.mkdir(parents=True, exist_ok=True)

archive_manifest = json.loads((ARCHIVE / "archive_manifest.json").read_text())
archive_records = {r["label"]: r for r in archive_manifest["records"]}
input_manifest = json.loads((ADD / "validation_distance_input_manifest.json").read_text())
input_records = {r["label"]: r for r in input_manifest["records"]}
if set(LABELS) - archive_records.keys() or set(LABELS) - input_records.keys():
    raise SystemExit("The independent DFT archive or input manifest lacks required v13 labels")

parent_train_path = V12 / "data/train.extxyz"
parent_valid_path = V12 / "data/valid.extxyz"
parent_test_path = V12 / "data/test.extxyz"
train_bytes = parent_train_path.read_bytes()
parent_train = read(parent_train_path, index=":")
valid = read(parent_valid_path, index=":")
test = read(parent_test_path, index=":")
if (len(parent_train), len(valid), len(test)) != (23, 2, 3):
    raise SystemExit("Unexpected v12 split sizes; refusing to build v13 from an unrecognized base")

new_frames = []
new_label_records = []
for label in LABELS:
    folder = ARCHIVE / label
    source = folder / f"{label}_PW_PBE.extxyz"
    summary_path = folder / "summary.json"
    summary = json.loads(summary_path.read_text())
    record = archive_records[label]
    input_record = input_records[label]
    if summary.get("scf_converged") is not True or summary.get("mpi_ranks") != 4:
        raise SystemExit(f"DFT reference not converged on four MPI ranks: {label}")
    if sha256(source.read_bytes()).hexdigest() != record["extxyz_sha256"]:
        raise SystemExit(f"DFT extxyz hash differs from the validated archive manifest: {label}")
    if sha256(summary_path.read_bytes()).hexdigest() != record["summary_sha256"]:
        raise SystemExit(f"DFT summary hash differs from the validated archive manifest: {label}")
    if not (
        summary.get("source_sha256")
        == input_record["input_sha256"]
        == record["input_sha256"]
    ):
        raise SystemExit(f"DFT input hash does not match the validated manifests: {label}")

    atoms = read(source)
    energy = float(atoms.info["PW_PBE_energy_eV"])
    forces = np.asarray(atoms.arrays["PW_PBE_forces"], dtype=float)
    if not np.isfinite(energy) or forces.shape != (len(atoms), 3) or not np.isfinite(forces).all():
        raise SystemExit(f"Invalid DFT labels for {label}")
    atoms.info["REF_energy"] = energy
    atoms.arrays["REF_forces"] = forces.copy()
    atoms.info["config_type"] = f"{label}_periodic_PW_PBE_energy_force"
    atoms.info["source_method"] = "GPAW 26.7.0 PW-PBE 500 eV Gamma single point"
    atoms.info.pop("PW_PBE_energy_eV", None)
    atoms.arrays.pop("PW_PBE_forces", None)
    atoms.calc = None
    new_frames.append(atoms)
    new_label_records.append({
        "label": label,
        "config_type": atoms.info["config_type"],
        "input_sha256": input_record["input_sha256"],
        "extxyz_sha256": record["extxyz_sha256"],
        "summary_sha256": record["summary_sha256"],
        "energy_eV_cell": energy,
        "atoms": len(atoms),
        "formula": atoms.get_chemical_formula(),
        "central_pair": record["central_pair"],
        "scf_iterations": summary["scf_iterations"],
        "previous_role": "independent_v12_validation; now assigned to v13 training",
    })

buffer = StringIO()
write(buffer, new_frames, format="extxyz")
if train_bytes and not train_bytes.endswith(b"\n"):
    train_bytes += b"\n"
(DATA / "train.extxyz").write_bytes(train_bytes + buffer.getvalue().encode())
shutil.copy2(parent_valid_path, DATA / "valid.extxyz")
shutil.copy2(parent_test_path, DATA / "test.extxyz")

combined_train = read(DATA / "train.extxyz", index=":")
if len(combined_train) != 26:
    raise SystemExit(f"Unexpected v13 training frame count: {len(combined_train)}")
for before, after in zip(parent_train, combined_train[:len(parent_train)]):
    if before.get_chemical_symbols() != after.get_chemical_symbols():
        raise SystemExit("A pre-existing v12 training frame changed element order")
    if not np.array_equal(before.arrays.get("lammps_id"), after.arrays.get("lammps_id")):
        raise SystemExit("A pre-existing v12 training frame changed atom IDs")
    if not np.allclose(before.positions, after.positions, atol=1e-12, rtol=0):
        raise SystemExit("A pre-existing v12 training frame changed coordinates")
    if not np.allclose(before.cell.array, after.cell.array, atol=1e-12, rtol=0):
        raise SystemExit("A pre-existing v12 training frame changed its cell")
    if not np.array_equal(before.pbc, after.pbc):
        raise SystemExit("A pre-existing v12 training frame changed periodic boundaries")
    if "REF_energy" in before.info and not np.isclose(before.info["REF_energy"], after.info["REF_energy"], atol=1e-12, rtol=0):
        raise SystemExit("A pre-existing v12 training energy changed")
    if "REF_forces" in before.arrays and not np.allclose(before.arrays["REF_forces"], after.arrays["REF_forces"], atol=1e-12, rtol=0):
        raise SystemExit("A pre-existing v12 training force changed")
if (DATA / "valid.extxyz").read_bytes() != parent_valid_path.read_bytes():
    raise SystemExit("The frozen v12 validation split changed")
if (DATA / "test.extxyz").read_bytes() != parent_test_path.read_bytes():
    raise SystemExit("The frozen v12 test split changed")

energy_frames = [atoms for atoms in combined_train if "REF_energy" in atoms.info]
composition = np.asarray(
    [[atoms.get_chemical_symbols().count(element) for element in ELEMENTS] for atoms in energy_frames],
    dtype=float,
)
rank = int(np.linalg.matrix_rank(composition))
if rank != 4:
    raise SystemExit(f"v13 energy-composition matrix rank is {rank}/4")

manifest = {
    "model_version": "MACE v13 Ag/Ti/Si/C with three independently evaluated PW-PBE distance labels",
    "parent_v12_model_sha256": sha256((V12 / "checkpoints/MACE_periodic_v12_interface_energy_run-45.model").read_bytes()).hexdigest(),
    "parent_v12_train_sha256": sha256(parent_train_path.read_bytes()).hexdigest(),
    "parent_v12_valid_sha256": sha256(parent_valid_path.read_bytes()).hexdigest(),
    "parent_v12_test_sha256": sha256(parent_test_path.read_bytes()).hexdigest(),
    "split_sizes": {"train": len(combined_train), "valid": len(valid), "test": len(test)},
    "energy_label_counts": {
        "train": len(energy_frames),
        "valid": sum("REF_energy" in atoms.info for atoms in valid),
        "test": sum("REF_energy" in atoms.info for atoms in test),
    },
    "force_label_counts": {
        "train": sum("REF_forces" in atoms.arrays for atoms in combined_train),
        "valid": sum("REF_forces" in atoms.arrays for atoms in valid),
        "test": sum("REF_forces" in atoms.arrays for atoms in test),
    },
    "energy_composition_matrix": {
        "elements": ELEMENTS,
        "rank": rank,
        "singular_values": np.linalg.svd(composition, compute_uv=False).tolist(),
    },
    "added_training_labels": new_label_records,
    "frozen_splits": {
        "valid_byte_identical_to_v12": True,
        "test_byte_identical_to_v12": True,
        "test_frames": [atoms.info.get("config_type", "") for atoms in test],
    },
    "validation_caveat": "The three added labels were independent checks of v12, but they share fixed-distance midpoint motifs with external structures. Do not reuse them as v13 holdouts; obtain additional independent Ag-C and Ag-Si DFT interface checks before production conclusions.",
    "output_sha256": {
        f"{name}.extxyz": sha256((DATA / f"{name}.extxyz").read_bytes()).hexdigest()
        for name in ("train", "valid", "test")
    },
}
(DATA / "dataset_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
(DATA / "README.md").write_text(
    "# MACE v13 data\n\n"
    "This training set preserves all v12 training frames and appends the three newly verified Ag-C, Ag-Si, and Ag-Ti PW-PBE distance labels. The v12 validation and test files are copied byte-for-byte. The three added points are training data in v13 and must not be scored as independent v13 holdouts. See `dataset_manifest.json` for hashes and caveats.\n"
)
print(json.dumps(manifest, indent=2))
