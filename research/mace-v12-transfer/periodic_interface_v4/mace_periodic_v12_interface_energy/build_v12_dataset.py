#!/usr/bin/env python3
"""Replace four LCAO-only AgSi/AgC training frames with their converged PW-PBE labels."""
from pathlib import Path
import hashlib
import json
import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
V11 = ROOT.parent / "mace_periodic_v11_energyfocus"
PWROOT = ROOT.parent / "pbe_interface_energy_additions_v12/calculations"
DATA = ROOT / "data"
if DATA.exists() and any(DATA.iterdir()):
    raise SystemExit(f"Refusing to overwrite existing dataset files in {DATA}")
DATA.mkdir(parents=True, exist_ok=True)
replacements = {
    "AgSi_2p2A_LCAO": "AgSi_d2p2",
    "AgSi_2p8A_LCAO": "AgSi_d2p8",
    "AgC_secondary_2.0A_LCAO": "AgC_d2p0",
    "AgC_secondary_2.4A_LCAO": "AgC_d2p4",
}
train = read(V11 / "data/train.extxyz", index=":")
valid = read(V11 / "data/valid.extxyz", index=":")
test = read(V11 / "data/test.extxyz", index=":")
new_train = []
seen = set()
for atoms in train:
    ctype = atoms.info.get("config_type")
    if ctype in replacements:
        if ctype in seen:
            raise SystemExit(f"Duplicate source config {ctype}")
        seen.add(ctype)
        continue
    new_train.append(atoms)
if seen != set(replacements):
    raise SystemExit(f"Did not find all replacement frames: missing={set(replacements)-seen}")

new_labels = []
for old_type, label in replacements.items():
    path = PWROOT / label / f"{label}_PW_PBE.extxyz"
    summary = json.loads((PWROOT / label / "summary.json").read_text())
    if summary.get("scf_converged") is not True:
        raise SystemExit(f"Non-converged DFT reference: {label}")
    atoms = read(path)
    if "PW_PBE_energy_eV" not in atoms.info or "PW_PBE_forces" not in atoms.arrays:
        raise SystemExit(f"PW energy/force properties absent in {path}")
    atoms.info["REF_energy"] = float(atoms.info["PW_PBE_energy_eV"])
    atoms.arrays["REF_forces"] = np.asarray(atoms.arrays["PW_PBE_forces"], dtype=float).copy()
    atoms.info["config_type"] = f"{label}_periodic_PW_PBE_energy_force"
    atoms.info["source_method"] = "GPAW 26.7.0 PW-PBE 500 eV Gamma single point"
    new_train.append(atoms)
    new_labels.append({
        "old_config_replaced": old_type,
        "new_config_type": atoms.info["config_type"],
        "label": label,
        "summary": str(PWROOT / label / "summary.json"),
        "extxyz_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "energy_eV_cell": float(atoms.info["REF_energy"]),
        "atoms": len(atoms),
        "formula": atoms.get_chemical_formula(),
    })

for split, frames in [("train", new_train), ("valid", valid), ("test", test)]:
    write(DATA / f"{split}.extxyz", frames, format="extxyz")

elements = ["Ag", "C", "Si", "Ti"]
energy_frames = [a for a in new_train if "REF_energy" in a.info]
composition = np.array([[a.get_chemical_symbols().count(e) for e in elements] for a in energy_frames], dtype=float)
rank = int(np.linalg.matrix_rank(composition))
if rank != 4:
    raise SystemExit(f"Training energy composition rank is {rank}/4")
source_hashes = {
    str(V11 / f"data/{name}.extxyz"): hashlib.sha256((V11 / f"data/{name}.extxyz").read_bytes()).hexdigest()
    for name in ["train", "valid", "test"]
}
manifest = {
    "model_version": "MACE v12 AgSi/AgC PW-PBE energy additions",
    "parent_v11_data": str(V11 / "data"),
    "split_sizes": {"train": len(new_train), "valid": len(valid), "test": len(test)},
    "energy_label_counts": {"train": len(energy_frames), "valid": sum("REF_energy" in a.info for a in valid), "test": sum("REF_energy" in a.info for a in test)},
    "force_label_counts": {"train": sum("REF_forces" in a.arrays for a in new_train), "valid": sum("REF_forces" in a.arrays for a in valid), "test": sum("REF_forces" in a.arrays for a in test)},
    "energy_composition_matrix": {"elements": elements, "rank": rank, "singular_values": np.linalg.svd(composition, compute_uv=False).tolist()},
    "replaced_lcao_force_frames": new_labels,
    "holdouts_preserved": ["v11 valid and test splits byte/structure-preserved", "external AgSi 2.4541 A PW-PBE", "external AgC 2.2026 A PW-PBE", "AgTi 2.70 A test"],
    "source_hashes": source_hashes,
    "caveat": "AgSi and AgC external structures are single midpoint/interface checks; they are not liquid or production-scale validation.",
}
manifest["output_sha256"] = {f"{name}.extxyz": hashlib.sha256((DATA / f"{name}.extxyz").read_bytes()).hexdigest() for name in ["train", "valid", "test"]}
(DATA / "dataset_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
(DATA / "README.md").write_text(
    "# MACE v12 data\n\n"
    "This split preserves v11 validation/test data and replaces four AgSi/AgC force-only LCAO training frames with converged periodic PW-PBE energy+force labels at the same geometries. The independent AgSi and AgC midpoint structures remain excluded. See `dataset_manifest.json` for hashes and SCF summaries.\n"
)
print(json.dumps(manifest, indent=2))
