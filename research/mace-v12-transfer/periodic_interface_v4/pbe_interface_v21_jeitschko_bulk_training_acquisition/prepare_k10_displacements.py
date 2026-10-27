#!/usr/bin/env python3
"""Register force-rich displacement labels at the higher k mesh selected by controls."""
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_17(path, atoms, label):
    lattice = " ".join(format(float(v), ".17g") for v in atoms.cell.array.ravel())
    header = f'Lattice="{lattice}" Properties=species:S:1:pos:R:3:local_new_id:I:1 pbc="T T T" config_type={label}'
    rows = []
    for symbol, pos, local_id in zip(atoms.get_chemical_symbols(), atoms.positions, atoms.arrays["local_new_id"]):
        rows.append(symbol + " " + " ".join(format(float(v), ".17g") for v in pos) + f" {int(local_id)}")
    path.write_text(f"{len(atoms)}\n{header}\n" + "\n".join(rows) + "\n")


def main():
    manifest_path = ROOT / "input_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest.get("revision") == 5
    baseline = next(r for r in manifest["records"] if r["label"] == "Ti3SiC2_JN1967_primitive_baseline_k10x10x2_mesh_control_v21")
    for record in manifest["records"]:
        if record["label"].startswith("Ti3SiC2_JN1967_primitive_") and record["method"].get("kpts") == [6, 6, 2]:
            record["status"] = "completed_underconverged_mesh_control"
            record["role"] = "method_diagnostic_not_training"
            record["launch_enabled"] = False

    new_records = []
    source_labels = [
        "Ti3SiC2_JN1967_primitive_Ti4f_id0003_z_m002A_k6x6x2_v21",
        "Ti3SiC2_JN1967_primitive_Ti4f_id0003_z_p002A_k6x6x2_v21",
        "Ti3SiC2_JN1967_primitive_C4f_id0009_z_m002A_k6x6x2_v21",
        "Ti3SiC2_JN1967_primitive_C4f_id0009_z_p002A_k6x6x2_v21",
    ]
    for source_label in source_labels:
        original = next(r for r in manifest["records"] if r["label"] == source_label)
        label = source_label.replace("_k6x6x2_v21", "_k10x10x2_v21")
        atoms = read(ROOT / original["input"])
        atoms.info["config_type"] = label
        path = ROOT / "inputs" / f"{label}.extxyz"
        write_17(path, atoms, label)
        reread = read(path)
        assert len(reread) == 12 and reread.get_chemical_formula() == "C4Si2Ti6"
        assert np.array_equal(reread.positions, atoms.positions)
        assert np.array_equal(reread.arrays["local_new_id"], atoms.arrays["local_new_id"])
        method = copy.deepcopy(baseline["method"])
        new_records.append({
            "label": label,
            "input": str(path.relative_to(ROOT)),
            "input_sha256": sha(path),
            "source_cif_sha256": baseline["source_cif_sha256"],
            "source_paper_sha256": baseline["source_paper_sha256"],
            "atom_count": 12,
            "formula": "C4Si2Ti6",
            "site_id": original.get("site_id"),
            "delta_z_A": original.get("delta_z_A"),
            "method": method,
            "role": "train_acquisition_candidate_pending_audit",
            "launch_enabled": True,
            "status": "queued",
            "mesh_selected_after_control": "10x10x2: 6x6x2-to-10x10x2 force-vector RMS difference was 0.06198 eV/A on identical geometry",
        })
    manifest["records"].extend(new_records)
    manifest["revision"] = 6
    manifest["current_batch"] = "four Ti/C fixed-cell displacement labels from the 1967 parent at PBE/PW500/10x10x2/Fermi0.10/TTT"
    manifest["mesh_revision_note"] = "The 6x6x2 force labels are retained only as method controls; the 10x10x2 mesh is selected for the new training-candidate displacement labels after identical-geometry comparison."
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"revision": 6, "registered": new_records}, indent=2))


if __name__ == "__main__":
    main()
