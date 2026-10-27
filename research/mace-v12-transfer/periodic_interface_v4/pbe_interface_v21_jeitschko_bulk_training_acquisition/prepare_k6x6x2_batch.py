#!/usr/bin/env python3
"""Register a denser-mesh set after the Gamma diagnostic proved too coarse."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

from ase.io import read


ROOT = Path(__file__).resolve().parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_17(path: Path, atoms, label: str) -> None:
    lattice = " ".join(format(float(value), ".17g") for value in atoms.cell.array.ravel())
    header = f'Lattice="{lattice}" Properties=species:S:1:pos:R:3:local_new_id:I:1 pbc="T T T" config_type={label}'
    rows = []
    for symbol, position, local_id in zip(atoms.get_chemical_symbols(), atoms.positions, atoms.arrays["local_new_id"]):
        rows.append(symbol + " " + " ".join(format(float(value), ".17g") for value in position) + f" {int(local_id)}")
    path.write_text(f"{len(atoms)}\n{header}\n" + "\n".join(rows) + "\n")


def main() -> None:
    manifest_path = ROOT / "input_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest.get("revision", 1) == 1, "do not regenerate a revised registration"
    assert len(manifest["records"]) == 5
    new_records = []
    for original in manifest["records"]:
        old_label = original["label"]
        if old_label == "Ti3SiC2_JN1967_baseline_gamma_v21":
            original["status"] = "completed_control_only"
            original["role"] = "coarse_Gamma_method_diagnostic_not_training"
            original["launch_enabled"] = False
        else:
            original["status"] = "not_launched_gamma_mesh_too_coarse"
            original["role"] = "coarse_Gamma_method_diagnostic_not_training"
            original["launch_enabled"] = False

        label = old_label.replace("_gamma_v21", "_k6x6x2_v21")
        assert label != old_label
        old_path = ROOT / original["input"]
        atoms = read(old_path)
        atoms.info["config_type"] = label
        input_path = ROOT / "inputs" / f"{label}.extxyz"
        write_17(input_path, atoms, label)
        reread = read(input_path)
        assert len(reread) == 48 and reread.get_chemical_formula() == "C16Si8Ti24"
        assert reread.get_chemical_symbols() == atoms.get_chemical_symbols()
        assert reread.pbc.tolist() == [True, True, True]

        method = copy.deepcopy(original["method"])
        method["kpts"] = [6, 6, 2]
        method["kpoint_label"] = "6x6x2"
        new_record = copy.deepcopy(original)
        new_record.update({
            "label": label,
            "input": str(input_path.relative_to(ROOT)),
            "input_sha256": sha(input_path),
            "method": method,
            "role": "train_acquisition_candidate_pending_audit",
            "status": "queued",
            "launch_enabled": True,
            "mesh_selection_basis": "existing Ti3SiC2 48-atom fixed-parent controls found 4x4x2 through 10x10x4 converged with max force below 0.10 eV/A; Gamma on this new source gave ~2 eV/A forces and is retained only as a coarse-mesh diagnostic",
        })
        new_records.append(new_record)
    manifest["revision"] = 2
    manifest["mesh_revision_note"] = "Gamma baseline converged numerically but was rejected for training because its Ti/C force magnitudes were about 2 eV/A; a 6x6x2 mesh is used for the new candidate labels."
    manifest["records"].extend(new_records)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"revision": 2, "new_records": new_records}, indent=2))


if __name__ == "__main__":
    main()
