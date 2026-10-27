#!/usr/bin/env python3
"""Register a denser-mesh fixed-geometry control for the 1967 primitive cell."""
import copy
import hashlib
import json
from pathlib import Path

from ase.io import read


ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest_path = ROOT / "input_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest.get("revision") == 3
    source_record = next(r for r in manifest["records"] if r["label"] == "Ti3SiC2_JN1967_primitive_baseline_k6x6x2_v21")
    label = "Ti3SiC2_JN1967_primitive_baseline_k10x10x2_mesh_control_v21"
    atoms = read(ROOT / source_record["input"])
    atoms.info["config_type"] = label
    path = ROOT / "inputs" / f"{label}.extxyz"
    lattice = " ".join(format(float(v), ".17g") for v in atoms.cell.array.ravel())
    header = f'Lattice="{lattice}" Properties=species:S:1:pos:R:3:local_new_id:I:1 pbc="T T T" config_type={label}'
    rows = []
    for symbol, pos, local_id in zip(atoms.get_chemical_symbols(), atoms.positions, atoms.arrays["local_new_id"]):
        rows.append(symbol + " " + " ".join(format(float(v), ".17g") for v in pos) + f" {int(local_id)}")
    path.write_text(f"{len(atoms)}\n{header}\n" + "\n".join(rows) + "\n")
    reread = read(path)
    assert len(reread) == 12 and reread.get_chemical_symbols() == atoms.get_chemical_symbols()
    assert reread.positions.tolist() == atoms.positions.tolist()
    method = copy.deepcopy(source_record["method"])
    method["kpts"] = [10, 10, 2]
    method["kpoint_label"] = "10x10x2"
    manifest["records"].append({
        "label": label,
        "input": str(path.relative_to(ROOT)),
        "input_sha256": sha(path),
        "source_cif_sha256": source_record["source_cif_sha256"],
        "source_paper_sha256": source_record["source_paper_sha256"],
        "atom_count": 12,
        "formula": "C4Si2Ti6",
        "method": method,
        "role": "fixed_geometry_kpoint_control_not_training",
        "launch_enabled": True,
        "status": "queued",
        "comparison_geometry": source_record["input_sha256"],
    })
    manifest["revision"] = 4
    manifest["current_control"] = "same Jeitschko primitive geometry at 10x10x2 to diagnose 6x6x2 force residual"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(label, sha(path))


if __name__ == "__main__":
    main()
