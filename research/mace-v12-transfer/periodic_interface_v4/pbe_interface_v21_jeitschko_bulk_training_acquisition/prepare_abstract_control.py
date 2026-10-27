#!/usr/bin/env python3
"""Register the paper-abstract coordinate variant as a fixed-geometry control."""
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
SOURCE_CIF = REPO / "research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_jeitschko1967_source_audit_20261027/Ti3SiC2_Jeitschko1967_abstract_reconstruction.cif"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest_path = ROOT / "input_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest.get("revision") == 4
    method_record = next(r for r in manifest["records"] if r["label"] == "Ti3SiC2_JN1967_primitive_baseline_k10x10x2_mesh_control_v21")
    atoms = read(SOURCE_CIF)
    assert len(atoms) == 12 and atoms.get_chemical_formula() == "C4Si2Ti6" and atoms.info["spacegroup"].no == 194
    atoms.pbc = True
    atoms.info.clear()
    for key in list(atoms.arrays):
        if key not in {"numbers", "positions"}:
            del atoms.arrays[key]
    atoms.new_array("local_new_id", np.arange(1, 13, dtype=int))
    label = "Ti3SiC2_JN1967_abstract_primitive_baseline_k10x10x2_mesh_control_v21"
    atoms.info["config_type"] = label
    path = ROOT / "inputs" / f"{label}.extxyz"
    lattice = " ".join(format(float(v), ".17g") for v in atoms.cell.array.ravel())
    header = f'Lattice="{lattice}" Properties=species:S:1:pos:R:3:local_new_id:I:1 pbc="T T T" config_type={label}'
    lines = [str(len(atoms)), header]
    for symbol, pos, local_id in zip(atoms.get_chemical_symbols(), atoms.positions, atoms.arrays["local_new_id"]):
        lines.append(symbol + " " + " ".join(format(float(v), ".17g") for v in pos) + f" {int(local_id)}")
    path.write_text("\n".join(lines) + "\n")
    reread = read(path)
    assert np.array_equal(reread.positions, atoms.positions)
    method = copy.deepcopy(method_record["method"])
    manifest["records"].append({
        "label": label,
        "input": str(path.relative_to(ROOT)),
        "input_sha256": sha(path),
        "source_cif_sha256": sha(SOURCE_CIF),
        "source_paper_sha256": manifest["source"]["uploaded_paper_sha256"],
        "atom_count": 12,
        "formula": "C4Si2Ti6",
        "method": method,
        "role": "fixed_geometry_coordinate_variant_control_not_training",
        "launch_enabled": True,
        "status": "queued",
        "compared_with": method_record["input_sha256"],
        "coordinate_variant": "publisher abstract zTi=0.135, zC=0.5675",
    })
    manifest["revision"] = 5
    manifest["current_control"] = "compare full-text Fourier zC=0.570(4) against abstract zC=0.5675 at identical 10x10x2 mesh"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(label, sha(path), sha(SOURCE_CIF))


if __name__ == "__main__":
    main()
