#!/usr/bin/env python3
"""Read-only geometry, identity, provenance, and config validation."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase.geometry import find_mic
from ase.io import read

ROOT = Path(__file__).resolve().parent
INTERFACE_PARENT = ROOT / "source/Ag_Ti3SiC2_12p5198deg_217atom_dft_lite.extxyz"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def min_pair(atoms, a_indices, b_indices):
    ai = np.asarray(a_indices, dtype=int)
    bi = np.asarray(b_indices, dtype=int)
    ii, jj = np.meshgrid(ai, bi, indexing="ij")
    keep = ii.ravel() != jj.ravel()
    ii, jj = ii.ravel()[keep], jj.ravel()[keep]
    dv, _ = find_mic(atoms.positions[jj] - atoms.positions[ii], atoms.cell, atoms.pbc)
    lengths = np.linalg.norm(dv, axis=1)
    k = int(np.argmin(lengths))
    return float(lengths[k]), int(ii[k]), int(jj[k])


manifest = json.loads((ROOT / "INPUT_MANIFEST.json").read_text())
parent = read(INTERFACE_PARENT)
assert sha(INTERFACE_PARENT) == "945c5e214c1ceaeba4125ff785f4092e14858c75a985061a7841da3e5ce6a9a2"
assert manifest["execution_authorized"] is False
assert manifest["dft_started"] is False
assert manifest["fitting_started"] is False
assert manifest["tscbop_modified"] is False
assert len(manifest["backbone_tasks"]) == 6
assert len(manifest["interface_scan_tasks"]) == 9
results = []
for task in manifest["backbone_tasks"] + manifest["interface_scan_tasks"]:
    structure = ROOT / task["structure_path"]
    config_path = ROOT / task["config_path"]
    assert sha(structure) == task["structure_sha256"]
    assert sha(config_path) == task["config_sha256"]
    atoms = read(structure)
    assert len(atoms) == int(task["natoms"])
    assert np.isfinite(atoms.positions).all() and np.isfinite(atoms.cell.array).all()
    assert np.array_equal(atoms.pbc, np.ones(3, dtype=bool))
    assert not any(k in atoms.info for k in ("energy", "REF_energy", "PW_PBE_energy_eV"))
    assert not any(k in atoms.arrays for k in ("forces", "REF_forces", "PW_PBE_forces"))
    symbols = np.asarray(atoms.get_chemical_symbols())
    comp = {s: int((symbols == s).sum()) for s in ("Ag", "Ti", "Si", "C") if (symbols == s).any()}
    assert comp == task["composition"]
    if task["parent_family"] == "TSC_KISI_COD9009647_BASELINE_DERIVED_FAMILY":
        assert len(atoms) == 48 and comp == {"Ti": 24, "Si": 8, "C": 16}
        assert task["method"]["kpoints"] == [10, 10, 4]
        min_d = min_pair(atoms, np.arange(len(atoms)), np.arange(len(atoms)))
        result = {"task_id": task["task_id"], "natoms": len(atoms), "composition": comp, "minimum_periodic_distance_A": min_d[0], "family": task["parent_family"], "geometry_sha256": task["structure_sha256"], "config_sha256": task["config_sha256"], "status": "PASS_GEOMETRY_ONLY"}
    else:
        assert len(atoms) == 217 and comp == {"Ag": 84, "Ti": 57, "Si": 38, "C": 38}
        assert task["method"]["kpoints"] == [2, 2, 1]
        assert np.array_equal(atoms.numbers, parent.numbers)
        assert np.array_equal(atoms.cell.array, parent.cell.array)
        assert np.array_equal(atoms.pbc, parent.pbc)
        for key in ("atom_id", "parent_atom_id", "region", "layer_id", "parent_layer_id"):
            assert np.array_equal(atoms.arrays[key], parent.arrays[key]), f"atom mapping changed: {key}"
        ag_mask = symbols == "Ag"
        expected_delta = float(task["geometry_operation"]["rigid_Ag_z_translation_A"])
        delta = atoms.positions - parent.positions
        assert np.allclose(delta[~ag_mask], 0.0, atol=1e-12)
        assert np.allclose(delta[ag_mask, :2], 0.0, atol=1e-12)
        assert np.allclose(delta[ag_mask, 2], expected_delta, atol=1e-12)
        ai = np.where(symbols == "Ag")[0]
        si = np.where(symbols == "Si")[0]
        # Surface gap is between top Si plane and bottom Ag plane, not an optimized d0.
        gap = float(atoms.positions[ai, 2].min() - atoms.positions[si, 2].max())
        target_gap = float(task["geometry_operation"]["target_gap_A"])
        assert abs(gap - target_gap) < 2e-6
        mins = {}
        for el in ("Si", "Ti", "C"):
            bi = np.where(symbols == el)[0]
            mins[f"Ag-{el}"] = min_pair(atoms, ai, bi)[0]
        assert min(mins.values()) >= 1.5
        all_min = min_pair(atoms, np.arange(len(atoms)), np.arange(len(atoms)))
        result = {"task_id": task["task_id"], "natoms": len(atoms), "composition": comp, "plane_gap_A": gap, "minimum_Ag_TSC_distance_A": min(mins.values()), "minimum_Ag_Si_A": mins["Ag-Si"], "minimum_Ag_Ti_A": mins["Ag-Ti"], "minimum_Ag_C_A": mins["Ag-C"], "minimum_periodic_distance_A": all_min[0], "family": task["parent_family"], "geometry_sha256": task["structure_sha256"], "config_sha256": task["config_sha256"], "status": "PASS_GEOMETRY_ONLY"}
    results.append(result)
out = {
    "status": "PASS_GEOMETRY_ONLY",
    "dft_started": False,
    "fitting_started": False,
    "tscbop_modified": False,
    "structures_checked": len(results),
    "backbone_count": 6,
    "interface_gap_scan_count": 9,
    "interface_minimum_plane_gap_A": min(x["plane_gap_A"] for x in results if "plane_gap_A" in x),
    "interface_maximum_plane_gap_A": max(x["plane_gap_A"] for x in results if "plane_gap_A" in x),
    "items": results,
}
(ROOT / "GEOMETRY_VALIDATION.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(f"PASS: {len(results)} geometries, no labels, no scientific calculations")
