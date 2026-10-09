#!/usr/bin/env python3
"""Generate two hash-pinned, label-free AgSi electrostatic control geometries."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
from ase.io import read, write


HERE = Path(__file__).resolve().parent
REPO = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], cwd=HERE, text=True).strip())
MANIFEST_REL = Path(
    "research/mace-v12-transfer/periodic_interface_v4/"
    "pbe_interface_v19_gap_scan/input_manifest.json")
SOURCE_LABEL = "AgSi_rigid_Ag_z_parent_k6x6x1_sigma0p10_dipoleXY_v19"
SOURCE_MANIFEST_SHA256 = "fa7a0508ab2c6d176f32244fae6cce399a5904395d78d791bc3ebc1123e39f76"
SOURCE_INPUT_SHA256 = "73eff17b1913b653017647e9af9c9681e340f0b8067936419f3377963a3bc02c"
SOURCE_PARENT_SHA256 = "8fe44bee35e17ebfc433e059925a7a85f0dbb4e775090a75ae1883b3eaa6c36c"
OFF_NAME = "AgSi_26A_vacuum_parent_dipole_OFF_PROPOSAL.extxyz"
ON_NAME = "AgSi_36A_vacuum_parent_dipole_ON_PROPOSAL.extxyz"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def check_label_free(atoms) -> bool:
    keys = list(atoms.info) + list(atoms.arrays)
    forbidden = ("energy", "force", "stress", "free_energy")
    return not any(any(term in key.lower() for term in forbidden) for key in keys)


def compare_identity(source, proposed, expected_shift_z: float) -> dict:
    assert np.array_equal(source.numbers, proposed.numbers)
    assert np.array_equal(source.pbc, proposed.pbc)
    assert source.info == proposed.info
    assert set(source.arrays) == set(proposed.arrays)
    for key in source.arrays:
        if key != "positions":
            assert np.array_equal(source.arrays[key], proposed.arrays[key]), key
    expected_positions = source.positions.copy()
    expected_positions[:, 2] += expected_shift_z
    assert np.allclose(proposed.positions, expected_positions, rtol=0.0, atol=1e-12)
    return {
        "atomic_numbers_equal": True,
        "atom_ids_equal": True,
        "pair_markers_equal": True,
        "all_nonposition_arrays_equal": True,
        "atom_order_preserved": True,
        "pbc_equal": True,
        "info_equal": True,
        "z_translation_A": expected_shift_z,
    }


def main() -> None:
    manifest_path = REPO / MANIFEST_REL
    assert sha256(manifest_path) == SOURCE_MANIFEST_SHA256, "Registered manifest changed"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["launch_enabled"] is False
    matches = [r for r in manifest["records"] if r.get("label") == SOURCE_LABEL]
    assert len(matches) == 1
    record = matches[0]
    assert record["launch_enabled"] is False
    assert record["input_sha256"] == SOURCE_INPUT_SHA256
    assert record["parent_geometry_sha256"] == SOURCE_PARENT_SHA256
    assert record["pbc"] == [True, True, False]
    assert record["poissonsolver"] == {"dipolelayer": "xy"}

    source_path = manifest_path.parent / record["input"]
    assert sha256(source_path) == SOURCE_INPUT_SHA256, "Registered source input changed"
    source = read(source_path, format="extxyz")
    assert len(source) == 52
    assert "lammps_id" in source.arrays and "central_pair" in source.arrays
    assert source.pbc.tolist() == [True, True, False]
    assert check_label_free(source), "Source unexpectedly contains labels"

    off_path = HERE / OFF_NAME
    on_path = HERE / ON_NAME
    manifest_out = HERE / "proposal_manifest.json"
    verification_out = HERE / "verification.json"
    for path in (off_path, on_path, manifest_out, verification_out):
        assert not path.exists(), f"Refusing to overwrite {path}"

    # OFF changes the proposed solver setting only; keep the source bytes exact.
    shutil.copyfile(source_path, off_path)
    off = read(off_path, format="extxyz")
    off_checks = compare_identity(source, off, 0.0)
    assert np.array_equal(source.cell.array, off.cell.array)
    assert check_label_free(off)

    # ON expands z by 10 A and translates the complete geometry by +5 A.
    expanded = source.copy()
    cell = source.cell.array.copy()
    cell[2, 2] += 10.0
    expanded.set_cell(cell, scale_atoms=False)
    expanded.positions[:, 2] += 5.0
    write(on_path, expanded, format="extxyz")
    on = read(on_path, format="extxyz")
    on_checks = compare_identity(source, on, 5.0)
    cell_delta = on.cell.array - source.cell.array
    assert np.allclose(cell_delta[:2], 0.0, rtol=0.0, atol=1e-12)
    assert np.allclose(cell_delta[2, :2], 0.0, rtol=0.0, atol=1e-12)
    assert np.isclose(cell_delta[2, 2], 10.0, rtol=0.0, atol=1e-12)
    assert check_label_free(on)

    source_span = float(np.ptp(source.positions[:, 2]))
    off_method = deepcopy(record["method"])
    off_method["poissonsolver"] = {}
    on_method = deepcopy(record["method"])
    method_without_solver = deepcopy(off_method)
    method_without_solver.pop("poissonsolver")
    registered_without_solver = deepcopy(on_method)
    registered_without_solver.pop("poissonsolver")
    assert method_without_solver == registered_without_solver
    assert on_method == record["method"]
    source_rel = str(source_path.relative_to(REPO))
    parent_span = source_span

    proposals = [
        {
            "proposal_id": "AgSi_26A_vacuum_parent_dipole_OFF_PROPOSAL",
            "input": OFF_NAME,
            "input_sha256": sha256(off_path),
            "source_label": SOURCE_LABEL,
            "source_input": source_rel,
            "source_input_sha256": SOURCE_INPUT_SHA256,
            "source_manifest": str(MANIFEST_REL),
            "source_manifest_sha256": SOURCE_MANIFEST_SHA256,
            "source_parent_geometry_sha256": SOURCE_PARENT_SHA256,
            "method": off_method,
            "pbc": [True, True, False],
            "atom_count": len(off),
            "formula": off.get_chemical_formula(),
            "cell_z_A": float(off.cell.lengths()[2]),
            "coordinate_span_z_A": float(np.ptp(off.positions[:, 2])),
            "coordinate_span_complement_A": float(off.cell.lengths()[2] - np.ptp(off.positions[:, 2])),
            "role": "boundary-control geometry proposal only; same-CIF training/numerical diagnostic, not independent validation",
            "launch_enabled": False,
            "owner": None,
            "identity_checks": off_checks,
            "purpose": "Same-geometry/same-cell contrast against the registered dipole-ON parent; only the proposed Poisson dipole correction is disabled."
        },
        {
            "proposal_id": "AgSi_36A_vacuum_parent_dipole_ON_PROPOSAL",
            "input": ON_NAME,
            "input_sha256": sha256(on_path),
            "source_label": SOURCE_LABEL,
            "source_input": source_rel,
            "source_input_sha256": SOURCE_INPUT_SHA256,
            "source_manifest": str(MANIFEST_REL),
            "source_manifest_sha256": SOURCE_MANIFEST_SHA256,
            "source_parent_geometry_sha256": SOURCE_PARENT_SHA256,
            "method": on_method,
            "pbc": [True, True, False],
            "atom_count": len(on),
            "formula": on.get_chemical_formula(),
            "cell_z_A": float(on.cell.lengths()[2]),
            "coordinate_span_z_A": float(np.ptp(on.positions[:, 2])),
            "coordinate_span_complement_A": float(on.cell.lengths()[2] - np.ptp(on.positions[:, 2])),
            "geometry_change": {"cell_z_delta_A": 10.0, "global_z_translation_A": 5.0},
            "role": "vacuum-size boundary-control geometry proposal only; same-CIF training/numerical diagnostic, not independent validation",
            "launch_enabled": False,
            "owner": None,
            "identity_checks": on_checks,
            "purpose": "Same dipole-ON method and same relative slab geometry as the registered 26-A parent, with only z-cell length increased by 10 A and all atoms shifted +5 A."
        }
    ]
    proposal_doc = {
        "task": "v19_gap_scan_dipole_vacuum_control_proposals_window_b",
        "status": "geometry_proposals_only",
        "launch_enabled": False,
        "owner": None,
        "source_record": record,
        "source_manifest_sha256": SOURCE_MANIFEST_SHA256,
        "source_input_sha256": SOURCE_INPUT_SHA256,
        "source_parent_geometry_sha256": SOURCE_PARENT_SHA256,
        "source_atom_count": len(source),
        "source_formula": source.get_chemical_formula(),
        "source_cell_z_A": float(source.cell.lengths()[2]),
        "source_coordinate_span_z_A": source_span,
        "source_coordinate_span_complement_A": float(source.cell.lengths()[2] - source_span),
        "source_has_energy_or_force_labels": False,
        "proposals": proposals
    }
    verification = {
        "status": "PASS",
        "checks": {
            "manifest_hash_pinned": True,
            "registered_parent_record_exactly_one": True,
            "registered_source_hash_matches": True,
            "source_is_label_free": True,
            "proposal_off_bytes_equal_source": sha256(off_path) == SOURCE_INPUT_SHA256,
            "off_geometry_cell_pbc_identity": True,
            "off_only_method_change_is_poissonsolver": True,
            "on_atom_order_ids_pair_markers_and_other_arrays_preserved": True,
            "on_cell_z_expanded_10_A": True,
            "on_all_atoms_translated_plus_5_A": True,
            "on_method_identical_to_registered_dipole_on_method": True,
            "proposals_label_free": True,
            "both_launch_disabled_and_unowned": True
        },
        "source": {"input": source_rel, "input_sha256": SOURCE_INPUT_SHA256, "manifest": str(MANIFEST_REL), "manifest_sha256": SOURCE_MANIFEST_SHA256},
        "proposal_hashes": {OFF_NAME: sha256(off_path), ON_NAME: sha256(on_path)},
        "vacuum_measure": "cell_z - (max atom z - min atom z); coordinate-span complement, not electron-density-defined vacuum"
    }
    write_json(manifest_out, proposal_doc)
    write_json(verification_out, verification)
    print(json.dumps({
        "source_sha256": SOURCE_INPUT_SHA256,
        "off_sha256": sha256(off_path),
        "on_sha256": sha256(on_path),
        "off_cell_z_A": float(off.cell.lengths()[2]),
        "on_cell_z_A": float(on.cell.lengths()[2]),
        "off_coordinate_span_complement_A": float(off.cell.lengths()[2] - np.ptp(off.positions[:, 2])),
        "on_coordinate_span_complement_A": float(on.cell.lengths()[2] - np.ptp(on.positions[:, 2])),
        "checks": "PASS"
    }, indent=2))


if __name__ == "__main__":
    main()
