#!/usr/bin/env python3
"""Register the 12-atom primitive-cell set after cancelling the costly 48-atom mesh run."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
SOURCE_CIF = REPO / "research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_jeitschko1967_source_audit_20261027/Ti3SiC2_Jeitschko1967_fulltext_Fourier_variant.cif"


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
    assert manifest.get("revision") == 2 and len(manifest["records"]) == 10
    for record in manifest["records"]:
        if record["label"] == "Ti3SiC2_JN1967_baseline_gamma_v21":
            assert record["status"] == "completed_control_only"
        if record["label"] == "Ti3SiC2_JN1967_baseline_k6x6x2_v21":
            record["status"] = "aborted_planned_primitive_cell_reduction"
            record["role"] = "interrupted_resource_estimate_not_label"
            record["launch_enabled"] = False
        elif record["label"].endswith("_k6x6x2_v21") and record["status"] == "queued":
            record["status"] = "not_launched_after_primitive_cell_replan"
            record["role"] = "interrupted_resource_estimate_not_label"
            record["launch_enabled"] = False

    primitive = read(SOURCE_CIF)
    assert len(primitive) == 12 and primitive.get_chemical_formula() == "C4Si2Ti6"
    primitive.pbc = True
    primitive.info.clear()
    for key in list(primitive.arrays):
        if key not in {"numbers", "positions"}:
            del primitive.arrays[key]
    primitive.new_array("local_new_id", np.arange(1, 13, dtype=int))
    baseline_positions = primitive.positions.copy()
    source_hash = sha(SOURCE_CIF)
    cases = [{"label": "Ti3SiC2_JN1967_primitive_baseline_k6x6x2_v21", "site_id": None, "delta_z_A": 0.0}]
    for element, local_id in [("Ti4f", 3), ("C4f", 9)]:
        for sign in (-1, 1):
            cases.append({
                "label": f"Ti3SiC2_JN1967_primitive_{element}_id{local_id:04d}_z_{'m' if sign < 0 else 'p'}002A_k6x6x2_v21",
                "site_id": local_id,
                "delta_z_A": sign * 0.02,
            })

    method_template = copy.deepcopy(manifest["method"])
    method_template["kpts"] = [6, 6, 2]
    method_template["kpoint_label"] = "6x6x2"
    new_records = []
    for case in cases:
        label = case["label"]
        frame = primitive.copy()
        if case["site_id"] is not None:
            idx = case["site_id"] - 1
            expected_symbol = "Ti" if case["site_id"] == 3 else "C"
            assert frame[idx].symbol == expected_symbol
            frame.positions[idx, 2] += case["delta_z_A"]
        frame.info.update({
            "config_type": label,
            "proposal_role": "train_acquisition_candidate_pending_audit",
            "source_family": "Jeitschko_Nowotny_1967_single_crystal",
            "source_doi": "10.1007/BF00899949",
            "source_cif_sha256": source_hash,
            "displaced_local_new_id": -1 if case["site_id"] is None else case["site_id"],
            "delta_z_A": case["delta_z_A"],
        })
        path = ROOT / "inputs" / f"{label}.extxyz"
        write_17(path, frame, label)
        reread = read(path)
        assert len(reread) == 12 and reread.get_chemical_formula() == "C4Si2Ti6"
        assert reread.get_chemical_symbols() == frame.get_chemical_symbols()
        assert np.array_equal(reread.positions, frame.positions)
        assert np.array_equal(reread.cell.array, frame.cell.array)
        assert np.array_equal(reread.pbc, frame.pbc)
        assert np.array_equal(reread.arrays["local_new_id"], frame.arrays["local_new_id"])
        new_records.append({
            **case,
            "input": str(path.relative_to(ROOT)),
            "input_sha256": sha(path),
            "source_cif_sha256": source_hash,
            "source_paper_sha256": manifest["source"]["uploaded_paper_sha256"],
            "atom_count": 12,
            "formula": reread.get_chemical_formula(),
            "method": copy.deepcopy(method_template),
            "role": "train_acquisition_candidate_pending_audit",
            "launch_enabled": True,
            "status": "queued",
        })
    for case in cases[1:]:
        changed = np.flatnonzero(np.any(np.abs(read(ROOT / "inputs" / f"{case['label']}.extxyz").positions - baseline_positions) > 1e-12, axis=1))
        assert changed.tolist() == [case["site_id"] - 1]

    manifest["revision"] = 3
    manifest["current_batch"] = "12-atom Jeitschko primitive cell; PBE/PW500/6x6x2/sigma0.10/TTT; one fixed-cell baseline plus four Ti/C z displacements"
    manifest["revision_note"] = "The 48-atom 6x6x2 preflight reached only SCF iteration 1 after about 4.8 min at 4.0 GiB. It was deliberately stopped before convergence; no OOM occurred. The production batch uses the equivalent 12-atom primitive cell to reduce wall time and retains dense 6x6x2 sampling."
    manifest["records"].extend(new_records)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"revision": 3, "new_records": new_records}, indent=2))


if __name__ == "__main__":
    main()
