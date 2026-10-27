#!/usr/bin/env python3
"""Read-only role/provenance audit for the four Jeitschko 1967 bulk labels.

Only the registered candidate geometries and geometry columns from approved
train36/dev3 are read.  No test or sealed path is opened.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

REPO = Path(__file__).resolve().parents[6]
TASK = REPO / "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_bulk_training_acquisition"
REPORT = Path(__file__).resolve().parent
HELPER = REPO / "research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set/generate_split_candidates.py"
DATA = REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v19_residual_response/data"
VERIFIER_PYTHON = Path("/workspace/.venvs/gpaw-mpi/bin/python")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_helper():
    spec = importlib.util.spec_from_file_location("geometry_helper", HELPER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> None:
    manifest = json.loads((TASK / "input_manifest.json").read_text())
    approved = json.loads((DATA / "dataset_manifest.json").read_text())
    dev_manifest = json.loads((REPO / "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_independent_dev/input_manifest.json").read_text())
    helper = load_helper()
    records = [r for r in manifest["records"] if r.get("role") == "train_acquisition_candidate_pending_audit"]
    assert len(records) == 4

    train_frames = helper.extxyz_geometry_frames(DATA / approved["train_path"])
    dev_frames = helper.extxyz_geometry_frames(DATA / approved["development_path"])
    assert len(train_frames) == 36 and len(dev_frames) == 3
    input_checks = []
    candidates = []
    for record in records:
        input_path = TASK / record["input"]
        result_dir = TASK / "calculations" / record["label"]
        result_path = result_dir / f"{record['label']}_PW_PBE.extxyz"
        input_digest = sha(input_path)
        assert input_digest == record["input_sha256"]
        verification = json.loads((result_dir / "verification.json").read_text())
        assert verification["status"] == "PASS"
        p = subprocess.run(["sha256sum", "-c", "SHA256SUMS.txt"], cwd=result_dir, capture_output=True, text=True)
        assert p.returncode == 0, p.stderr
        # The project verifier rechecks exact positions, cell, species order,
        # source hashes, method, convergence and finite values.
        py = str(VERIFIER_PYTHON if VERIFIER_PYTHON.is_file() else Path(sys.executable))
        q = subprocess.run([py, str(TASK / "verify_result.py"), record["label"]], cwd=REPO, capture_output=True, text=True)
        assert q.returncode == 0, q.stderr
        verified = json.loads(q.stdout)
        assert verified["status"] == "PASS"
        atoms = helper.extxyz_geometry_frames(input_path)[0]
        expected_formula = {"C": 4, "Si": 2, "Ti": 6}
        counts = {symbol: atoms.get_chemical_symbols().count(symbol) for symbol in sorted(set(atoms.get_chemical_symbols()))}
        assert len(atoms) == 12 and counts == expected_formula
        summary = json.loads((result_dir / "summary.json").read_text())
        assert summary["native_extrapolated_energy_eV_cell"] is not None
        assert summary["free_energy_eV_cell"] is not None
        assert summary["forces_method"].endswith("preserve free_energy separately")
        assert summary["training_integration"].startswith("not performed")

        comparisons = []
        for split_name, frames in (("train36", train_frames), ("dev3", dev_frames)):
            comparable = [helper.compare_geometry(atoms, other) for other in frames]
            comparable = [value for value in comparable if value is not None]
            comparisons.append({
                "split": split_name,
                "frames_read_geometry_only": len(frames),
                "same_formula_cell_pbc_comparable_frames": len(comparable),
                "exact_hits": sum(bool(value["exact_duplicate"]) for value in comparable),
                "near_hits": sum(bool(value["near_duplicate"]) for value in comparable),
                "note": "The registered train36/dev3 frames are 52 atoms and contain Ag; candidate is 12-atom Ag-free Ti3SiC2. This audit screen does not claim supercell-equivalence checking.",
            })
        input_checks.append({
            "label": record["label"], "input_sha256": input_digest,
            "source_cif_sha256": record["source_cif_sha256"],
            "role": record["role"], "status": record["status"],
            "formula": atoms.get_chemical_formula(), "atom_count": len(atoms),
            "method": {"xc": record["method"]["xc"], "cutoff_eV": record["method"]["cutoff_eV"], "kpts": record["method"]["kpts"], "smearing_eV": record["method"]["smearing_eV"]},
            "result_verifier": "PASS", "per_result_SHA256": "PASS",
            "native_and_free_energy_separate": True,
            "force_energy_convention": summary["forces_method"],
            "training_integration": "not performed",
            "official_split_geometry_comparisons": comparisons,
        })
        candidates.append(atoms)

    pairwise = []
    for i, left in enumerate(candidates):
        for j, right in enumerate(candidates[i + 1:], start=i + 1):
            value = helper.compare_geometry(left, right)
            pairwise.append({"left": records[i]["label"], "right": records[j]["label"], "comparison": value})

    # The independent family claim is paper/source lineage, not a CIF-byte
    # claim: the candidate is a reconstruction of a 1967 single-crystal paper,
    # while the frozen V19 dev substrate is the 1998 COD/Kisi source.
    assert manifest["source"]["doi"] == "10.1007/BF00899949"
    assert dev_manifest["source_CIF_sha256"] != manifest["source"]["source_cif_sha256"]
    result = {
        "status": "PASS_PROVENANCE_AND_LABEL_INTEGRITY; external_role_audit_pending",
        "candidate_task": manifest["task_id"],
        "candidate_source": {
            "source_family": "Jeitschko_Nowotny_1967_single_crystal",
            "doi": "10.1007/BF00899949",
            "source_cif_sha256": manifest["source"]["source_cif_sha256"],
            "uploaded_paper_sha256": manifest["source"]["uploaded_paper_sha256"],
            "geometry_assumption": "Full occupancy; the original paper discusses possible nonstoichiometry.",
        },
        "frozen_v19_development_source": {
            "source_cif_sha256": dev_manifest["source_CIF_sha256"],
            "source_evidence_limit": dev_manifest["source_evidence_limit"],
            "independent_crystallographic_family": True,
        },
        "official_split": {
            "status": approved["status"], "train_rows": approved["roles"]["train"],
            "dev_rows": approved["roles"]["development_validation"], "test_rows": approved["roles"]["test"],
            "train_sha256": approved["train_sha256"], "dev_sha256": approved["development_sha256"],
            "test_or_sealed_data_opened": False,
            "candidate_training_data_added": False,
        },
        "record_audits": input_checks,
        "candidate_pairwise_geometry": pairwise,
        "interpretation": [
            "The 1967 paper/source family is distinct from the COD 9009647/Kisi 1998 family used by V19 dev3.",
            "All four result archives re-pass the registered geometry/method/convergence/finite-value verifier and individual checksums.",
            "The four labels are a single source-family displacement panel, not four independent source families.",
            "Native extrapolated energy and finite-width free energy are separately retained; any training target still requires an explicit energy-convention decision.",
            "The scope does not perform full supercell-equivalence screening against every historical geometry or grant train/dev assignment.",
        ],
        "decision": "Keep all four as train_acquisition_candidate_pending_audit; no official dataset or model changed.",
    }
    (REPORT / "candidate_role_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "records": len(input_checks), "pairwise_comparisons": len(pairwise), "train_rows": len(train_frames), "dev_rows": len(dev_frames)}, indent=2))


if __name__ == "__main__":
    main()
