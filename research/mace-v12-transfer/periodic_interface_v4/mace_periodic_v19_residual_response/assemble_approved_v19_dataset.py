#!/usr/bin/env python3
"""Assemble the reviewed V19 train36/dev3 dataset; never reads test/sealed data."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
DATA = ROOT / "data"
INDEPENDENT = ROOT.parent / "pbe_interface_v19_independent_dev"
REPORT = REPO / "research/mace-v12-transfer/coordination/reports/window-b/v19_dataset_entry_independent_audit"
AUDIT = REPORT / "audit.json"
REPORT_FILES = ["README.md", "audit.json", "audit.py", "SHA256SUMS.txt"]
FOUNDATION = ROOT.parent / "mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model"
OWNER = "9c31b999-b2fe-4f0a-86f0-52c605d705ad"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(REPO.resolve()))


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def verify_report() -> dict:
    subprocess.run(["sha256sum", "-c", "SHA256SUMS.txt"], cwd=REPORT, check=True)
    completed = subprocess.run(
        [sys.executable, str(REPORT / "audit.py")], cwd=REPO, check=True,
        text=True, capture_output=True,
    )
    reproduced = json.loads(completed.stdout)
    assert reproduced["status"] == "PASS_WITH_SHARED_SOURCE_LIMITATION"
    report = load_json(AUDIT)
    assert report["status"] == reproduced["status"]
    assert all(row["status"] == "PASS" for row in report["checks"])
    assert report["decision_boundary"]["dataset_approved"] is False
    assert report["decision_boundary"]["test_or_sealed_data_opened"] is False
    return report


def registered_frames(report: dict, input_manifest: dict) -> tuple[list, list[dict], dict[str, str]]:
    records = {row["label"]: row for row in input_manifest["records"]}
    ordered = ["AgC_v19_independent_dev_01", "AgSi_v19_independent_dev_01", "AgTi_v19_independent_dev_01"]
    frames, rows, pins = [], [], {}
    source_group = "Ti3SiC2_COD_9009647_shared_crystallographic_source_family"
    for index, label in enumerate(ordered):
        entry = records[label]
        inp = REPO / entry["input"]
        calc = REPO / entry["result_directory"] / "calculations" / label
        output = calc / f"{label}_PW_PBE.extxyz"
        summary_path = calc / "summary.json"
        verification_path = calc / "verification.json"
        sums_path = calc / "SHA256SUMS.txt"
        summary = load_json(summary_path)
        verification = load_json(verification_path)
        archive_evidence = report["candidate_entries"][label]

        assert entry["role"] == "independent_development"
        assert entry["launch_enabled"] is True
        assert digest(inp) == entry["input_sha256"]
        assert verification["status"] == "PASS"
        assert archive_evidence["archive"]["scf_converged"] is True
        assert archive_evidence["fixed_geometry"]["positions_max_abs_delta_A"] == 0.0
        assert archive_evidence["fixed_geometry"]["cell_max_abs_delta_A"] == 0.0
        assert archive_evidence["fixed_geometry"]["species_and_order_equal"] is True
        assert archive_evidence["fixed_geometry"]["lammps_ids_and_order_equal"] is True

        atoms = read(output)
        energy = float(summary["energy_eV_cell"])
        free_energy = float(summary["free_energy_eV_cell"])
        forces = np.asarray(atoms.arrays["PW_PBE_forces"], dtype=float)
        assert np.isfinite(energy) and np.isfinite(free_energy)
        assert forces.shape == (len(atoms), 3) and np.isfinite(forces).all()
        assert np.isclose(float(atoms.info["PW_PBE_energy_eV"]), energy, rtol=0.0, atol=1e-12)
        atoms.arrays["REF_forces"] = forces.copy()
        atoms.info["REF_energy"] = energy
        atoms.info["free_energy_eV"] = free_energy
        atoms.info["dataset_role"] = "independent_development"
        atoms.info["source_group"] = source_group
        atoms.info["source_method"] = "GPAW 26.7.0 PBE/PW 500 eV, 6x6x1, Fermi-Dirac sigma 0.10 eV, xy dipole"
        atoms.info["energy_convention"] = "native extrapolated energy; free energy separately retained"
        atoms.info["source_input_sha256"] = entry["input_sha256"]
        atoms.info["source_label_sha256"] = digest(output)
        atoms.info["B_audit_status"] = "PASS_WITH_SHARED_SOURCE_LIMITATION"
        atoms.info["A_lineage_review"] = "APPROVED_NOT_REUSED_V17_DEVELOPMENT_FAMILY"
        frames.append(atoms)

        geometry = archive_evidence["geometry_sha256"]
        rows.append({
            "frame_index": index,
            "config_type": label,
            "label": label,
            "interface": entry["interface"],
            "role": "independent_development",
            "A_lineage_review": "APPROVED_NOT_REUSED_V17_DEVELOPMENT_FAMILY",
            "source_group": source_group,
            "source_CIF": input_manifest["source_CIF"],
            "source_CIF_sha256": input_manifest["source_CIF_sha256"],
            "input_path": entry["input"],
            "input_sha256": digest(inp),
            "label_path": rel(output),
            "label_sha256": digest(output),
            "geometry_sha256": geometry,
            "method": entry["method"],
            "energy_eV_cell_extrapolated": energy,
            "free_energy_eV_cell": free_energy,
            "scf_iterations": archive_evidence["archive"]["scf_iterations"],
            "archive_verification": "PASS",
            "geometry_screen": "PASS_NO_EXACT_OR_NEAR_HITS_IN_AUDITED_SCOPE",
            "source_family_independent": False,
        })
        for path in (inp, output, summary_path, verification_path, sums_path):
            pins[rel(path)] = digest(path)

    return frames, rows, pins


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--approve-and-assemble", action="store_true", help="record A's approved train36/dev3/test0 dataset")
    args = parser.parse_args()
    if not args.approve_and_assemble:
        parser.error("explicit --approve-and-assemble is required")
    assert not DATA.exists(), f"Refusing to overwrite existing dataset directory: {DATA}"

    audit = verify_report()
    input_manifest_path = INDEPENDENT / "input_manifest.json"
    input_manifest = load_json(input_manifest_path)
    assert input_manifest["status"] == "READY_LAUNCH_ENABLED"
    assert input_manifest["source_evidence_limit"].find("remote bytes/checksum remain unverified") >= 0
    train_source = ROOT / "training_candidate/train.extxyz"
    candidate_manifest = load_json(ROOT / "training_candidate/manifest.json")
    assert candidate_manifest["roles"]["train_candidate"] == 36
    assert candidate_manifest["withheld_labels_opened"] is False
    assert digest(train_source) == candidate_manifest["train_sha256"]
    subprocess.run([sys.executable, str(ROOT / "build_training_candidate.py"), "--check"], cwd=REPO, check=True)

    frames, dev_rows, pins = registered_frames(audit, input_manifest)
    pins[rel(train_source)] = digest(train_source)
    pins[rel(ROOT / "training_candidate/manifest.json")] = digest(ROOT / "training_candidate/manifest.json")
    pins[rel(input_manifest_path)] = digest(input_manifest_path)
    geometry_audit = INDEPENDENT / "geometry_audit.json"
    pins[rel(geometry_audit)] = digest(geometry_audit)
    cif_path = REPO / input_manifest["source_CIF"]
    pins[rel(cif_path)] = digest(cif_path)
    for name in REPORT_FILES:
        path = REPORT / name
        pins[rel(path)] = digest(path)

    DATA.mkdir(parents=True)
    shutil.copyfile(train_source, DATA / "train.extxyz")
    write(DATA / "valid.extxyz", frames, format="extxyz")
    train_sha = digest(DATA / "train.extxyz")
    dev_sha = digest(DATA / "valid.extxyz")
    assert train_sha == digest(train_source)

    (DATA / "README.md").write_text(
        "# A-approved V19 dataset\n\n"
        "This package contains the frozen train36/dev3 split approved by window A after B's read-only audit. "
        "`train.extxyz` preserves the existing 36 training rows byte-for-byte; `valid.extxyz` contains the three "
        "fixed-geometry AgC, AgSi, and AgTi development labels with native extrapolated energies and forces. "
        "There is no test file, and no sealed/test labels were opened.\n\n"
        "## Independence boundary\n\n"
        "The B audit passed its geometry, hash, method, and archive checks, with a shared-source limitation: all "
        "three development structures inherit the same uploaded Ti3SiC2 COD 9009647 CIF. They are distinct from "
        "the screened train36, six V19 residual geometries, and inspected V17/V18 development geometries, but "
        "they are one crystallographic source family, not three independent source families. The upstream COD "
        "bytes remain externally unverified. This split supports development screening, not source-diverse or "
        "blind-test generalization claims. Slab k-point, smearing, and vacuum convergence remain unresolved.\n\n"
        "See `dataset_manifest.json`, `../entry_review_approval.json`, and `SHA256SUMS.txt` for roles, lineage, "
        "approval evidence, and reproducibility hashes.\n"
    )

    source_files = dict(sorted(pins.items()))
    manifest = {
        "schema_version": 1,
        "status": "A_APPROVED_V19_DATASET",
        "approved_by": "window-a",
        "approved_owner_instance": OWNER,
        "roles": {"train": 36, "development_validation": 3, "test": 0},
        "train_path": "train.extxyz",
        "development_path": "valid.extxyz",
        "train_sha256": train_sha,
        "development_sha256": dev_sha,
        "source_files_sha256": source_files,
        "development_rows": dev_rows,
        "B_audit": {
            "status": audit["status"],
            "commit": "a13259dc663d793ab78f61763284ca557b3496a9",
            "report_dir": rel(REPORT),
            "audit_json_sha256": digest(AUDIT),
            "all_checks_pass": True,
            "shared_source_group_count": 1,
        },
        "source_independence_boundary": "Three geometrically separated development rows share the COD 9009647 Ti3SiC2 crystallographic source family; not three independent source families.",
        "withheld_labels_opened": False,
        "training_authorized": True,
        "energy_convention": "GPAW native extrapolated total energy (force_consistent=False); free energy separately recorded",
        "foundation_model_sha256": digest(FOUNDATION),
        "limitations": [
            "All three development rows share one crystallographic source family; source-diverse generalization is not established.",
            "The uploaded COD 9009647 CIF bytes have internal lineage checks but were not externally verified against COD.",
            "Slab k-point, smearing, and vacuum convergence remain unresolved.",
            "Development labels do not constitute blind test evidence or establish model accuracy.",
        ],
    }
    (DATA / "dataset_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    review_files = {
        rel(REPORT / name): digest(REPORT / name) for name in REPORT_FILES
    }
    review_files.update({
        rel(input_manifest_path): digest(input_manifest_path),
        rel(geometry_audit): digest(geometry_audit),
        rel(train_source): digest(train_source),
        rel(ROOT / "training_candidate/manifest.json"): digest(ROOT / "training_candidate/manifest.json"),
    })
    approval = {
        "status": "A_APPROVED_AFTER_B_AUDIT",
        "approved_utc": datetime.now(timezone.utc).isoformat(),
        "approved_by": "window-a",
        "owner_instance": OWNER,
        "dataset_status": "A_APPROVED_V19_DATASET",
        "dataset_manifest_sha256": digest(DATA / "dataset_manifest.json"),
        "roles": {"train": 36, "development_validation": 3, "test": 0},
        "decision": "Approve the audited train36/dev3 split for the existing fixed V19 screening fit, with the shared-source and slab-convergence limitations retained.",
        "B_audit_status": audit["status"],
        "B_audit_commit": "a13259dc663d793ab78f61763284ca557b3496a9",
        "shared_source_limitation_acknowledged": True,
        "review_files_sha256": dict(sorted(review_files.items())),
        "withheld_labels_opened": False,
        "preflight_required_before_training": True,
        "training_authorized_only_if_preflight_passes": True,
    }
    (ROOT / "entry_review_approval.json").write_text(json.dumps(approval, indent=2, sort_keys=True) + "\n")
    files = ["README.md", "train.extxyz", "valid.extxyz", "dataset_manifest.json"]
    (DATA / "SHA256SUMS.txt").write_text("".join(f"{digest(DATA / name)}  {name}\n" for name in files))

    assert len(read(DATA / "train.extxyz", ":")) == 36
    valid = read(DATA / "valid.extxyz", ":")
    assert len(valid) == 3
    assert [a.info["config_type"] for a in valid] == [r["config_type"] for r in dev_rows]
    for atoms in valid:
        assert atoms.info["dataset_role"] == "independent_development"
        assert np.isfinite(atoms.info["REF_energy"])
        assert np.isfinite(atoms.arrays["REF_forces"]).all()
    print(json.dumps({
        "status": manifest["status"],
        "roles": manifest["roles"],
        "train_sha256": train_sha,
        "development_sha256": dev_sha,
        "B_audit": audit["status"],
        "source_group_count": 1,
        "withheld_labels_opened": False,
    }, indent=2))


if __name__ == "__main__":
    main()
