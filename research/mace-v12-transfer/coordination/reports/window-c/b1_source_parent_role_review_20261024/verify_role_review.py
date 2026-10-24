#!/usr/bin/env python3
"""Read-only integrity and parent-role audit for the B1 geometry proposals."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from ase.io import read


ROOT = Path(__file__).resolve().parents[6]
PROJECT = ROOT / "research/mace-v12-transfer"
REPORTS = PROJECT / "coordination/reports/window-c"
OUT = Path(__file__).resolve().parent / "role_review.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_package(folder: str, manifest_name: str, *, subdir: str = "") -> list[dict]:
    base = REPORTS / folder
    manifest = json.loads((base / manifest_name).read_text())
    checked = []
    for record in manifest["records"]:
        rel = Path(subdir) / record["file"] if subdir else Path(record["file"])
        path = base / rel
        actual_hash = sha256(path)
        expected_hash = record["file_sha256"]
        assert actual_hash == expected_hash, f"file hash mismatch: {path}"
        atoms = read(path)
        expected_n = record.get("atoms", record.get("natoms"))
        expected_formula = record.get("formula")
        assert len(atoms) == expected_n, f"atom count mismatch: {path}"
        assert atoms.get_chemical_formula() == expected_formula, f"formula mismatch: {path}"
        assert atoms.calc is None or not atoms.calc.results, f"unexpected calculator labels: {path}"
        info_keys = {str(k).lower() for k in atoms.info}
        array_keys = {str(k).lower() for k in atoms.arrays}
        label_keys = {"energy", "energy_ev", "free_energy", "free_energy_ev", "forces", "stress"}
        assert not (label_keys & (info_keys | array_keys)), f"label-like extxyz key found: {path}"
        assert record.get("DFT_launch_enabled", record.get("dft_launch_enabled")) is False
        label_free = record.get("no_DFT_label") is True or record.get("dft_label_present") is False
        assert label_free, f"manifest does not mark candidate label-free: {path}"
        checked.append({
            "file": str(path.relative_to(ROOT)),
            "file_sha256": actual_hash,
            "formula": expected_formula,
            "atoms": expected_n,
            "info_keys": sorted(info_keys),
            "arrays": sorted(array_keys),
        })
    assert len(checked) == manifest.get("candidate_count", len(manifest["records"]))
    return checked


def main() -> None:
    solid = verify_package(
        "ag_solid_cod9013045_reference_candidates_20261024", "candidate_manifest.json"
    )
    bai = verify_package(
        "independent_ti3sic2_bai2010_candidate_20261024",
        "supercell_geometry_manifest.json",
        subdir="supercell_candidates",
    )
    liquid = verify_package(
        "abdullaev2020_liquid_ag_density_rescaled_20261024", "rescaled_candidate_manifest.json"
    )

    liquid_manifest = json.loads(
        (REPORTS / "abdullaev2020_liquid_ag_density_rescaled_20261024/rescaled_candidate_manifest.json").read_text()
    )
    source_base = REPORTS / "continuation_geometry_prep_20261024"
    for rec in liquid_manifest["records"]:
        source = source_base / rec["source_file"]
        assert sha256(source) == rec["source_sha256"], f"liquid parent hash mismatch: {source}"

    legacy_csv = REPORTS / "continuation_geometry_prep_20261024/B1_ROLE_FREEZE_DRAFT.csv"
    with legacy_csv.open(newline="") as f:
        legacy_rows = list(csv.DictReader(f))
    assert len(legacy_rows) == 30
    assert all(r["DFT_label_present"] == "false" for r in legacy_rows)
    assert all(r["DFT_launch_enabled"] == "false" for r in legacy_rows)
    assert all(r["dev_or_test_role"] == "unassigned" for r in legacy_rows)
    legacy_by_batch = dict(Counter(r["material_batch"] for r in legacy_rows))
    legacy_liquid_temperatures = dict(Counter(r["temperature_K"] for r in legacy_rows if r["material_batch"] == "B1-LIQ"))

    gate = json.loads((REPORTS / "b1_gate_reassessment_20261024/gate_status.json").read_text())
    task_a = json.loads((PROJECT / "coordination/tasks/window-a.json").read_text())
    v20 = {
        "task_snapshot_utc": task_a.get("updated_utc"),
        "owner_instance": task_a.get("owner_instance"),
        "sigma_0p05": task_a.get("current_jobs", {}).get("AgSi_slab_k6x6x1_sigma0p05_v20"),
        "sigma_0p20": task_a.get("current_jobs", {}).get("AgSi_slab_k6x6x1_sigma0p20_v20"),
        "live_telemetry_available": False,
    }

    result = {
        "report_date": "2026-10-24",
        "status": "INTEGRITY_PASS_SOURCE_PARENT_ROLE_REVIEW_COMPLETE_ROLES_NOT_FROZEN",
        "verification": {
            "reviewed_geometry_candidate_count": len(solid) + len(bai) + len(liquid),
            "geometry_source_verified_candidate_count": len(solid) + len(bai),
            "density_reference_only_candidate_count": len(liquid),
            "conditionally_held_ti3sic2_geometry_count": 0,
            "solid_ag_cod9013045_count": len(solid),
            "ti3sic2_bai2010_count": len(bai),
            "liquid_ag_abdullaev2020_count": len(liquid),
            "all_candidate_file_hashes_match_manifests": True,
            "all_formula_and_atom_counts_match": True,
            "all_candidate_extxyz_are_label_free": True,
            "all_candidate_files_disable_dft_launch": True,
            "all_liquid_source_parent_hashes_match": True,
            "legacy_role_draft_count": len(legacy_rows),
            "legacy_role_draft_by_batch": legacy_by_batch,
            "legacy_liquid_temperature_counts": legacy_liquid_temperatures,
        },
        "role_recommendation": {
            "dev_assignments": 0,
            "test_assignments": 0,
            "training_roles_frozen": False,
            "dft_authorized": False,
            "energy_target_frozen": gate.get("energy_target_frozen"),
            "dft_method_frozen": gate.get("dft_method_frozen"),
            "groups": [
                {"group": "COD9013045_fcc_Ag_293K", "count": len(solid), "recommendation": "keep all COD-derived supercells and displacement variants together; geometry candidates only"},
                {"group": "Bai2010_Ti3SiC2_alpha_phase", "count": len(bai), "recommendation": "one source family across both supercells; source parameters verified; geometry candidates only"},
                {"group": "Ag_liquid_EAM_seed501", "count": 6, "recommendation": "keep all temperatures and density-ratio copies together; geometry lineage is EAM-derived"},
                {"group": "Ag_liquid_EAM_seed902", "count": 6, "recommendation": "keep all temperatures and density-ratio copies together; geometry lineage is EAM-derived"},
            ],
            "split_protection": [
                {"group": "Ag_fcc_COD9013045", "members": "all six COD supercells and displacement variants", "rule": "one lineage; never split across train/dev/test"},
                {"group": "Ti3SiC2_alpha_phase", "members": "Bai2010 set and any retained COD9009647 set", "rule": "phase/source-correlated structures; no independent dev/test claim"},
                {"group": "Ag_liquid_EAM_seed501", "members": "all six reference-scaled candidates from seed 501", "rule": "keep temperature and density-ratio variants together"},
                {"group": "Ag_liquid_EAM_seed902", "members": "all six reference-scaled candidates from seed 902", "rule": "keep temperature and density-ratio variants together"},
            ],
        },
        "legacy_draft_disposition": {
            "superseded_solid_ag_eam_candidates": 6,
            "same_cod_family_ti3sic2_candidates_not_independent": 6,
            "liquid_ag_1050K_excluded_from_liquid_reference_pool": 6,
            "liquid_ag_1250K_1600K_model_anchor_versions_replaced_by_reference_scaled_copies": 12,
        },
        "v20_owner_snapshot": v20,
        "candidate_files": {
            "solid_ag": solid,
            "ti3sic2_bai2010": bai,
            "liquid_ag": liquid,
        },
    }
    assert result["verification"]["reviewed_geometry_candidate_count"] == 24
    assert result["verification"]["geometry_source_verified_candidate_count"] == 12
    assert result["verification"]["density_reference_only_candidate_count"] == 12
    assert result["verification"]["conditionally_held_ti3sic2_geometry_count"] == 0
    assert result["role_recommendation"]["energy_target_frozen"] is False
    assert result["role_recommendation"]["dft_method_frozen"] is False
    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "candidate_files"}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
