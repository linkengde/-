#!/usr/bin/env python3
"""Build and verify the frozen B1 training-candidate-only geometry ledger."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

from ase.io import read


ROOT = Path(__file__).resolve().parents[6]
REPORTS = ROOT / "research/mace-v12-transfer/coordination/reports/window-c"
OUT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def collect(folder: str, manifest_name: str, group_kind: str, subdir: str = "") -> list[dict]:
    base = REPORTS / folder
    manifest = json.loads((base / manifest_name).read_text())
    rows = []
    for record in manifest["records"]:
        rel = Path(subdir) / record["file"] if subdir else Path(record["file"])
        path = base / rel
        digest = sha256(path)
        assert digest == record["file_sha256"], f"source manifest hash mismatch: {path}"
        atoms = read(path)
        assert len(atoms) == record.get("atoms", record.get("natoms")), f"atom count mismatch: {path}"
        assert atoms.get_chemical_formula() == record["formula"], f"formula mismatch: {path}"
        assert atoms.calc is None or not atoms.calc.results, f"unexpected calculator labels: {path}"
        labels = {"energy", "energy_ev", "free_energy", "free_energy_ev", "forces", "stress"}
        keys = {str(k).lower() for k in atoms.info} | {str(k).lower() for k in atoms.arrays}
        assert not (labels & keys), f"label-like fields present: {path}"
        assert atoms.info.get("dft_launch_enabled") is False, f"DFT launch not disabled: {path}"

        if group_kind == "solid_ag":
            lineage = "B1_SOLID_AG_COD9013045_SUH1988"
            source_scope = "direct_structure_source_same_experiment_as_COD"
            temperature = 293
            density_ratio = ""
            seed = record.get("seed", "")
        elif group_kind == "ti3sic2":
            lineage = "B1_Ti3SiC2_Bai2010_ALPHA"
            source_scope = "independent_computational_structure_parameters"
            temperature = ""
            density_ratio = ""
            seed = record.get("seed", "")
        else:
            match = re.search(r"_seed(\d+)_", record.get("config_type", path.stem))
            if match is None:
                match = re.search(r"_seed(\d+)_", path.stem)
            assert match is not None, f"missing liquid source seed: {path}"
            seed = int(match.group(1))
            assert seed in (501, 902), f"unexpected liquid seed {seed}: {path}"
            lineage = f"B1_LIQUID_AG_EAM_SEED_{seed}"
            source_scope = "density_reference_only_EAM_generated_coordinates"
            temperature = record["temperature_K"]
            density_ratio = record["density_ratio"]

        rows.append({
            "config_type": record.get("config_type", path.stem),
            "candidate_file": str(path.relative_to(ROOT)),
            "file_sha256": digest,
            "formula": record["formula"],
            "natoms": record.get("atoms", record.get("natoms")),
            "lineage_group": lineage,
            "source_scope": source_scope,
            "temperature_K": temperature,
            "density_ratio": density_ratio,
            "seed": seed,
            "candidate_role": "TRAINING_CANDIDATE_ONLY",
            "dev_test_role": "NONE",
            "dft_label_present": False,
            "dft_launch_enabled": False,
            "dataset_integrated": False,
        })
    return rows


def main() -> None:
    rows = []
    rows += collect("ag_solid_cod9013045_reference_candidates_20261024", "candidate_manifest.json", "solid_ag")
    rows += collect("independent_ti3sic2_bai2010_candidate_20261024", "supercell_geometry_manifest.json", "ti3sic2", "supercell_candidates")
    rows += collect("abdullaev2020_liquid_ag_density_rescaled_20261024", "rescaled_candidate_manifest.json", "liquid_ag")
    rows.sort(key=lambda row: (row["lineage_group"], row["config_type"]))

    counts: dict[str, int] = {}
    for row in rows:
        counts[row["lineage_group"]] = counts.get(row["lineage_group"], 0) + 1
    assert len(rows) == 24
    assert len({row["file_sha256"] for row in rows}) == 24
    assert counts == {
        "B1_SOLID_AG_COD9013045_SUH1988": 6,
        "B1_Ti3SiC2_Bai2010_ALPHA": 6,
        "B1_LIQUID_AG_EAM_SEED_501": 6,
        "B1_LIQUID_AG_EAM_SEED_902": 6,
    }, counts
    assert all(row["candidate_role"] == "TRAINING_CANDIDATE_ONLY" for row in rows)
    assert all(row["dev_test_role"] == "NONE" for row in rows)
    assert all(not row["dft_label_present"] and not row["dft_launch_enabled"] for row in rows)
    assert all(not row["dataset_integrated"] for row in rows)

    columns = list(rows[0])
    with (OUT / "role_ledger.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    audit = {
        "status": "PASS_FROZEN_TRAINING_CANDIDATE_POOL_ONLY",
        "candidate_count": len(rows),
        "unique_candidate_hash_count": len({row["file_sha256"] for row in rows}),
        "group_counts": counts,
        "all_label_free": all(not row["dft_label_present"] for row in rows),
        "all_dft_launch_disabled": all(not row["dft_launch_enabled"] for row in rows),
        "all_dev_test_roles_none": all(row["dev_test_role"] == "NONE" for row in rows),
        "all_unintegrated": all(not row["dataset_integrated"] for row in rows),
        "energy_target_for_future_b1_finite_smearing_arm": "free_energy_eV",
        "native_energy_retained_separately": True,
        "dft_authorized": False,
        "dataset_integration_authorized": False,
        "dft_method_frozen": False,
        "v20_owner_controls_complete": False,
        "v20_task_snapshot_utc": "2026-10-24T03:28:13.980616+00:00",
        "v20_live_telemetry_available": False,
    }
    (OUT / "ledger_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
