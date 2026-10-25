#!/usr/bin/env python3
"""Read-only audit of the smearing width in approved V19 train36/dev3 evidence.

Reads only named, approved lineage manifests, archive summaries, and existing
train36/dev3 audit records. It does not open dataset extxyz files, test/sealed
data, or run GPAW, MACE, or any other calculation.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
V19 = REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v19_residual_response"
V16 = REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy"
RECON = REPO / "research/mace-v12-transfer/coordination/reports/window-a/v19_inherited_train36_archive_reconciliation_20261023/reconciliation.json"
LINEAGE = REPO / "research/mace-v12-transfer/coordination/reports/window-a/v19_train36_energy_method_lineage_audit_20261023/audit.json"
V16_MANIFEST = V16 / "data/dataset_manifest.json"
V19_CANDIDATE = V19 / "training_candidate/manifest.json"
V19_DATASET = V19 / "data/dataset_manifest.json"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def row(label: str, group: str, method: dict, evidence: Path | str) -> dict:
    assert method.get("xc") == "PBE", (label, method)
    assert method.get("basis") in ("PW", "plane wave"), (label, method)
    assert method.get("cutoff_eV") == 500, (label, method)
    assert method.get("kpts") == [1, 1, 1], (label, method)
    assert method.get("smearing_eV") == 0.1, (label, method)
    return {
        "label": label,
        "group": group,
        "method": method,
        "evidence": str(evidence),
    }


def main() -> None:
    reconciliation = read_json(RECON)
    lineage = read_json(LINEAGE)
    v16 = read_json(V16_MANIFEST)
    candidate = read_json(V19_CANDIDATE)
    dataset = read_json(V19_DATASET)
    rows: list[dict] = []

    # Twelve formerly inherited rows, now reconciled to their archived summaries.
    assert reconciliation["status"] == "PASS" and reconciliation["row_count"] == 12
    for item in reconciliation["rows"]:
        assert item["status"] == "PASS"
        rows.append(row(item["label"], "inherited_reconciled_12", item["source_method"], RECON))

    # Three V15 rows carry exact method metadata in the V16 per-frame ledger.
    v15 = [item for item in v16["per_frame_provenance"]
           if item.get("config_type", "").startswith((
               "AgC_residual_shell_v15_01", "AgSi_residual_shell_v15_01", "AgTi_residual_shell_v15_01"))]
    assert len(v15) == 3
    for item in v15:
        rows.append(row(item["config_type"], "v15_archived_3", item["method"], V16_MANIFEST))

    # Six V16 acquisitions are checked against their exact archived summaries.
    additions = v16["V16_additions"]
    assert len(additions) == 6
    for item in additions:
        summary_path = REPO / item["source"] / "summary.json"
        summary = read_json(summary_path)
        rows.append(row(item["label"], "v16_archived_6", summary["method"], summary_path))

    # Eight V17 acquisitions plus the separately caveated rank-recovery row.
    v17 = [item for item in lineage["frames"] if item["source_category"].startswith("V17_")]
    assert len(v17) == 9
    for item in v17:
        output_path = REPO / item["lineage_evidence"]["source_path"]
        summary_path = output_path.parent / "summary.json"
        summary = read_json(summary_path)
        rows.append(row(item["config_type"], "v17_archived_9", summary["method"], summary_path))

    # Six latest V19 acquisitions use the method recorded in their source manifest.
    assert len(candidate["new_rows"]) == 6
    for item in candidate["new_rows"]:
        rows.append(row(item["config_type"], "v19_archived_6", item["source_method"], V19_CANDIDATE))

    assert len(rows) == 36
    assert len({item["label"] for item in rows}) == 36
    train_frames = lineage["frames"]
    assert len(train_frames) == 36
    assert all(item["pbc"] == [True, True, True] for item in train_frames)

    # The approved development manifest lists three sigma=0.10, 6x6x1 slab labels.
    dev = dataset["development_rows"]
    assert len(dev) == 3
    assert all(item["method"]["smearing_eV"] == 0.1 for item in dev)
    assert all(item["method"]["kpts"] == [6, 6, 1] for item in dev)
    assert all(item["method"]["xc"] == "PBE" and item["method"]["cutoff_eV"] == 500 for item in dev)
    assert lineage["dev3_summary"]["pbc_counts"] == {"(True, True, False)": 3}

    counts = {}
    for item in rows:
        counts[item["group"]] = counts.get(item["group"], 0) + 1
    checks = {
        "train36_row_coverage": len(rows) == 36,
        "train36_all_sigma_0p10": all(item["method"]["smearing_eV"] == 0.1 for item in rows),
        "train36_all_gamma": all(item["method"]["kpts"] == [1, 1, 1] for item in rows),
        "dev3_all_sigma_0p10": all(item["method"]["smearing_eV"] == 0.1 for item in dev),
        "dev3_all_6x6x1": all(item["method"]["kpts"] == [6, 6, 1] for item in dev),
        "train_vs_dev_mesh_differs": all(item["method"]["kpts"] == [1, 1, 1] for item in rows)
        and all(item["method"]["kpts"] == [6, 6, 1] for item in dev),
        "train_vs_dev_pbc_differs": all(item["pbc"] == [True, True, True] for item in train_frames)
        and lineage["dev3_summary"]["pbc_counts"] == {"(True, True, False)": 3},
    }
    assert all(checks.values()), checks

    inputs = [RECON, LINEAGE, V16_MANIFEST, V19_CANDIDATE, V19_DATASET]
    inputs.extend(REPO / item["source"] / "summary.json" for item in additions)
    inputs.extend((REPO / item["lineage_evidence"]["source_path"]).parent / "summary.json" for item in v17)
    inputs = sorted(set(inputs))
    result = {
        "status": "PASS",
        "scope": "Read-only method audit of approved V19 train36/dev3 metadata and exact source summaries. No dataset extxyz, test/sealed data, GPAW, MACE, or calculation was opened/run.",
        "counts": {"train36": len(rows), "dev3": len(dev), **counts},
        "train_method": {"xc": "PBE", "cutoff_eV": 500, "kpts": [1, 1, 1], "smearing_eV": 0.1, "pbc": [True, True, True]},
        "development_method": {"xc": "PBE", "cutoff_eV": 500, "kpts": [6, 6, 1], "smearing_eV": 0.1, "pbc": [True, True, False]},
        "checks": checks,
        "training_rows": rows,
        "development_rows": [{"label": item["label"], "method": item["method"]} for item in dev],
        "source_sha256": {path.relative_to(REPO).as_posix(): sha256(path) for path in inputs},
        "interpretation": "The approved V19 train36 and dev3 both use 0.10 eV smearing. The training labels use Gamma and full periodicity, while dev3 uses 6x6x1 and a nonperiodic-z slab. The V20 smearing comparison supports keeping 0.10 eV as the current width reference but does not resolve the remaining mesh/PBC mismatch.",
        "dataset_changed": False,
    }
    out = HERE / "approved_v19_smearing_consistency.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "counts": result["counts"], "checks": checks, "output": out.name}, indent=2))


if __name__ == "__main__":
    main()
