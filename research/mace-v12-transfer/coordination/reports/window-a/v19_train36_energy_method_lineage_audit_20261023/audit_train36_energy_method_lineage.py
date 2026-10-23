#!/usr/bin/env python3
"""Read-only lineage audit of approved V19 train36 and dev3 metadata only."""
from __future__ import annotations

import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from ase.io import read

REPO = Path(__file__).resolve().parents[6]
PKG = REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v19_residual_response"
TRAIN_PATH = PKG / "data/train.extxyz"
DEV_PATH = PKG / "data/valid.extxyz"
OUTPUT = Path(__file__).resolve().parent / "audit.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    train = read(TRAIN_PATH, index=":")
    dev = read(DEV_PATH, index=":")
    candidate = json.loads((PKG / "training_candidate/manifest.json").read_text())
    v16 = json.loads((REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy/data/dataset_manifest.json").read_text())
    v18 = json.loads((REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/data/dataset_manifest.json").read_text())
    v19 = json.loads((PKG / "data/dataset_manifest.json").read_text())
    stdout = (PKG / "runs/v19_train36_dev3_20261023_A/train.stdout").read_text()

    v16_provenance = {row.get("config_type"): row for row in v16.get("per_frame_provenance", []) if row.get("config_type")}
    v16_additions = {row["label"]: row for row in v16.get("V16_additions", [])}
    base_rows = {
        row.get("config", row.get("config_type")): row
        for row in candidate.get("base_rows", [])
    }
    new_rows = {
        row.get("config_type"): row
        for row in candidate.get("new_rows", [])
    }

    categories: dict[str, list[dict[str, object]]] = defaultdict(list)
    frame_rows = []
    for index, atoms in enumerate(train):
        label = str(atoms.info.get("config_type", ""))
        base_source = base_rows.get(label, {})
        v16_addition_label = next((key for key in v16_additions if label.startswith(key)), None)
        if label in new_rows:
            source = new_rows[label]
            category = "V19_archive_verified_training_acquisition"
            lineage = {
                "evidence_tier": source.get("lineage_role"),
                "energy_convention": source.get("energy_convention"),
                "source_method": source.get("source_method"),
                "input_sha256": source.get("input_sha256"),
                "label_sha256": source.get("label_sha256"),
            }
        elif v16_addition_label is not None:
            source = v16_additions[v16_addition_label]
            category = "V16_archive_verified_acquisition"
            lineage = {
                "evidence_tier": source.get("provenance_status"),
                "energy_convention": source.get("energy_convention"),
                "source_path": source.get("source"),
                "source_sha256": source.get("source_sha256"),
            }
        elif base_source.get("source") == "V16" and label in v16_provenance:
            evidence = v16_provenance[label]
            tier = evidence.get("evidence_tier")
            if tier == "inherited_declaration_not_reverified":
                category = "V16_inherited_source_not_independently_verified"
            elif tier == "verified_residual_archive":
                category = "V15_residual_archive_verified"
            else:
                category = f"V16_other_{tier or 'unclassified'}"
            method_evidence = evidence.get("energy_convention_evidence") or {}
            original_evidence = evidence.get("original_run_evidence") or {}
            lineage = {
                "evidence_tier": tier,
                "selected_energy_field": method_evidence.get("selected_field"),
                "original_getter_verified": method_evidence.get("independent_original_getter_verified"),
                "source_method_verified": original_evidence.get("proof_level") == "archive_verified",
                "source_frame": evidence.get("frame"),
            }
        else:
            source = base_source
            if "Ti_gap_2p3_rank_recovery" in label:
                category = "V17_rank_recovery_kpoint_convergence_unresolved"
                lineage = {
                    "evidence_tier": source.get("evidence"),
                    "source_path": source.get("source"),
                    "caveat": source.get("caveat"),
                }
            else:
                category = "V17_archive_verified_targeted_acquisition"
                lineage = {"evidence_tier": source.get("evidence"), "source_path": source.get("source"), "source_hash": source.get("source_hash")}

        method = str(atoms.info.get("source_method", "<missing>"))
        frame = {
            "frame": index,
            "config_type": label,
            "formula": atoms.get_chemical_formula(),
            "n_atoms": len(atoms),
            "source_category": category,
            "source_method_as_serialized": method,
            "pbc": [bool(value) for value in atoms.pbc],
            "has_REF_energy": "REF_energy" in atoms.info,
            "has_PW_PBE_energy_eV_alias": "PW_PBE_energy_eV" in atoms.info,
            "has_free_energy_eV": "free_energy_eV" in atoms.info,
            "has_energy_info_key": "energy" in atoms.info,
            "lineage_evidence": lineage,
        }
        categories[category].append(frame)
        frame_rows.append(frame)

    elements = ["Ag", "Ti", "Si", "C"]
    composition = np.array(
        [[frame.get_chemical_symbols().count(element) for element in elements] for frame in train],
        dtype=float,
    )
    singular = np.linalg.svd(composition, compute_uv=False)
    e0 = {
        int(z): {"foundation_e0_eV": float(base), "estimated_correction_eV": float(corr), "estimated_e0_eV": float(new)}
        for z, base, corr, new in re.findall(
            r"Element (\d+): foundation E0 =\s*([-+0-9.]+) eV, correction =\s*([-+0-9.]+) eV, new E0 =\s*([-+0-9.]+) eV",
            stdout,
        )
    }

    result = {
        "status": "TRAIN36_READ_ONLY_LINEAGE_AUDIT",
        "scope": "Approved train36 and dev3 metadata only. No test or sealed labels read.",
        "input_sha256": {"train36": sha256(TRAIN_PATH), "dev3": sha256(DEV_PATH)},
        "roles": v19["roles"],
        "train36_summary": {
            "n_frames": len(train),
            "source_category_counts": {key: len(value) for key, value in sorted(categories.items())},
            "source_method_serialized_counts": dict(Counter(str(frame.info.get("source_method", "<missing>")) for frame in train)),
            "pbc_counts": dict(Counter(str(tuple(bool(x) for x in frame.pbc)) for frame in train)),
            "energy_info_key_counts": {
                "REF_energy": sum("REF_energy" in frame.info for frame in train),
                "PW_PBE_energy_eV_alias": sum("PW_PBE_energy_eV" in frame.info for frame in train),
                "free_energy_eV": sum("free_energy_eV" in frame.info for frame in train),
                "energy": sum("energy" in frame.info for frame in train),
            },
            "formula_counts": dict(Counter(frame.get_chemical_formula() for frame in train)),
            "composition_matrix_elements": elements,
            "composition_matrix_rank": int(np.linalg.matrix_rank(composition)),
            "composition_matrix_singular_values": [float(value) for value in singular],
            "e0_estimates_from_training_log": e0,
            "v16_declared_physical_method_consistency": v16.get("physical_method_consistency"),
            "v18_energy_convention": v18.get("energy_convention"),
            "v18_energy_convention_limitations": v18.get("limitations"),
        },
        "dev3_summary": {
            "n_frames": len(dev),
            "formula_counts": dict(Counter(frame.get_chemical_formula() for frame in dev)),
            "source_method_counts": dict(Counter(str(frame.info.get("source_method", "<missing>")) for frame in dev)),
            "pbc_counts": dict(Counter(str(tuple(bool(x) for x in frame.pbc)) for frame in dev)),
        },
        "comparison": {
            "recorded_training_method": "Gamma-point PW-PBE 500 eV in all train36 source_method strings; original source settings are not independently confirmed for the 12 inherited rows.",
            "development_method": "PW-PBE 500 eV, 6x6x1, sigma 0.10 eV, xy dipole; recorded consistently in the approved dev manifest.",
            "training_vs_development_pbc": "all train36 TTT versus all dev3 TTF",
            "unresolved_energy_reference_rows": [frame["config_type"] for frame in frame_rows if frame["source_category"] == "V16_inherited_source_not_independently_verified"],
            "rank_recovery_convergence_unresolved_rows": [frame["config_type"] for frame in frame_rows if frame["source_category"] == "V17_rank_recovery_kpoint_convergence_unresolved"],
        },
        "frames": frame_rows,
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["train36_summary"], indent=2))
    print("Wrote", OUTPUT.relative_to(REPO))


if __name__ == "__main__":
    main()
