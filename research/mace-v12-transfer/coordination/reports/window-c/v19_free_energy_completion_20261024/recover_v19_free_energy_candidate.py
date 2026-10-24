#!/usr/bin/env python3
"""Build a separate 36-frame V19 free-energy sidecar from verified archives."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.io import read, write


ROOT = Path(__file__).resolve().parents[6]
HERE = Path(__file__).resolve().parent
V19_DIR = ROOT / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v19_residual_response"
V19_TRAIN = V19_DIR / "training_candidate/train.extxyz"
V19_MANIFEST = V19_DIR / "training_candidate/manifest.json"
V18_DIR = ROOT / "research/mace-v12-transfer/coordination/reports/window-c/continuation_readiness_20261024"
V18_SIDECAR = V18_DIR / "v18_train30_free_energy_candidate.extxyz"
V18_AUDIT = V18_DIR / "v18_free_energy_recovery_audit.json"
OUT = HERE / "v19_train36_free_energy_candidate.extxyz"
MANIFEST_OUT = HERE / "recovery_audit.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare_geometry(a, b, *, tolerance: float = 1.0e-12) -> dict:
    assert a.get_chemical_symbols() == b.get_chemical_symbols()
    assert np.array_equal(a.pbc, b.pbc)
    pos_error = float(np.max(np.abs(a.positions - b.positions)))
    cell_error = float(np.max(np.abs(np.asarray(a.cell) - np.asarray(b.cell))))
    assert pos_error <= tolerance
    assert cell_error <= tolerance
    return {"max_position_error_A": pos_error, "max_cell_error_A": cell_error}


def find_source_file(source_hashes: dict[str, str], config: str, suffix: str) -> Path:
    matches = [
        ROOT / path for path in source_hashes
        if config in path and path.endswith(suffix)
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one {suffix} source for {config}, got {len(matches)}")
    return matches[0]


def logged_value(log_text: str, label: str) -> float:
    matches = re.findall(rf"^\s*{re.escape(label)}:\s*([-+0-9.eE]+)", log_text, re.MULTILINE)
    if not matches:
        raise ValueError(f"GPAW log has no {label!r} line")
    return float(matches[-1])


def main() -> None:
    training_manifest = json.loads(V19_MANIFEST.read_text())
    v18_audit = json.loads(V18_AUDIT.read_text())
    if training_manifest.get("training_authorized") is not False:
        raise ValueError("V19 candidate no longer records training as unauthorized")
    if v18_audit.get("status") != "PASS_ALL_30_FREE_ENERGIES_RECOVERED":
        raise ValueError("V18 free-energy source audit is not PASS_ALL_30")
    if sha256(V19_TRAIN) != training_manifest["train_sha256"]:
        raise ValueError("V19 candidate input hash differs from its manifest")

    v19 = read(V19_TRAIN, index=":")
    v18 = read(V18_SIDECAR, index=":")
    if len(v19) != 36 or len(v18) != 30 or len(v18_audit.get("records", [])) != 30:
        raise ValueError(f"Unexpected frame counts V19={len(v19)}, V18={len(v18)}")

    audit_rows: list[dict] = []
    for i in range(30):
        v19_frame, v18_frame = v19[i], v18[i]
        geometry = compare_geometry(v19_frame, v18_frame)
        assert v19_frame.info["config_type"] == v18_frame.info["config_type"]
        assert v19_frame.get_chemical_formula() == v18_frame.get_chemical_formula()
        assert v19_frame.info["REF_energy"] == v18_frame.info["REF_energy"]
        force_error = float(np.max(np.abs(v19_frame.arrays["REF_forces"] - v18_frame.arrays["REF_forces"])))
        assert force_error <= 1.0e-8
        free_energy = float(v18_frame.info["free_energy_eV"])
        if "free_energy_eV" in v19_frame.info:
            assert float(v19_frame.info["free_energy_eV"]) == free_energy
        else:
            v19_frame.info["free_energy_eV"] = free_energy
        v19_frame.info["free_energy_provenance"] = "verified V18 exact-archive sidecar; V19 base30 bytes/geometry/forces matched"
        v19_frame.info["free_energy_source_audit_sha256"] = sha256(V18_AUDIT)
        audit_rows.append({
            "frame_index": i,
            "config_type": v19_frame.info["config_type"],
            "source": "V18 exact-archive recovery sidecar",
            "native_energy_eV": float(v19_frame.info["REF_energy"]),
            "free_energy_eV": free_energy,
            "native_minus_free_energy_eV": float(v19_frame.info["REF_energy"]) - free_energy,
            "geometry_comparison": geometry,
            "force_max_abs_error_eV_A": force_error,
            "native_energy_exact_match": True,
        })

    source_hashes = training_manifest["source_files_sha256"]
    for i in range(30, 36):
        frame = v19[i]
        config = str(frame.info["config_type"])
        manifest_row = next(row for row in training_manifest["new_rows"] if row["config_type"] == config)
        result_path = find_source_file(source_hashes, config, "_PW_PBE.extxyz")
        summary_path = find_source_file(source_hashes, config, "summary.json")
        verification_path = find_source_file(source_hashes, config, "verification.json")
        summary = json.loads(summary_path.read_text())
        verification = json.loads(verification_path.read_text())
        result = read(result_path)
        log_path = result_path.parent / "gpaw.log"

        if sha256(result_path) != manifest_row["label_sha256"]:
            raise ValueError(f"Label hash mismatch for {config}")
        if sha256(summary_path) != source_hashes[str(summary_path.relative_to(ROOT))]:
            raise ValueError(f"Summary hash mismatch for {config}")
        if sha256(verification_path) != source_hashes[str(verification_path.relative_to(ROOT))]:
            raise ValueError(f"Verification hash mismatch for {config}")
        if verification.get("status") != "PASS" or verification.get("owner_task") not in {"window-a", "window-b"}:
            raise ValueError(f"Archive verification is not PASS for {config}")
        log_hash = sha256(log_path)
        if log_hash != verification["sha256"]["gpaw.log"]:
            raise ValueError(f"GPAW log hash mismatch for {config}")
        if summary["source_sha256"] != manifest_row["input_sha256"]:
            raise ValueError(f"Source input hash mismatch for {config}")
        if not summary.get("scf_converged") or summary.get("label") != config:
            raise ValueError(f"SCF or label metadata mismatch for {config}")

        geometry = compare_geometry(frame, result)
        assert frame.get_chemical_formula() == result.get_chemical_formula()
        assert np.array_equal(frame.arrays["REF_forces"], result.arrays["PW_PBE_forces"])
        native = float(summary["energy_eV_cell"])
        free = float(summary["free_energy_eV_cell"])
        assert abs(float(frame.info["REF_energy"]) - native) <= 1.0e-12
        if "PW_PBE_energy_eV" in frame.info:
            assert abs(float(frame.info["PW_PBE_energy_eV"]) - native) <= 1.0e-12
        assert abs(float(frame.info["free_energy_eV"]) - free) <= 1.0e-12

        log_text = log_path.read_text(errors="replace")
        logged_native = logged_value(log_text, "Extrapolated")
        logged_free = logged_value(log_text, "Free energy")
        assert abs(logged_native - native) <= 5.1e-7
        assert abs(logged_free - free) <= 5.1e-7
        frame.info["free_energy_provenance"] = "V19 owner-verified archive summary and GPAW log"
        frame.info["free_energy_source_summary_sha256"] = sha256(summary_path)
        audit_rows.append({
            "frame_index": i,
            "config_type": config,
            "source": "V19 owner-verified summary/log/output archive",
            "result_sha256": sha256(result_path),
            "summary_sha256": sha256(summary_path),
            "verification_sha256": sha256(verification_path),
            "gpaw_log_sha256": log_hash,
            "native_energy_eV": native,
            "free_energy_eV": free,
            "logged_native_rounding_error_eV": logged_native - native,
            "logged_free_rounding_error_eV": logged_free - free,
            "geometry_comparison": geometry,
            "force_max_abs_error_eV_A": 0.0,
            "archive_verification": "PASS",
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    write(OUT, v19, format="extxyz")
    recovered = read(OUT, index=":")
    if len(recovered) != 36 or any("free_energy_eV" not in a.info for a in recovered):
        raise ValueError("Serialized sidecar is missing free-energy values")
    max_roundtrip_force_error = 0.0
    max_roundtrip_position_error = 0.0
    for source, saved in zip(v19, recovered):
        assert source.get_chemical_symbols() == saved.get_chemical_symbols()
        assert np.array_equal(source.pbc, saved.pbc)
        position_error = float(np.max(np.abs(source.positions - saved.positions)))
        force_error = float(np.max(np.abs(source.arrays["REF_forces"] - saved.arrays["REF_forces"])))
        assert position_error <= 1.0e-8 and force_error <= 1.0e-8
        assert np.array_equal(np.asarray(source.cell), np.asarray(saved.cell))
        assert float(source.info["REF_energy"]) == float(saved.info["REF_energy"])
        assert abs(float(source.info["free_energy_eV"]) - float(saved.info["free_energy_eV"])) <= 1.0e-8
        max_roundtrip_position_error = max(max_roundtrip_position_error, position_error)
        max_roundtrip_force_error = max(max_roundtrip_force_error, force_error)
    native_free_differences = [float(a.info["REF_energy"]) - float(a.info["free_energy_eV"]) for a in recovered]
    if any(a.info.get("training_authorized") is True for a in recovered):
        raise ValueError("A candidate sidecar must not authorize training")

    manifest = {
        "status": "PASS_ALL_36_FREE_ENERGIES_RECOVERED_CANDIDATE_ONLY",
        "frame_count": len(recovered),
        "input_training_candidate_sha256": sha256(V19_TRAIN),
        "input_manifest_sha256": sha256(V19_MANIFEST),
        "v18_sidecar_sha256": sha256(V18_SIDECAR),
        "v18_audit_sha256": sha256(V18_AUDIT),
        "output_sha256": sha256(OUT),
        "native_energy_key_preserved": "REF_energy",
        "free_energy_key": "free_energy_eV",
        "all_frames_have_both_energy_keys": True,
        "v18_inherited_frames_reconciled": 30,
        "v19_new_frames_reconciled_to_owner_archives": 6,
        "training_authorized": False,
        "dataset_modified": False,
        "development_or_test_data_read": False,
        "dft_or_training_run": False,
        "native_minus_free_energy_eV_cell_range": [min(native_free_differences), max(native_free_differences)],
        "serialized_max_position_error_A": max_roundtrip_position_error,
        "serialized_max_force_error_eV_A": max_roundtrip_force_error,
        "interpretation": "Complete alternative force-consistent energy-key sidecar is available for all V19 train-candidate frames. It does not choose the project energy convention, authorize V19 training, replace the existing train.extxyz, or provide free energies for validation frames.",
        "records": audit_rows,
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "status": manifest["status"],
        "frames": len(recovered),
        "V18_inherited": 30,
        "V19_new_archives": 6,
        "all_have_free_energy": all("free_energy_eV" in a.info for a in recovered),
        "training_authorized": False,
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
