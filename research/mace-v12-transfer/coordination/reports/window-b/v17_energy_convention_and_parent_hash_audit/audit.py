#!/usr/bin/env python3
"""Reproduce the scoped V14-parent/V17-diagnostic provenance audit.

This script reads only the three assigned V14 training rows and their source
archives, the eight V17 paired training-diagnostic inputs/archives, and the
named provenance reports/manifests. It does not inspect V17 development or
withheld-test outputs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.geometry import find_mic
from ase.io import read


EXPECTED = {
    "AgC_registry_strain_acq_v14_01": {
        "interface": "AgC", "frame": 29, "formula": "C13Ag2Si10Ti23",
        "declared_hash": "72d649542158aade39a3cf50dc73759ad25df8a7ff7b2b117f1aa0c6ffc513a3",
    },
    "AgSi_registry_strain_acq_v14_01": {
        "interface": "AgSi", "frame": 30, "formula": "C13Ag17Si8Ti13",
        "declared_hash": "4a5540c7fd824a01df8fcb6155e62c7ff045ec36d637ac5f405247c9b5dbe4e5",
    },
    "AgTi_registry_probe_v13_01": {
        "interface": "AgTi", "frame": 28, "formula": "C11Ag28SiTi9",
        "declared_hash": "18dc1be733a506b748386175293d2f094e720e3e1c9dc96e259623088c005b7c",
    },
}

EXPECTED_V17_LABELS = {
    "AgC_framework_neighbor_p_acq_v17_01",
    "AgC_framework_neighbor_m_acq_v17_01",
    "AgSi_transverse_Ag_p_acq_v17_01",
    "AgSi_transverse_Ag_m_acq_v17_01",
    "AgSi_framework_neighbor_p_acq_v17_01",
    "AgSi_framework_neighbor_m_acq_v17_01",
    "AgTi_framework_neighbor_p_acq_v17_01",
    "AgTi_framework_neighbor_m_acq_v17_01",
}

EXPECTED_METHOD = {
    "xc": "PBE",
    "basis": "plane wave",
    "cutoff_eV": 500,
    "kpts": [1, 1, 1],
    "smearing_eV": 0.1,
    "scf_thresholds": "energy/density/eigenstates 1e-5",
    "mixer": "Pulay beta=0.05, nmaxold=8, weight=100",
    "GPAW_version": "26.7.0",
}


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def geometry_sha_v17(atoms) -> str:
    """Exact V17 canonical geometry serialization used by the V17 generator."""
    payload = {
        "symbols": atoms.get_chemical_symbols(),
        "positions_A": np.round(np.asarray(atoms.positions, dtype=float), 10).tolist(),
        "cell_A": np.round(np.asarray(atoms.cell.array, dtype=float), 10).tolist(),
        "pbc": np.asarray(atoms.pbc, dtype=bool).tolist(),
        "lammps_id": np.asarray(
            atoms.arrays.get("lammps_id", np.array([], dtype=int)), dtype=int
        ).tolist(),
    }
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256_bytes(serialized)


def exact_geometry_identity(a, b) -> bool:
    return (
        a.get_chemical_symbols() == b.get_chemical_symbols()
        and np.array_equal(a.arrays.get("lammps_id"), b.arrays.get("lammps_id"))
        and np.array_equal(np.asarray(a.positions), np.asarray(b.positions))
        and np.array_equal(np.asarray(a.cell.array), np.asarray(b.cell.array))
        and np.array_equal(np.asarray(a.pbc), np.asarray(b.pbc))
    )


def assert_method(method: dict, context: str) -> None:
    for key, expected in EXPECTED_METHOD.items():
        if method.get(key) != expected:
            raise AssertionError(f"{context}: method field {key}={method.get(key)!r}, expected {expected!r}")


def log_energy_values(log_text: str) -> dict:
    free = re.search(r"^\s*Free energy:\s*(-?\d+(?:\.\d+)?)\s*$", log_text, re.M)
    extrap = re.search(r"^\s*Extrapolated:\s*(-?\d+(?:\.\d+)?)\s*$", log_text, re.M)
    width = re.search(
        r"occupations=Occupations\(name=['\"]fermi-dirac['\"],\s*width=([0-9.]+)",
        log_text,
    )
    converged = re.search(r"Converged in\s+(\d+)\s+steps", log_text)
    return {
        "free_energy_log_eV_rounded": float(free.group(1)) if free else None,
        "extrapolated_log_eV_rounded": float(extrap.group(1)) if extrap else None,
        "fermi_dirac_width_eV": float(width.group(1)) if width else None,
        "converged_steps": int(converged.group(1)) if converged else None,
    }


def assert_rounding_matches(exact: float, rounded: float | None, context: str) -> None:
    if rounded is None or abs(float(exact) - rounded) > 0.00000051:
        raise AssertionError(f"{context}: log value {rounded!r} does not match {exact!r} to six decimals")


def output_equal_input_geometry(output, input_atoms, context: str) -> bool:
    same = (
        output.get_chemical_symbols() == input_atoms.get_chemical_symbols()
        and np.array_equal(output.arrays.get("lammps_id"), input_atoms.arrays.get("lammps_id"))
        and np.array_equal(output.positions, input_atoms.positions)
        and np.array_equal(output.cell.array, input_atoms.cell.array)
        and np.array_equal(output.pbc, input_atoms.pbc)
    )
    if not same:
        raise AssertionError(f"{context}: archived output geometry differs from its input")
    return same


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[6])
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    project = repo / "research/mace-v12-transfer"
    outdir = Path(__file__).resolve().parent

    v14_root = project / "periodic_interface_v4/pbe_interface_v14_parallel_acquisition"
    dataset_path = project / "periodic_interface_v4/mace_periodic_v14_interface_energy/data/train.extxyz"
    dataset_manifest_path = dataset_path.with_name("dataset_manifest.json")
    v14_input_manifest_path = v14_root / "input_manifest.json"
    parent_hash_manifest_path = project / "coordination/reports/window-b/v15_registry_candidates/candidate_screening_manifest.json"
    v17_prior_path = project / "coordination/reports/window-b/v17_paired_response_analysis/analysis_final/analysis.json"
    v16_audit_path = project / "coordination/reports/window-b/v16_dataset_integration_audit/audit.json"
    v14_runner_path = v14_root / "run_pw_reference.py"
    v14_builder_path = project / "periodic_interface_v4/mace_periodic_v14_interface_energy/build_v14_dataset.py"

    dataset_sha = sha256_file(dataset_path)
    expected_dataset_sha = "3323cdbd3081fc8d91673662c5f83c8da93ba3b4d7a23428f6b70d0e7a4c9faf"
    if dataset_sha != expected_dataset_sha:
        raise AssertionError(f"V14 dataset SHA changed: {dataset_sha}")
    dataset_manifest = read_json(dataset_manifest_path)
    v14_input_manifest = read_json(v14_input_manifest_path)
    parent_hash_manifest = read_json(parent_hash_manifest_path)
    v17_prior = read_json(v17_prior_path)
    v16_audit = read_json(v16_audit_path)

    source_rows = read(dataset_path, format="extxyz", index=":")
    labels_by_name = {x["label"]: x for x in dataset_manifest["added_training_labels"]}
    v14_inputs = {x["label"]: x for x in v14_input_manifest["records"]}
    parent_configs = parent_hash_manifest["source_dataset"]["parent_configs"]
    if v17_prior.get("all_eight_archives_verified") is not True:
        raise AssertionError("Prior V17 paired audit did not record all eight archives as verified")
    if set(v17_prior["archive_records"]) != EXPECTED_V17_LABELS:
        raise AssertionError("V17 paired-analysis scope differs from the eight assigned diagnostic labels")

    runner_source = v14_runner_path.read_text(encoding="utf-8")
    builder_source = v14_builder_path.read_text(encoding="utf-8")
    if "energy = float(atoms.get_potential_energy())" not in runner_source:
        raise AssertionError("V14 runner no longer shows default atoms.get_potential_energy() extraction")
    if "forces = np.asarray(atoms.get_forces(), dtype=float)" not in runner_source:
        raise AssertionError("V14 runner no longer shows direct atoms.get_forces() extraction")
    if "energy=float(atoms.info['PW_PBE_energy_eV'])" not in builder_source:
        raise AssertionError("V14 dataset builder does not copy PW_PBE_energy_eV to REF_energy as expected")
    if "forces=np.asarray(atoms.arrays['PW_PBE_forces'],dtype=float)" not in builder_source:
        raise AssertionError("V14 dataset builder does not copy PW_PBE_forces to REF_forces as expected")

    parent_rows = []
    parent_by_config = {}
    for label, spec in EXPECTED.items():
        interface = spec["interface"]
        frame_index = spec["frame"]
        if frame_index >= len(source_rows):
            raise AssertionError(f"Frame {frame_index} is absent from the V14 data")
        row = source_rows[frame_index]
        record = v14_inputs[label]
        manifest_row = labels_by_name[label]
        source_parent = parent_configs[interface]
        if row.info.get("config_type") != source_parent["source_config_type"]:
            raise AssertionError(f"{label}: config_type does not identify the expected V14 row")
        if row.get_chemical_formula() != spec["formula"] or len(row) != manifest_row["atoms"]:
            raise AssertionError(f"{label}: formula or atom count differs")
        if source_parent["source_dataset_sha256"] != dataset_sha:
            raise AssertionError(f"{label}: candidate manifest points to another V14 dataset SHA")
        if source_parent["source_frame_geometry_sha256"] != spec["declared_hash"]:
            raise AssertionError(f"{label}: declared legacy hash differs from the recorded source value")
        parent_canonical_sha = geometry_sha_v17(row)

        archive_dir = v14_root / "calculations" / label
        summary_path = archive_dir / "summary.json"
        verification_path = archive_dir / "verification.json"
        log_path = archive_dir / "gpaw.log"
        output_path = archive_dir / f"{label}_PW_PBE.extxyz"
        input_path = v14_root / record["input"]
        verification = read_json(verification_path)
        summary = read_json(summary_path)
        archived = read(output_path, format="extxyz")
        input_atoms = read(input_path, format="extxyz")
        log_text = log_path.read_text(encoding="utf-8", errors="replace")
        log_values = log_energy_values(log_text)
        actual_hashes = {
            "input": sha256_file(input_path),
            "output": sha256_file(output_path),
            "summary": sha256_file(summary_path),
            "log": sha256_file(log_path),
        }
        expected_hashes = {
            "input": record["input_sha256"],
            "output": manifest_row["output_sha256"],
            "summary": manifest_row["summary_sha256"],
            "verification_input": verification["input_sha256"],
            "verification_output": verification["sha256"][output_path.name],
            "verification_summary": verification["sha256"][summary_path.name],
            "verification_log": verification["sha256"][log_path.name],
        }
        if actual_hashes["input"] != expected_hashes["input"]:
            raise AssertionError(f"{label}: input SHA mismatch")
        if actual_hashes["output"] != expected_hashes["output"] or actual_hashes["output"] != expected_hashes["verification_output"]:
            raise AssertionError(f"{label}: output SHA mismatch")
        if actual_hashes["summary"] != expected_hashes["summary"] or actual_hashes["summary"] != expected_hashes["verification_summary"]:
            raise AssertionError(f"{label}: summary SHA mismatch")
        if actual_hashes["log"] != expected_hashes["verification_log"]:
            raise AssertionError(f"{label}: GPAW log SHA mismatch")
        if verification["status"] != "PASS" or not all(verification["checks"].values()):
            raise AssertionError(f"{label}: source archive verification failed")
        assert_method(summary["method"], label)
        if summary["mpi_ranks"] != 4 or summary["scf_converged"] is not True:
            raise AssertionError(f"{label}: archive rank/convergence evidence failed")
        if log_values["fermi_dirac_width_eV"] != 0.1 or log_values["converged_steps"] != summary["scf_iterations"]:
            raise AssertionError(f"{label}: GPAW log smearing or convergence evidence differs")
        assert_rounding_matches(summary["energy_eV_cell"], log_values["extrapolated_log_eV_rounded"], label + " native energy")
        if abs(float(summary["energy_eV_cell"]) - float(manifest_row["energy_eV_cell"])) > 0.0:
            raise AssertionError(f"{label}: dataset manifest energy differs from source summary")
        if not exact_geometry_identity(row, archived):
            raise AssertionError(f"{label}: V14 dataset row is not geometry/ID-identical to source archive")
        output_equal_input_geometry(archived, input_atoms, label)
        if row.info.get("REF_energy") != archived.info.get("PW_PBE_energy_eV"):
            raise AssertionError(f"{label}: REF_energy differs from original GPAW output energy")
        if float(row.info["REF_energy"]) != float(summary["energy_eV_cell"]):
            raise AssertionError(f"{label}: REF_energy differs from source summary")
        if not np.array_equal(row.arrays["REF_forces"], archived.arrays["PW_PBE_forces"]):
            raise AssertionError(f"{label}: REF_forces are not an exact array copy of source forces")
        if row.info.get("source_method") != "GPAW 26.7.0 PW-PBE 500 eV Gamma single point":
            raise AssertionError(f"{label}: source_method declaration differs")

        row_ids = np.asarray(row.arrays["lammps_id"], dtype=int)
        id_index = {int(atom_id): i for i, atom_id in enumerate(row_ids)}
        central_ids = source_parent["central_pair_lammps_ids"]
        if any(int(atom_id) not in id_index for atom_id in central_ids):
            raise AssertionError(f"{label}: candidate manifest control IDs are not in the parent row")
        declared_meta = source_parent["source_frame_geometry_sha256"]
        parent_record = {
            "label": label,
            "interface": interface,
            "frame_index_zero_based": frame_index,
            "config_type": row.info["config_type"],
            "formula": row.get_chemical_formula(),
            "atoms": len(row),
            "dataset_sha256": dataset_sha,
            "config_identity_matches": True,
            "control_atom_ids": [int(x) for x in central_ids],
            "declared_parent_geometry_sha256": declared_meta,
            "recomputed_v17_canonical_geometry_sha256": parent_canonical_sha,
            "legacy_hash_exact_match": parent_canonical_sha == declared_meta,
            "energy": {
                "REF_energy_eV_cell": float(row.info["REF_energy"]),
                "source_output_PW_PBE_energy_eV": float(archived.info["PW_PBE_energy_eV"]),
                "source_summary_energy_eV_cell": float(summary["energy_eV_cell"]),
                "dataset_manifest_energy_eV_cell": float(manifest_row["energy_eV_cell"]),
                "log_extrapolated_energy_eV_rounded": log_values["extrapolated_log_eV_rounded"],
                "log_free_energy_eV_rounded": log_values["free_energy_log_eV_rounded"],
                "ref_equals_native_extrapolated": True,
                "ref_equals_free_energy": False,
                "energy_difference_free_minus_native_eV_cell_rounded": float(
                    log_values["free_energy_log_eV_rounded"] - log_values["extrapolated_log_eV_rounded"]
                ),
                "free_energy_precision": "GPAW log prints six decimal places; original high-precision free-energy value is not stored in the V14 output/summary.",
            },
            "forces": {
                "ref_forces_exact_array_copy_of_source_PW_PBE_forces": True,
                "source_call": "atoms.get_forces() at Fermi-Dirac width 0.1 eV",
                "force_consistent_energy_flag_persisted": False,
            },
            "method": {
                **summary["method"],
                "smearing_direct_log_width_eV": log_values["fermi_dirac_width_eV"],
                "mpi_ranks": int(summary["mpi_ranks"]),
                "scf_iterations": int(summary["scf_iterations"]),
                "converged": bool(summary["scf_converged"]),
                "source_method_metadata": row.info["source_method"],
                "source_method_omits_smearing_and_energy_convention": True,
            },
            "archive": {
                "verification_status": verification["status"],
                "input_sha256": actual_hashes["input"],
                "output_sha256": actual_hashes["output"],
                "summary_sha256": actual_hashes["summary"],
                "log_sha256": actual_hashes["log"],
                "row_geometry_matches_source_output": True,
                "verification_checks_all_true": True,
            },
            "child_inputs": [],
        }
        parent_rows.append(parent_record)
        parent_by_config[source_parent["source_config_type"]] = (row, parent_record, interface)

    # Reconcile the six additions explicitly named in the V16 provenance review.
    v16_new = v16_audit["v16_additions"]
    v16_conventions = sorted({x.get("energy_convention") for x in v16_new})
    if v16_conventions != ["GPAW native extrapolated energy, force_consistent=False; free_energy recorded separately"]:
        raise AssertionError("V16 provenance review does not show the expected explicit convention")

    v17_manifest_paths = [
        project / "periodic_interface_v4/pbe_interface_v17_main_targeted_acquisition/input_manifest.json",
        project / "periodic_interface_v4/pbe_interface_v17_parallel_targeted_acquisition/input_manifest.json",
    ]
    v17_records = {}
    v17_label_rows = []
    for manifest_path in v17_manifest_paths:
        manifest = read_json(manifest_path)
        base = manifest_path.parent
        for record in manifest["records"]:
            label = record["label"]
            if label not in EXPECTED_V17_LABELS:
                continue
            if label in v17_records:
                raise AssertionError(f"V17 diagnostic label occurs in multiple manifests: {label}")
            v17_records[label] = (manifest_path, base, record)
    if set(v17_records) != EXPECTED_V17_LABELS:
        raise AssertionError("Could not find exactly the eight assigned V17 training diagnostics")

    for label in sorted(EXPECTED_V17_LABELS):
        manifest_path, base, record = v17_records[label]
        source_parent = record["source_parent"]
        if source_parent["dataset_sha256"] != dataset_sha:
            raise AssertionError(f"{label}: source parent dataset SHA mismatch")
        if source_parent["frame_index_zero_based"] >= len(source_rows):
            raise AssertionError(f"{label}: source parent frame is out of bounds")
        parent = source_rows[source_parent["frame_index_zero_based"]]
        parent_config = source_parent["config_type"]
        if parent_config != parent.info.get("config_type"):
            raise AssertionError(f"{label}: source parent config type does not match V14 row")
        _, parent_audit, interface = parent_by_config[parent_config]
        if source_parent["frame_geometry_sha256"] != parent_audit["declared_parent_geometry_sha256"]:
            raise AssertionError(f"{label}: child manifest did not preserve the upstream legacy hash")
        input_path = base / record["input"]
        input_atoms = read(input_path, format="extxyz")
        input_sha = sha256_file(input_path)
        if input_sha != record["input_sha256"] or input_sha != v17_prior["archive_records"][label]["input_sha256"]:
            raise AssertionError(f"{label}: diagnostic input hash mismatch")
        if len(input_atoms) != len(parent) or not np.array_equal(input_atoms.arrays["lammps_id"], parent.arrays["lammps_id"]):
            raise AssertionError(f"{label}: child atom count/order/IDs do not match parent")
        if input_atoms.get_chemical_symbols() != parent.get_chemical_symbols():
            raise AssertionError(f"{label}: child element order differs from parent")
        if not np.array_equal(input_atoms.cell.array, parent.cell.array) or not np.array_equal(input_atoms.pbc, parent.pbc):
            raise AssertionError(f"{label}: child cell or PBC differs from parent")
        parent_ids = np.asarray(parent.arrays["lammps_id"], dtype=int)
        id_to_idx = {int(atom_id): i for i, atom_id in enumerate(parent_ids)}
        declared_changes = record["design"]["changed_atoms"]
        expected_change_ids = {int(x["id"]) for x in declared_changes}
        if not expected_change_ids:
            raise AssertionError(f"{label}: no declared changed atoms")
        max_changed_delta_error = 0.0
        max_unlisted_displacement = 0.0
        actual_moved_ids = set()
        for atom_id, idx in id_to_idx.items():
            raw = np.asarray(input_atoms.positions[idx] - parent.positions[idx], dtype=float)
            mic = np.asarray(find_mic(raw, cell=parent.cell, pbc=parent.pbc)[0], dtype=float)
            magnitude = float(np.linalg.norm(mic))
            if magnitude > 1e-8:
                actual_moved_ids.add(atom_id)
            if atom_id in expected_change_ids:
                expected = next(np.asarray(x["delta_xyz_A"], dtype=float) for x in declared_changes if int(x["id"]) == atom_id)
                error = float(np.linalg.norm(mic - expected))
                max_changed_delta_error = max(max_changed_delta_error, error)
            else:
                max_unlisted_displacement = max(max_unlisted_displacement, magnitude)
        if actual_moved_ids != expected_change_ids:
            raise AssertionError(f"{label}: actual moved ID set differs from manifest: {actual_moved_ids ^ expected_change_ids}")
        if max_changed_delta_error > 5e-8 or max_unlisted_displacement > 1e-8:
            raise AssertionError(f"{label}: child displacement differs from declared design")

        archive_dir = base / "calculations" / label
        output_path = archive_dir / f"{label}_PW_PBE.extxyz"
        summary_path = archive_dir / "summary.json"
        verification_path = archive_dir / "verification.json"
        log_path = archive_dir / "gpaw.log"
        progress_path = archive_dir / "progress.json"
        summary = read_json(summary_path)
        verification = read_json(verification_path)
        output = read(output_path, format="extxyz")
        log_values = log_energy_values(log_path.read_text(encoding="utf-8", errors="replace"))
        actual_archive_hashes = {
            "verification.json": sha256_file(verification_path),
            output_path.name: sha256_file(output_path),
            "summary.json": sha256_file(summary_path),
            "gpaw.log": sha256_file(log_path),
            "progress.json": sha256_file(progress_path),
        }
        expected_archive_hashes = v17_prior["archive_records"][label]["archive_files_sha256"]
        if actual_archive_hashes != expected_archive_hashes:
            raise AssertionError(f"{label}: archive hashes differ from paired-analysis record")
        if verification.get("status") != "PASS" or not all(verification.get("checks", {}).values()):
            raise AssertionError(f"{label}: archive verification failed")
        if verification.get("sha256") != {k: v for k, v in actual_archive_hashes.items() if k != "verification.json"}:
            raise AssertionError(f"{label}: verification archive hash map differs")
        if verification.get("input_sha256") != input_sha:
            raise AssertionError(f"{label}: archive verification points to another input")
        assert_method(summary["method"], label)
        if summary["mpi_ranks"] != 4 or summary["scf_converged"] is not True:
            raise AssertionError(f"{label}: archive rank/convergence status failed")
        if summary.get("energy_convention") != "GPAW native extrapolated energy, force_consistent=False; free_energy recorded separately":
            raise AssertionError(f"{label}: explicit V17 energy convention field absent")
        if summary["source_sha256"] != input_sha:
            raise AssertionError(f"{label}: summary source hash differs")
        if float(summary["energy_eV_cell"]) != float(output.info["PW_PBE_energy_eV"]):
            raise AssertionError(f"{label}: summary native energy differs from output")
        prior = v17_prior["archive_records"][label]
        if float(summary["energy_eV_cell"]) != float(prior["native_energy_eV"]):
            raise AssertionError(f"{label}: summary native energy differs from paired-analysis record")
        if float(summary["free_energy_eV_cell"]) != float(prior["free_energy_eV"]):
            raise AssertionError(f"{label}: summary free energy differs from paired-analysis record")
        assert_rounding_matches(summary["energy_eV_cell"], log_values["extrapolated_log_eV_rounded"], label + " native energy")
        assert_rounding_matches(summary["free_energy_eV_cell"], log_values["free_energy_log_eV_rounded"], label + " free energy")
        if log_values["fermi_dirac_width_eV"] != 0.1:
            raise AssertionError(f"{label}: GPAW log does not evidence Fermi-Dirac 0.1 eV")
        output_equal_input_geometry(output, input_atoms, label)
        force_array = np.asarray(output.arrays["PW_PBE_forces"], dtype=float)
        if force_array.shape != (len(input_atoms), 3) or not np.isfinite(force_array).all():
            raise AssertionError(f"{label}: output force array is missing/invalid")
        child_record = {
            "label": label,
            "interface": interface,
            "owner": prior["owner"],
            "input_sha256": input_sha,
            "archive_verification": "PASS",
            "source_parent_frame_index_zero_based": int(source_parent["frame_index_zero_based"]),
            "source_parent_declared_geometry_sha256": source_parent["frame_geometry_sha256"],
            "source_parent_recomputed_v17_canonical_geometry_sha256": parent_audit["recomputed_v17_canonical_geometry_sha256"],
            "parent_row_identity": "PASS; dataset SHA, frame, config, formula, atom count, atom ID/order, symbols, cell and PBC verified",
            "declared_changed_atom_ids": sorted(expected_change_ids),
            "max_changed_atom_delta_error_A": max_changed_delta_error,
            "max_undeclared_atom_displacement_A": max_unlisted_displacement,
            "input_output_geometry_identity": "PASS",
            "native_energy_eV_cell": float(summary["energy_eV_cell"]),
            "free_energy_eV_cell": float(summary["free_energy_eV_cell"]),
            "energy_convention": summary["energy_convention"],
            "method": summary["method"],
            "mpi_ranks": int(summary["mpi_ranks"]),
            "scf_iterations": int(summary["scf_iterations"]),
            "log_extrapolated_eV_rounded": log_values["extrapolated_log_eV_rounded"],
            "log_free_eV_rounded": log_values["free_energy_log_eV_rounded"],
            "archive_sha256": actual_archive_hashes,
        }
        parent_audit["child_inputs"].append(child_record)
        v17_label_rows.append(child_record)

    # The previous matched-response report is the quantitative force/energy cross-check.
    paired_response_context = [
        {"pair": "AgC framework neighbor", "native_residual_eV_A": 0.0006, "free_residual_eV_A": 0.0018},
        {"pair": "AgSi transverse Ag", "native_residual_eV_A": -0.0716, "free_residual_eV_A": -0.0273},
        {"pair": "AgSi framework neighbor", "native_residual_eV_A": 0.0325, "free_residual_eV_A": 0.0014},
        {"pair": "AgTi framework neighbor", "native_residual_eV_A": 0.0583, "free_residual_eV_A": -0.0015},
    ]

    result = {
        "audit": "V17 energy-convention and V14 parent-geometry-hash trace",
        "scope": {
            "v14_parent_rows": [x["frame_index_zero_based"] for x in parent_rows],
            "v17_diagnostic_labels": sorted(EXPECTED_V17_LABELS),
            "v17_dev_validation_outputs_read": False,
            "v17_withheld_test_labels_or_outputs_read": False,
            "DFT_or_MACE_run": False,
            "production_metadata_modified": False,
        },
        "v14_dataset": {"path": str(dataset_path.relative_to(repo)), "sha256": dataset_sha, "frame_count": len(source_rows)},
        "canonical_hash_spec": {
            "source_implementation": "coordination/reports/window-b/v17_split_geometry_design/final_set/generate_split_candidates.py::geometry_sha",
            "payload_keys": ["symbols", "positions_A", "cell_A", "pbc", "lammps_id"],
            "symbols_and_ids": "ordered as stored in the extxyz row; atom identity is not sorted before serialization",
            "positions_and_cell": "ASE Cartesian positions/cell in angstrom, rounded with numpy.round(..., 10)",
            "pbc": "boolean list in stored order",
            "json": "json.dumps(sort_keys=True, separators=(',', ':')).encode('utf-8')",
            "hash": "SHA256 of serialized JSON bytes",
        },
        "legacy_hash_origin": {
            "first_tracked_artifact": "coordination/reports/window-b/v15_registry_candidates/candidate_screening_manifest.json",
            "introducing_commit": "a3b1299",
            "introduced_with_generator_or_hash_spec": False,
            "v14_builder_records_this_field": False,
            "v16_and_v17_copy_the_field_without_recomputing": True,
            "reproduction_status": "unresolved: the committed manifest gives the three digest strings but no payload schema, rounding, unit, ID ordering, serializer, or generating script",
            "canonical_v17_recomputation_is_a_replacement": False,
        },
        "parents": parent_rows,
        "v16_provenance_revision_evidence": {
            "source": "coordination/reports/window-b/v16_dataset_integration_audit/audit.json",
            "v16_additions": len(v16_new),
            "new_addition_energy_conventions": v16_conventions,
            "inherited_method_consistency": v16_audit.get("inherited_physical_method_consistency"),
            "scope_note": "The V16 report explicitly records native extrapolated energy/force_consistent=False and free energy separately for its six verified additions. Its inherited method-consistency result remains UNKNOWN; it does not retroactively define the legacy hash algorithm or certify all inherited frames.",
        },
        "paired_response_energy_force_context_from_existing_report": {
            "source": "coordination/reports/window-b/v17_paired_response_analysis/report.md, matched-response results",
            "interpretation": "The free-energy secant has a smaller force-consistency residual for three of four finite-displacement pairs; this motivates retaining both energy fields but is not proof that any DFT archive or force label is wrong.",
            "pairs": paired_response_context,
        },
        "recommendation": {
            "parent_relabeling_required_for_native_energy_target": False,
            "v17_energy_target_for_these_parents": "Use the existing REF_energy values as native extrapolated energy; new V17 REF_energy values should use summary.energy_eV_cell and preserve free_energy_eV_cell separately.",
            "forces": "Keep the archived force arrays unchanged. Record the GPAW get_forces() extraction and 0.1 eV Fermi smearing; do not imply REF_forces are force-consistent with native extrapolated energy because the legacy row lacks an explicit force_consistent field.",
            "parent_hashes": "Preserve each historical digest as legacy provenance; add the reproducible V17 canonical digest in a separate field and mark the legacy algorithm unknown. Row identity and signed child-displacement checks support the mapping, so the mismatch alone does not invalidate the DFT labels.",
            "whole_dataset_limit": "This resolves only the three assigned V14 parents. The V16 provenance audit still reports inherited physical-method consistency UNKNOWN; one cannot claim all inherited training rows are uniformly sourced from this method on this evidence alone.",
            "if_free_energy_is_chosen_instead": "Do not silently replace REF_energy. The V14 archived free-energy log values are rounded to six decimals, and a corpus-wide free-energy conversion would require a deliberate label-precision and provenance plan.",
        },
        "checks": {
            "v14_parent_archives_and_ref_copy": "PASS for all three rows",
            "v17_diagnostic_archives": "PASS for all eight labels",
            "parent_rows_and_child_displacements": "PASS for all eight inputs",
            "legacy_geometry_hash_reproduction": "UNRESOLVED because original algorithm/serialization evidence is absent",
            "expected_legacy_vs_v17_canonical_hash_equality": "FAIL for all three parents (documented metadata discrepancy)",
        },
    }
    (outdir / "audit.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    csv_path = outdir / "parent_evidence.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        fields = [
            "interface", "label", "frame_index_zero_based", "formula", "atoms", "dataset_sha256",
            "declared_legacy_geometry_sha256", "recomputed_v17_canonical_geometry_sha256", "hash_match",
            "REF_energy_eV_cell", "native_energy_eV_cell", "free_energy_log_eV_rounded",
            "REF_forces_exact_copy", "scf_iterations", "mpi_ranks", "smearing_eV",
            "child_diagnostic_count", "child_delta_max_error_A", "child_undeclared_move_max_A", "archive_status",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for parent in parent_rows:
            children = parent["child_inputs"]
            writer.writerow({
                "interface": parent["interface"],
                "label": parent["label"],
                "frame_index_zero_based": parent["frame_index_zero_based"],
                "formula": parent["formula"],
                "atoms": parent["atoms"],
                "dataset_sha256": parent["dataset_sha256"],
                "declared_legacy_geometry_sha256": parent["declared_parent_geometry_sha256"],
                "recomputed_v17_canonical_geometry_sha256": parent["recomputed_v17_canonical_geometry_sha256"],
                "hash_match": parent["legacy_hash_exact_match"],
                "REF_energy_eV_cell": parent["energy"]["REF_energy_eV_cell"],
                "native_energy_eV_cell": parent["energy"]["source_summary_energy_eV_cell"],
                "free_energy_log_eV_rounded": parent["energy"]["log_free_energy_eV_rounded"],
                "REF_forces_exact_copy": parent["forces"]["ref_forces_exact_array_copy_of_source_PW_PBE_forces"],
                "scf_iterations": parent["method"]["scf_iterations"],
                "mpi_ranks": parent["method"]["mpi_ranks"],
                "smearing_eV": parent["method"]["smearing_eV"],
                "child_diagnostic_count": len(children),
                "child_delta_max_error_A": max(x["max_changed_atom_delta_error_A"] for x in children),
                "child_undeclared_move_max_A": max(x["max_undeclared_atom_displacement_A"] for x in children),
                "archive_status": parent["archive"]["verification_status"],
            })
    print(f"Wrote {outdir / 'audit.json'}")
    print(f"Wrote {csv_path}")
    print("V14 archive lineage, V17 diagnostic archives, and parent/child geometry checks PASS; legacy hash algorithm unresolved as expected.")


if __name__ == "__main__":
    main()
