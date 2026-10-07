#!/usr/bin/env python3
"""Build V17's provisional 43-frame training pool and 3-frame dev validation set."""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

PROJECT = Path("research/mace-v12-transfer/periodic_interface_v4")
ENTRY = PROJECT / "mace_periodic_v17_interface_energy"
ELEMENTS = ("Ag", "C", "Si", "Ti")
ENERGY_CONVENTION = "GPAW native extrapolated energy, force_consistent=False; free_energy recorded separately"
EXPECTED_POSITIVE = {
    "AgC_framework_neighbor_p_acq_v17_01",
    "AgSi_transverse_Ag_p_acq_v17_01",
    "AgSi_framework_neighbor_p_acq_v17_01",
    "AgTi_framework_neighbor_p_acq_v17_01",
}
EXPECTED_NEGATIVE = {
    "AgC_framework_neighbor_m_acq_v17_01",
    "AgSi_transverse_Ag_m_acq_v17_01",
    "AgSi_framework_neighbor_m_acq_v17_01",
    "AgTi_framework_neighbor_m_acq_v17_01",
}
EXPECTED_DEV = {
    "AgC_v17_dev_validation_01",
    "AgSi_v17_dev_validation_01",
    "AgTi_v17_dev_validation_01",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load(path: Path):
    return json.loads(path.read_text())


def safe_path(repo: Path, relative: str) -> Path:
    path = (repo / relative).resolve()
    require(path.is_relative_to(repo.resolve()), f"Path escapes repository: {relative}")
    return path


def valid_labels(atoms, prefix: str = "REF"):
    energy_key = f"{prefix}_energy"
    force_key = f"{prefix}_forces"
    require(energy_key in atoms.info and force_key in atoms.arrays, f"Missing {energy_key}/{force_key}")
    energy = float(atoms.info[energy_key])
    forces = np.asarray(atoms.arrays[force_key], dtype=float)
    require(np.isfinite(energy), f"Nonfinite energy in {atoms.info.get('config_type', '?')}")
    require(forces.shape == (len(atoms), 3) and np.isfinite(forces).all(), "Invalid/nonfinite forces")
    require(np.isfinite(atoms.positions).all() and np.isfinite(atoms.cell.array).all(), "Nonfinite geometry")
    if "lammps_id" in atoms.arrays:
        require(len(set(map(int, atoms.arrays["lammps_id"]))) == len(atoms), "Duplicate atom IDs")
    return energy, forces


def exact_integer_rank(matrix: np.ndarray) -> int:
    """Exact rational row rank for the small integer composition matrix."""
    basis: list[tuple[int, list[Fraction]]] = []
    for source in matrix.tolist():
        row = [Fraction(int(x)) for x in source]
        for pivot, existing in basis:
            if row[pivot]:
                factor = row[pivot]
                row = [x - factor * y for x, y in zip(row, existing)]
        pivot = next((i for i, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        scale = row[pivot]
        row = [x / scale for x in row]
        basis.append((pivot, row))
        basis.sort(key=lambda item: item[0])
    return len(basis)


def frame_geometry_matches(source, result) -> None:
    require(np.array_equal(source.numbers, result.numbers), "Element/order mismatch")
    require(np.array_equal(source.pbc, result.pbc) and bool(np.all(source.pbc)), "PBC mismatch")
    require("lammps_id" in source.arrays and "lammps_id" in result.arrays, "Atom ID missing")
    require(np.array_equal(source.arrays["lammps_id"], result.arrays["lammps_id"]), "Atom ID/order mismatch")
    if "central_pair" in source.arrays or "central_pair" in result.arrays:
        require("central_pair" in source.arrays and "central_pair" in result.arrays, "Marked pair missing")
        require(np.array_equal(source.arrays["central_pair"], result.arrays["central_pair"]), "Marked pair mismatch")
    # GPAW extxyz results serialize positions at limited decimal precision; these
    # match the established archive verifier tolerances for that serialization.
    require(np.allclose(source.positions, result.positions, atol=1e-7, rtol=0), "Coordinates mismatch")
    require(np.allclose(source.cell.array, result.cell.array, atol=1e-8, rtol=0), "Cell mismatch")


def verify_training_archive(repo: Path, folder_rel: str, manifest_rel: str, expected_labels: set[str], owner: str, sources: dict[str, str]):
    folder = safe_path(repo, folder_rel)
    manifest_path = safe_path(repo, manifest_rel)
    manifest = load(manifest_path)
    sources[str(manifest_path.relative_to(repo))] = sha(manifest_path)
    require(manifest.get("owner_task") == owner, f"Wrong owner in {manifest_rel}")
    records = manifest.get("records", [])
    require({r.get("label") for r in records} == expected_labels and len(records) == len(expected_labels), f"Unexpected records in {manifest_rel}")
    frames = []
    for record in records:
        label = record["label"]
        require(record.get("owner_task") == owner and record.get("role") == "v17_targeted_diagnostic_training_acquisition", f"Wrong role/owner for {label}")
        input_path = safe_path(folder, record["input"])
        require(input_path.parent == folder / "inputs" and sha(input_path) == record["input_sha256"], f"Input hash mismatch: {label}")
        sources[str(input_path.relative_to(repo))] = sha(input_path)
        archive = folder / "calculations" / label
        verification = load(archive / "verification.json")
        summary = load(archive / "summary.json")
        progress = load(archive / "progress.json")
        sources[str((archive / "verification.json").relative_to(repo))] = sha(archive / "verification.json")
        require(verification.get("status") == "PASS" and verification.get("label") == label, f"Archive not verified: {label}")
        require(verification.get("checks") and all(v is True for v in verification["checks"].values()), f"Archive check failed: {label}")
        require(verification.get("input_sha256") == record["input_sha256"] == summary.get("source_sha256") == progress.get("source_sha256"), f"Archive/input linkage mismatch: {label}")
        require(summary.get("scf_converged") is True and progress.get("scf_converged") is True and progress.get("status") == "complete", f"SCF incomplete: {label}")
        require(summary.get("mpi_ranks") == progress.get("mpi_ranks") == 4, f"Expected four MPI ranks: {label}")
        method = summary.get("method", {})
        require(method.get("xc") == "PBE" and method.get("basis") == "plane wave" and method.get("cutoff_eV") == 500 and method.get("kpts") == [1, 1, 1] and method.get("smearing_eV") == 0.1, f"Unexpected DFT settings: {label}")
        require(summary.get("energy_convention") == ENERGY_CONVENTION, f"Unexpected energy convention: {label}")
        for name, digest in verification.get("sha256", {}).items():
            member = archive / name
            require(member.parent == archive and member.is_file() and sha(member) == digest, f"Archive hash mismatch: {label}/{name}")
            sources[str(member.relative_to(repo))] = digest
        result_path = archive / f"{label}_PW_PBE.extxyz"
        source = read(input_path)
        result = read(result_path)
        frame_geometry_matches(source, result)
        energy = float(result.info["PW_PBE_energy_eV"])
        forces = np.asarray(result.arrays["PW_PBE_forces"], dtype=float)
        require(np.isfinite(energy) and forces.shape == (len(result), 3) and np.isfinite(forces).all(), f"Invalid labels: {label}")
        require(abs(energy - float(summary["energy_eV_cell"])) < 1e-10, f"Summary energy mismatch: {label}")
        result.calc = None
        result.info["REF_energy"] = energy
        result.arrays["REF_forces"] = forces.copy()
        result.info.pop("PW_PBE_energy_eV", None)
        result.arrays.pop("PW_PBE_forces", None)
        result.info["config_type"] = label + "_PW_PBE_v17_train"
        result.info["v17_role"] = "paired_diagnostic_training"
        result.info["v17_provenance_status"] = "directly_archive_verified"
        result.info["v17_energy_convention"] = "native_extrapolated; free energy retained only in source archive"
        result.info["v17_reference_sha256"] = sha(result_path)
        frames.append((label, result))
    return frames


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--refresh-generated-data", action="store_true", help="replace only a prior hash-valid V17 provisional dataset produced by this builder")
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    entry = repo / ENTRY
    data = entry / "data"
    if data.exists() and any(data.iterdir()):
        require(args.refresh_generated_data, "Refusing to overwrite existing V17 dataset without --refresh-generated-data")
        old = load(data / "dataset_manifest.json")
        require(old.get("version") == "V17 provisional data-only numerical screening" and old.get("status", "").startswith("PROVISIONAL"), "Refresh refused: existing data are not a prior V17 provisional build")
        for line in (data / "SHA256SUMS.txt").read_text().splitlines():
            expected, name = line.split(maxsplit=1)
            require(sha(data / name.strip()) == expected, f"Refresh refused: previous generated file changed: {name}")
    sources: dict[str, str] = {}

    v16 = repo / PROJECT / "mace_periodic_v16_interface_energy"
    v16_data = v16 / "data"
    v16_manifest = load(v16_data / "dataset_manifest.json")
    v16_pins = load(v16 / "input_pins.json")
    sources[str((v16_data / "dataset_manifest.json").relative_to(repo))] = sha(v16_data / "dataset_manifest.json")
    sources[str((v16 / "input_pins.json").relative_to(repo))] = sha(v16 / "input_pins.json")
    provenance_path = repo / "research/mace-v12-transfer/coordination/reports/window-b/v17_training_pool_provenance_audit/ledger.json"
    provenance_ledger = load(provenance_path)
    sources[str(provenance_path.relative_to(repo))] = sha(provenance_path)
    train_rel = str((v16_data / "train.extxyz").relative_to(repo))
    expected_train_sha = v16_pins["sha256"][train_rel]
    require(sha(v16_data / "train.extxyz") == expected_train_sha == v16_manifest["output_sha256"]["train"], "V16 training split hash mismatch")
    sources[train_rel] = expected_train_sha
    train = read(v16_data / "train.extxyz", index=":")
    frame_map = {int(r["frame"]): r for r in provenance_ledger.get("rows", [])}
    require(len(train) == 35 and len(frame_map) == 35, "Expected 35 V16 training frames and complete B-audited provenance ledger")
    for i, atoms in enumerate(train):
        valid_labels(atoms)
        provenance = frame_map[i]
        require(provenance.get("config_type") == atoms.info.get("config_type"), f"V16 provenance/frame mismatch at index {i}")
        require(provenance.get("evidence_status") in ("directly_archive_verified", "partially_supported"), f"Missing audited evidence status at index {i}")
        atoms.calc = None
        atoms.info["v17_role"] = "inherited_v16_training"
        atoms.info["v17_source_frame"] = i
        atoms.info["v17_provenance_status"] = str(provenance["evidence_status"])
        atoms.info["v17_evidence_tier"] = str(provenance["evidence_status"])
        atoms.info["v17_source_energy_convention_status"] = "inherited; see per-frame provenance; not upgraded by V17"

    positive_rel = str(PROJECT / "pbe_interface_v17_main_targeted_acquisition")
    negative_rel = str(PROJECT / "pbe_interface_v17_parallel_targeted_acquisition")
    positive = verify_training_archive(repo, positive_rel, positive_rel + "/input_manifest.json", EXPECTED_POSITIVE, "window-a", sources)
    negative = verify_training_archive(repo, negative_rel, negative_rel + "/input_manifest.json", EXPECTED_NEGATIVE, "window-b", sources)
    additions = [atoms for _, atoms in positive + negative]
    train.extend(additions)
    require(len(train) == 43, f"Expected 43 training frames, found {len(train)}")

    # B's aggregate archive summary predates the fourth negative label. Each per-label
    # verification above must pass; record this bookkeeping discrepancy transparently.
    b_aggregate = load(repo / PROJECT / "pbe_interface_v17_parallel_targeted_acquisition" / "partial_archive_manifest.json")
    b_aggregate_stale = b_aggregate.get("complete") is not True or bool(b_aggregate.get("missing_labels"))
    sources[str((PROJECT / "pbe_interface_v17_parallel_targeted_acquisition" / "partial_archive_manifest.json"))] = sha(repo / PROJECT / "pbe_interface_v17_parallel_targeted_acquisition" / "partial_archive_manifest.json")

    role_root = repo / "research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set"
    role_path = role_root / "role_manifest.json"
    role_manifest = load(role_path)
    sources[str(role_path.relative_to(repo))] = sha(role_path)
    dev_candidates = [r for r in role_manifest["candidates"] if r.get("proposed_role") == "development_validation"]
    require({r["candidate_id"] for r in dev_candidates} == EXPECTED_DEV and len(dev_candidates) == 3, "Unexpected V17 dev geometry manifest")
    dev_root = repo / PROJECT / "pbe_interface_v17_main_targeted_acquisition" / "development_validation_archives"
    validation = []
    for candidate in dev_candidates:
        label = candidate["candidate_id"]
        geometry_path = role_root / candidate["file"]
        require(sha(geometry_path) == candidate["file_sha256"], f"Frozen dev input hash mismatch: {label}")
        sources[str(geometry_path.relative_to(repo))] = sha(geometry_path)
        execution_input = repo / "research/mace-v12-transfer/coordination/reports/window-b/v17_holdout_dft_execution_pack/inputs" / (label + ".extxyz")
        require(execution_input.is_file() and sha(execution_input) == candidate["file_sha256"], f"DFT input mismatch: {label}")
        sources[str(execution_input.relative_to(repo))] = sha(execution_input)
        archive = dev_root / label
        verification = load(archive / "verification.json")
        summary = load(archive / "summary.json")
        sources[str((archive / "verification.json").relative_to(repo))] = sha(archive / "verification.json")
        require(verification.get("status") == "PASS" and verification.get("role") == "development_validation", f"Dev archive not verified: {label}")
        require(verification.get("checks") and all(v is True for v in verification["checks"].values()), f"Dev archive check failed: {label}")
        require(verification.get("input_sha256") == summary.get("source_sha256") == candidate["file_sha256"], f"Dev input linkage mismatch: {label}")
        require(summary.get("scf_converged") is True and summary.get("finite") is True and summary.get("mpi_ranks") == 4, f"Dev DFT incomplete: {label}")
        require(str(summary.get("energy_convention", "")).startswith("GPAW native extrapolated energy, force_consistent=False"), f"Dev convention mismatch: {label}")
        dev_method = summary.get("method", {})
        require(dev_method.get("xc") == "PBE" and dev_method.get("basis") == "plane wave" and dev_method.get("cutoff_eV") == 500 and dev_method.get("kpts") == [1, 1, 1] and dev_method.get("smearing_eV") == 0.1, f"Dev DFT settings mismatch: {label}")
        if summary.get("free_energy_supported") is True:
            require(np.isfinite(float(summary["free_energy_eV_cell"])), f"Nonfinite free energy in dev archive: {label}")
        for name, digest in verification.get("compact_file_sha256", {}).items():
            member = archive / name
            require(member.parent == archive and member.is_file() and sha(member) == digest, f"Dev archive hash mismatch: {label}/{name}")
            sources[str(member.relative_to(repo))] = digest
        ref_path = archive / f"{label}_PW_PBE.extxyz"
        source = read(execution_input)
        result = read(ref_path)
        frame_geometry_matches(source, result)
        energy = float(result.info["PW_PBE_energy_eV"])
        forces = np.asarray(result.arrays["PW_PBE_forces"], dtype=float)
        require(np.isfinite(energy) and forces.shape == (len(result), 3) and np.isfinite(forces).all(), f"Dev labels invalid: {label}")
        require(abs(energy - float(summary["energy_eV_cell"])) < 1e-10, f"Dev summary energy mismatch: {label}")
        result.calc = None
        result.info["REF_energy"] = energy
        result.arrays["REF_forces"] = forces.copy()
        result.info.pop("PW_PBE_energy_eV", None)
        result.arrays.pop("PW_PBE_forces", None)
        result.info["v17_role"] = "development_validation"
        result.info["v17_energy_convention"] = "native_extrapolated; free energy retained only in source archive"
        result.info["v17_reference_sha256"] = sha(ref_path)
        validation.append((label, result))
    validation.sort(key=lambda item: item[0])

    matrix = np.asarray([[atoms.get_chemical_symbols().count(el) for el in ELEMENTS] for atoms in train], dtype=int)
    rank = exact_integer_rank(matrix)
    require(rank == 4, f"Composition rank is {rank}, expected 4")
    for atoms in train:
        valid_labels(atoms)
    for _, atoms in validation:
        valid_labels(atoms)
        require(len(np.flatnonzero(np.asarray(atoms.arrays.get("central_pair", []), dtype=bool))) == 2, "Dev marked pair must contain two atoms")

    for name in ("README.md", "build_v17_dataset.py", "preflight_v17.py", "run_v17_training.sh", "evaluate_v17_dev.py"):
        path = repo / ENTRY / name
        sources[str(path.relative_to(repo))] = sha(path)

    # A geometry-only split audit exists for these frozen roles. This builder never
    # opens any V17 withheld DFT output or label.
    data.mkdir(parents=True, exist_ok=True)
    write(data / "train.extxyz", train, format="extxyz")
    write(data / "valid.extxyz", [atoms for _, atoms in validation], format="extxyz")
    frame_status_counts = {}
    for atoms in train:
        key = str(atoms.info.get("v17_provenance_status", "unknown"))
        frame_status_counts[key] = frame_status_counts.get(key, 0) + 1
    manifest = {
        "version": "V17 provisional data-only numerical screening",
        "status": "PROVISIONAL; inherited source/method evidence unresolved",
        "roles": {"train": 43, "development_validation": 3, "test": 0},
        "energy_label_counts": {"train": 43, "development_validation": 3, "test": 0},
        "force_label_counts": {"train": 43, "development_validation": 3, "test": 0},
        "train_input_sha256": sha(data / "train.extxyz"),
        "development_validation_input_sha256": sha(data / "valid.extxyz"),
        "composition_matrix_elements": list(ELEMENTS),
        "composition_matrix_rank_exact": rank,
        "composition_matrix": matrix.tolist(),
        "train_provenance_status_counts": frame_status_counts,
        "V16_manifest_per_frame_provenance_records": len(v16_manifest.get("per_frame_provenance", [])),
        "V16_manifest_frame_provenance_count_gap": 35 - len(v16_manifest.get("per_frame_provenance", [])),
        "frame_provenance_source": "B's 35-row hash-pinned v17_training_pool_provenance_audit/ledger.json",
        "inherited_V16_frame_count": 35,
        "newly_verified_V17_pair_label_count": 8,
        "V17_development_validation_count": 3,
        "verified_V17_pair_labels": [label for label, _ in positive + negative],
        "development_validation_labels": [label for label, _ in validation],
        "training_pool_policy": "Preserve all 35 V16 frames with their frame-level provenance metadata; add eight verified V17 +/- diagnostics. Do not represent provisional inherited labels as newly verified.",
        "energy_policy": "Use the stored inherited REF_energy unchanged; map new verified GPAW native extrapolated energies to REF_energy; keep free energy separate in source archives.",
        "validation_policy": "Use the three frozen V17 development structures only. They have been scored against V16 and are not independent final blind tests.",
        "withheld_test_policy": "No V17 withheld output, label, summary, log, or model score was opened. No test split is passed to MACE.",
        "B_aggregate_archive_manifest_was_stale": b_aggregate_stale,
        "B_aggregate_missing_labels_as_recorded": b_aggregate.get("missing_labels", []),
        "per_label_archive_checks": "All eight signed diagnostic archives were independently checked for verification PASS, input/result identity, SHA256, SCF convergence, finite labels, summary energy, and common PW-PBE500/Gamma/0.1 eV settings.",
        "source_files_sha256": dict(sorted(sources.items())),
        "foundation_model_sha256": sha(repo / PROJECT / "mace_periodic_v12_interface_energy/foundation_models/mace-mp-0b3-medium.model"),
        "provenance_audit_sha256": sha(provenance_path),
        "split_leakage_audit_sha256": sha(repo / "research/mace-v12-transfer/coordination/reports/window-b/v17_split_leakage_audit/report.md"),
    }
    (data / "dataset_manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    files = [data / "train.extxyz", data / "valid.extxyz", data / "dataset_manifest.json"]
    (data / "SHA256SUMS.txt").write_text("".join(f"{sha(path)}  {path.name}\n" for path in files))
    print(json.dumps({"status": "PROVISIONAL_BUILD_PASS", "train": len(train), "development_validation": len(validation), "test": 0, "composition_rank_exact": rank, "B_aggregate_manifest_stale": b_aggregate_stale, "provenance_status_counts": frame_status_counts}, indent=2))


if __name__ == "__main__":
    main()
