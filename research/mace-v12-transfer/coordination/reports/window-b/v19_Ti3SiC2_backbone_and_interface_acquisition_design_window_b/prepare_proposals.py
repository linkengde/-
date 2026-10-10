#!/usr/bin/env python3
"""Build launch-disabled geometry proposals for the assigned B review."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from ase.data import atomic_numbers, covalent_radii
from ase.geometry import find_mic
from ase.io import read, write
from ase.neighborlist import neighbor_list


REPORT_DIR = Path(__file__).resolve().parent
ROOT = next(p for p in (REPORT_DIR, *REPORT_DIR.parents) if (p / ".git").exists())
PURE_DIR = ROOT / "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls"
SOURCE_DIR = ROOT / "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures"
V18_AUDIT_DIR = ROOT / "research/mace-v12-transfer/coordination/reports/window-b/v19_V18_Ti3SiC2_training_coverage_audit_window_b"
INTERFACE_PREP_DIR = ROOT / "research/mace-v12-transfer/coordination/reports/window-b/v19_interface_normal_separation_entry_preparation"
TASK_ID = "v19_Ti3SiC2_backbone_and_interface_acquisition_design_window_b"
COD_ID = "9009647"
COD_URL = "https://www.crystallography.net/cod/9009647.cif"
COD_DOI = "10.1016/S0022-3697(98)00226-1"
OVERLAP_SCALE = 0.80
DUPLICATE_TOL_A = 1.0e-8


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_hash(obj: object) -> str:
    data = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return sha256_bytes(data)


def verify_sums(directory: Path) -> dict[str, str]:
    sums_path = directory / "SHA256SUMS.txt"
    expected: dict[str, str] = {}
    for raw in sums_path.read_text(encoding="utf-8").splitlines():
        raw = raw.strip()
        if not raw:
            continue
        digest, rel = raw.split(maxsplit=1)
        rel = rel.lstrip("* ")
        target = directory / rel
        if not target.is_file() or sha256_file(target) != digest:
            raise ValueError(f"checksum mismatch in {directory}: {rel}")
        expected[rel] = digest
    return expected


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def atom_ids(atoms) -> list[int]:
    key = "local_new_id" if "local_new_id" in atoms.arrays else "lammps_id" if "lammps_id" in atoms.arrays else None
    if key is None:
        raise ValueError("proposal input is missing its declared atom ID array")
    values = [int(x) for x in atoms.arrays[key]]
    if len(set(values)) != len(values):
        raise ValueError("duplicate local_new_id")
    return values


def assert_unlabelled(atoms, label: str) -> None:
    label_arrays = {"forces", "energies", "energy", "free_energy"} & set(atoms.arrays)
    label_info = {"energy", "free_energy", "forces"} & set(atoms.info)
    calc_results = getattr(getattr(atoms, "calc", None), "results", {}) or {}
    label_calc = {"energy", "free_energy", "forces"} & set(calc_results)
    if label_arrays or label_info or label_calc:
        raise ValueError(f"label-bearing input is not allowed for proposal {label}")


def formula_counts(atoms) -> dict[str, int]:
    return dict(sorted(Counter(atoms.get_chemical_symbols()).items()))


def ordered_identity_sha(atoms) -> str:
    payload = [{"id": i, "symbol": s} for i, s in zip(atom_ids(atoms), atoms.get_chemical_symbols())]
    return canonical_hash(payload)


def periodic_pair_minima(atoms, cutoff_A: float = 10.0) -> dict[str, float]:
    ii, jj, shifts, distances = neighbor_list("ijSd", atoms, cutoff_A, self_interaction=True)
    minima: dict[str, float] = {}
    for i, j, shift, distance in zip(ii, jj, shifts, distances):
        i, j = int(i), int(j)
        shift_tuple = tuple(int(x) for x in shift)
        if i == j:
            if shift_tuple == (0, 0, 0) or shift_tuple <= (0, 0, 0):
                continue
        elif i > j:
            continue
        pair = "-".join(sorted((atoms[i].symbol, atoms[j].symbol)))
        minima[pair] = min(minima.get(pair, float("inf")), float(distance))
    return {k: round(v, 10) for k, v in sorted(minima.items())}


def overlap_audit(atoms) -> dict:
    distances = periodic_pair_minima(atoms)
    results = {}
    for pair, distance in distances.items():
        first, second = pair.split("-")
        floor = OVERLAP_SCALE * (
            covalent_radii[atomic_numbers[first]] + covalent_radii[atomic_numbers[second]]
        )
        results[pair] = {
            "minimum_periodic_distance_A": distance,
            "screen_floor_A": round(float(floor), 10),
            "margin_A": round(distance - float(floor), 10),
            "pass": bool(distance >= floor - 1.0e-10),
        }
    return {"scale": OVERLAP_SCALE, "pairs": results, "pass": all(x["pass"] for x in results.values())}


def geometry_fingerprint(atoms) -> str:
    scaled = atoms.get_scaled_positions(wrap=True)
    payload = {
        "symbols": atoms.get_chemical_symbols(),
        "local_new_id": atom_ids(atoms),
        "pbc": [bool(x) for x in atoms.pbc],
        "cell_A": np.round(np.asarray(atoms.cell), 10).tolist(),
        "scaled_positions": np.round(scaled, 10).tolist(),
    }
    return canonical_hash(payload)


def same_geometry(a, b) -> bool:
    if (
        a.get_chemical_symbols() != b.get_chemical_symbols()
        or atom_ids(a) != atom_ids(b)
        or list(map(bool, a.pbc)) != list(map(bool, b.pbc))
        or not np.allclose(np.asarray(a.cell), np.asarray(b.cell), atol=1.0e-10, rtol=0.0)
    ):
        return False
    mic, _ = find_mic(a.positions - b.positions, a.cell, a.pbc)
    return bool(np.max(np.linalg.norm(mic, axis=1), initial=0.0) <= DUPLICATE_TOL_A)


def decorate(atoms, *, label: str, family: str, source_sha: str, method_hash: str | None) -> None:
    atoms.calc = None
    atoms.info["proposal_label"] = label
    atoms.info["proposal_role"] = "unlabelled_geometry_only_pending_A_review"
    atoms.info["launch_enabled"] = False
    atoms.info["source_family"] = family
    atoms.info["source_geometry_sha256"] = source_sha
    atoms.info["method_reference_sha256"] = method_hash or "pending_A_reference_method_review"


def verify_identities_unchanged(before, after, label: str) -> None:
    if (
        before.get_chemical_symbols() != after.get_chemical_symbols()
        or atom_ids(before) != atom_ids(after)
        or not np.allclose(np.asarray(before.cell), np.asarray(after.cell), atol=1e-12, rtol=0.0)
        or not np.array_equal(before.pbc, after.pbc)
    ):
        raise ValueError(f"identity/order/cell/PBC changed in {label}")


def add_record(records: list[dict], atoms_by_label: dict[str, object], *, label: str, path: Path,
               parent, source_family: str, source_description: dict, operation: dict,
               method_hash: str | None, method_reference: dict | None,
               external_source_verified: bool) -> None:
    atoms = atoms_by_label[label]
    assert_unlabelled(atoms, label)
    overlap = overlap_audit(atoms)
    if not overlap["pass"]:
        raise ValueError(f"periodic species-aware overlap screen failed for {label}")
    verify_identities_unchanged(parent, atoms, label)
    write(path, atoms, format="extxyz")
    reread = read(path, format="extxyz")
    assert_unlabelled(reread, label)
    verify_identities_unchanged(parent, reread, label)
    if geometry_fingerprint(reread) != geometry_fingerprint(atoms):
        raise ValueError(f"extxyz round-trip changed geometry for {label}")
    rec = {
        "label": label,
        "path": path.relative_to(REPORT_DIR).as_posix(),
        "sha256": sha256_file(path),
        "atom_count": len(atoms),
        "formula_counts": formula_counts(atoms),
        "ordered_identity_sha256": ordered_identity_sha(atoms),
        "geometry_fingerprint_sha256": geometry_fingerprint(atoms),
        "source_family": source_family,
        "source": source_description,
        "parent_input_sha256": sha256_file_from_record(source_description),
        "operation": operation,
        "pbc": [bool(x) for x in atoms.pbc],
        "cell_A": np.round(np.asarray(atoms.cell), 12).tolist(),
        "periodic_species_pair_overlap_screen": overlap,
        "method_reference_sha256": method_hash,
        "method_reference": method_reference,
        "external_source_bytes_verified": external_source_verified,
        "energy_label_present": False,
        "force_label_present": False,
        "launch_enabled": False,
        "owner": None,
    }
    records.append(rec)


def sha256_file_from_record(source: dict) -> str:
    return str(source.get("parent_input_sha256") or source.get("input_sha256") or source.get("sha256"))


def main() -> None:
    REPORT_DIR.joinpath("inputs").mkdir(parents=True, exist_ok=True)
    # These files were generated before a steering update clarified that an
    # externally unverifiable public source must remain blocked.
    for interface in ("AgC", "AgSi", "AgTi"):
        for stale in (REPORT_DIR / "inputs").glob(f"{interface}_COD{COD_ID}_*.extxyz"):
            stale.unlink()
    pure_manifest = read_json(PURE_DIR / "input_manifest.json")
    pure_record = next(r for r in pure_manifest["records"] if r["label"] == "Ti3SiC2_baseline_k10x10x4_v19")
    bulk_input = PURE_DIR / pure_record["input"]
    bulk_sha = sha256_file(bulk_input)
    if bulk_sha != pure_record["input_sha256"]:
        raise ValueError("bulk baseline source hash does not match registered manifest")
    method = pure_record["method"]
    method_hash = canonical_hash(method)
    bulk_parent = read(bulk_input, format="extxyz")
    assert_unlabelled(bulk_parent, "Ti3SiC2 baseline")
    if len(bulk_parent) != 48 or formula_counts(bulk_parent) != {"C": 16, "Si": 8, "Ti": 24}:
        raise ValueError("unexpected Ti3SiC2 baseline composition")
    if not np.array_equal(bulk_parent.pbc, np.array([True, True, True])):
        raise ValueError("unexpected Ti3SiC2 baseline PBC")
    bulk_ids = atom_ids(bulk_parent)
    if bulk_ids != list(range(1, 49)):
        raise ValueError("unexpected Ti3SiC2 baseline local ID order")
    if "spacegroup_kinds" not in bulk_parent.arrays:
        raise ValueError("baseline lacks spacegroup_kinds needed to document representative sites")
    bulk_by_id = {int(atom_id): i for i, atom_id in enumerate(bulk_ids)}
    ti_id, c_id = 3, 9
    ti_i, c_i = bulk_by_id[ti_id], bulk_by_id[c_id]
    if bulk_parent[ti_i].symbol != "Ti" or int(bulk_parent.arrays["spacegroup_kinds"][ti_i]) != 1:
        raise ValueError("selected Ti is not the intended 4f representative")
    if bulk_parent[c_i].symbol != "C" or int(bulk_parent.arrays["spacegroup_kinds"][c_i]) != 3:
        raise ValueError("selected C is not the intended 4f representative")
    bulk_overlap = overlap_audit(bulk_parent)
    if not bulk_overlap["pass"]:
        raise ValueError("baseline fails periodic species-aware overlap screen")

    source_manifest_path = SOURCE_DIR / "manifest.json"
    source_manifest_hash = sha256_file(source_manifest_path)
    source_manifest = read_json(source_manifest_path)
    cif_path = SOURCE_DIR / "source/Ti3SiC2_COD_9009647.cif"
    cif_hash = sha256_file(cif_path)
    if cif_hash != source_manifest["source_cif_sha256"]:
        raise ValueError("source CIF hash does not match source manifest")
    source_package_sums = verify_sums(SOURCE_DIR)
    source_records = {r["interface"]: r for r in source_manifest["records"]}
    audit_sums = verify_sums(V18_AUDIT_DIR)
    train_csv = V18_AUDIT_DIR / "training_provenance.csv"
    training_rows = list(csv.DictReader(train_csv.open(newline="", encoding="utf-8")))
    if len(training_rows) != 30 or any(r.get("proposal_role") != "train" for r in training_rows):
        raise ValueError("train-only provenance report has unexpected row count or roles")
    hash_columns = [
        "source_structure_sha256", "source_geometry_sha256", "source_frame_geometry_sha256",
        "parent_file_sha256", "parent_geometry_sha256",
    ]
    training_geometry_hashes = {
        r[col].strip().lower()
        for r in training_rows for col in hash_columns
        if r.get(col, "").strip()
    }
    training_text = " ".join(" ".join(r.values()) for r in training_rows)
    train_has_cod = COD_ID in training_text
    explicit_source_count = sum(bool(r.get("source_label") or r.get("source_config_type")) for r in training_rows)
    family_counts = dict(sorted(Counter(r.get("interface_family", "") for r in training_rows).items()))

    proposal_records: list[dict] = []
    atoms_by_label: dict[str, object] = {}
    parents_by_label: dict[str, object] = {}
    source_by_label: dict[str, dict] = {}

    bulk_family = f"Ti3SiC2 bulk; COD {COD_ID}"
    bulk_source = {
        "type": "verified_existing_pure_phase_baseline_geometry",
        "input_path": str(bulk_input.relative_to(ROOT)),
        "input_sha256": bulk_sha,
        "source_cif_sha256": cif_hash,
        "source_family": bulk_family,
        "source_manifest_sha256": source_manifest_hash,
        "parent_input_sha256": bulk_sha,
    }
    bulk_operations = [
        ("Ti", ti_id, "Ti4f", -0.02), ("Ti", ti_id, "Ti4f", 0.02),
        ("C", c_id, "C4f", -0.02), ("C", c_id, "C4f", 0.02),
    ]
    for species, atom_id, site, dz in bulk_operations:
        sign = "m" if dz < 0 else "p"
        label = f"Ti3SiC2_{site}_id{atom_id:04d}_z_{sign}020A_PROPOSAL"
        atoms = bulk_parent.copy()
        atoms.calc = None
        atoms.positions[bulk_by_id[atom_id], 2] += dz
        expected_delta = np.zeros_like(atoms.positions)
        expected_delta[bulk_by_id[atom_id], 2] = dz
        if not np.allclose(atoms.positions - bulk_parent.positions, expected_delta, atol=1e-12, rtol=0.0):
            raise ValueError(f"bulk proposal changed atoms outside the requested displacement: {label}")
        decorate(atoms, label=label, family=bulk_family, source_sha=bulk_sha, method_hash=method_hash)
        atoms_by_label[label] = atoms
        parents_by_label[label] = bulk_parent
        source_by_label[label] = bulk_source
        add_record(
            proposal_records, atoms_by_label, label=label, path=REPORT_DIR / "inputs" / f"{label}.extxyz",
            parent=bulk_parent, source_family=bulk_family, source_description=bulk_source,
            operation={"kind": "single_atom_z_displacement", "species": species, "local_new_id": atom_id,
                       "representative_site": site, "delta_z_A": dz, "fixed_cell": True},
            method_hash=method_hash, method_reference=method, external_source_verified=False,
        )

    interface_specs = {"AgC": "C", "AgSi": "Si", "AgTi": "Ti"}
    blocked_interfaces = []
    for interface, termination in interface_specs.items():
        source_record = source_records[interface]
        parent_path = ROOT / source_record["path"]
        parent_sha = sha256_file(parent_path)
        if parent_sha != source_record["sha256"]:
            raise ValueError(f"{interface} parent source hash does not match manifest")
        blocked_interfaces.append({
            "chemistry": interface,
            "termination_species": termination,
            "status": "BLOCKED_SOURCE_BYTES_NOT_EXTERNALLY_VERIFIED",
            "proposal_generated": False,
            "proposal_geometry_count": 0,
            "public_record": COD_ID,
            "doi": COD_DOI,
            "source_url": COD_URL,
            "source_cif_sha256": cif_hash,
            "source_manifest_sha256": source_manifest_hash,
            "parent_label": source_record["label"],
            "parent_path": str(parent_path.relative_to(ROOT)),
            "parent_input_sha256": parent_sha,
            "repository_parent_hash_matches_manifest": True,
            "repository_CIF_package_hash_checks_pass": True,
            "external_source_bytes_verified": False,
            "blocker": "The declared COD URL returned HTTP 403 twice. Per task boundary, do not generate parent or gap-scan geometries until the public source bytes can be verified.",
            "shared_source_family": f"COD {COD_ID} Ti3SiC2 CIF; shared among all three terminations if later verified",
            "V18_train_chemistry_present": True,
            "exact_COD_or_parent_hash_in_train_only_provenance": False,
            "train_provenance_limitation": "No matching accession or exact recorded hash; some V18 train rows lack explicit source hashes.",
        })

    # Every candidate must be distinct under identity-preserving minimum-image comparison.
    seen = []
    for label, atoms in atoms_by_label.items():
        for previous_label, previous_atoms in seen:
            if same_geometry(atoms, previous_atoms):
                raise ValueError(f"duplicate geometry: {label} == {previous_label}")
        seen.append((label, atoms))

    for rec in proposal_records:
        rec["source_family_novelty_vs_v18_train_provenance"] = {
            "cod_accession_text_present": train_has_cod,
            "exact_source_geometry_hash_match": rec["source"]["source_cif_sha256"].lower() in training_geometry_hashes
                or rec["parent_input_sha256"].lower() in training_geometry_hashes,
            "train_family_counts": family_counts,
            "train_frames": len(training_rows),
            "train_rows_with_explicit_source_label_or_type": explicit_source_count,
            "interpretation": "No exact COD accession/source/parent geometry hash appears in the selected train-only provenance report; incomplete legacy hashes limit a stronger family-independence claim.",
        }

    manifest = {
        "task_id": TASK_ID,
        "status": "BULK_PROPOSALS_READY_INTERFACE_SOURCE_BLOCKED",
        "proposal_count": len(proposal_records),
        "launch_enabled": False,
        "owner": None,
        "energy_labels_present": False,
        "force_labels_present": False,
        "source_verification": {
            "COD_record": COD_ID,
            "DOI": COD_DOI,
            "source_url": COD_URL,
            "repository_source_CIF_path": str(cif_path.relative_to(ROOT)),
            "source_CIF_sha256": cif_hash,
            "source_manifest_sha256": source_manifest_hash,
            "repository_source_package_sha256_checks_pass": True,
            "live_public_bytes_verified": False,
            "live_fetch_result": "HTTP 403 on initial and approved retry; remote COD bytes unavailable for hash comparison",
            "scope_note": "The COD accession/DOI and repository source-manifest lineage are recorded; this bundle does not claim a live remote-byte checksum match.",
        },
        "train_only_novelty_basis": {
            "report_path": str((V18_AUDIT_DIR / "training_provenance.csv").relative_to(ROOT)),
            "report_sha256": sha256_file(train_csv),
            "report_package_checks_pass": True,
            "rows": len(training_rows),
            "rows_with_explicit_source_label_or_type": explicit_source_count,
            "family_counts": family_counts,
            "COD_accession_present": train_has_cod,
            "exact_source_hashes_considered": sorted(training_geometry_hashes),
            "limitation": "Several V18 train rows have no explicit source hashes; source-family novelty is based only on this selected train-only provenance report, not asserted as statistical independence.",
        },
        "bulk_reference": {
            "baseline_input_path": str(bulk_input.relative_to(ROOT)),
            "baseline_input_sha256": bulk_sha,
            "baseline_registered_label": pure_record["label"],
            "baseline_method": method,
            "method_manifest_canonical_sha256": method_hash,
            "hash_scheme": "SHA256 of UTF-8 JSON serialized with sorted keys and compact separators",
            "method_status": "reference only; proposals are not registered jobs",
            "representative_sites": {
                "Ti": {"local_new_id": ti_id, "spacegroup_kind": int(bulk_parent.arrays["spacegroup_kinds"][ti_i]), "Wyckoff_label": "4f"},
                "C": {"local_new_id": c_id, "spacegroup_kind": int(bulk_parent.arrays["spacegroup_kinds"][c_i]), "Wyckoff_label": "4f"},
            },
            "geometry_source_CIF_sha256": cif_hash,
        },
        "interface_settings": "pending A slab mesh/smearing/vacuum/PBC/dipole reference review; no method assigned",
        "shared_interface_source_family_count": 1,
        "blocked_interface_proposals": blocked_interfaces,
        "proposals": proposal_records,
    }
    manifest_path = REPORT_DIR / "proposal_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    ledger_path = REPORT_DIR / "source_novelty_ledger.csv"
    ledger_rows = []
    ledger_rows.append({
        "proposal_family": "Ti3SiC2_bulk", "proposal_status": "PROPOSED_GEOMETRY_ONLY",
        "proposal_geometry_count": 4, "chemistry_present_in_v18_train": False,
        "public_source_id": f"COD {COD_ID}", "source_cif_sha256": cif_hash,
        "parent_input_sha256": bulk_sha,
        "exact_source_or_parent_hash_in_v18_train_provenance": False,
        "COD_accession_in_v18_train_provenance": train_has_cod,
        "source_record_verification": "existing baseline input hash matches its registered manifest; local CIF/package hash matches source manifest; live COD fetch HTTP 403",
        "source_family_count_across_this_bundle": 1,
        "novelty_assessment": "no exact source/COD match in selected train-only provenance; no pure Ti3SiC2 training frame per completed audit",
        "limitation": "train-only provenance contains incomplete source hashes; report does not prove statistical independence",
    })
    for blocked in blocked_interfaces:
        interface = blocked["chemistry"]
        parent_sha = blocked["parent_input_sha256"]
        hash_match = parent_sha.lower() in training_geometry_hashes or cif_hash.lower() in training_geometry_hashes
        ledger_rows.append({
            "proposal_family": interface,
            "proposal_status": "BLOCKED_SOURCE_BYTES_NOT_EXTERNALLY_VERIFIED",
            "proposal_geometry_count": 0,
            "chemistry_present_in_v18_train": True,
            "public_source_id": f"COD {COD_ID}",
            "source_cif_sha256": cif_hash,
            "parent_input_sha256": parent_sha,
            "exact_source_or_parent_hash_in_v18_train_provenance": hash_match,
            "COD_accession_in_v18_train_provenance": train_has_cod,
            "source_record_verification": "local CIF/package and parent hashes match repository manifest; live COD fetch HTTP 403",
            "source_family_count_across_this_bundle": 1,
            "novelty_assessment": "same chemistry occurs in V18 train; this COD parent hash/accession is not recorded in selected train-only provenance; candidate withheld pending public-source verification",
            "limitation": "train-only provenance has incomplete source hashes; three interfaces share one possible CIF family and no geometry was proposed",
        })
    with ledger_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(dict.fromkeys(key for row in ledger_rows for key in row))
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(ledger_rows)

    verify = {
        "proposal_count": len(proposal_records),
        "blocked_interface_count": len(blocked_interfaces),
        "blocked_interface_labels": [x["chemistry"] for x in blocked_interfaces],
        "interface_geometry_generated": False,
        "all_launch_disabled": all(not r["launch_enabled"] for r in proposal_records),
        "all_unowned": all(r["owner"] is None for r in proposal_records),
        "all_unlabelled": all(not r["energy_label_present"] and not r["force_label_present"] for r in proposal_records),
        "all_overlap_screens_pass": all(r["periodic_species_pair_overlap_screen"]["pass"] for r in proposal_records),
        "all_geometry_fingerprints_unique": len({r["geometry_fingerprint_sha256"] for r in proposal_records}) == len(proposal_records),
        "all_ID_order_preserved": True,
        "all_parent_and_input_hashes_match": (
            all(x["repository_parent_hash_matches_manifest"] for x in blocked_interfaces)
            and all(r["sha256"] == sha256_file(REPORT_DIR / r["path"]) for r in proposal_records)
        ),
        "CIF_hash": cif_hash,
        "source_manifest_hash": source_manifest_hash,
        "source_package_entries_verified": len(source_package_sums),
        "V18_train_provenance_sha256": sha256_file(train_csv),
        "V18_train_report_entries_verified": len(audit_sums),
        "remote_COD_hash_verified": False,
        "remote_COD_fetch_result": "HTTP 403",
        "prohibited_data_access": "No dev/test/holdout/sealed geometry or labels; no raw V18 train file opened by this generator.",
        "prohibited_compute": "No DFT/MPI, model inference/training, or dataset integration invoked.",
    }
    if not all(verify[k] for k in (
        "all_launch_disabled", "all_unowned", "all_unlabelled", "all_overlap_screens_pass",
        "all_geometry_fingerprints_unique", "all_ID_order_preserved", "all_parent_and_input_hashes_match",
    )):
        raise ValueError(f"proposal verification failed: {verify}")
    (REPORT_DIR / "validation.json").write_text(json.dumps(verify, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    checksum_lines = []
    for path in sorted(p for p in REPORT_DIR.rglob("*") if p.is_file() and p.name != "SHA256SUMS.txt"):
        checksum_lines.append(f"{sha256_file(path)}  {path.relative_to(REPORT_DIR).as_posix()}")
    (REPORT_DIR / "SHA256SUMS.txt").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
