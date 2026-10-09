#!/usr/bin/env python3
"""Recompute the registered read-only AgSi slab reference-method audit."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read


TASK = "v19_slab_reference_method_readiness_review_window_b"
REPORT_DIR = Path(__file__).resolve().parent
ARCHIVE_ROOT = Path("research/mace-v12-transfer/periodic_interface_v4")
ARCHIVES = {
    "AgSi_COD9009647_pilot_gamma": "pbe_interface_v19_cif_convergence_pilot",
    "AgSi_COD9009647_pilot_k2x2": "pbe_interface_v19_cif_convergence_pilot",
    "AgSi_COD9009647_pilot_k3x3": "pbe_interface_v19_cif_convergence_pilot",
    "AgSi_COD9009647_pilot_k4x4": "pbe_interface_v19_cif_convergence_pilot",
    "AgSi_COD9009647_pilot_k5x5": "pbe_interface_v19_cif_convergence_pilot",
    "AgSi_COD9009647_pilot_k6x6": "pbe_interface_v19_cif_convergence_pilot",
    "AgSi_COD9009647_pilot_k8x8": "pbe_interface_v19_cif_evenmesh_pilot",
    "AgSi_COD9009647_pilot_k5x5_sigma0p05": "pbe_interface_v19_cif_smearing_pilot",
    "AgSi_COD9009647_pilot_k6x6_sigma0p05": "pbe_interface_v19_cif_smearing_pilot",
    "AgSi_COD9009647_pilot_k5x5_sigma0p20": "pbe_interface_v19_cif_smearing020_pilot",
    "AgSi_COD9009647_pilot_k6x6_sigma0p20": "pbe_interface_v19_cif_smearing020_pilot",
    "AgSi_originalvac_slab_dipole_off_v19": "pbe_interface_v19_cif_slab_controls",
    "AgSi_originalvac_slab_dipole_on_v19": "pbe_interface_v19_cif_slab_controls",
}
MESH_LADDER = [
    "AgSi_COD9009647_pilot_gamma",
    "AgSi_COD9009647_pilot_k2x2",
    "AgSi_COD9009647_pilot_k3x3",
    "AgSi_COD9009647_pilot_k4x4",
    "AgSi_COD9009647_pilot_k5x5",
    "AgSi_COD9009647_pilot_k6x6",
    "AgSi_COD9009647_pilot_k8x8",
]
FORCE_KEY = "PW_PBE_forces"
ID_KEY = "lammps_id"
PAIR_KEY = "central_pair"
THRESHOLDS = {
    "absolute_native_energy_delta_meV_atom": 2.0,
    "absolute_free_energy_delta_meV_atom": 2.0,
    "all_atom_force_vector_rms_eV_A": 0.01,
}


def find_repo_root() -> Path:
    for parent in REPORT_DIR.parents:
        if (parent / ".git").exists():
            return parent
    raise RuntimeError("Cannot locate repository root above report directory")


REPO = find_repo_root()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_sha256_manifest(folder: Path) -> dict:
    manifest = folder / "SHA256SUMS.txt"
    if not manifest.is_file():
        raise RuntimeError(f"Missing checksum inventory: {manifest}")
    records = []
    for line in manifest.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        expected, name = line.split(maxsplit=1)
        name = name.lstrip("*")
        path = (folder / name).resolve()
        if folder.resolve() not in path.parents:
            raise RuntimeError(f"Unsafe checksum path in {manifest}: {name}")
        if not path.is_file():
            records.append({"file": name, "expected": expected, "actual": None, "pass": False})
            continue
        actual = sha256(path)
        records.append({"file": name, "expected": expected, "actual": actual, "pass": actual == expected})
    if not records or not all(record["pass"] for record in records):
        raise RuntimeError(f"Checksum verification failed: {manifest}")
    return {"status": "PASS", "files_checked": len(records), "files": records}


def source_path(repo: Path, recorded: str) -> Path:
    path = Path(recorded)
    if path.is_file():
        return path
    marker = "research/mace-v12-transfer/"
    value = recorded.replace("\\", "/")
    index = value.find(marker)
    if index >= 0:
        candidate = repo / value[index:]
        if candidate.is_file():
            return candidate
    if not path.is_absolute():
        candidate = repo / path
        if candidate.is_file():
            return candidate
    raise RuntimeError(f"Cannot resolve archived source input: {recorded}")


def required_arrays(atoms, context: str) -> None:
    for key in (ID_KEY, PAIR_KEY):
        if key not in atoms.arrays:
            raise RuntimeError(f"{context} is missing the {key} identity array")


def geometry_equal(left, right) -> bool:
    required_arrays(left, "left structure")
    required_arrays(right, "right structure")
    return bool(
        np.array_equal(left.numbers, right.numbers)
        and np.array_equal(left.arrays[ID_KEY], right.arrays[ID_KEY])
        and np.array_equal(left.arrays[PAIR_KEY], right.arrays[PAIR_KEY])
        and np.array_equal(left.positions, right.positions)
        and np.array_equal(left.cell.array, right.cell.array)
        and np.array_equal(left.pbc, right.pbc)
    )


def physical_geometry_hash(atoms) -> str:
    """Hash atom identity and geometry, excluding PBC to expose that difference separately."""
    required_arrays(atoms, "structure")
    payload = b"".join(
        (
            np.asarray(atoms.numbers, dtype="<i4").tobytes(),
            np.asarray(atoms.arrays[ID_KEY], dtype="<i8").tobytes(),
            np.asarray(atoms.arrays[PAIR_KEY], dtype="<i4").tobytes(),
            np.asarray(atoms.positions, dtype="<f8").tobytes(),
            np.asarray(atoms.cell.array, dtype="<f8").tobytes(),
        )
    )
    return hashlib.sha256(payload).hexdigest()


def archive_record(repo: Path, label: str, parent: str) -> dict:
    folder = repo / ARCHIVE_ROOT / parent / "calculations" / label
    checksum = verify_sha256_manifest(folder)
    summary_path = folder / "summary.json"
    verification_path = folder / "verification.json"
    summary = json.loads(summary_path.read_text())
    verification = json.loads(verification_path.read_text())
    if summary.get("label") != label or summary.get("scf_converged") is not True:
        raise RuntimeError(f"Archive summary identity/convergence failed for {label}")
    checks = verification.get("checks", {})
    if verification.get("status") != "PASS" or not checks or not all(checks.values()):
        raise RuntimeError(f"Archive verifier did not pass all checks for {label}")
    recorded_source = summary.get("source")
    src = source_path(repo, recorded_source)
    source_hash = sha256(src)
    if source_hash != summary.get("source_sha256"):
        raise RuntimeError(f"Archived source hash mismatch for {label}")
    output = folder / f"{label}_PW_PBE.extxyz"
    if not output.is_file():
        raise RuntimeError(f"Missing result structure for {label}")
    atoms = read(output)
    original = read(src)
    if not geometry_equal(atoms, original):
        raise RuntimeError(f"Archived result structure differs from its source for {label}")
    if FORCE_KEY not in atoms.arrays or atoms.arrays[FORCE_KEY].shape != (len(atoms), 3):
        raise RuntimeError(f"Missing or malformed force vectors for {label}")
    if len(set(np.asarray(atoms.arrays[ID_KEY], dtype=int).tolist())) != len(atoms):
        raise RuntimeError(f"Atom IDs are not unique for {label}")
    method = summary.get("method", {})
    cell_z = float(atoms.cell[2, 2])
    z_span = float(np.ptp(atoms.positions[:, 2]))
    return {
        "label": label,
        "relative_directory": str(folder.relative_to(repo)),
        "sha256_inventory": checksum,
        "archive_verification": {
            "status": verification["status"],
            "checks_passed": len(checks),
            "all_checks_pass": True,
        },
        "source_input": str(src.relative_to(repo)),
        "source_sha256": source_hash,
        "result_extxyz_sha256": sha256(output),
        "scf_converged": summary["scf_converged"],
        "scf_iterations": summary.get("scf_iterations"),
        "atom_count": len(atoms),
        "formula": summary.get("formula"),
        "energy_eV_cell": summary.get("energy_eV_cell"),
        "free_energy_eV_cell": summary.get("free_energy_eV_cell"),
        "energy_convention": summary.get("energy_convention"),
        "method": method,
        "pbc": atoms.pbc.tolist(),
        "cell_z_A": cell_z,
        "coordinate_span_z_A": z_span,
        "geometric_empty_length_A": cell_z - z_span,
        "physical_geometry_sha256_excluding_pbc": physical_geometry_hash(atoms),
        "atom_ids_unique": True,
        "source_and_result_geometry_exact": True,
        "_summary": summary,
        "_atoms": atoms,
    }


def method_except(left: dict, right: dict, varied_key: str) -> dict:
    lhs, rhs = dict(left), dict(right)
    lhs.pop(varied_key, None)
    rhs.pop(varied_key, None)
    keys = sorted(set(lhs) | set(rhs))
    differences = {key: {"left": lhs.get(key), "right": rhs.get(key)} for key in keys if lhs.get(key) != rhs.get(key)}
    return {"only_expected_method_change": not differences, "unexpected_method_differences": differences}


def compare_pair(records: dict, left_label: str, right_label: str, varied_key: str) -> dict:
    left, right = records[left_label], records[right_label]
    la, ra = left["_atoms"], right["_atoms"]
    same_physical = physical_geometry_hash(la) == physical_geometry_hash(ra)
    same_pbc = np.array_equal(la.pbc, ra.pbc)
    same_source = left["source_sha256"] == right["source_sha256"]
    method_check = method_except(left["method"], right["method"], varied_key)
    if not same_physical or not same_pbc:
        raise RuntimeError(f"Geometry/PBC mismatch in matched pair {left_label} -> {right_label}")
    if not method_check["only_expected_method_change"]:
        raise RuntimeError(f"Unexpected method differences in {left_label} -> {right_label}: {method_check}")
    if not same_source:
        raise RuntimeError(f"Source input hashes differ in {left_label} -> {right_label}")
    delta_force = np.asarray(ra.arrays[FORCE_KEY]) - np.asarray(la.arrays[FORCE_KEY])
    atom_vector_delta = np.linalg.norm(delta_force, axis=1)
    n_atoms = left["atom_count"]
    delta_native = (right["energy_eV_cell"] - left["energy_eV_cell"]) * 1000.0 / n_atoms
    delta_free = (right["free_energy_eV_cell"] - left["free_energy_eV_cell"]) * 1000.0 / n_atoms
    force_vector_rms = float(np.sqrt(np.mean(atom_vector_delta**2)))
    return {
        "left": left_label,
        "right": right_label,
        "delta_direction": "right minus left",
        "varied_method_key": varied_key,
        "same_source_input_hash": same_source,
        "same_atom_identity_geometry_cell_and_pbc": same_physical and same_pbc,
        "method_differs_only_in_expected_key": method_check["only_expected_method_change"],
        "kpts_left": left["method"].get("kpts"),
        "kpts_right": right["method"].get("kpts"),
        "smearing_eV_left": left["method"].get("smearing_eV"),
        "smearing_eV_right": right["method"].get("smearing_eV"),
        "pbc": la.pbc.tolist(),
        "native_energy_delta_meV_atom": float(delta_native),
        "free_energy_delta_meV_atom": float(delta_free),
        "all_atom_force_vector_rms_eV_A": force_vector_rms,
        "max_atom_force_vector_delta_eV_A": float(np.max(atom_vector_delta)),
        "componentwise_force_rms_eV_A": float(np.sqrt(np.mean(delta_force**2))),
        "pass": {
            "native_energy": abs(delta_native) <= THRESHOLDS["absolute_native_energy_delta_meV_atom"],
            "free_energy": abs(delta_free) <= THRESHOLDS["absolute_free_energy_delta_meV_atom"],
            "all_atom_force_vector_rms": force_vector_rms <= THRESHOLDS["all_atom_force_vector_rms_eV_A"],
        },
    }


def proposal_audit(repo: Path) -> dict:
    folder = repo / "research/mace-v12-transfer/coordination/reports/window-b/v19_gap_scan_dipole_vacuum_control_proposals_window_b"
    inventory = verify_sha256_manifest(folder)
    manifest = json.loads((folder / "proposal_manifest.json").read_text())
    verification = json.loads((folder / "verification.json").read_text())
    if verification.get("status") != "PASS" or not all(verification.get("checks", {}).values()):
        raise RuntimeError("The published 26A/36A proposal verification is not PASS")
    if manifest.get("launch_enabled") is not False or manifest.get("owner") is not None:
        raise RuntimeError("A vacuum proposal is unexpectedly owned or launch-enabled")
    off_name = "AgSi_26A_vacuum_parent_dipole_OFF_PROPOSAL.extxyz"
    on_name = "AgSi_36A_vacuum_parent_dipole_ON_PROPOSAL.extxyz"
    off_path, on_path = folder / off_name, folder / on_name
    off, on = read(off_path), read(on_path)
    source = source_path(repo, manifest["proposals"][0]["source_input"])
    parent_geometry = source_path(repo, manifest["source_record"]["parent_geometry_path"])
    if off_path.read_bytes() != source.read_bytes():
        raise RuntimeError("The 26A OFF proposal is not byte-identical to its registered parent")
    if sha256(parent_geometry) != manifest["source_record"]["parent_geometry_sha256"]:
        raise RuntimeError("The proposal's shared CIF-parent geometry hash does not match its manifest")
    if not (
        np.array_equal(off.numbers, on.numbers)
        and np.array_equal(off.arrays[ID_KEY], on.arrays[ID_KEY])
        and np.array_equal(off.arrays[PAIR_KEY], on.arrays[PAIR_KEY])
        and np.array_equal(off.cell.array[:2], on.cell.array[:2])
        and np.array_equal(off.pbc, on.pbc)
        and np.allclose(on.positions - off.positions, np.array([0.0, 0.0, 5.0]), rtol=0, atol=1e-12)
        and abs(float(on.cell[2, 2] - off.cell[2, 2]) - 10.0) <= 1e-12
    ):
        raise RuntimeError("The 26A/36A geometry identity or expansion checks failed")
    results = []
    for name, atoms, method in (
        (off_name, off, manifest["proposals"][0]["method"]),
        (on_name, on, manifest["proposals"][1]["method"]),
    ):
        span = float(np.ptp(atoms.positions[:, 2]))
        results.append(
            {
                "file": name,
                "sha256": sha256(folder / name),
                "pbc": atoms.pbc.tolist(),
                "cell_z_A": float(atoms.cell[2, 2]),
                "coordinate_span_z_A": span,
                "geometric_empty_length_A": float(atoms.cell[2, 2]) - span,
                "poissonsolver": method.get("poissonsolver"),
                "has_energy_or_force_labels": any(key in atoms.info for key in ("energy", "free_energy", "forces")) or FORCE_KEY in atoms.arrays,
            }
        )
    return {
        "sha256_inventory": inventory,
        "verification_status": verification["status"],
        "launch_enabled": manifest["launch_enabled"],
        "owner": manifest["owner"],
        "source_parent_hash": manifest["source_record"]["parent_geometry_sha256"],
        "source_parent_geometry_sha256_verified": sha256(parent_geometry),
        "source_input_sha256": sha256(source),
        "proposals": results,
        "interpretation": "The proposals change only geometry/method as recorded, but neither has energy/force results. The 26A OFF versus 36A ON comparison changes vacuum and dipole together and is therefore confounded.",
    }


def ti_bulk_report(repo: Path) -> dict:
    folder = repo / "research/mace-v12-transfer/coordination/reports/window-a/v19_Ti3SiC2_bulk_mesh_comparison"
    inventory = verify_sha256_manifest(folder)
    matrix = json.loads((folder / "mesh_matrix_kz2_kz4.json").read_text())
    selected = [
        pair
        for pair in matrix.get("pairs", [])
        if (pair.get("left"), pair.get("right"))
        in {
            ("Ti3SiC2_baseline_k8x8x2_v19", "Ti3SiC2_baseline_k10x10x2_v19"),
            ("Ti3SiC2_baseline_k8x8x4_v19", "Ti3SiC2_baseline_k10x10x4_v19"),
        }
    ]
    if len(selected) != 2:
        raise RuntimeError("A's Ti3SiC2 fixed-geometry report lacks the 8x8 -> 10x10 matched pairs")
    return {
        "sha256_inventory": inventory,
        "status": matrix.get("status"),
        "scope": matrix.get("scope"),
        "thresholds": matrix.get("thresholds"),
        "8x8_to_10x10_pairs": [
            {
                "left": pair["left"],
                "right": pair["right"],
                "native_energy_delta_meV_atom": pair["native_energy_delta_meV_atom"],
                "free_energy_delta_meV_atom": pair["free_energy_delta_meV_atom"],
                "force_vector_rms_eV_A": pair["force_vector_rms_eV_A"],
                "pass": pair["pass"],
            }
            for pair in selected
        ],
        "slab_transferable": False,
        "scope_note": "Useful bulk-only control; it does not prove AgSi slab/interface mesh convergence.",
    }


def clean_record(record: dict) -> dict:
    return {key: value for key, value in record.items() if not key.startswith("_")}


def main() -> None:
    records = {}
    for label, parent in ARCHIVES.items():
        records[label] = archive_record(REPO, label, parent)

    mesh_pairs = [
        compare_pair(records, left, right, "kpts")
        for left, right in zip(MESH_LADDER, MESH_LADDER[1:])
    ]
    smearing_pairs = []
    for mesh in ("k5x5", "k6x6"):
        low = f"AgSi_COD9009647_pilot_{mesh}_sigma0p05"
        mid = f"AgSi_COD9009647_pilot_{mesh}"
        high = f"AgSi_COD9009647_pilot_{mesh}_sigma0p20"
        smearing_pairs.append(compare_pair(records, low, mid, "smearing_eV"))
        smearing_pairs.append(compare_pair(records, mid, high, "smearing_eV"))
    dipole_pair = compare_pair(
        records,
        "AgSi_originalvac_slab_dipole_off_v19",
        "AgSi_originalvac_slab_dipole_on_v19",
        "poissonsolver",
    )

    mesh_pbcs = {tuple(records[label]["pbc"]) for label in MESH_LADDER}
    width_pbcs = {tuple(records[label]["pbc"]) for label in ARCHIVES if "sigma" in label}
    slab_pbcs = {
        tuple(records[label]["pbc"])
        for label in ("AgSi_originalvac_slab_dipole_off_v19", "AgSi_originalvac_slab_dipole_on_v19")
    }
    mesh_geometry_hashes = {records[label]["physical_geometry_sha256_excluding_pbc"] for label in MESH_LADDER}
    slab_geometry_hashes = {
        records[label]["physical_geometry_sha256_excluding_pbc"]
        for label in ("AgSi_originalvac_slab_dipole_off_v19", "AgSi_originalvac_slab_dipole_on_v19")
    }

    proposals = proposal_audit(REPO)
    bulk = ti_bulk_report(REPO)
    source_hashes = {records[label]["source_sha256"] for label in MESH_LADDER}

    metrics = {
        "task": TASK,
        "status": "PASS_READ_ONLY_AUDIT",
        "scope": "Existing published AgSi CIF slab archives/proposals and A's published Ti3SiC2 fixed-geometry mesh matrix; no DFT/MPI or input changes.",
        "recomputed_with": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "ase": __import__("ase").__version__,
            "force_array": FORCE_KEY,
            "force_vector_rms_definition": "sqrt(mean_i(sum_xyz((F_right-F_left)^2))) over atoms; atom vector magnitudes are RMS-averaged, not componentwise RMS.",
        },
        "thresholds": THRESHOLDS,
        "archive_count": len(records),
        "all_archive_checksums_and_verifiers_pass": True,
        "archives": [clean_record(records[label]) for label in ARCHIVES],
        "geometry_and_boundary_audit": {
            "mesh_ladder_single_physical_geometry": len(mesh_geometry_hashes) == 1,
            "smearing_archives_share_mesh_geometry": all(
                records[label]["physical_geometry_sha256_excluding_pbc"] in mesh_geometry_hashes
                for label in ARCHIVES
                if "sigma" in label
            ),
            "dipole_pair_same_physical_geometry": len(slab_geometry_hashes) == 1,
            "mesh_ladder_source_hash_count": len(source_hashes),
            "mesh_ladder_pbc_values": [list(value) for value in sorted(mesh_pbcs)],
            "smearing_pbc_values": [list(value) for value in sorted(width_pbcs)],
            "dipole_pair_pbc_values": [list(value) for value in sorted(slab_pbcs)],
            "mesh_and_smearing_pbc_match_target_nonperiodic_z_slab": mesh_pbcs == slab_pbcs and width_pbcs == slab_pbcs,
            "physical_geometry_hash_excluding_pbc": next(iter(mesh_geometry_hashes)) if len(mesh_geometry_hashes) == 1 else None,
            "slab_geometry_hash_excluding_pbc": next(iter(slab_geometry_hashes)) if len(slab_geometry_hashes) == 1 else None,
            "interpretation": "Mesh and smearing archives are fully periodic in z; the off/on slab controls are nonperiodic in z. Same atomic geometry does not remove this method difference.",
        },
        "mesh_ladder_pairs": mesh_pairs,
        "smearing_pairs": smearing_pairs,
        "dipole_off_to_on_pair": dipole_pair,
        "vacuum_proposals": proposals,
        "ti3sic2_bulk_mesh_report": bulk,
        "recommendation": {
            "status": "PROVISIONAL_REFERENCE_RECIPE; NOT slab-mesh/smearing converged",
            "recipe": {
                "xc": "PBE",
                "basis": "plane wave",
                "cutoff_eV": 500,
                "kpts": [6, 6, 1],
                "smearing_eV": 0.1,
                "pbc": [True, True, False],
                "poissonsolver": {"dipolelayer": "xy"},
                "GPAW_version": "26.7.0",
                "PAW_data_version": "1.2.1",
                "energy_convention": "GPAW native extrapolated energy (force_consistent=False); report free_energy separately",
            },
            "basis": "This is the closest already archived nonperiodic-z reference method: the 6x6x1, 0.10 eV, 500 eV, PBE dipole OFF/ON pair uses identical geometry/cell and differs only by the xy dipole correction. The 26A/36A proposal settings retain the same general method family. This is a starting recipe, not a claim of slab convergence.",
            "controls_before_enabling_gap_scan": [
                "Run a matched 4x4x1/6x6x1/8x8x1 in-plane mesh ladder on the target nonperiodic-z slab, keeping cell, dipole setting, smearing, cutoff, PBE, and geometry fixed; compare native/free energies and all-atom force-vector RMS.",
                "Run sigma 0.05/0.10/0.20 eV at fixed target slab geometry, boundary, mesh, and cutoff; compare both energy conventions and force vectors.",
                "Isolate vacuum from dipole by matching 26A and 36A cells under the same dipole setting. To compare dipole treatment at both vacuum widths, complete the 2x2 vacuum-by-dipole matrix; the current 26A-OFF/36A-ON proposal pair is confounded.",
                "After freezing reference settings, evaluate forces and separation-force projections on the intended AgSi gap geometries before enabling the gap scan.",
            ],
            "not_supported_by_current_evidence": [
                "The fully periodic-z mesh or smearing ladder proves convergence for the nonperiodic-z slab.",
                "The 26A or 36A vacuum proposals have converged energies or forces.",
                "Ti3SiC2 bulk mesh convergence establishes AgSi slab/interface force accuracy.",
            ],
        },
    }
    output = REPORT_DIR / "metrics.json"
    output.write_text(json.dumps(metrics, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    print(f"Wrote {output.relative_to(REPO)}")
    print(f"Verified archives: {len(records)}; proposal verification: {proposals['verification_status']}; bulk matrix: {bulk['status']}")
    print(f"Mesh comparisons: {len(mesh_pairs)}; smearing comparisons: {len(smearing_pairs)}; dipole pair: {dipole_pair['pass']}")
    print("Mesh ladder PBC:", metrics["geometry_and_boundary_audit"]["mesh_ladder_pbc_values"])
    print("Dipole-control PBC:", metrics["geometry_and_boundary_audit"]["dipole_pair_pbc_values"])
    print("Recommendation:", metrics["recommendation"]["status"])


if __name__ == "__main__":
    main()
