#!/usr/bin/env python3
"""Audit and compare the published Ti3SiC2 2x2 k-mesh matrix.

Run with the GPAW environment's Python from any working directory. With the
8x8x2 archive absent, this writes verified pairwise comparisons and marks the
factorial decomposition pending; it never estimates the missing point.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
import subprocess
from pathlib import Path

import numpy as np
from ase.io import read


LABELS = {
    (6, 6, 2): "Ti3SiC2_baseline_k6x6x2_v19",
    (6, 6, 4): "Ti3SiC2_baseline_k6x6x4_v19",
    (8, 8, 2): "Ti3SiC2_baseline_k8x8x2_v19",
    (8, 8, 4): "Ti3SiC2_baseline_k8x8x4_v19",
}
THRESHOLDS = {
    "absolute_native_energy_delta_meV_atom": 2.0,
    "absolute_free_energy_delta_meV_atom": 2.0,
    "force_vector_rms_eV_A": 0.01,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_sha_manifest(folder: Path) -> dict[str, str]:
    manifest = folder / "SHA256SUMS.txt"
    checked: dict[str, str] = {}
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        expected, name = line.split(maxsplit=1)
        name = name.lstrip(" *")
        path = folder / name
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"SHA256 mismatch: {path}")
        checked[name] = actual
    return checked


def rms_and_max(force_delta: np.ndarray, symbols: np.ndarray) -> dict:
    vector_norms = np.linalg.norm(force_delta, axis=1)
    component_rms = np.sqrt(np.mean(force_delta**2, axis=0))
    result = {
        "force_vector_rms_eV_A": float(np.sqrt(np.mean(vector_norms**2))),
        "maximum_atom_force_vector_difference_eV_A": float(vector_norms.max()),
        "force_component_rms_eV_A_xyz": component_rms.tolist(),
        "force_component_max_abs_eV_A_xyz": np.max(np.abs(force_delta), axis=0).tolist(),
        "species": {},
    }
    for symbol in sorted(set(symbols.tolist())):
        mask = symbols == symbol
        local = force_delta[mask]
        local_vectors = np.linalg.norm(local, axis=1)
        result["species"][symbol] = {
            "atom_count": int(mask.sum()),
            "force_vector_rms_eV_A": float(np.sqrt(np.mean(local_vectors**2))),
            "maximum_atom_force_vector_difference_eV_A": float(local_vectors.max()),
            "force_component_rms_eV_A_xyz": np.sqrt(np.mean(local**2, axis=0)).tolist(),
            "force_component_max_abs_eV_A_xyz": np.max(np.abs(local), axis=0).tolist(),
        }
    return result


def load_point(root: Path, calc_root: Path, manifest_records: dict, point: tuple[int, int, int]):
    label = LABELS[point]
    record = manifest_records.get(label)
    if record is None:
        return None
    folder = calc_root / label
    if not folder.is_dir():
        return None
    checked_hashes = verify_sha_manifest(folder)
    summary = json.loads((folder / "summary.json").read_text())
    verification = json.loads((folder / "verification.json").read_text())
    if verification.get("status") != "PASS" or not all(verification.get("checks", {}).values()):
        raise ValueError(f"Archive verification failed for {label}")
    if not summary.get("scf_converged"):
        raise ValueError(f"SCF is not converged for {label}")

    source = root / "research/mace-v12-transfer/periodic_interface_v4/" \
        "pbe_interface_v19_pure_phase_controls" / record["input"]
    source_hash = sha256(source)
    if source_hash != record["input_sha256"]:
        raise ValueError(f"Input manifest hash does not match source file for {label}")
    if summary.get("source_sha256") != source_hash:
        raise ValueError(f"Summary source hash mismatch for {label}")
    if verification.get("input_sha256") != source_hash:
        raise ValueError(f"Verification input hash mismatch for {label}")
    if list(record["kpts"]) != list(point) or list(summary["method"]["kpts"]) != list(point):
        raise ValueError(f"Mesh mismatch for {label}")
    if record["method"] != summary["method"]:
        raise ValueError(f"Input manifest method differs from completed summary for {label}")
    if record.get("launch_status") not in ("completed_archive_verified", "registered_ready"):
        raise ValueError(f"Unexpected manifest status for {label}: {record.get('launch_status')}")

    atoms = read(folder / f"{label}_PW_PBE.extxyz", format="extxyz")
    if len(atoms) != int(summary["atoms"]) or atoms.get_chemical_formula() != summary["formula"]:
        raise ValueError(f"Atom count or formula mismatch for {label}")
    if atoms.info.get("PW_PBE_energy_eV") != summary["energy_eV_cell"]:
        raise ValueError(f"Extxyz and summary native energies differ for {label}")
    if len(atoms) != 48 or atoms.get_chemical_formula() != "C16Si8Ti24":
        raise ValueError(f"Unexpected fixed baseline geometry for {label}")

    return {
        "label": label,
        "point": list(point),
        "folder": folder,
        "record": record,
        "summary": summary,
        "verification": verification,
        "atoms": atoms,
        "archive_sha256": checked_hashes,
    }


def assert_same_geometry(points: dict[tuple[int, int, int], dict]) -> None:
    values = list(points.values())
    first = values[0]["atoms"]
    first_summary = values[0]["summary"]
    first_method = copy.deepcopy(first_summary["method"])
    first_method.pop("kpts")
    first_ids = first.arrays["local_new_id"]
    for value in values[1:]:
        atoms = value["atoms"]
        summary = value["summary"]
        method = copy.deepcopy(summary["method"])
        method.pop("kpts")
        if method != first_method:
            raise ValueError("Completed calculations differ by more than k-points")
        if summary["source_sha256"] != first_summary["source_sha256"]:
            raise ValueError("Completed calculations use different input hashes")
        if not np.array_equal(atoms.numbers, first.numbers):
            raise ValueError("Atom symbols/order differ between archives")
        if not np.array_equal(atoms.positions, first.positions):
            raise ValueError("Frozen positions differ between archives")
        if not np.array_equal(atoms.cell.array, first.cell.array):
            raise ValueError("Cells differ between archives")
        if not np.array_equal(atoms.pbc, first.pbc):
            raise ValueError("PBC differ between archives")
        if not np.array_equal(atoms.arrays["local_new_id"], first_ids):
            raise ValueError("Atom IDs/order differ between archives")


def compare_pair(first: dict, second: dict) -> dict:
    a0, a1 = first["atoms"], second["atoms"]
    s0, s1 = first["summary"], second["summary"]
    n = len(a0)
    native_delta = (s1["energy_eV_cell"] - s0["energy_eV_cell"]) * 1000.0 / n
    free_delta = (s1["free_energy_eV_cell"] - s0["free_energy_eV_cell"]) * 1000.0 / n
    force_delta = a1.arrays["PW_PBE_forces"] - a0.arrays["PW_PBE_forces"]
    metrics = rms_and_max(force_delta, np.asarray(a0.get_chemical_symbols()))
    checks = {
        "native_energy": abs(native_delta) <= THRESHOLDS["absolute_native_energy_delta_meV_atom"],
        "free_energy": abs(free_delta) <= THRESHOLDS["absolute_free_energy_delta_meV_atom"],
        "force_vector_rms": metrics["force_vector_rms_eV_A"] <= THRESHOLDS["force_vector_rms_eV_A"],
    }
    return {
        "from": first["label"],
        "to": second["label"],
        "delta_direction": "to minus from",
        "native_energy_delta_meV_atom": float(native_delta),
        "free_energy_delta_meV_atom": float(free_delta),
        "force": metrics,
        "threshold_checks": checks,
    }


def scalar_delta(points: dict, to_point: tuple, from_point: tuple, field: str, scale: float) -> float:
    return (points[to_point]["summary"][field] - points[from_point]["summary"][field]) * scale


def factorial_decomposition(points: dict[tuple[int, int, int], dict]) -> dict:
    a = (6, 6, 2)
    b = (6, 6, 4)
    c = (8, 8, 2)
    d = (8, 8, 4)
    n = len(points[a]["atoms"])
    in_plane_kz2 = compare_pair(points[a], points[c])
    in_plane_kz4 = compare_pair(points[b], points[d])
    c_axis_k6 = compare_pair(points[a], points[b])
    c_axis_k8 = compare_pair(points[c], points[d])
    force = lambda p: points[p]["atoms"].arrays["PW_PBE_forces"]
    symbols = np.asarray(points[a]["atoms"].get_chemical_symbols())
    interaction_forces = (force(d) - force(c)) - (force(b) - force(a))
    interaction_force_stats = rms_and_max(interaction_forces, symbols)
    native_interaction = (
        scalar_delta(points, d, c, "energy_eV_cell", 1000 / n)
        - scalar_delta(points, b, a, "energy_eV_cell", 1000 / n)
    )
    free_interaction = (
        scalar_delta(points, d, c, "free_energy_eV_cell", 1000 / n)
        - scalar_delta(points, b, a, "free_energy_eV_cell", 1000 / n)
    )
    return {
        "directional_contrasts": {
            "in_plane_6_to_8_at_kz2": in_plane_kz2,
            "in_plane_6_to_8_at_kz4": in_plane_kz4,
            "kz_2_to_4_at_in_plane6": c_axis_k6,
            "kz_2_to_4_at_in_plane8": c_axis_k8,
        },
        "difference_in_differences_interaction": {
            "definition": "(in-plane 6->8 at kz=4) - (in-plane 6->8 at kz=2); equivalently (kz 2->4 at 8x8) - (kz 2->4 at 6x6)",
            "native_energy_meV_atom": float(native_interaction),
            "free_energy_meV_atom": float(free_interaction),
            "force": interaction_force_stats,
            "thresholds_are_not_defined_for_interaction": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("report.json"))
    args = parser.parse_args()
    root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
    base = root / "research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls"
    calc_root = base / "calculations"
    manifest = json.loads((base / "input_manifest.json").read_text())
    records = {record["label"]: record for record in manifest["records"]}

    points = {}
    missing = []
    for point, label in LABELS.items():
        loaded = load_point(root, calc_root, records, point)
        if loaded is None:
            missing.append(label)
        else:
            points[point] = loaded
    if not points:
        raise SystemExit("No published mesh archives found")
    assert_same_geometry(points)

    pairs = []
    for first_point, second_point in itertools.combinations(sorted(points), 2):
        pairs.append(compare_pair(points[first_point], points[second_point]))

    report = {
        "task": "v19_Ti3SiC2_mesh_2x2_factorial_review_window_b",
        "status": "pending_missing_point" if missing else "complete",
        "scope": "Frozen 48-atom Ti3SiC2 pure-phase numerical control only; no validation, stability, or interface claim.",
        "grid": {
            "verified_points": {"x".join(map(str, point)): value["label"] for point, value in sorted(points.items())},
            "missing_points": missing,
            "missing_values_inferred": False,
        },
        "identity_checks": {
            "source_input_sha256": next(iter(points.values()))["summary"]["source_sha256"],
            "method_equal_except_kpts": True,
            "same_formula_count_ordered_ids_positions_cell_pbc": True,
            "all_present_archive_manifests_and_verifiers_pass": True,
        },
        "points": {
            "x".join(map(str, point)): {
                "label": value["label"],
                "native_energy_eV_cell": value["summary"]["energy_eV_cell"],
                "free_energy_eV_cell": value["summary"]["free_energy_eV_cell"],
                "scf_iterations": value["summary"]["scf_iterations"],
                "elapsed_s": value["summary"]["elapsed_s"],
                "archive_verified": value["verification"]["status"] == "PASS",
                "archive_sha256": value["archive_sha256"],
            }
            for point, value in sorted(points.items())
        },
        "available_pairwise_comparisons": pairs,
        "thresholds": THRESHOLDS,
        "factorial_decomposition": factorial_decomposition(points) if not missing else None,
        "interpretation": (
            "The fourth point is missing; only available pairwise comparisons are reported. "
            "The in-plane, c-axis, and interaction decomposition is withheld until all four "
            "archives are verified."
            if missing else
            "All four archives passed identity checks; directional effects and the interaction "
            "are reported for this frozen geometry only."
        ),
    }
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
