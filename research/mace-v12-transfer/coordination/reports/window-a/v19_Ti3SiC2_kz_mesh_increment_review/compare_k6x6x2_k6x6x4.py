#!/usr/bin/env python3
"""Recompute the fixed-geometry Ti3SiC2 k_z mesh comparison from verified archives."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np
from ase.io import read

LABELS = ["Ti3SiC2_baseline_k6x6x2_v19", "Ti3SiC2_baseline_k6x6x4_v19"]
THRESHOLDS = {
    "absolute_native_energy_delta_meV_atom": 2.0,
    "absolute_free_energy_delta_meV_atom": 2.0,
    "force_vector_rms_eV_A": 0.01,
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_manifest(folder):
    manifest = folder / "SHA256SUMS.txt"
    checked = {}
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        expected, filename = line.split(maxsplit=1)
        filename = filename.lstrip(" *")
        actual = sha256(folder / filename)
        assert actual == expected, f"SHA256 mismatch: {folder / filename}"
        checked[filename] = actual
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("report.json"))
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"Refusing to overwrite {output}")

    root = Path(subprocess.check_output(
        ["git", "rev-parse", "--show-toplevel"], text=True).strip())
    calc = root / "research/mace-v12-transfer/periodic_interface_v4/" \
        "pbe_interface_v19_pure_phase_controls/calculations"
    rows = []
    archive_checks = {}
    for label in LABELS:
        folder = calc / label
        archive_checks[label] = verify_manifest(folder)
        summary = json.loads((folder / "summary.json").read_text())
        verification = json.loads((folder / "verification.json").read_text())
        assert verification["status"] == "PASS"
        assert all(verification["checks"].values())
        assert summary["scf_converged"] is True
        assert summary["source_sha256"] == verification["input_sha256"]
        atoms = read(folder / f"{label}_PW_PBE.extxyz", format="extxyz")
        rows.append((label, folder, summary, verification, atoms))

    (label0, folder0, summary0, verification0, atoms0), \
    (label1, folder1, summary1, verification1, atoms1) = rows
    assert len(atoms0) == len(atoms1) == 48
    assert atoms0.get_chemical_formula() == atoms1.get_chemical_formula() == "C16Si8Ti24"
    assert np.array_equal(atoms0.numbers, atoms1.numbers)
    assert np.array_equal(atoms0.positions, atoms1.positions)
    assert np.array_equal(atoms0.cell.array, atoms1.cell.array)
    assert np.array_equal(atoms0.pbc, atoms1.pbc)
    assert np.array_equal(atoms0.arrays["local_new_id"], atoms1.arrays["local_new_id"])
    assert summary0["source_sha256"] == summary1["source_sha256"]
    assert verification0["input_sha256"] == verification1["input_sha256"]

    method0 = copy.deepcopy(summary0["method"])
    method1 = copy.deepcopy(summary1["method"])
    kpts0 = method0.pop("kpts")
    kpts1 = method1.pop("kpts")
    assert method0 == method1
    assert kpts0 == [6, 6, 2] and kpts1 == [6, 6, 4]

    force_delta = atoms1.arrays["PW_PBE_forces"] - atoms0.arrays["PW_PBE_forces"]
    atom_force_norms = np.linalg.norm(force_delta, axis=1)
    n_atoms = len(atoms0)
    metrics = {
        "native_energy_delta_meV_atom":
            (summary1["energy_eV_cell"] - summary0["energy_eV_cell"]) * 1000 / n_atoms,
        "free_energy_delta_meV_atom":
            (summary1["free_energy_eV_cell"] - summary0["free_energy_eV_cell"]) * 1000 / n_atoms,
        "force_vector_rms_eV_A":
            float(np.sqrt(np.mean(np.sum(force_delta ** 2, axis=1)))),
        "maximum_atom_force_vector_difference_eV_A": float(atom_force_norms.max()),
        "force_component_rms_eV_A": np.sqrt(np.mean(force_delta ** 2, axis=0)).tolist(),
        "species_force_vector_rms_eV_A": {},
    }
    symbols = np.array(atoms0.get_chemical_symbols())
    for symbol in sorted(set(symbols.tolist())):
        mask = symbols == symbol
        local = force_delta[mask]
        metrics["species_force_vector_rms_eV_A"][symbol] = {
            "atom_count": int(mask.sum()),
            "force_vector_rms_eV_A": float(np.sqrt(np.mean(np.sum(local ** 2, axis=1)))),
            "max_atom_force_vector_difference_eV_A": float(np.linalg.norm(local, axis=1).max()),
        }

    checks = {
        "native_energy": abs(metrics["native_energy_delta_meV_atom"])
            <= THRESHOLDS["absolute_native_energy_delta_meV_atom"],
        "free_energy": abs(metrics["free_energy_delta_meV_atom"])
            <= THRESHOLDS["absolute_free_energy_delta_meV_atom"],
        "force_vector": metrics["force_vector_rms_eV_A"]
            <= THRESHOLDS["force_vector_rms_eV_A"],
    }
    report = {
        "task": "v19_Ti3SiC2_kz_mesh_increment_review",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "comparison": "same frozen 48-atom Ti3SiC2 baseline; k6x6x2 versus k6x6x4",
        "delta_direction": "k6x6x4 minus k6x6x2",
        "geometry_checks": {
            "symbols_exact": True, "positions_exact": True, "cell_exact": True,
            "pbc_exact": True, "atom_ids_exact": True, "source_and_input_hash_equal": True,
            "both_archives_verified": True,
        },
        "meshes": [kpts0, kpts1],
        "method_equal_except_kpts": True,
        "source_input_sha256": summary0["source_sha256"],
        "archive_sha256_manifests_verified": archive_checks,
        "scf": {
            label: {"iterations": summary["scf_iterations"],
                    "elapsed_s": summary["elapsed_s"], "converged": summary["scf_converged"]}
            for label, _, summary, _, _ in rows
        },
        "metrics": metrics,
        "thresholds": THRESHOLDS,
        "checks": checks,
        "interpretation": (
            "At this fixed geometry and fixed 6x6 in-plane mesh, increasing kz from 2 to 4 "
            "passes the stated native-energy, free-energy, and force-vector RMS budgets. "
            "This is a numerical mesh control only; it does not establish general k-mesh "
            "convergence, material stability, or model accuracy. The registered 8x8x2 result "
            "is still needed to isolate in-plane sensitivity and complete the 2x2 mesh comparison."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
