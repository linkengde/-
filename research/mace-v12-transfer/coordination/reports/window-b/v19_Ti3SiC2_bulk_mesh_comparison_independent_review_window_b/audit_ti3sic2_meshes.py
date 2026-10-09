#!/usr/bin/env python3
"""Independently audit published Ti3SiC2 k-mesh archives without DFT."""
import argparse
import csv
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np
from ase.io import read


HERE = Path(__file__).resolve().parent
REPO = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], cwd=HERE, text=True).strip())
A_REPORT_REL = Path(
    "research/mace-v12-transfer/coordination/reports/window-a/"
    "v19_Ti3SiC2_bulk_mesh_comparison")
PURE_CONTROLS_REL = Path(
    "research/mace-v12-transfer/periodic_interface_v4/"
    "pbe_interface_v19_pure_phase_controls")
LABELS = ["Ti3SiC2_baseline_k6x6x2_v19", "Ti3SiC2_baseline_k8x8x4_v19"]
CONTROL_LABEL = "Ti3SiC2_baseline_k8x8x2_v19"
TOL = 2e-12


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_sha256sums(folder: Path) -> list[str]:
    manifest = folder / "SHA256SUMS.txt"
    checked = []
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        expected, filename = line.split(maxsplit=1)
        filename = filename.lstrip(" *")
        path = folder / filename
        actual = sha256(path)
        assert actual == expected, f"SHA256 mismatch: {path}"
        checked.append(filename)
    return checked


def close(actual, expected, name):
    assert math.isclose(float(actual), float(expected), rel_tol=0.0, abs_tol=TOL), (
        name, actual, expected)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    out_names = ["audit.json", "audit.csv", "verification.json", "SHA256SUMS.txt"]
    for name in out_names:
        assert not (output / name).exists(), f"Refusing to overwrite {output / name}"
    if output != HERE:
        assert not (output / "audit_ti3sic2_meshes.py").exists(), (
            f"Refusing to overwrite {output / 'audit_ti3sic2_meshes.py'}")

    report_dir = REPO / A_REPORT_REL
    checked_report_files = verify_sha256sums(report_dir)
    ref_path = report_dir / "k6_k8x8x4_comparison.json"
    ref = json.loads(ref_path.read_text())
    assert ref["labels"] == LABELS

    calc_root = REPO / PURE_CONTROLS_REL / "calculations"
    rows = []
    checked_archives = {}
    for label in LABELS:
        folder = calc_root / label
        checked_archives[label] = verify_sha256sums(folder)
        summary_path = folder / "summary.json"
        verification_path = folder / "verification.json"
        extxyz_path = folder / f"{label}_PW_PBE.extxyz"
        summary = json.loads(summary_path.read_text())
        verification = json.loads(verification_path.read_text())
        atoms = read(extxyz_path, format="extxyz")
        assert verification["status"] == "PASS"
        assert all(verification["checks"].values())
        assert summary["scf_converged"] is True
        assert ref["source_sha256"][label][extxyz_path.name] == sha256(extxyz_path)
        assert ref["source_sha256"][label][summary_path.name] == sha256(summary_path)
        assert ref["source_sha256"][label][verification_path.name] == sha256(verification_path)
        assert verification["input_sha256"] == summary["source_sha256"]
        rows.append((label, folder, summary, verification, atoms))

    (label0, folder0, summary0, verification0, atoms0), \
    (label1, folder1, summary1, verification1, atoms1) = rows
    assert ref["geometry_checks"]["both_archives_verified"] is True
    assert len(atoms0) == len(atoms1) == 48
    id_key = "local_new_id"
    assert id_key in atoms0.arrays and id_key in atoms1.arrays
    assert np.array_equal(atoms0.numbers, atoms1.numbers)
    assert np.array_equal(atoms0.positions, atoms1.positions)
    assert np.array_equal(atoms0.cell.array, atoms1.cell.array)
    assert np.array_equal(atoms0.pbc, atoms1.pbc)
    assert np.array_equal(atoms0.arrays[id_key], atoms1.arrays[id_key])

    method0 = copy.deepcopy(summary0["method"])
    method1 = copy.deepcopy(summary1["method"])
    kpts0 = method0.pop("kpts")
    kpts1 = method1.pop("kpts")
    assert method0 == method1
    assert kpts0 == [6, 6, 2] and kpts1 == [8, 8, 4]
    assert verification0["input_sha256"] == verification1["input_sha256"]
    assert summary0["source_sha256"] == summary1["source_sha256"]

    n_atoms = len(atoms0)
    delta_native = (summary1["energy_eV_cell"] - summary0["energy_eV_cell"]) * 1000 / n_atoms
    delta_free = (summary1["free_energy_eV_cell"] - summary0["free_energy_eV_cell"]) * 1000 / n_atoms
    force_delta = atoms1.arrays["PW_PBE_forces"] - atoms0.arrays["PW_PBE_forces"]
    force_norm = np.linalg.norm(force_delta, axis=1)
    all_metrics = {
        "native_energy_delta_meV_atom": float(delta_native),
        "free_energy_delta_meV_atom": float(delta_free),
        "force_vector_rms_eV_A": float(np.sqrt(np.mean(np.sum(force_delta ** 2, axis=1)))),
        "maximum_atom_force_vector_difference_eV_A": float(force_norm.max()),
        "force_component_rms_eV_A": np.sqrt(np.mean(force_delta ** 2, axis=0)).tolist(),
    }
    species_metrics = {}
    symbols = np.array(atoms0.get_chemical_symbols())
    for symbol in sorted(set(symbols.tolist())):
        mask = symbols == symbol
        species_delta = force_delta[mask]
        species_norm = np.linalg.norm(species_delta, axis=1)
        species_metrics[symbol] = {
            "atom_count": int(mask.sum()),
            "force_vector_rms_eV_A": float(np.sqrt(np.mean(np.sum(species_delta ** 2, axis=1)))),
            "force_component_rms_eV_A": np.sqrt(np.mean(species_delta ** 2, axis=0)).tolist(),
            "max_atom_force_vector_difference_eV_A": float(species_norm.max()),
        }

    for key, value in all_metrics.items():
        expected = ref["metrics"][key]
        if isinstance(value, list):
            assert np.allclose(value, expected, rtol=0.0, atol=TOL), (key, value, expected)
        else:
            close(value, expected, key)
    for symbol, metrics in species_metrics.items():
        expected = ref["metrics"]["species_force_vector_rms_eV_A"][symbol]
        for key, value in metrics.items():
            target = expected[key]
            if isinstance(value, list):
                assert np.allclose(value, target, rtol=0.0, atol=TOL), (symbol, key, value, target)
            elif isinstance(value, (int, float)):
                close(value, target, f"{symbol}.{key}")

    thresholds = ref["thresholds"]
    gates = {
        "native_energy": abs(delta_native) <= thresholds["absolute_native_energy_delta_meV_atom"],
        "free_energy": abs(delta_free) <= thresholds["absolute_free_energy_delta_meV_atom"],
        "force_vector": all_metrics["force_vector_rms_eV_A"] <= thresholds["force_vector_rms_eV_A"],
    }
    assert gates == ref["checks"]

    controls_root = REPO / PURE_CONTROLS_REL
    input_manifest_path = controls_root / "input_manifest.json"
    input_manifest = json.loads(input_manifest_path.read_text())
    manifest_rows = {r["label"]: r for r in input_manifest["records"]
                     if r.get("label") in [LABELS[0], CONTROL_LABEL]}
    assert set(manifest_rows) == {LABELS[0], CONTROL_LABEL}
    low_record = manifest_rows[LABELS[0]]
    same_kz_record = manifest_rows[CONTROL_LABEL]
    low_method = copy.deepcopy(low_record["method"])
    control_method = copy.deepcopy(same_kz_record["method"])
    low_kpts = low_method.pop("kpts")
    control_kpts = control_method.pop("kpts")
    same_kz_registration = {
        "label": CONTROL_LABEL,
        "input_sha256": same_kz_record["input_sha256"],
        "matches_k6_input_sha256": same_kz_record["input_sha256"] == low_record["input_sha256"],
        "method_equal_to_k6_except_kpts": control_method == low_method,
        "kpts": control_kpts,
        "kz_same_as_k6": control_kpts[2] == low_kpts[2],
        "launch_enabled": same_kz_record["launch_enabled"],
        "owner_task": same_kz_record["owner_task"],
    }
    assert same_kz_registration["matches_k6_input_sha256"]
    assert same_kz_registration["method_equal_to_k6_except_kpts"]
    assert same_kz_registration["kz_same_as_k6"]
    task_a = json.loads((REPO / "research/mace-v12-transfer/coordination/tasks/window-a.json").read_text())
    same_kz_registration["label_status_snapshot"] = task_a.get("labels", {}).get(CONTROL_LABEL)
    same_kz_registration["progress_snapshot"] = task_a.get("current_jobs", {}).get(CONTROL_LABEL)

    source_hashes = {
        label: {
            "extxyz_sha256": sha256(folder / f"{label}_PW_PBE.extxyz"),
            "summary_sha256": sha256(folder / "summary.json"),
            "verification_sha256": sha256(folder / "verification.json"),
            "input_sha256": verification["input_sha256"],
            "source_sha256": summary["source_sha256"],
        }
        for label, folder, summary, verification, _ in rows
    }
    audit = {
        "task": "v19_Ti3SiC2_bulk_mesh_comparison_independent_review_window_b",
        "source_report": str(A_REPORT_REL / "k6_k8x8x4_comparison.json"),
        "source_report_sha256": sha256(ref_path),
        "source_report_sha256_manifest_files_verified": checked_report_files,
        "source_archives": source_hashes,
        "archive_sha256_manifests_verified": checked_archives,
        "pair": {
            "labels": LABELS,
            "atom_count": n_atoms,
            "formula": atoms0.get_chemical_formula(),
            "identity_field": id_key,
            "positions_exact": True,
            "cell_exact": True,
            "pbc": atoms0.pbc.tolist(),
            "atom_ids_exact": True,
            "source_and_input_hash_equal": True,
            "method_equal_except_kpts": True,
            "kpts": [kpts0, kpts1],
        },
        "delta_direction": "k8x8x4 minus k6x6x2",
        "metrics": all_metrics,
        "species_metrics": species_metrics,
        "thresholds": thresholds,
        "gates": gates,
        "registered_same_kz_control": same_kz_registration,
        "interpretation": {
            "direct_pair_conclusion": "Energy gates pass; all-atom force-vector RMS gate fails.",
            "confound": "The direct pair changes in-plane mesh and kz together.",
            "same_kz_control": "If archived and verified with the registered identical input/method, k6x6x2 versus k8x8x2 isolates the in-plane increment at kz=2. k8x8x2 versus k8x8x4 then isolates the z-mesh increment at in-plane 8x8. This path does not measure a full 2x2 mesh interaction (no k6x6x4 point) and the current control has no result at this audit snapshot.",
        },
        "no_calculation_performed": True,
    }

    csv_path = output / "audit.csv"
    columns = ["row_type", "species", "atom_count", "native_energy_delta_meV_atom",
               "free_energy_delta_meV_atom", "force_vector_rms_eV_A",
               "max_force_vector_difference_eV_A", "force_component_rms_x_eV_A",
               "force_component_rms_y_eV_A", "force_component_rms_z_eV_A"]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerow({
            "row_type": "all", "species": "ALL", "atom_count": n_atoms,
            "native_energy_delta_meV_atom": all_metrics["native_energy_delta_meV_atom"],
            "free_energy_delta_meV_atom": all_metrics["free_energy_delta_meV_atom"],
            "force_vector_rms_eV_A": all_metrics["force_vector_rms_eV_A"],
            "max_force_vector_difference_eV_A": all_metrics["maximum_atom_force_vector_difference_eV_A"],
            "force_component_rms_x_eV_A": all_metrics["force_component_rms_eV_A"][0],
            "force_component_rms_y_eV_A": all_metrics["force_component_rms_eV_A"][1],
            "force_component_rms_z_eV_A": all_metrics["force_component_rms_eV_A"][2],
        })
        for symbol, metrics in species_metrics.items():
            writer.writerow({
                "row_type": "species", "species": symbol, "atom_count": metrics["atom_count"],
                "force_vector_rms_eV_A": metrics["force_vector_rms_eV_A"],
                "max_force_vector_difference_eV_A": metrics["max_atom_force_vector_difference_eV_A"],
                "force_component_rms_x_eV_A": metrics["force_component_rms_eV_A"][0],
                "force_component_rms_y_eV_A": metrics["force_component_rms_eV_A"][1],
                "force_component_rms_z_eV_A": metrics["force_component_rms_eV_A"][2],
            })

    (output / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
    verification = {
        "status": "PASS",
        "checks": {
            "A_report_sha_manifest": True,
            "both_archive_sha_manifests": True,
            "both_archive_verifications_pass": True,
            "exact_geometry_and_ids": True,
            "method_differs_only_by_kpts": True,
            "A_report_arithmetic_matches": True,
            "energy_and_force_gate_recomputation": True,
            "same_kz_control_registration_matches_input_and_method": True,
            "no_DFT_or_MPI": True,
        },
        "A_report_files_checked": checked_report_files,
        "archive_files_checked": checked_archives,
        "metrics": all_metrics,
        "gates": gates,
    }
    (output / "verification.json").write_text(json.dumps(verification, ensure_ascii=False, indent=2) + "\n")

    if output == HERE:
        checksum_names = ["README.md", "audit_ti3sic2_meshes.py", "audit.json", "audit.csv", "verification.json"]
    else:
        shutil.copyfile(Path(__file__).resolve(), output / "audit_ti3sic2_meshes.py")
        checksum_names = ["audit_ti3sic2_meshes.py", "audit.json", "audit.csv", "verification.json"]
    (output / "SHA256SUMS.txt").write_text("".join(
        f"{sha256(output / name)}  {name}\n" for name in checksum_names))
    print(json.dumps({"status": "PASS", "metrics": all_metrics, "gates": gates,
                      "output_dir": str(output)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
