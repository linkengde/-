#!/usr/bin/env python3
"""Future compact archive verifier; withheld-test checks are integrity-only."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

import numpy as np
from ase.io import read


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--workspace-root", required=True, type=Path)
    p.add_argument("--pack-root", required=True, type=Path)
    p.add_argument("--input-root", required=True, type=Path)
    p.add_argument("--run-root", required=True, type=Path)
    p.add_argument("--archive-root", required=True, type=Path)
    p.add_argument("--authorization", required=True, type=Path)
    p.add_argument("--owner", required=True, choices=("window-a", "window-b"))
    p.add_argument("--label", required=True)
    args = p.parse_args()

    workspace = args.workspace_root.resolve()
    pack = args.pack_root.resolve()
    input_root = args.input_root.resolve()
    run_root = args.run_root.resolve()
    archive_root = args.archive_root.resolve()
    auth = json.loads(args.authorization.resolve().read_text())
    manifest = json.loads((pack / "input_manifest.json").read_text())
    rec = next((row for row in manifest["records"] if row["label"] == args.label), None)
    if rec is None:
        raise SystemExit(f"Unknown label: {args.label}")
    if auth.get("status") != "authorized" or auth.get("reviewed_by") != "window-a":
        raise SystemExit("A authorization is required for archive verification")
    if auth.get("owner_assignments", {}).get(args.label) != args.owner:
        raise SystemExit("Archive owner does not match A assignment")
    if auth.get("role_assignments", {}).get(args.label) != rec["role"]:
        raise SystemExit("Archive role does not match A role assignment")
    if auth.get("role_manifest_sha256") != manifest["source_geometry_bundle"]["role_manifest_sha256"]:
        raise SystemExit("Role-manifest hash mismatch")
    if auth.get("input_sha256", {}).get(args.label) != rec["input_sha256"]:
        raise SystemExit("Authorization input hash mismatch")
    if rec["role"] == "withheld_test":
        frozen = auth.get("frozen_selection", {})
        if frozen.get("frozen_before_test_label_read") is not True:
            raise SystemExit("Frozen model selection record required before withheld-test label access")

    source = input_root / Path(rec["input"]).name
    run = run_root / args.label
    output = run / f"{args.label}_PW_PBE.extxyz"
    summary_path = run / "summary.json"
    progress_path = run / "progress.json"
    gpaw_log = run / "gpaw.log"
    launcher_log = run_root / f"{args.label}.launcher.log"
    needed = (source, output, summary_path, progress_path, gpaw_log, launcher_log)
    missing = [str(path) for path in needed if not path.is_file()]
    if missing:
        raise SystemExit("Required archive source is missing: " + ", ".join(missing))
    if digest(source) != rec["input_sha256"]:
        raise SystemExit("Source geometry SHA256 mismatch")
    original = read(source, format="extxyz")
    atoms = read(output, format="extxyz")
    summary = json.loads(summary_path.read_text())
    progress = json.loads(progress_path.read_text())
    energy = float(atoms.info["PW_PBE_energy_eV"])
    forces = np.asarray(atoms.arrays["PW_PBE_forces"], dtype=float)
    free_energy = summary.get("free_energy_eV_cell")
    log_text = gpaw_log.read_text(errors="replace")
    launcher_text = launcher_log.read_text(errors="replace")
    log_converged = re.search(r"Converged in\s+\d+\s+steps", log_text) is not None
    checks = {
        "source_hash": digest(source) == rec["input_sha256"] == summary.get("source_sha256"),
        "authorization_hash": auth["input_sha256"][args.label] == rec["input_sha256"],
        "role_owner": summary.get("role") == rec["role"] and summary.get("owner") == args.owner,
        "scf_convergence": summary.get("scf_converged") is True and progress.get("status") == "complete" and log_converged,
        "mpi_ranks": summary.get("mpi_ranks") == 4 and progress.get("mpi_ranks") == 4,
        "atom_count": len(atoms) == len(original) == rec["atom_count"] == summary.get("atom_count"),
        "formula": atoms.get_chemical_formula() == original.get_chemical_formula() == rec["formula"],
        "element_order": atoms.get_chemical_symbols() == original.get_chemical_symbols() == rec["symbols_in_file_order"],
        "atom_id_order": np.array_equal(atoms.arrays["lammps_id"], original.arrays["lammps_id"]) and np.asarray(atoms.arrays["lammps_id"], dtype=int).tolist() == rec["atom_ids_in_file_order"],
        "positions": np.allclose(atoms.positions, original.positions, atol=5e-8, rtol=0),
        "cell": np.allclose(atoms.cell.array, original.cell.array, atol=1e-10, rtol=0),
        "pbc": np.array_equal(atoms.pbc, original.pbc) and np.asarray(atoms.pbc, dtype=bool).tolist() == rec["pbc"],
        "finite_energy_forces": np.isfinite(energy) and forces.shape == (len(atoms), 3) and np.isfinite(forces).all() and summary.get("finite") is True,
        "summary_energy_matches_extxyz": np.isclose(energy, float(summary["energy_eV_cell"]), atol=1e-10, rtol=0),
        "finite_free_energy_when_supported": (free_energy is None and summary.get("free_energy_supported") is False) or (free_energy is not None and np.isfinite(float(free_energy)) and summary.get("free_energy_supported") is True),
        "method": summary.get("method", {}).get("GPAW_version") == "26.7.0" and summary.get("method", {}).get("cutoff_eV") == 500 and summary.get("method", {}).get("kpts") == [1, 1, 1] and summary.get("method", {}).get("smearing_eV") == 0.1,
        "log_evidence": log_converged and len(launcher_text) > 0,
        "no_state_gpw_in_archive_plan": True,
    }
    checks = {key: bool(value) for key, value in checks.items()}
    if not all(checks.values()):
        raise SystemExit("Archive validation failed: " + json.dumps(checks, sort_keys=True))

    target = archive_root / args.label
    if target.exists():
        raise SystemExit(f"Refusing to overwrite existing archive: {target}")
    archive_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{args.label}.staging-", dir=archive_root))
    try:
        copies = {
            output.name: output,
            "summary.json": summary_path,
            "progress.json": progress_path,
            "gpaw.log": gpaw_log,
            "launcher.log": launcher_log,
        }
        for name, source_path in copies.items():
            shutil.copy2(source_path, staging / name)
        checksums = {name: digest(staging / name) for name in sorted(copies)}
        verification = {
            "status": "PASS",
            "label": args.label,
            "role": rec["role"],
            "owner": args.owner,
            "checks": checks,
            "input_sha256": rec["input_sha256"],
            "compact_file_sha256": checksums,
            "state_gpw_archived": False,
            "withheld_test_handling": "integrity-only; do not model-score or interpret; A controls label unsealing",
        }
        (staging / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
        checksum_path = staging / "SHA256SUMS.txt"
        checksum_path.write_text("".join(f"{checksums[name]}  {name}\n" for name in sorted(checksums)))
        staging.rename(target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    # Do not print energy, forces, or any scored metric, especially for withheld-test roles.
    print(f"Archive verification PASS: {args.label}; role={rec['role']}; values not displayed.")


if __name__ == "__main__":
    main()
