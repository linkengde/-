#!/usr/bin/env python3
"""Guarded Ag-Si phase-2 runner. Default invocation is read-only preflight."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import traceback
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REGISTRY_ROOT = HERE.parent
REPO_ROOT = HERE.parents[6]
LOCK_PATH = HERE / "PHASE2_PARAMETER_LOCK.json"
AUTH_REQUEST_PATH = HERE / "PHASE2_EXEC_AUTH_REQUEST.json"
LOCK = json.loads(LOCK_PATH.read_text(encoding="utf-8"))
REGISTRY_IDS = tuple(LOCK["allowed_structure_ids"])
FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stop(message: str) -> None:
    raise SystemExit("HOLD: " + message)


def bound_files() -> dict[str, Any]:
    entries: dict[str, Any] = {}
    for registry_id in REGISTRY_IDS:
        spec = LOCK["structures"][registry_id]
        input_path = REPO_ROOT / spec["input_path"]
        structure_path = REPO_ROOT / spec["structure_path"]
        actual_input_sha = sha256(input_path) if input_path.is_file() else None
        actual_structure_sha = sha256(structure_path) if structure_path.is_file() else None
        entries[registry_id] = {
            "input_path": str(input_path),
            "input_sha256": actual_input_sha,
            "input_sha_match": actual_input_sha == spec["input_sha256"],
            "structure_path": str(structure_path),
            "structure_sha256": actual_structure_sha,
            "structure_sha_match": actual_structure_sha == spec["structure_sha256"],
        }
    return entries


def requested_hash_bindings() -> dict[str, Any]:
    return {
        "parameter_lock_sha256": sha256(LOCK_PATH),
        "runner_sha256": sha256(Path(__file__).resolve()),
        "collector_sha256": sha256(HERE / "collect_registry_phase2_results.py"),
        "structures": {
            registry_id: {
                "input_sha256": LOCK["structures"][registry_id]["input_sha256"],
                "structure_sha256": LOCK["structures"][registry_id]["structure_sha256"],
            }
            for registry_id in REGISTRY_IDS
        },
    }


def preflight() -> dict[str, Any]:
    observed = bound_files()
    hashes_ok = all(
        item["input_sha_match"] and item["structure_sha_match"]
        for item in observed.values()
    )
    request = json.loads(AUTH_REQUEST_PATH.read_text(encoding="utf-8")) if AUTH_REQUEST_PATH.exists() else {}
    return {
        "status": "PREPARED_NOT_ISSUED",
        "execution_authorized": False,
        "dft_started": False,
        "registry_ids": list(REGISTRY_IDS),
        "kpoint_mesh": LOCK["kpoint_mesh"],
        "fixed_geometry": LOCK["geometry"]["fixed_geometry"],
        "relax_atoms": LOCK["geometry"]["relax_atoms"],
        "relax_cell": LOCK["geometry"]["relax_cell"],
        "mpi_ranks": LOCK["resources"]["mpi_ranks"],
        "omp_threads": LOCK["resources"]["omp_threads"],
        "role": LOCK["role"],
        "source_hashes_match": hashes_ok,
        "source_bindings": observed,
        "requested_bindings": requested_hash_bindings(),
        "current_auth_status": request.get("authorization_status", "PREPARED_NOT_ISSUED"),
    }


def verify_auth(path: Path, expected_sha: str) -> tuple[dict[str, Any], str]:
    actual_sha = sha256(path)
    if actual_sha.lower() != expected_sha.lower():
        stop("authorization SHA-256 mismatch")
    auth = json.loads(path.read_text(encoding="utf-8"))
    if auth.get("authorization_status") != "ISSUED" or auth.get("execution_authorized") is not True:
        stop("C authorization is not ISSUED")
    approval = auth.get("C_approval", {})
    if approval.get("signed_by") != "C" or not approval.get("approval_id"):
        stop("C approval signature/ID is absent")
    if auth.get("allowed_structure_ids") != list(REGISTRY_IDS):
        stop("authorization registry allowlist mismatch")
    bindings = auth.get("bindings", {})
    current = requested_hash_bindings()
    if bindings.get("parameter_lock_sha256") != current["parameter_lock_sha256"]:
        stop("parameter-lock SHA binding mismatch")
    if bindings.get("runner_sha256") != current["runner_sha256"]:
        stop("runner SHA binding mismatch")
    if bindings.get("collector_sha256") != current["collector_sha256"]:
        stop("collector SHA binding mismatch")
    if bindings.get("structures") != current["structures"]:
        stop("R1/R2 input or structure SHA binding mismatch")
    if auth.get("kpoint_mesh") != LOCK["kpoint_mesh"] or auth.get("role") != LOCK["role"]:
        stop("scientific scope mismatch")
    return auth, actual_sha


def parse_scf_log(path: Path) -> tuple[bool | None, int | None]:
    text = path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""
    converged = list(re.finditer(r"Converged in\s+(\d+)\s+steps", text, re.I))
    if converged:
        return True, int(converged[-1].group(1))
    iterations = [int(item) for item in re.findall(r"\|iter:\s*(\d+)\|", text)]
    return (False, max(iterations)) if iterations else (None, None)


def check_runtime_environment() -> tuple[Any, Any, Any, Any]:
    import ase
    import gpaw
    import gpaw_data
    from gpaw.mpi import world
    from importlib.metadata import version

    if ase.__version__ != LOCK["engine_versions"]["ase"]:
        stop("ASE version mismatch")
    if gpaw.__version__ != LOCK["engine_versions"]["gpaw"]:
        stop("GPAW version mismatch")
    if version("gpaw-data") != LOCK["engine_versions"]["gpaw_data"]:
        stop("gpaw-data version mismatch")
    if world.size != LOCK["resources"]["mpi_ranks"]:
        stop("MPI world size must be exactly four")
    if os.environ.get("OMP_NUM_THREADS") != "1":
        stop("OMP_NUM_THREADS must be 1")
    setup_path = os.environ.get("GPAW_SETUP_PATH")
    if not setup_path or Path(setup_path).resolve() != Path(gpaw_data.datapath()).resolve():
        stop("GPAW_SETUP_PATH does not match gpaw_data.datapath()")
    return ase, gpaw, gpaw_data, world


def verify_paw_files() -> dict[str, dict[str, str]]:
    from gpaw.setup_data import search_for_file

    observed: dict[str, dict[str, str]] = {}
    for element, expected_sha in LOCK["paw_sha256_by_element"].items():
        filename, _ = search_for_file(f"{element}.PBE")
        actual_sha = sha256(Path(filename))
        if actual_sha != expected_sha:
            stop(f"PAW SHA mismatch for {element}")
        observed[element] = {"path": str(filename), "sha256": actual_sha}
    return observed


def run_registry(registry_id: str, work_root: Path, auth: dict[str, Any], auth_sha: str, paw_files: dict[str, Any], ase: Any, gpaw: Any, world: Any) -> dict[str, Any]:
    from ase.io import read
    from gpaw import FermiDirac, GPAW, PW
    from gpaw.mixer import Mixer

    spec = LOCK["structures"][registry_id]
    structure_path = REPO_ROOT / spec["structure_path"]
    input_path = REPO_ROOT / spec["input_path"]
    case_dir = work_root / registry_id
    if world.rank == 0:
        case_dir.mkdir(parents=False, exist_ok=False)
    world.barrier()
    log_path = case_dir / "gpaw.log"
    result_path = case_dir / "result.json"
    atoms = read(structure_path, format="extxyz")
    before_positions = atoms.positions.copy()
    before_cell = atoms.cell.array.copy()
    before_pbc = atoms.pbc.copy()

    record = {
        "registry_id": registry_id,
        "authorization_sha256": auth_sha,
        "authorization_id": auth["C_approval"]["approval_id"],
        "parameter_lock_sha256": sha256(LOCK_PATH),
        "runner_sha256": sha256(Path(__file__).resolve()),
        "collector_sha256": sha256(HERE / "collect_registry_phase2_results.py"),
        "input_sha256": sha256(input_path),
        "structure_sha256": sha256(structure_path),
        "kpoint_mesh": LOCK["kpoint_mesh"],
        "role": LOCK["role"],
        "mpi_ranks": LOCK["resources"]["mpi_ranks"],
        "omp_threads": LOCK["resources"]["omp_threads"],
        "fixed_geometry": True,
        "relax_atoms": False,
        "relax_cell": False,
        "C_approval": auth["C_approval"],
    }
    if world.rank == 0:
        (case_dir / "C_EXECUTION_RELEASE_RECORD.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    world.barrier()

    labels: dict[str, Any] = {
        "total_energy_eV": None,
        "force_consistent_free_energy_eV": None,
        "forces_eV_A": None,
        "stress_eV_A3": None,
        "stress_status": "NOT_AVAILABLE",
    }
    run_error = None
    try:
        method = LOCK["method"]
        calc = GPAW(
            mode=PW(method["cutoff_eV"]),
            xc=method["xc"],
            kpts=tuple(LOCK["kpoint_mesh"]),
            occupations=FermiDirac(method["occupations"]["width_eV"]),
            convergence=method["convergence"],
            maxiter=method["max_scf_iterations"],
            mixer=Mixer(method["mixer"]["beta"], method["mixer"]["nmaxold"], method["mixer"]["weight"]),
            symmetry={"point_group": True, "time_reversal": True},
            poissonsolver={},
            parallel={"sl_auto": True},
            txt=str(log_path),
        )
        atoms.calc = calc
        labels["total_energy_eV"] = float(atoms.get_potential_energy(force_consistent=False))
        labels["force_consistent_free_energy_eV"] = float(atoms.get_potential_energy(force_consistent=True))
        labels["forces_eV_A"] = atoms.get_forces(apply_constraint=False).tolist()
        labels["stress_eV_A3"] = atoms.get_stress(voigt=False, apply_constraint=False).tolist()
        labels["stress_status"] = "AVAILABLE"
    except Exception as exc:  # retain log and continue with the other registry
        run_error = f"{type(exc).__name__}: {exc}"
        if world.rank == 0:
            (case_dir / "stderr.log").write_text(traceback.format_exc(), encoding="utf-8")
    world.barrier()
    scf_converged, scf_steps = parse_scf_log(log_path)
    geometry_unchanged = (
        (atoms.positions == before_positions).all()
        and (atoms.cell.array == before_cell).all()
        and (atoms.pbc == before_pbc).all()
    )
    if not geometry_unchanged:
        stop(f"geometry changed unexpectedly for {registry_id}")
    result = {
        "structure_id": registry_id,
        "input_sha256": record["input_sha256"],
        "structure_sha256": record["structure_sha256"],
        "configuration_sha256": record["parameter_lock_sha256"],
        "kpoint_mesh": LOCK["kpoint_mesh"],
        "total_energy_eV": labels["total_energy_eV"],
        "force_consistent_free_energy_eV": labels["force_consistent_free_energy_eV"],
        "forces_eV_A": labels["forces_eV_A"],
        "forces_status": "AVAILABLE" if labels["forces_eV_A"] is not None else "NOT_AVAILABLE",
        "stress_eV_A3": labels["stress_eV_A3"],
        "stress_status": labels["stress_status"],
        "scf_converged": scf_converged,
        "scf_steps": scf_steps,
        "run_status": "SUCCESS" if scf_converged is True and run_error is None else "FAILED_SCF_OR_RUNTIME",
        "failure_reason": run_error,
        "geometry_mutated": False,
        "role": LOCK["role"],
        "paw_files": paw_files,
    }
    if world.rank == 0:
        result_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    world.barrier()
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--preflight", action="store_true", help="read-only binding check (default)")
    mode.add_argument("--execute", action="store_true", help="run only with a C-issued, SHA-bound authorization")
    parser.add_argument("--authorization", type=Path)
    parser.add_argument("--authorization-sha256")
    parser.add_argument("--work-root", type=Path)
    args = parser.parse_args()
    if not args.execute:
        report = preflight()
        print(json.dumps(report, indent=2))
        return 0 if report["source_hashes_match"] else 2

    if args.authorization is None or args.authorization_sha256 is None or args.work_root is None:
        stop("--execute requires --authorization, --authorization-sha256, and --work-root")
    auth, auth_sha = verify_auth(args.authorization, args.authorization_sha256)
    observed = bound_files()
    if not all(item["input_sha_match"] and item["structure_sha_match"] for item in observed.values()):
        stop("source input/structure SHA check failed")
    if not args.work_root.is_dir():
        stop("work-root must already exist as an isolated empty directory")
    if any(args.work_root.iterdir()):
        stop("work-root must be empty; refusing to overwrite existing data")
    ase, gpaw, _gpaw_data, world = check_runtime_environment()
    paw_files = verify_paw_files()
    for registry_id in REGISTRY_IDS:
        run_registry(registry_id, args.work_root, auth, auth_sha, paw_files, ase, gpaw, world)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
