#!/usr/bin/env python3
"""Guarded continue-on-error scheduler for the authorized R1-R4 6x6x1 batch.

Default mode is dry-run. DFT starts only with --execute after all frozen source,
authorization, environment, and MPI checks pass. One registry failure never
prevents subsequent registries from being attempted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from importlib.metadata import version
from pathlib import Path
from typing import Callable

AUTH_SHA256 = "86b5058f12745e1484680bfc926659b3a58fe426a873076514ef614467fc6b9f"
CONFIG_SHA256 = "515cd81f16e491d8d521426ba63f17946cb47dd5fa0361fa9a8d946e07f0d035"
REGISTRIES = ("R1", "R2", "R3", "R4")
ROOT = Path(__file__).resolve().parent


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def write_json(path: Path, data: object) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_authority(auth_path: Path, config_path: Path) -> tuple[dict, dict]:
    auth_raw = auth_path.read_bytes()
    config_raw = config_path.read_bytes()
    if sha256_bytes(auth_raw) != AUTH_SHA256:
        raise RuntimeError("authorization SHA-256 mismatch")
    if sha256_bytes(config_raw) != CONFIG_SHA256:
        raise RuntimeError("configuration SHA-256 mismatch")
    auth = json.loads(auth_raw)
    config = json.loads(config_raw)
    if auth.get("execution_authorized") is not True:
        raise RuntimeError("C authorization is not enabled")
    if auth.get("C_approval", {}).get("signed_by") != "C":
        raise RuntimeError("C approval identity mismatch")
    if auth.get("ROLE") != "REGISTRY_SCREENING_ONLY":
        raise RuntimeError("authorization role mismatch")
    if auth.get("CONFIG_SHA256") != CONFIG_SHA256:
        raise RuntimeError("authorization is bound to another configuration")
    if auth.get("scope", {}).get("authorized_registries") != list(REGISTRIES):
        raise RuntimeError("authorized registry set/order mismatch")
    if auth.get("kpoint_mesh") != [6, 6, 1]:
        raise RuntimeError("authorized k-point mesh mismatch")
    runtime = auth.get("runtime", {})
    if runtime.get("mpi_ranks") != 4 or runtime.get("omp_threads") != 1 or runtime.get("OMP_NUM_THREADS") != "1":
        raise RuntimeError("authorized MPI/OMP allocation mismatch")
    method = config.get("method", {})
    if (method.get("engine") != "GPAW 26.7.0" or method.get("xc") != "PBE"
            or method.get("pseudopotential") != "PAW" or method.get("cutoff_eV") != 500
            or method.get("smearing_eV") != 0.1):
        raise RuntimeError("frozen scientific parameters mismatch")
    phase1 = config.get("phase1_registry_screening", {})
    if phase1.get("kpoint_mesh") != [6, 6, 1]:
        raise RuntimeError("configuration k-point mesh mismatch")
    convention = phase1.get("kpoint_convention", "").lower()
    if not all(token in convention for token in ("monkhorst-pack", "no explicit gamma override or extra offset", "slab kz=1")):
        raise RuntimeError("configuration k-point convention mismatch")
    if config.get("geometry", {}).get("fixed") is not True:
        raise RuntimeError("fixed-geometry setting mismatch")
    if config.get("geometry", {}).get("relax_atoms") is not False or config.get("geometry", {}).get("relax_cell") is not False:
        raise RuntimeError("geometry relaxation is not disabled")
    if config.get("data_role", {}).get("role") != "REGISTRY_SCREENING_ONLY":
        raise RuntimeError("data role mismatch")
    return auth, config


def validate_candidate_files(work_root: Path, auth: dict, *, execute: bool) -> list[dict]:
    records = {r["registry_id"]: r for r in auth["registries"]}
    if set(records) != set(REGISTRIES):
        raise RuntimeError("authorization registry records mismatch")
    plans = []
    for registry_id in REGISTRIES:
        authorized = records[registry_id]
        directory = work_root / registry_id
        input_path = directory / "input.py"
        structure_path = directory / "structure.xyz"
        release_path = directory / "C_EXECUTION_RELEASE_RECORD.json"
        for path in (input_path, structure_path, release_path):
            if not path.is_file():
                raise RuntimeError(f"{registry_id}: required frozen file missing: {path.name}")
        input_sha = sha256_file(input_path)
        structure_sha = sha256_file(structure_path)
        release_raw = release_path.read_bytes()
        release_sha = sha256_bytes(release_raw)
        release = json.loads(release_raw)
        if input_sha != authorized["input_sha256"]:
            raise RuntimeError(f"{registry_id}: input.py SHA-256 mismatch")
        if structure_sha != authorized["structure_sha256"]:
            raise RuntimeError(f"{registry_id}: structure SHA-256 mismatch")
        runner_fields = authorized["runner_release_fields"]
        if runner_fields.get("kpoint_convention") != "GPAW default Monkhorst-Pack integer mesh; no explicit gamma override or extra offset":
            raise RuntimeError(f"{registry_id}: authorized k-point convention mismatch")
        for key, value in runner_fields.items():
            if release.get(key) != value:
                raise RuntimeError(f"{registry_id}: release field mismatch: {key}")
        if release.get("C_approval") != auth["C_approval"]:
            raise RuntimeError(f"{registry_id}: release C_approval mismatch")
        if release.get("authorization_sha256") != AUTH_SHA256:
            raise RuntimeError(f"{registry_id}: release authorization binding mismatch")
        if release.get("configuration_sha256") != CONFIG_SHA256:
            raise RuntimeError(f"{registry_id}: release configuration binding mismatch")
        output_dir = directory / "results" / "coarse"
        prior_artifacts = [directory / name for name in ("stdout.log", "stderr.log", "exit_code.txt", "SCF_STATUS.json", "RUN_STATUS.json")]
        occupied = output_dir.exists() or any(p.exists() for p in prior_artifacts)
        if execute and occupied:
            raise RuntimeError(f"{registry_id}: output already exists; use a fresh run package")
        plans.append({
            "registry_id": registry_id,
            "input_sha256": input_sha,
            "structure_sha256": structure_sha,
            "release_sha256": release_sha,
            "output_directory": str(output_dir),
            "execution_blocked_by_existing_output": occupied,
        })
    return plans


def validate_environment(gpaw_exe: Path, mpiexec: str, launcher_dir: Path, work_root: Path) -> dict:
    if not gpaw_exe.is_file() or not os.access(gpaw_exe, os.X_OK):
        raise RuntimeError(f"GPAW executable unavailable: {gpaw_exe}")
    if not shutil.which(mpiexec):
        raise RuntimeError(f"MPI launcher unavailable: {mpiexec}")
    if not (launcher_dir / "sitecustomize.py").is_file() or not (launcher_dir / "mpi_guard_smoke.py").is_file():
        raise RuntimeError("rank-0 directory guard helper is incomplete")
    import ase
    import gpaw
    if gpaw.__version__ != "26.7.0" or ase.__version__ != "3.29.0" or version("gpaw-data") != "1.2.1":
        raise RuntimeError("GPAW/ASE/PAW runtime versions do not match the frozen configuration")
    from gpaw.setup_data import search_for_file
    expected_paw = {
        "Ag": "bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97",
        "Ti": "bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f",
        "Si": "2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0",
        "C": "fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983",
    }
    paw_hashes = {}
    for symbol, expected in expected_paw.items():
        filename, _ = search_for_file(f"{symbol}.PBE")
        actual = sha256_file(Path(filename))
        if actual != expected:
            raise RuntimeError(f"PAW SHA-256 mismatch for {symbol}")
        paw_hashes[symbol] = actual
    cpu_quota = None
    cpu_max = Path("/sys/fs/cgroup/cpu.max")
    if cpu_max.exists():
        quota, period = cpu_max.read_text().split()
        if quota != "max":
            cpu_quota = int(quota) / int(period)
            if cpu_quota < 4:
                raise RuntimeError("CPU quota is below the authorized 4 MPI ranks")
    env = os.environ.copy()
    env.update({
        "GPAW_MPI_OUTPUT_CREATE_GUARD": "1",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "PYTHONPATH": str(launcher_dir) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""),
    })
    with tempfile.TemporaryDirectory(prefix="gpaw-mpi-dir-guard-") as scratch:
        env["GPAW_MPI_GUARD_TEST_TARGET"] = str(Path(scratch) / "results" / "coarse")
        smoke = subprocess.run(
            [mpiexec, "--oversubscribe", "-n", "4", str(gpaw_exe), "python", str(launcher_dir / "mpi_guard_smoke.py")],
            cwd=work_root, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        lines = [line for line in smoke.stdout.splitlines() if line.strip()]
        if smoke.returncode != 0 or len(lines) != 4:
            raise RuntimeError("four-rank MPI/output-directory smoke test failed: " + smoke.stderr[-2000:])
        ranks = sorted(json.loads(line)["rank"] for line in lines)
        sizes = {json.loads(line)["size"] for line in lines}
        if ranks != [0, 1, 2, 3] or sizes != {4}:
            raise RuntimeError("MPI smoke test did not report exactly ranks 0-3 of 4")
    disk = shutil.disk_usage(work_root)
    if disk.free < 5 * 1024**3:
        raise RuntimeError("less than 5 GiB free storage")
    return {
        "python": sys.version.split()[0], "gpaw": gpaw.__version__, "ase": ase.__version__,
        "gpaw_data": version("gpaw-data"), "mpi_smoke_ranks": 4,
        "omp_num_threads": "1", "cpu_quota_cores": cpu_quota,
        "free_disk_bytes": disk.free, "paw_sha256": paw_hashes,
    }


def failure_reason(stderr: str, stdout: str, scf_log: str, return_code: int | None) -> tuple[str, str]:
    combined = "\n".join((stderr, stdout, scf_log))
    # Only classify an SCF failure from an actual exception, not source-code
    # context printed in a traceback (e.g. the line `require(..., "SCF did not
    # converge")` can appear even when the real exception is AttributeError).
    if ("KohnShamConvergenceError" in combined
            or re.search(r"(?:SystemExit|RuntimeError|AssertionError):\s+(?:HOLD:\s+)?SCF did not converge", combined)
            or "Maximum number of iterations reached without convergence" in combined):
        return "FAILED_SCF", "SCF reached the frozen iteration limit without convergence; full SCF log and partial output were retained."
    if "AttributeError: 'ASECalculator' object has no attribute 'scf'" in combined:
        return "FAILED_POSTPROCESS", "GPAW single-point quantities were calculated, but the frozen input's convergence-status check accessed a missing ASECalculator.scf attribute; convergence is unverified and no result.json was written."
    if return_code is None:
        return "FAILED", "process could not be started"
    exception_lines = [line.replace("\x00", "").strip() for line in stderr.splitlines()
                       if re.search(r"(?:AttributeError|RuntimeError|SystemExit|ValueError|TypeError|ImportError):", line)]
    if exception_lines:
        return "FAILED", exception_lines[-1][-1000:]
    nonempty = [line.replace("\x00", "").strip() for line in stderr.splitlines() if line.replace("\x00", "").strip()]
    if not nonempty:
        nonempty = [line.replace("\x00", "").strip() for line in stdout.splitlines() if line.replace("\x00", "").strip()]
    detail = nonempty[-1][-1000:] if nonempty else f"process exited with code {return_code}"
    return "FAILED", detail


def run_one(registry_id: str, work_root: Path, *, gpaw_exe: Path, mpiexec: str, launcher_dir: Path) -> dict:
    directory = work_root / registry_id
    release_path = directory / "C_EXECUTION_RELEASE_RECORD.json"
    release_sha = sha256_file(release_path)
    command = [
        mpiexec, "--oversubscribe", "-n", "4", str(gpaw_exe), "python", "input.py", "--",
        "--run", "--mesh", "coarse", "--release-config", release_path.name,
        "--release-sha256", release_sha,
    ]
    env = os.environ.copy()
    env.update({
        "GPAW_MPI_OUTPUT_CREATE_GUARD": "1", "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
        "PYTHONPATH": str(launcher_dir) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""),
    })
    start = time.time()
    return_code: int | None = None
    launch_error = None
    with (directory / "stdout.log").open("wb") as stdout_file, (directory / "stderr.log").open("wb") as stderr_file:
        try:
            completed = subprocess.run(command, cwd=directory, env=env, stdout=stdout_file, stderr=stderr_file, check=False)
            return_code = completed.returncode
        except Exception as exc:  # Continue to the next registry after any launch failure.
            launch_error = f"{type(exc).__name__}: {exc}"
    (directory / "exit_code.txt").write_text((str(return_code) if return_code is not None else "NOT_STARTED") + "\n", encoding="utf-8")
    result_path = directory / "results" / "coarse" / "result.json"
    log_path = directory / "results" / "coarse" / "gpaw.log"
    stdout = (directory / "stdout.log").read_text(errors="replace")
    stderr = (directory / "stderr.log").read_text(errors="replace")
    scf_log = log_path.read_text(errors="replace") if log_path.exists() else ""
    result = None
    if return_code == 0 and result_path.is_file():
        try:
            result = json.loads(result_path.read_text(encoding="utf-8"))
        except Exception:
            result = None
    success = bool(result and result.get("scf_converged") is True and result.get("kpoint_mesh") == [6, 6, 1])
    if success:
        scf_status = "CONVERGED"
        record = {
            "registry_id": registry_id, "status": "SUCCESS", "scf_status": scf_status,
            "exit_code": return_code, "runtime_seconds": round(time.time() - start, 3),
            "total_energy_eV": result["total_energy_eV"],
            "force_consistent_free_energy_eV": result["free_energy_eV"],
            "forces_eV_A": result["forces_eV_A"], "stress_eV_A3": result["stress_eV_A3"],
            "result_path": str(result_path.relative_to(work_root)),
        }
    else:
        failure_status, reason = failure_reason(stderr, stdout, scf_log, return_code)
        if launch_error:
            reason = launch_error
        if failure_status == "FAILED_SCF":
            scf_status = "FAILED_SCF"
        elif failure_status == "FAILED_POSTPROCESS":
            scf_status = "CONVERGENCE_UNVERIFIED"
        else:
            scf_status = "NOT_REACHED" if not scf_log else "UNVERIFIED"
        partial = sorted(str(p.relative_to(directory)) for p in (directory / "results" / "coarse").glob("*") if p.is_file()) if (directory / "results" / "coarse").is_dir() else []
        record = {
            "registry_id": registry_id, "status": failure_status, "scf_status": scf_status,
            "exit_code": return_code, "runtime_seconds": round(time.time() - start, 3),
            "failure_reason": reason, "partial_output_files": partial,
        }
    write_json(directory / "SCF_STATUS.json", {"registry_id": registry_id, "scf_status": scf_status, "exit_code": return_code, "result_written": result_path.is_file(), "partial_output_preserved": bool(scf_log), "last_scf_log_lines": scf_log.splitlines()[-6:]})
    write_json(directory / "RUN_STATUS.json", record)
    return record


def dispatch(registry_ids: tuple[str, ...], execute_one: Callable[[str], dict]) -> list[dict]:
    """Invoke each registry independently and keep going after returned errors/exceptions."""
    outcomes = []
    for registry_id in registry_ids:
        try:
            outcomes.append(execute_one(registry_id))
        except Exception as exc:
            outcomes.append({"registry_id": registry_id, "status": "FAILED", "scf_status": "NOT_REACHED", "exit_code": None, "failure_reason": f"{type(exc).__name__}: {exc}"})
    return outcomes


def build_ranking(outcomes: list[dict], authorization_sha256: str, configuration_sha256: str) -> dict:
    successes = sorted((r for r in outcomes if r.get("status") == "SUCCESS"), key=lambda r: (r["total_energy_eV"], r["registry_id"]))
    failures = []
    for r in outcomes:
        if r.get("status") != "SUCCESS":
            failures.append({
                "registry_id": r["registry_id"], "status": r.get("status", "FAILED"),
                "scf_status": r.get("scf_status", "NOT_REACHED"),
                "exit_code": r.get("exit_code"),
                "failure_reason": r.get("failure_reason", "unspecified failure"),
            })
    return {
        "schema_version": 1,
        "configuration_sha256": configuration_sha256,
        "authorization_sha256": authorization_sha256,
        "role": "REGISTRY_SCREENING_ONLY",
        "kpoint_mesh": [6, 6, 1],
        "SUCCESS": successes,
        "FAILED": failures,
        "status": "READY_FOR_8X8_VALIDATION" if len(successes) >= 2 else "INSUFFICIENT_FOR_8X8",
        "auto_launch_8x8": False,
    }


def write_checksums(root: Path) -> None:
    lines = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p.name != "SHA256SUMS.txt"):
        lines.append(f"{sha256_file(path)}  {path.relative_to(root).as_posix()}")
    (root / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_self_test() -> dict:
    seen = []
    def mock_execute(registry_id: str) -> dict:
        seen.append(registry_id)
        if registry_id in ("R1", "R4"):
            raise RuntimeError(f"synthetic {registry_id} failure")
        return {"registry_id": registry_id, "status": "SUCCESS", "scf_status": "CONVERGED", "exit_code": 0, "total_energy_eV": -10.0 if registry_id == "R2" else -9.0, "forces_eV_A": [], "stress_eV_A3": []}
    outcomes = dispatch(REGISTRIES, mock_execute)
    ranking = build_ranking(outcomes, AUTH_SHA256, CONFIG_SHA256)
    assert seen == list(REGISTRIES), seen
    assert [r["registry_id"] for r in ranking["SUCCESS"]] == ["R2", "R3"]
    assert [r["registry_id"] for r in ranking["FAILED"]] == ["R1", "R4"]
    assert all(r["failure_reason"] for r in ranking["FAILED"])
    assert ranking["status"] == "READY_FOR_8X8_VALIDATION" and ranking["auto_launch_8x8"] is False
    failed_scf, reason = failure_reason("gpaw.KohnShamConvergenceError", "", "iteration 160", 42)
    assert failed_scf == "FAILED_SCF" and "partial output" in reason
    postprocess, reason = failure_reason("Traceback\nrequire(bool(calc.scf.converged), \"SCF did not converge\")\nAttributeError: 'ASECalculator' object has no attribute 'scf'", "", "iteration 124", 42)
    assert postprocess == "FAILED_POSTPROCESS" and "convergence is unverified" in reason
    fewer = build_ranking([outcomes[0], outcomes[1]], AUTH_SHA256, CONFIG_SHA256)
    assert fewer["status"] == "INSUFFICIENT_FOR_8X8"
    return {"status": "PASS", "test": "continue_on_error_failure_classification_and_ranking_thresholds", "registries_attempted": seen, "synthetic_failures_did_not_stop_dispatch": True, "max_iteration_failure_classified": failed_scf, "calculator_attribute_error_not_misclassified_as_scf": postprocess, "dft_started": False}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true", help="test scheduler logic with synthetic outcomes; no DFT")
    parser.add_argument("--execute", action="store_true", help="execute authorized R1-R4 jobs; otherwise only validate and print commands")
    parser.add_argument("--work-root", type=Path, help="fresh root containing R1-R4 input packages")
    parser.add_argument("--authorization", type=Path, default=ROOT / "REGISTRY_SCREENING_EXEC_AUTH_v2.json")
    parser.add_argument("--configuration", type=Path, default=ROOT / "AG_SI_REGISTRY_DFT_RELEASE_v2.json")
    parser.add_argument("--gpaw", type=Path, default=Path(shutil.which("gpaw") or "gpaw"))
    parser.add_argument("--mpiexec", default="mpiexec")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(run_self_test(), indent=2))
        return 0
    if args.work_root is None:
        parser.error("--work-root is required outside --self-test")
    work_root = args.work_root.resolve()
    auth, _ = validate_authority(args.authorization, args.configuration)
    plans = validate_candidate_files(work_root, auth, execute=args.execute)
    launcher_dir = ROOT / "launcher"
    preflight = validate_environment(args.gpaw.resolve(), args.mpiexec, launcher_dir, work_root)
    if not args.execute:
        print(json.dumps({"status": "DRY_RUN_ONLY", "dft_started": False, "scheduler_policy": {"fail_fast": False, "continue_on_error": True}, "preflight": preflight, "plans": plans, "note": "R1-R4 were not launched. Use --execute only for an authorized fresh work package."}, indent=2))
        return 0
    outcomes = dispatch(REGISTRIES, lambda rid: run_one(rid, work_root, gpaw_exe=args.gpaw.resolve(), mpiexec=args.mpiexec, launcher_dir=launcher_dir))
    ranking = build_ranking(outcomes, AUTH_SHA256, CONFIG_SHA256)
    ranking["preflight"] = preflight
    write_json(work_root / "RANKING_6x6.json", ranking)
    write_checksums(work_root)
    print(json.dumps({"status": ranking["status"], "dft_started": True, "successful": [r["registry_id"] for r in ranking["SUCCESS"]], "failed": [r["registry_id"] for r in ranking["FAILED"]], "ranking": str(work_root / "RANKING_6x6.json")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
