#!/usr/bin/env python3
"""Archive logs/manifests/checksums after an authorized Phase2 runner exits."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REGISTRY_ROOT = HERE.parent
REGISTRY_IDS = ("R1", "R2")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def write_manifest(case_dir: Path, root: Path, auth_sha: str) -> dict[str, Any]:
    result_path = case_dir / "result.json"
    record_path = case_dir / "C_EXECUTION_RELEASE_RECORD.json"
    record = read_json(record_path) or {}
    result = read_json(result_path) or {}
    artifacts = {}
    for name in ("C_EXECUTION_RELEASE_RECORD.json", "gpaw.log", "result.json", "stdout.log", "stderr.log"):
        path = case_dir / name
        if path.is_file():
            artifacts[name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    manifest = {
        "registry_id": case_dir.name,
        "authorization_sha256": auth_sha,
        "parameter_lock_sha256": record.get("parameter_lock_sha256"),
        "runner_sha256": record.get("runner_sha256"),
        "collector_sha256": record.get("collector_sha256"),
        "input_sha256": record.get("input_sha256"),
        "structure_sha256": record.get("structure_sha256"),
        "kpoint_mesh": record.get("kpoint_mesh"),
        "role": record.get("role"),
        "run_status": result.get("run_status", "RESULT_NOT_WRITTEN"),
        "scf_converged": result.get("scf_converged"),
        "scf_steps": result.get("scf_steps"),
        "artifacts": artifacts,
        "stdout_scope": "Shared MPI runner stdout captured for ordered R1 then R2 invocation.",
        "stderr_scope": "Per-registry runner traceback when available; otherwise shared MPI launcher stderr.",
    }
    (case_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    sums = []
    for path in sorted(case_dir.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            sums.append(f"{sha256(path)}  {path.name}")
    (case_dir / "SHA256SUMS.txt").write_text("\n".join(sums) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--authorization", type=Path, required=True)
    parser.add_argument("--launcher-stdout", type=Path, required=True)
    parser.add_argument("--launcher-stderr", type=Path, required=True)
    args = parser.parse_args()
    root = args.work_root
    if not root.is_dir():
        raise SystemExit("work root is missing")
    auth_sha = sha256(args.authorization)
    archive_auth = root / "REGISTRY_PHASE2_EXEC_AUTH_v2.json"
    if not archive_auth.exists():
        shutil.copyfile(args.authorization, archive_auth)
    shutil.copyfile(args.launcher_stdout, root / "launcher.stdout.log")
    shutil.copyfile(args.launcher_stderr, root / "launcher.stderr.log")

    registry_manifests = []
    for registry_id in REGISTRY_IDS:
        case_dir = root / registry_id
        if not case_dir.is_dir():
            continue
        stdout = case_dir / "stdout.log"
        stderr = case_dir / "stderr.log"
        if not stdout.exists():
            shutil.copyfile(root / "launcher.stdout.log", stdout)
        if not stderr.exists():
            shutil.copyfile(root / "launcher.stderr.log", stderr)
        registry_manifests.append(write_manifest(case_dir, root, auth_sha))

    execution = {
        "document_id": "AG_SI_PHASE2_EXECUTION_MANIFEST",
        "authorization_file": archive_auth.name,
        "authorization_sha256": auth_sha,
        "source_authorization_file_sha256": auth_sha,
        "allowed_registry_order": list(REGISTRY_IDS),
        "kpoint_mesh": [8, 8, 1],
        "calculation": "PBE/PAW fixed-geometry single point",
        "mpi_ranks": 4,
        "omp_threads": 1,
        "role": "REGISTRY_SCREENING_ONLY",
        "dft_started": True,
        "registries": registry_manifests,
    }
    (root / "EXECUTION_MANIFEST.json").write_text(json.dumps(execution, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Ag–Si Phase2 Execution Report",
        "",
        f"- Authorization SHA-256: `{auth_sha}`",
        "- Scope: R1 then R2; 8×8×1; PBE/PAW; fixed geometry; MPI 4; OMP 1.",
        "- Role: `REGISTRY_SCREENING_ONLY`.",
        "- Failures are recorded per registry; the authorized runner proceeds to the next registry.",
        "",
        "| Registry | Status | Energy (eV) | Free energy (eV) | SCF converged | SCF steps |",
        "|---|---|---:|---:|---|---:|",
    ]
    for registry_id in REGISTRY_IDS:
        result = read_json(root / registry_id / "result.json") or {}
        lines.append(
            f"| {registry_id} | {result.get('run_status', 'NOT_STARTED_OR_NO_RESULT')} | "
            f"{result.get('total_energy_eV')} | {result.get('force_consistent_free_energy_eV')} | "
            f"{result.get('scf_converged')} | {result.get('scf_steps')} |"
        )
    (root / "PHASE2_EXECUTION_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    root_sums = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            root_sums.append(f"{sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256SUMS.txt").write_text("\n".join(root_sums) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
