#!/usr/bin/env python3
"""Collect labels only from a completed, converged, hash-bound R2 run."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

from runner import CONFIG_PATH, REPO_ROOT, load_config, sha256_file, write_json_atomic


CONVERGED_RE = re.compile(r"Converged in\s+(\d+)\s+steps")


def read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def log_converged_steps(path: Path) -> int | None:
    try:
        content = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    values = [int(match.group(1)) for match in CONVERGED_RE.finditer(content)]
    return values[-1] if values else None


def finite_number(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def valid_array(value: Any, shape: tuple[int, ...]) -> bool:
    if not isinstance(value, list):
        return False
    try:
        import numpy as np
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        return False
    return array.shape == shape and bool(np.isfinite(array).all())


def collect(run_dir: Path, stdout_path: Path | None = None,
            stderr_path: Path | None = None) -> dict[str, Any]:
    config, config_sha = load_config()
    source = config["source_bindings"]
    status = read_json(run_dir / "RUN_STATUS.json")
    raw = read_json(run_dir / "RESULTS_RAW.json")
    convergence_steps = log_converged_steps(run_dir / "gpaw.log")
    failure_reasons: list[str] = []

    if status is None:
        failure_reasons.append("RUN_STATUS_MISSING_OR_INVALID")
    if raw is None:
        failure_reasons.append("RESULTS_RAW_MISSING_OR_INVALID")
    if convergence_steps is None:
        failure_reasons.append("GPAW_CONVERGENCE_MARKER_MISSING")

    successful = bool(
        status is not None
        and raw is not None
        and status.get("run_status") == "SUCCESS"
        and status.get("scf_converged") is True
        and raw.get("status") == "SUCCESS"
        and raw.get("scf_converged") is True
        and convergence_steps is not None
        and status.get("scf_steps") == convergence_steps
        and raw.get("scf_steps") == convergence_steps
        and raw.get("registry_id") == "R2"
        and raw.get("role") == "REGISTRY_SCREENING_ONLY_DIAGNOSTIC_COPY"
        and raw.get("input_sha256") == source["input_sha256"]
        and raw.get("structure_sha256") == source["structure_sha256"]
        and raw.get("parameter_lock_sha256") == source["parameter_lock_sha256"]
        and raw.get("config_sha256") == config_sha
        and bool(raw.get("authorization_sha256"))
        and raw.get("authorization_sha256") == status.get("authorization_sha256")
        and bool(raw.get("authorization_id"))
        and raw.get("authorization_id") == status.get("authorization_id")
        and finite_number(raw.get("total_energy_eV"))
        and finite_number(raw.get("force_consistent_free_energy_eV"))
        and valid_array(raw.get("forces_eV_A"), (20, 3))
        and valid_array(raw.get("stress_eV_A3"), (3, 3))
    )

    if not successful:
        if status and status.get("failure_reason"):
            failure_reasons.append(str(status["failure_reason"]))
        # Deliberately omit all partial energy/force/stress values.
        result = {
            "status": "FAILED_SCF",
            "registry_id": "R2",
            "role": "REGISTRY_SCREENING_ONLY_DIAGNOSTIC_COPY",
            "scf_converged": False,
            "scf_steps": status.get("scf_steps") if status else None,
            "total_energy_eV": None,
            "force_consistent_free_energy_eV": None,
            "forces_eV_A": None,
            "stress_eV_A3": None,
            "forces_status": "NOT_A_LABEL",
            "stress_status": "NOT_A_LABEL",
            "failure_reason": "; ".join(dict.fromkeys(failure_reasons)) or "SCF_NOT_CONVERGED",
            "authorization_sha256": status.get("authorization_sha256") if status else None,
            "authorization_id": status.get("authorization_id") if status else None,
        }
    else:
        result = {
            "status": "SUCCESS",
            "registry_id": "R2",
            "role": raw["role"],
            "scf_converged": True,
            "scf_steps": convergence_steps,
            "total_energy_eV": float(raw["total_energy_eV"]),
            "force_consistent_free_energy_eV": float(raw["force_consistent_free_energy_eV"]),
            "forces_eV_A": raw["forces_eV_A"],
            "stress_eV_A3": raw["stress_eV_A3"],
            "forces_status": "AVAILABLE",
            "stress_status": "AVAILABLE",
            "failure_reason": None,
            "authorization_sha256": raw["authorization_sha256"],
            "authorization_id": raw["authorization_id"],
        }

    result.update({
        "input_sha256": source["input_sha256"],
        "structure_sha256": source["structure_sha256"],
        "parameter_lock_sha256": source["parameter_lock_sha256"],
        "config_sha256": config_sha,
        "kpoint_mesh": config["method"]["kpoint_mesh"],
        "fixed_geometry": True,
        "raw_artifacts": {},
    })
    candidates = [
        run_dir / "RUN_STATUS.json",
        run_dir / "RESULTS_RAW.json",
        run_dir / "gpaw.log",
        stdout_path,
        stderr_path,
    ]
    for path in candidates:
        if path is not None and path.is_file():
            try:
                result["raw_artifacts"][str(path)] = {
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            except OSError:
                pass
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--launcher-stdout", type=Path)
    parser.add_argument("--launcher-stderr", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = collect(args.run_dir, args.launcher_stdout, args.launcher_stderr)
    output = args.output or args.run_dir / "RESULTS_COLLECTED.json"
    if args.run_dir.exists():
        write_json_atomic(output, result)
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "SUCCESS" else 5


if __name__ == "__main__":
    raise SystemExit(main())
