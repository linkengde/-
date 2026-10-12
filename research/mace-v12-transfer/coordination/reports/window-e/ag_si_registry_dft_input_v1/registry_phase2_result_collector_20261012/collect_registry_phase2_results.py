#!/usr/bin/env python3
"""Collect existing Ag-Si phase-2 GPAW outputs without starting calculations.

The collector is deliberately read-only with respect to input/result trees. It
only writes the requested summary JSON. Phase 2 means the 8x8x1 refinement;
6x6x1 phase-1 output is rejected as a mesh mismatch and is never substituted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

ALLOWED_IDS = ("R2", "R1")
EXPECTED_KPOINT = [8, 8, 1]
PACKAGE_DIR = Path(__file__).resolve().parent
REGISTRY_ROOT = PACKAGE_DIR.parent
PHASE1_ARCHIVE = REGISTRY_ROOT / "phase1_6x6x1_20261011_exec04"
FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def first_value(objects: list[dict[str, Any]], keys: tuple[str, ...]) -> Any:
    for obj in objects:
        for key in keys:
            value = obj.get(key)
            if value is not None:
                return value
    return None


def numeric(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def numeric_matrix(value: Any) -> list[list[float]] | None:
    if not isinstance(value, list) or not value:
        return None
    rows: list[list[float]] = []
    for row in value:
        if not isinstance(row, list) or not row:
            return None
        converted = [numeric(item) for item in row]
        if any(item is None for item in converted):
            return None
        rows.append([float(item) for item in converted])
    return rows


def parse_kpoint(value: Any) -> list[int] | None:
    if isinstance(value, dict):
        value = value.get("size") or value.get("mesh")
    if isinstance(value, (list, tuple)) and len(value) == 3:
        try:
            mesh = [int(item) for item in value]
        except (TypeError, ValueError):
            return None
        return mesh if all(item > 0 for item in mesh) else None
    return None


def is_sha256(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-fA-F]{64}", value) is not None


def parse_log(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace")
    result: dict[str, Any] = {"source_path": str(path)}

    mesh_match = re.search(
        r"kpts\s*=\s*MonkhorstPack\(size=array\(\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\]\)\)",
        text,
    )
    if mesh_match:
        result["kpoint"] = [int(item) for item in mesh_match.groups()]

    convergence_matches = list(re.finditer(r"Converged in\s+(\d+)\s+steps", text, re.I))
    if convergence_matches:
        result["scf_converged"] = True
        result["scf_steps"] = int(convergence_matches[-1].group(1))
    else:
        iteration_numbers = [int(item) for item in re.findall(r"\|iter:\s*(\d+)\|", text)]
        if iteration_numbers:
            result["scf_steps"] = max(iteration_numbers)
        result["scf_converged"] = False if iteration_numbers else None

    free_matches = re.findall(r"^\s*Free energy:\s*([-+0-9.eE]+)", text, re.M)
    energy_matches = re.findall(r"^\s*Extrapolated:\s*([-+0-9.eE]+)", text, re.M)
    if free_matches:
        result["free_energy"] = numeric(float(free_matches[-1]))
    if energy_matches:
        result["energy"] = numeric(float(energy_matches[-1]))

    force_heading = re.search(r"^Forces\s+#\s*eV/Ang\s*$", text, re.M)
    if force_heading:
        force_rows: list[list[float]] = []
        number = rf"({FLOAT_PATTERN})"
        for line in text[force_heading.end():].splitlines():
            match = re.match(
                rf"^\|\s*{number}\s*\|\s*{number}\s*\|\s*{number}\s*\|",
                line,
            )
            if match:
                force_rows.append([float(item) for item in match.groups()])
            elif force_rows:
                break
        if force_rows:
            result["forces"] = force_rows

    stress_match = re.search(
        r"Stress(?: tensor)?(?:\s+#\s*eV(?:/Ang(?:strom)?\^?3|/\u00c5\^?3))?\s*:?\s*\n((?:\s*[-+0-9.eE]+(?:\s+|\s*\n)){6,})",
        text,
        re.I,
    )
    if stress_match:
        values = [float(item) for item in re.findall(r"[-+0-9.eE]+", stress_match.group(1))]
        if len(values) >= 9:
            result["stress"] = [values[i:i + 3] for i in range(0, 9, 3)]
    return result


def result_files(case_dir: Path) -> list[Path]:
    if not case_dir.is_dir():
        return []
    preferred = [
        case_dir / "results" / "medium" / "result.json",
        case_dir / "results" / "8x8x1" / "result.json",
        case_dir / "results" / "phase2" / "result.json",
        case_dir / "result.json",
    ]
    selected = [path for path in preferred if path.is_file()]
    return selected if selected else sorted(case_dir.rglob("result.json"))


def log_files(case_dir: Path) -> list[Path]:
    if not case_dir.is_dir():
        return []
    preferred = [
        case_dir / "results" / "medium" / "gpaw.log",
        case_dir / "results" / "8x8x1" / "gpaw.log",
        case_dir / "results" / "phase2" / "gpaw.log",
    ]
    selected = [path for path in preferred if path.is_file()]
    return selected if selected else sorted(case_dir.rglob("gpaw.log"))


def sidecar_objects(case_dir: Path, result_obj: dict[str, Any] | None) -> list[dict[str, Any]]:
    objects: list[dict[str, Any]] = []
    for name in ("C_EXECUTION_RELEASE_RECORD.json", "SCF_STATUS.json", "RUN_STATUS.json"):
        obj = read_json(case_dir / name)
        if obj:
            objects.append(obj)
    if result_obj:
        objects.insert(0, result_obj)
        nested = result_obj.get("results")
        if isinstance(nested, dict):
            objects.insert(0, nested)
    return objects


def baseline_paths(registry_id: str) -> tuple[Path, Path]:
    input_path = REGISTRY_ROOT / registry_id / "input.py"
    structure_path = PHASE1_ARCHIVE / registry_id / "structure.xyz"
    return input_path, structure_path


def collect_one(registry_id: str, phase2_root: Path) -> dict[str, Any]:
    case_dir = phase2_root / registry_id
    input_path = case_dir / "input.py"
    structure_path = case_dir / "structure.xyz"
    if not input_path.is_file() or not structure_path.is_file():
        default_input, default_structure = baseline_paths(registry_id)
        input_path = input_path if input_path.is_file() else default_input
        structure_path = structure_path if structure_path.is_file() else default_structure

    input_sha = sha256(input_path) if input_path.is_file() else None
    structure_sha = sha256(structure_path) if structure_path.is_file() else None

    json_paths = result_files(case_dir)
    log_paths = log_files(case_dir)
    record = read_json(case_dir / "C_EXECUTION_RELEASE_RECORD.json") or {}
    result_obj: dict[str, Any] | None = None
    log_obj: dict[str, Any] = {}
    data_status = "RESULTS_NOT_FOUND"
    issues: list[str] = []

    if len(json_paths) > 1:
        issues.append("Multiple result.json files found; refusing ambiguous selection")
        data_status = "AMBIGUOUS_RESULT_FILES"
    elif json_paths:
        result_obj = read_json(json_paths[0])
        if result_obj is None:
            issues.append("result.json is malformed or not a JSON object")
            data_status = "MALFORMED_RESULT_JSON"

    if len(log_paths) > 1:
        issues.append("Multiple GPAW logs found; refusing ambiguous selection")
        data_status = "AMBIGUOUS_GPAW_LOGS"
    elif log_paths:
        log_obj = parse_log(log_paths[0])

    objects = sidecar_objects(case_dir, result_obj)
    # Metadata from the result or its execution record takes precedence. The
    # phase-1 6x6 config hash is intentionally never used for phase 2.
    config_sha = first_value(objects, ("config_sha", "config_sha256", "configuration_sha256", "configuration_sha"))
    kpoint = parse_kpoint(first_value(objects, ("kpoint", "kpoint_mesh", "kpts")))
    if kpoint is None:
        kpoint = parse_kpoint(log_obj.get("kpoint"))

    energy = numeric(first_value(objects, ("energy", "energy_eV", "total_energy_eV")))
    free_energy = numeric(first_value(objects, ("free_energy", "free_energy_eV", "force_consistent_free_energy_eV")))
    forces = numeric_matrix(first_value(objects, ("forces", "forces_eV_A", "atomic_forces_eV_A")))
    stress = numeric_matrix(first_value(objects, ("stress", "stress_eV_A3", "stress_tensor_eV_A3")))
    scf_steps_val = first_value(objects, ("scf_steps", "iterations", "n_scf_steps"))
    scf_steps = int(scf_steps_val) if isinstance(scf_steps_val, int) and not isinstance(scf_steps_val, bool) else None
    scf_converged = first_value(objects, ("scf_converged", "converged"))
    scf_status = first_value(objects, ("scf_status", "SCF_status", "status"))

    if energy is None:
        energy = numeric(log_obj.get("energy"))
    if free_energy is None:
        free_energy = numeric(log_obj.get("free_energy"))
    if forces is None:
        forces = numeric_matrix(log_obj.get("forces"))
    if stress is None:
        stress = numeric_matrix(log_obj.get("stress"))
    if scf_steps is None and isinstance(log_obj.get("scf_steps"), int):
        scf_steps = log_obj["scf_steps"]
    if scf_converged is None:
        scf_converged = log_obj.get("scf_converged")

    recorded_input_sha = first_value(objects, ("input_sha", "input_sha256"))
    recorded_structure_sha = first_value(objects, ("structure_sha", "structure_sha256"))
    if recorded_input_sha is not None and input_sha and (
        not is_sha256(recorded_input_sha) or recorded_input_sha.lower() != input_sha.lower()
    ):
        data_status = "PROVENANCE_MISMATCH"
        issues.append("Recorded input SHA does not match the available input file")
    if recorded_structure_sha is not None and structure_sha and (
        not is_sha256(recorded_structure_sha) or recorded_structure_sha.lower() != structure_sha.lower()
    ):
        data_status = "PROVENANCE_MISMATCH"
        issues.append("Recorded structure SHA does not match the available structure file")
    if config_sha is not None and not is_sha256(config_sha):
        issues.append("Phase-2 config SHA is not a 64-character SHA-256")
        config_sha = None

    if not any(issue.startswith("Recorded ") for issue in issues) and data_status not in (
        "AMBIGUOUS_RESULT_FILES", "AMBIGUOUS_GPAW_LOGS", "MALFORMED_RESULT_JSON"
    ) and (result_obj is not None or log_obj):
        if kpoint != EXPECTED_KPOINT:
            data_status = "KPOINT_MISMATCH" if kpoint is not None else "KPOINT_UNVERIFIED"
            issues.append(f"Expected phase-2 kpoint {EXPECTED_KPOINT}; observed {kpoint}")
            # Do not expose numerical labels from a result that is not proven
            # to belong to the requested phase-2 mesh.
            energy = free_energy = None
            forces = stress = None
            scf_steps = None
            scf_converged = None
        else:
            data_status = "RESULTS_COLLECTED"
    elif data_status == "RESULTS_NOT_FOUND":
        issues.append("No phase-2 result.json or gpaw.log was found for this registry")
    if config_sha is None:
        issues.append("No phase-2-bound C configuration SHA is available")

    if data_status not in ("RESULTS_COLLECTED", "RESULTS_NOT_FOUND"):
        energy = free_energy = None
        forces = stress = None
        scf_steps = None
        scf_converged = None

    if scf_converged is None and scf_status:
        status_text = str(scf_status).upper()
        if "CONVERGED" in status_text and "NOT_CONVERGED" not in status_text:
            scf_converged = True
        elif any(token in status_text for token in ("FAILED", "UNCONVERGED", "NOT_CONVERGED")):
            scf_converged = False

    forces_status = "AVAILABLE" if forces is not None else "NOT_AVAILABLE"
    stress_status = "AVAILABLE" if stress is not None else "NOT_AVAILABLE"
    if result_obj and isinstance(result_obj.get("forces_status"), str) and forces is None:
        forces_status = result_obj["forces_status"]
    if result_obj and isinstance(result_obj.get("stress_status"), str) and stress is None:
        stress_status = result_obj["stress_status"]

    return {
        "structure_id": registry_id,
        "input_sha": recorded_input_sha or input_sha,
        "structure_sha": recorded_structure_sha or structure_sha,
        "config_sha": config_sha,
        "kpoint": kpoint or EXPECTED_KPOINT,
        "energy": energy,
        "free_energy": free_energy,
        "forces": forces,
        "forces_status": forces_status,
        "stress": stress,
        "stress_status": stress_status,
        "scf_converged": scf_converged,
        "scf_steps": scf_steps,
        "data_status": data_status,
        "issues": issues,
        "result_source": str((json_paths[0] if len(json_paths) == 1 else log_paths[0]) if (len(json_paths) == 1 or len(log_paths) == 1) else "") or None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase2-root",
        type=Path,
        default=REGISTRY_ROOT / "phase2_8x8x1_20261012",
        help="read-only root containing only phase-2 R1/R2 result directories",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PACKAGE_DIR / "REGISTRY_PHASE2_RESULTS.json",
        help="summary JSON to write (defaults beside this script)",
    )
    args = parser.parse_args()

    results = [collect_one(registry_id, args.phase2_root) for registry_id in ALLOWED_IDS]
    all_collected = all(
        item["data_status"] == "RESULTS_COLLECTED"
        and is_sha256(item["config_sha"])
        and is_sha256(item["input_sha"])
        and is_sha256(item["structure_sha"])
        for item in results
    )
    document = {
        "document_id": "REGISTRY_PHASE2_RESULTS",
        "phase": "phase2_8x8x1",
        "authority_basis": "REGISTRY_PHASE2_RELEASE_REVIEW; allowed IDs supplied by C review: R2, R1",
        "allowed_structure_ids": list(ALLOWED_IDS),
        "expected_kpoint": EXPECTED_KPOINT,
        "phase2_root": str(args.phase2_root),
        "collection_status": "COMPLETE" if all_collected else "PARTIAL_OR_PENDING",
        "dft_started_by_collector": False,
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"collection_status": document["collection_status"], "output": str(args.output), "items": {item["structure_id"]: item["data_status"] for item in results}}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
