#!/usr/bin/env python3
"""C-lock-gated Ag64 rho100 continuation runner (not run in this task).

The only current configuration has execution_authorized=false.  The --run
path additionally requires a C release blob that binds every executable and
input hash.  This script never creates or reinitializes an Ag structure.
"""
from __future__ import annotations

import argparse
import csv
import copy
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import platform
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
import time

import ase
from ase import Atoms
from ase.calculators.eam import EAM
from ase.constraints import FixCom
from ase.io import Trajectory, write
from ase.md.langevin import Langevin
from ase.units import fs
import numpy as np

import archive_backend
import liquid_criteria as criteria


PACKAGE_DIR = Path(__file__).resolve().parent
CONFIG_NAME = "FOLLOWUP_EXECUTION_CONFIG.json"
SCIENTIFIC_LOCK_SHA256 = "e1adaf12d542184c709d2eb2a5363cbfe9578c168ebc5e7ad260ee39a4faf09b"
STRATEGY_LOCK_SHA256 = "054f36dcda113e37b23b112e5919ba1c8201c669626db37b79bbd2281b020e1c"
EXPECTED_CHECKPOINT_SHA256 = "c820437a2512c4b3d7dc9e60bb60ac3ced6588aa1f2a629d58cfdabb24f34046"
EXPECTED_POTENTIAL_SHA256 = "958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447"
EXPECTED_ENVIRONMENT_SHA256 = "f3f6d27f7e8b387bb5cc7b91048dba7fa26280f9ed54aa1572c73879946e867f"
TASK_ID = "C-D-AGL-1250K-rho100-followup01"
SOURCE_FAMILY_ID = "PF-WMH2006-Ag64-rho100-Dpilot-20261012-seed003"
RNG_NAMES = ("packing", "velocity_initialization", "thermostat")


class GateError(RuntimeError):
    pass


_STOP_REQUESTED = False


def _request_graceful_stop(_signum, _frame):
    global _STOP_REQUESTED
    _STOP_REQUESTED = True


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_sha(value: dict) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(payload).hexdigest()


def cooling_target_K(completed_ramp_steps: int) -> float:
    """C-frozen ramp: update 3.5 K per 100 completed steps over 10,000."""
    if completed_ramp_steps < 0 or completed_ramp_steps > 10000:
        raise GateError("cooling-ramp position is outside its frozen 10,000-step boundary")
    updates = min(100, completed_ramp_steps // 100)
    return 1600.0 - 3.5 * updates


def next_chunk_end(current_step: int, phase_end_step: int) -> int:
    if current_step < 1000 or phase_end_step < current_step:
        raise GateError("invalid segment boundary")
    return min(60000, phase_end_step, current_step + 1000)


def advance_completed_phase(state: dict, current_step: int) -> None:
    """Apply only the frozen B-to-C boundary after a completed ramp chunk."""
    if (
        state.get("stage") == "B_COOL_1600_TO_1250K"
        and current_step >= int(state["stage_b_start_step"]) + 10000
    ):
        state["stage"] = "C_HOLD_1250K"
        state["stage_c_start_step"] = current_step


def apply_monitor_transition(state: dict, current_step: int, monitor: dict, candidate_spacing: dict | None = None) -> str:
    """Apply C-frozen phase transitions from an already computed QA result."""
    status = "CHUNK_COMPLETE"
    if current_step % 1000 or state.get("stage") not in {
        "A_HOT_1600K", "B_COOL_1600_TO_1250K", "C_HOLD_1250K", "D_CANDIDATE_1250K"
    }:
        return status
    gate = monitor.get("gate_status")
    if gate == "PASS" and state["stage"] == "A_HOT_1600K":
        state["liquid_onset_step"] = current_step
        state["stage_b_start_step"] = current_step
        state["stage"] = "B_COOL_1600_TO_1250K"
    elif gate == "PASS" and state["stage"] == "C_HOLD_1250K":
        if not isinstance(candidate_spacing, dict) or "snapshot_spacing_steps" not in candidate_spacing:
            raise GateError("passing phase-C gate requires measured candidate snapshot spacing")
        state["candidate_snapshot_spacing"] = candidate_spacing
        state["candidate_snapshot_spacing_steps"] = candidate_spacing["snapshot_spacing_steps"]
        state["stage_d_start_step"] = current_step
        state["stage"] = "D_CANDIDATE_1250K"
    elif state["stage"] == "D_CANDIDATE_1250K" and gate == "FAIL_WINDOW":
        status = "HOLD_CANDIDATE_WINDOW_QA_FAILED"
        state["stage"] = status
    return status


def apply_completed_chunk_limits(state: dict, current_step: int, status: str, maximum_global: int) -> str:
    """Apply frozen stage/global stop gates after QA and before archiving."""
    if status != "CHUNK_COMPLETE":
        return status
    stage = state["stage"]
    if stage == "A_HOT_1600K" and current_step >= 20000:
        status = "HOLD_STAGE_A_MAXIMUM_WITHOUT_LIQUID_GATE"
    elif stage == "C_HOLD_1250K" and current_step >= min(60000, int(state["stage_c_start_step"]) + 20000):
        status = "HOLD_STAGE_C_MAXIMUM_WITHOUT_LIQUID_GATE"
    elif stage == "D_CANDIDATE_1250K" and current_step >= int(state["stage_d_start_step"]) + 10000:
        status = (
            "CANDIDATE_WINDOW_COMPLETE_FINAL_QA_PASS"
            if state.get("last_monitor", {}).get("gate_status") == "PASS"
            else "HOLD_CANDIDATE_FINAL_QA_FAILED"
        )
    if status == "CHUNK_COMPLETE" and current_step >= maximum_global:
        status = "HOLD_GLOBAL_STEP_LIMIT"
    if status != "CHUNK_COMPLETE":
        state["stage"] = status
    return status


def _verify_release_blob(repo: Path, commit: str, repo_path: str, expected_sha: str, auth_path: Path) -> None:
    if len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit):
        raise GateError("C release commit must be a full hexadecimal commit SHA")
    pure = PurePosixPath(repo_path)
    if pure.is_absolute() or ".." in pure.parts:
        raise GateError("C release repository path is unsafe")
    result = subprocess.run(
        ["git", "-C", str(repo), "show", f"{commit}:{pure.as_posix()}"],
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if result.returncode or hashlib.sha256(result.stdout).hexdigest() != expected_sha:
        raise GateError("authorization file is not the SHA-matching blob from the supplied C commit")
    if sha256_file(auth_path) != expected_sha or auth_path.read_bytes() != result.stdout:
        raise GateError("local authorization copy differs from the C commit blob")


def _git_root(path: Path) -> Path:
    result = subprocess.run(["git", "-C", str(path), "rev-parse", "--show-toplevel"], check=False, capture_output=True, text=True)
    if result.returncode:
        raise GateError("package must run inside the D Git checkout")
    return Path(result.stdout.strip()).resolve()


def _environment_snapshot() -> dict:
    import gpaw
    import gpaw_data
    runtime = Path("/workspace/onboarding/runtime.sh")
    if not runtime.is_file():
        raise GateError("runtime.sh is missing")
    mpi = subprocess.run(["mpiexec", "--version"], check=True, capture_output=True, text=True, timeout=10).stdout.splitlines()[0].strip()
    return {
        "python_version": ".".join(map(str, sys.version_info[:3])),
        "python_executable": sys.executable,
        "ase_version": ase.__version__,
        "numpy_version": np.__version__,
        "scipy_version": __import__("scipy").__version__,
        "runtime_sh": str(runtime),
        "runtime_sh_sha256": sha256_file(runtime),
        "GPAW_SETUP_PATH": os.environ.get("GPAW_SETUP_PATH"),
        "GPAW_MPI_BACKEND": os.environ.get("GPAW_MPI_BACKEND"),
        "openmpi_version": mpi,
        "host_os": platform.platform(),
        "gpaw_version": gpaw.__version__,
        "gpaw_data_path": str(gpaw_data.datapath()),
    }


def _resource_snapshot(path: Path) -> dict:
    def read_text(p: Path):
        try:
            return p.read_text(encoding="ascii").strip()
        except OSError:
            return None

    cpu_raw = read_text(Path("/sys/fs/cgroup/cpu.max"))
    cpu_quota = cpu_period = None
    if cpu_raw and cpu_raw != "max":
        parts = cpu_raw.split()
        if len(parts) == 2:
            cpu_quota, cpu_period = int(parts[0]), int(parts[1])
    memory_raw = read_text(Path("/sys/fs/cgroup/memory.max"))
    memory_current_raw = read_text(Path("/sys/fs/cgroup/memory.current"))
    memory_limit = None if memory_raw in (None, "max") else int(memory_raw)
    memory_current = None if memory_current_raw is None else int(memory_current_raw)
    usage = shutil.disk_usage(path)
    vfs = os.statvfs(path)
    rss_kib = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    rss_bytes = rss_kib * (1024 if sys.platform != "darwin" else 1)
    cpu_count = os.cpu_count() or 1
    cpu_quota_count = (cpu_quota / cpu_period) if cpu_quota and cpu_period else float(cpu_count)
    return {
        "cpu_max_raw": cpu_raw,
        "cpu_quota_cores": cpu_quota_count,
        "cpu_count_visible": cpu_count,
        "memory_limit_bytes": memory_limit,
        "memory_current_bytes": memory_current,
        "process_max_rss_bytes": rss_bytes,
        "disk_available_bytes": usage.free,
        "disk_total_bytes": usage.total,
        "free_inodes": vfs.f_favail,
        "resource_path": str(path),
    }


def _resource_gate(snapshot: dict, limits: dict) -> bool:
    return bool(
        snapshot["cpu_quota_cores"] >= 1.0
        and (snapshot["memory_limit_bytes"] is None or snapshot["memory_limit_bytes"] >= limits["maximum_rss_bytes"])
        and snapshot["process_max_rss_bytes"] <= limits["maximum_rss_bytes"]
        and snapshot["disk_available_bytes"] >= limits["minimum_free_bytes"] + limits["checkpoint_reserve_bytes"]
        and snapshot["free_inodes"] > 0
    )


def _rng_restore(states: dict) -> dict[str, np.random.Generator]:
    if set(states) != set(RNG_NAMES):
        raise GateError("step1000 checkpoint lacks exactly the three frozen PCG64 child streams")
    result = {}
    for name in RNG_NAMES:
        bitgen = np.random.PCG64()
        bitgen.state = copy.deepcopy(states[name])
        result[name] = np.random.Generator(bitgen)
    return result


def _restore_atoms(checkpoint: dict) -> tuple[Atoms, dict[str, np.random.Generator]]:
    if checkpoint.get("global_step") != 1000:
        raise GateError("initial continuation checkpoint is not exactly step 1000")
    if checkpoint.get("symbols") != ["Ag"] * 64 or checkpoint.get("atom_ids") != list(range(1, 65)):
        raise GateError("checkpoint Ag64 composition or atom ordering differs from C's lock")
    positions = np.asarray(checkpoint["positions_A"], dtype=float)
    momenta = np.asarray(checkpoint["momenta_amu_A_per_fs"], dtype=float)
    if positions.shape != (64, 3) or momenta.shape != (64, 3) or not np.isfinite(positions).all() or not np.isfinite(momenta).all():
        raise GateError("checkpoint positions or momenta have invalid shape or values")
    atoms = Atoms(
        symbols=checkpoint["symbols"], positions=positions,
        cell=np.asarray(checkpoint["cell_A"], dtype=float), pbc=checkpoint["pbc"],
    )
    atoms.set_masses(np.asarray(checkpoint["masses_amu"], dtype=float))
    atoms.new_array("atom_id", np.asarray(checkpoint["atom_ids"], dtype=int))
    atoms.set_constraint(FixCom())
    # ASE stores momenta in its internal time unit; the existing checkpoint
    # stores amu*Angstrom/fs, so divide by ASE's fs conversion exactly once.
    atoms.set_momenta(momenta / fs, apply_constraint=False)
    restored_rngs = _rng_restore(checkpoint["rng_states"])
    if not np.array_equal(atoms.get_positions(), positions):
        raise GateError("checkpoint position restore was not exact")
    if not np.array_equal(atoms.get_momenta() * fs, momenta):
        raise GateError("checkpoint momentum restore was not exact")
    if not np.array_equal(
        np.asarray(atoms.cell), np.asarray(checkpoint["cell_A"], dtype=float)
    ) or not np.array_equal(np.asarray(atoms.pbc), np.asarray(checkpoint["pbc"], dtype=bool)):
        raise GateError("checkpoint cell or PBC restore was not exact")
    return atoms, restored_rngs


def _verify_json_sha(path: Path, expected: str, label: str) -> dict:
    if not path.is_file():
        raise GateError(f"missing {label}: {path}")
    actual = sha256_file(path)
    if actual != expected:
        raise GateError(f"{label} SHA mismatch: expected {expected}, got {actual}")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise GateError(f"{label} JSON parse failed: {type(exc).__name__}") from exc


def run_preflight(config_path: Path, expected_config_sha: str) -> dict:
    if not config_path.is_file():
        raise GateError(f"configuration missing: {config_path}")
    config_sha = sha256_file(config_path)
    if config_sha != expected_config_sha:
        raise GateError("execution config SHA mismatch")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("execution_authorized") is not False:
        raise GateError("packaged preparation config must stay execution_authorized=false")
    if config.get("task_id") != TASK_ID or config.get("owner") != "D":
        raise GateError("configuration owner/task identity mismatch")
    strategy = _verify_json_sha(PACKAGE_DIR / "C_STRATEGY_LOCK.json", STRATEGY_LOCK_SHA256, "C strategy lock")
    scientific = _verify_json_sha(PACKAGE_DIR / "C_SCIENTIFIC_LOCK.json", SCIENTIFIC_LOCK_SHA256, "C scientific lock")
    if strategy.get("execution_authorized") is not False:
        raise GateError("C frozen strategy unexpectedly grants execution")
    if scientific.get("potential", {}).get("sha256") != EXPECTED_POTENTIAL_SHA256:
        raise GateError("scientific lock potential identity mismatch")
    repo = _git_root(PACKAGE_DIR)
    branch = subprocess.run(["git", "-C", str(repo), "branch", "--show-current"], check=True, capture_output=True, text=True).stdout.strip()
    if branch != config["archive"]["branch"]:
        raise GateError("current checkout is not the D-owned archive branch")
    source_commit = config["source_checkpoint"]["D_source_commit"]
    if subprocess.run(["git", "-C", str(repo), "merge-base", "--is-ancestor", source_commit, "HEAD"], check=False).returncode != 0:
        raise GateError("D source archive commit is not an ancestor of the current D branch")
    if subprocess.run(["git", "-C", str(repo), "cat-file", "-e", config["strategy"]["C_commit"]], check=False).returncode != 0:
        raise GateError("fixed C strategy commit is unavailable in this checkout")
    source_rel = config["source_checkpoint"]["archive_path"]
    source_dir = repo / source_rel
    source_shas = config["source_checkpoint"]["member_sha256"]
    archive_backend.verify_sha256sums(source_dir, filename="SHA256SUMS.txt")
    actual_members = {}
    for rel, expected in source_shas.items():
        p = source_dir / rel
        actual = sha256_file(p) if p.is_file() else None
        actual_members[rel] = actual
        if actual != expected:
            raise GateError(f"archived source member SHA mismatch: {rel}")
    checkpoint_path = source_dir / "checkpoint.json"
    checkpoint = _verify_json_sha(checkpoint_path, EXPECTED_CHECKPOINT_SHA256, "step1000 checkpoint")
    sidecar = (source_dir / "checkpoint.sha256").read_text(encoding="ascii").split()
    if not sidecar or sidecar[0] != EXPECTED_CHECKPOINT_SHA256:
        raise GateError("step1000 checkpoint sidecar SHA mismatch")
    if checkpoint.get("lock_sha256") != SCIENTIFIC_LOCK_SHA256 or checkpoint.get("potential_sha256") != EXPECTED_POTENTIAL_SHA256:
        raise GateError("step1000 checkpoint is not bound to the frozen scientific lock/potential")
    if checkpoint.get("task_id") != config["source_checkpoint"]["task_id"]:
        raise GateError("source checkpoint task identity mismatch")
    if checkpoint.get("environment_sha256") != EXPECTED_ENVIRONMENT_SHA256:
        raise GateError("step1000 checkpoint environment identity differs from pilot evidence")
    atoms, restored_rngs = _restore_atoms(checkpoint)
    if len(atoms) != 64 or atoms.get_chemical_formula() != "Ag64":
        raise GateError("restored checkpoint composition mismatch")
    if len({int(x) for x in atoms.arrays["atom_id"]}) != 64:
        raise GateError("atom ID uniqueness check failed")
    expected_cell = np.asarray(scientific["density"]["cell"]["matrix_A"], dtype=float)
    if not np.allclose(np.asarray(atoms.cell), expected_cell, rtol=0.0, atol=1e-12) or not np.array_equal(np.asarray(atoms.pbc), np.asarray([True, True, True])):
        raise GateError("step1000 checkpoint cell or PBC differs from the original scientific lock")
    if checkpoint.get("time_fs") != 500.0 or checkpoint.get("parameter_sha256") != config["source_checkpoint"]["parameter_sha256"]:
        raise GateError("step1000 checkpoint time/parameter-lock identity mismatch")
    schedule = checkpoint.get("temperature_schedule_position", {})
    if schedule.get("last_completed_step") != 1000 or schedule.get("next_integrator_step_target_K") != 1600.0:
        raise GateError("step1000 checkpoint temperature schedule position mismatch")
    if checkpoint.get("runner_sha256") != config["source_checkpoint"]["pilot_runner_sha256"]:
        raise GateError("step1000 checkpoint is not bound to the archived pilot runner")
    if checkpoint.get("rng_provenance", {}).get("bit_generator") != "PCG64" or checkpoint.get("rng_provenance", {}).get("child_order") != list(RNG_NAMES):
        raise GateError("checkpoint RNG provenance differs from the frozen three-stream rule")
    actual_code_sha = {name: sha256_file(PACKAGE_DIR / name) for name in config["code_sha256"]}
    code_pass = actual_code_sha == config["code_sha256"]
    if not code_pass:
        raise GateError("runner, criteria, archive backend, or test code SHA differs from config")
    env = _environment_snapshot()
    env_sha = canonical_sha(env)
    expected_env = config["runtime"]["expected_environment"]
    env_compare = {k: {"expected": expected_env[k], "actual": env.get(k), "pass": env.get(k) == expected_env[k]} for k in expected_env}
    env_pass = all(v["pass"] for v in env_compare.values()) and env_sha == EXPECTED_ENVIRONMENT_SHA256
    potential_path = Path(config["potential"]["path"])
    potential_sha = sha256_file(potential_path) if potential_path.is_file() else None
    potential_pass = potential_sha == EXPECTED_POTENTIAL_SHA256 and potential_path.stat().st_size == config["potential"]["bytes"] if potential_sha else False
    resources = _resource_snapshot(Path("/workspace"))
    resource_pass = _resource_gate(resources, config["limits"])
    return {
        "schema_version": 1,
        "owner": "D",
        "task_id": TASK_ID,
        "config_sha256": config_sha,
        "code_sha256": actual_code_sha,
        "code_sha_pass": code_pass,
        "strategy_lock_sha256": STRATEGY_LOCK_SHA256,
        "scientific_lock_sha256": SCIENTIFIC_LOCK_SHA256,
        "checkpoint_step": checkpoint["global_step"],
        "checkpoint_sha256": sha256_file(checkpoint_path),
        "checkpoint_position_momentum_rng_restore_pass": True,
        "checkpoint_member_sha256_pass": True,
        "potential_path": str(potential_path),
        "potential_sha256": potential_sha,
        "potential_sha_pass": potential_pass,
        "resources": resources,
        "resource_preflight_pass": resource_pass,
        "environment": env,
        "environment_sha256": env_sha,
        "environment_expected_sha256": EXPECTED_ENVIRONMENT_SHA256,
        "environment_sha_pass": env_pass,
        "environment_field_comparison": env_compare,
        "execution_authorized": False,
        "preflight_status": "PREPARED_HOLD" if potential_pass and env_pass and resource_pass else "BLOCKED",
        "formal_md_started": False,
        "formal_steps": 0,
        "restored_rng_streams": sorted(restored_rngs),
        "current_task_only": True,
        "notes": [
            "Read-only check only; no calculator, force, energy, dynamics, or output directory was opened.",
            "A C execution release bound to this config/code/environment is still required.",
            "The frozen strategy's qualitative RDF shape rule lacks a numerical pass threshold; this pack fails closed pending C operationalization.",
        ],
    }


def _authorization_check(auth_path: Path, auth_sha: str, config_sha: str, config: dict, repo: Path, release_commit: str, release_repo_path: str) -> dict:
    _verify_release_blob(repo, release_commit, release_repo_path, auth_sha, auth_path)
    auth = _verify_json_sha(auth_path, auth_sha, "C execution release")
    required = {
        "execution_authorized": True,
        "owner": "D",
        "task_id": TASK_ID,
        "strategy_lock_sha256": STRATEGY_LOCK_SHA256,
        "scientific_lock_sha256": SCIENTIFIC_LOCK_SHA256,
        "config_sha256": config_sha,
        "runner_sha256": config["code_sha256"]["followup_runner.py"],
        "criteria_sha256": config["code_sha256"]["liquid_criteria.py"],
        "archive_backend_sha256": config["code_sha256"]["archive_backend.py"],
        "test_suite_sha256": config["code_sha256"]["test_execution_pack.py"],
        "checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
        "potential_sha256": EXPECTED_POTENTIAL_SHA256,
        "environment_sha256": EXPECTED_ENVIRONMENT_SHA256,
        "maximum_global_step": 60000,
        "maximum_added_steps": 59000,
        "maximum_chunk_steps": 1000,
        "maximum_chunk_wall_seconds": 3600,
        "d_branch_head": subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip(),
    }
    for key, value in required.items():
        if auth.get(key) != value:
            raise GateError(f"C execution release binding mismatch: {key}")
    policy = auth.get("rdf_shape_policy")
    if not isinstance(policy, dict) or canonical_sha(policy) != auth.get("rdf_shape_policy_sha256"):
        raise GateError("C release must operationalize and hash the frozen qualitative RDF-shape rule")
    if not {
        "minimum_first_peak_fwhm_A",
        "secondary_search_start_A",
        "maximum_secondary_peak_prominence_fraction",
    }.issubset(policy):
        raise GateError("C-signed RDF shape policy is missing required thresholds")
    return auth


def _atomic_json(path: Path, value: dict) -> None:
    temp = path.with_name(path.name + f".tmp-{os.getpid()}")
    with temp.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temp, path)


def _csv_append(path: Path, header: list[str], row: list[object]) -> None:
    exists = path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        if not exists:
            writer.writerow(header)
        writer.writerow(row)
        stream.flush()
        os.fsync(stream.fileno())


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _load_frame_metrics(directory: Path) -> list[criteria.FrameMetrics]:
    result = []
    path = directory / "frame_metrics.jsonl"
    if not path.is_file():
        return result
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        item = json.loads(line)
        for key in ("rdf_r_A", "rdf_g", "stress_eV_A3"):
            item[key] = tuple(item[key])
        result.append(criteria.FrameMetrics(**item))
    return result


def _minimum_distance(atoms: Atoms) -> float:
    from ase.neighborlist import neighbor_list
    distances = neighbor_list("d", atoms, 1.8, self_interaction=False)
    return float(np.min(distances)) if len(distances) else math.inf


def _save_checkpoint(path: Path, atoms: Atoms, rngs: dict, state: dict, config: dict, start_wall: float) -> str:
    payload = {
        "schema_version": 2,
        "task_id": TASK_ID,
        "source_step1000_checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
        "global_step": int(state["global_step"]),
        "time_fs": float(state["global_step"] * 0.5),
        "symbols": atoms.get_chemical_symbols(),
        "atom_ids": atoms.arrays["atom_id"].astype(int).tolist(),
        "positions_A": atoms.get_positions().tolist(),
        "momenta_amu_A_per_fs": (atoms.get_momenta() * fs).tolist(),
        "velocities_A_per_fs": (atoms.get_velocities() * fs).tolist(),
        "velocities_ASE_internal_A_per_time_unit": atoms.get_velocities().tolist(),
        "masses_amu": atoms.get_masses().tolist(),
        "cell_A": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(),
        "rng_states": {name: copy.deepcopy(rngs[name].bit_generator.state) for name in RNG_NAMES},
        "rng_provenance": {
            "algorithm": "numpy.SeedSequence(master_seed).spawn(3)",
            "bit_generator": "PCG64",
            "master_seed": 125010003,
            "child_order": list(RNG_NAMES),
        },
        "temperature_schedule_position": {
            "last_completed_step": int(state["global_step"]),
            "stage": state["stage"],
            "next_integrator_step_target_K": (
                1600.0 if state["stage"] == "A_HOT_1600K" else
                cooling_target_K(int(state["global_step"]) - int(state["stage_b_start_step"])) if state["stage"] == "B_COOL_1600_TO_1250K" else
                1250.0
            ),
        },
        "lock_sha256": SCIENTIFIC_LOCK_SHA256,
        "parameter_sha256": "86ea1bb9ad269d0cd3aa626a9b01e39c3bd0cd5ac5c64c951c5d44e2deea22e4",
        "strategy_lock_sha256": STRATEGY_LOCK_SHA256,
        "potential_sha256": EXPECTED_POTENTIAL_SHA256,
        "runner_sha256": config["code_sha256"]["followup_runner.py"],
        "criteria_sha256": config["code_sha256"]["liquid_criteria.py"],
        "environment_sha256": EXPECTED_ENVIRONMENT_SHA256,
        "ase_version": ase.__version__,
        "python_version": ".".join(map(str, sys.version_info[:3])),
        "stage_state": state,
        "trajectory": {"path": "trajectory.traj", "sha256": sha256_file(path.parent / "trajectory.traj")},
        "thermal_csv": {"path": "thermal.csv", "sha256": sha256_file(path.parent / "thermal.csv")},
        "energy_csv": {"path": "energy.csv", "sha256": sha256_file(path.parent / "energy.csv")},
        "stress_csv": {"path": "stress.csv", "sha256": sha256_file(path.parent / "stress.csv")},
        "unwrapped_positions": {"path": "unwrapped_positions.jsonl", "sha256": sha256_file(path.parent / "unwrapped_positions.jsonl")},
        "elapsed_wall_s": time.monotonic() - start_wall,
    }
    _atomic_json(path, payload)
    digest = sha256_file(path)
    path.with_suffix(path.suffix + ".sha256").write_text(f"{digest}  {path.name}\n", encoding="ascii")
    return digest


def _bundle_and_archive(repo: Path, out: Path, global_step: int, config: dict, auth: dict) -> dict:
    chunk_start = max(1000, global_step - 1000)
    with tempfile.TemporaryDirectory(prefix="ag64-rho100-bundle-") as temp_name:
        bundle = Path(temp_name) / "bundle"
        shutil.copytree(out, bundle)
        bundle_files = [p for p in bundle.rglob("*") if p.is_file()]
        if any(p.stat().st_size > int(config["limits"]["maximum_single_file_bytes"]) for p in bundle_files):
            raise GateError("archive single-file size cap exceeded; retained local checkpoint")
        if sum(p.stat().st_size for p in bundle_files) > int(config["limits"]["maximum_bundle_bytes"]):
            raise GateError("archive bundle size cap exceeded; retained local checkpoint")
        (bundle / "SEGMENT_MANIFEST.json").write_text(json.dumps({
            "task_id": TASK_ID,
            "source_family_id": SOURCE_FAMILY_ID,
            "step_start_exclusive": chunk_start,
            "step_end_inclusive": global_step,
            "strategy_lock_sha256": STRATEGY_LOCK_SHA256,
            "scientific_lock_sha256": SCIENTIFIC_LOCK_SHA256,
            "execution_release_sha256": auth["_release_sha256"],
            "checkpoint_sha256": sha256_file(bundle / "checkpoint.json"),
            "potential_sha256": EXPECTED_POTENTIAL_SHA256,
            "runner_sha256": config["code_sha256"]["followup_runner.py"],
            "criteria_sha256": config["code_sha256"]["liquid_criteria.py"],
        }, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        archive_relpath = f"{config['archive']['root']}/chunks/step_{global_step:05d}"
        restore_dir = Path(tempfile.gettempdir()) / f"ag64-rho100-restored-{global_step:05d}-{os.getpid()}"
        head = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip()
        try:
            return archive_backend.archive_push_fresh_restore(
                repo_root=repo,
                source_dir=bundle,
                archive_relpath=archive_relpath,
                branch=config["archive"]["branch"],
                remote=config["archive"]["remote"],
                expected_remote_head=head,
                restore_dir=restore_dir,
                commit_message=f"Archive Ag64 rho100 verified segment through step {global_step}",
            )
        finally:
            if restore_dir.exists():
                shutil.rmtree(restore_dir)


def _run_loop(config_path: Path, config_sha: str, config: dict, auth: dict, atoms: Atoms, rngs: dict, state: dict, output_dir: Path, repo: Path) -> dict:
    global _STOP_REQUESTED
    _STOP_REQUESTED = False
    potential = Path(config["potential"]["path"])
    if sha256_file(potential) != EXPECTED_POTENTIAL_SHA256:
        raise GateError("potential changed after preflight")
    if output_dir.exists():
        raise GateError("run output directory already exists; use only a separately verified resume archive")
    output_dir.mkdir(parents=True)
    source = repo / config["source_checkpoint"]["archive_path"]
    for name in ("checkpoint.json", "trajectory.traj", "thermal.csv", "energy.csv"):
        shutil.copy2(source / name, output_dir / name)
    (output_dir / "stress.csv").write_text("step,time_fs,xx_eV_A3,yy_eV_A3,zz_eV_A3,yz_eV_A3,xz_eV_A3,xy_eV_A3\n", encoding="utf-8")
    (output_dir / "unwrapped_positions.jsonl").write_text("", encoding="utf-8")
    (output_dir / "frame_metrics.jsonl").write_text("", encoding="utf-8")
    calculator = EAM(potential=str(potential), form="alloy")
    atoms.calc = calculator
    dynamics = Langevin(
        atoms, timestep=0.5 * fs, temperature_K=1600.0,
        friction=0.01 / fs, rng=rngs["thermostat"], fixcm=False,
    )
    start_wall = time.monotonic()
    last_checkpoint_sha = sha256_file(output_dir / "checkpoint.json")
    candidate_start = None
    current_step = int(state["global_step"])
    target_K = 1600.0
    sample_interval = 20
    trajectory_interval = 100
    chunk_steps = 1000
    maximum_global = min(60000, current_step + 59000)
    status = "RUNNING"
    previous_handlers = {
        signal.SIGTERM: signal.signal(signal.SIGTERM, _request_graceful_stop),
        signal.SIGINT: signal.signal(signal.SIGINT, _request_graceful_stop),
    }
    trajectory = Trajectory(output_dir / "trajectory.traj", "a", properties=[])
    archive_result = None
    try:
        while current_step < maximum_global:
            chunk_start_wall = time.monotonic()
            try:
                live_preflight = run_preflight(config_path, config_sha)
                if live_preflight["preflight_status"] != "PREPARED_HOLD":
                    stop_reason = "LIVE_SEGMENT_PREFLIGHT_BLOCKED"
                    state["live_preflight"] = live_preflight
                else:
                    stop_reason = None
            except BaseException as exc:
                stop_reason = "LIVE_SEGMENT_PREFLIGHT_ERROR"
                state["live_preflight_error"] = f"{type(exc).__name__}: {exc}"
            chunk_resources = _resource_snapshot(output_dir)
            if stop_reason is None and not _resource_gate(chunk_resources, config["limits"]):
                stop_reason = "RESOURCE_PREFLIGHT_BLOCKED"
                state["resource_preflight"] = chunk_resources
                segment_end = current_step
            if state["stage"] == "A_HOT_1600K":
                target_K = 1600.0
                phase_boundary = min(20000, current_step + chunk_steps)
            elif state["stage"] == "B_COOL_1600_TO_1250K":
                ramp_start = int(state["stage_b_start_step"])
                ramp_end = ramp_start + 10000
                completed = current_step - ramp_start
                target_K = cooling_target_K(completed)
                dynamics.set_temperature(temperature_K=target_K)
                phase_boundary = min(ramp_end, current_step + chunk_steps)
            elif state["stage"] in {"C_HOLD_1250K", "D_CANDIDATE_1250K"}:
                target_K = 1250.0
                dynamics.set_temperature(temperature_K=target_K)
                phase_start = int(state["stage_c_start_step"] if state["stage"] == "C_HOLD_1250K" else state["stage_d_start_step"])
                phase_max = min(60000, phase_start + (20000 if state["stage"] == "C_HOLD_1250K" else 10000))
                phase_boundary = min(phase_max, current_step + chunk_steps)
            else:
                status = str(state["stage"])
                break

            # Finish one currently started integration step, then check stop
            # conditions before starting the next one.
            if stop_reason is None:
                segment_end = min(next_chunk_end(current_step, phase_boundary), maximum_global)
            while current_step < segment_end:
                if _STOP_REQUESTED:
                    stop_reason = "GRACEFUL_SIGNAL_STOP"
                    break
                elapsed = time.monotonic() - chunk_start_wall
                if elapsed >= 3600.0:
                    stop_reason = "CHUNK_WALLTIME_STOP"
                    break
                if current_step >= 60000 or current_step >= maximum_global:
                    stop_reason = "GLOBAL_STEP_LIMIT"
                    break
                if state["stage"] == "B_COOL_1600_TO_1250K":
                    ramp_offset = current_step - int(state["stage_b_start_step"])
                    if ramp_offset > 0 and ramp_offset % 100 == 0:
                        target_K = cooling_target_K(ramp_offset)
                        dynamics.set_temperature(temperature_K=target_K)
                old_positions = atoms.get_positions().copy()
                old_momenta = atoms.get_momenta().copy()
                old_rng_state = copy.deepcopy(rngs["thermostat"].bit_generator.state)
                try:
                    dynamics.run(1)
                except BaseException as exc:
                    atoms.set_positions(old_positions, apply_constraint=False)
                    atoms.set_momenta(old_momenta, apply_constraint=False)
                    rngs["thermostat"].bit_generator.state = old_rng_state
                    stop_reason = f"INTEGRATOR_EXCEPTION_{type(exc).__name__}"
                    state["failure"] = f"{type(exc).__name__}: {exc}"
                    break
                current_step += 1
                state["global_step"] = current_step
                # ASE positions are not wrapped by this runner; carry explicit
                # unwrapped coordinates for multiple-origin MSD after restart.
                delta = atoms.get_positions() - old_positions
                cell = np.asarray(atoms.cell)
                frac = delta @ np.linalg.inv(cell)
                frac -= np.rint(frac)
                state["unwrapped_positions_A"] = (np.asarray(state["unwrapped_positions_A"]) + frac @ cell).tolist()
                if current_step % sample_interval == 0:
                    temperature = float(atoms.get_temperature())
                    potential_e = float(atoms.get_potential_energy())
                    kinetic_e = float(atoms.get_kinetic_energy())
                    stress = np.asarray(atoms.get_stress(include_ideal_gas=True), dtype=float)
                    forces = np.asarray(atoms.get_forces(), dtype=float)
                    if not all(math.isfinite(x) for x in (temperature, potential_e, kinetic_e)) or not np.isfinite(stress).all() or not np.isfinite(forces).all():
                        atoms.set_positions(old_positions, apply_constraint=False)
                        atoms.set_momenta(old_momenta, apply_constraint=False)
                        rngs["thermostat"].bit_generator.state = old_rng_state
                        current_step -= 1
                        state["global_step"] = current_step
                        stop_reason = "NONFINITE_STATE"
                        break
                    if _minimum_distance(atoms) < 1.8:
                        stop_reason = "MINIMUM_PAIR_DISTANCE_BELOW_1P8_A"
                        break
                    _csv_append(output_dir / "thermal.csv", ["step", "time_fs", "target_temperature_K", "temperature_K", "potential_energy_eV", "kinetic_energy_eV", "total_energy_eV"], [current_step, current_step * 0.5, target_K, temperature, potential_e, kinetic_e, potential_e + kinetic_e])
                    _csv_append(output_dir / "energy.csv", ["step", "time_fs", "potential_energy_eV", "kinetic_energy_eV", "total_energy_eV"], [current_step, current_step * 0.5, potential_e, kinetic_e, potential_e + kinetic_e])
                    _csv_append(output_dir / "stress.csv", ["step", "time_fs", "xx_eV_A3", "yy_eV_A3", "zz_eV_A3", "yz_eV_A3", "xz_eV_A3", "xy_eV_A3"], [current_step, current_step * 0.5, *stress.tolist()])
                if current_step % trajectory_interval == 0:
                    trajectory.write(atoms)
                    with (output_dir / "unwrapped_positions.jsonl").open("a", encoding="utf-8") as stream:
                        stream.write(json.dumps({"step": current_step, "positions_A": state["unwrapped_positions_A"]}, allow_nan=False) + "\n")
                    frame = criteria.frame_metrics(
                        atoms, current_step, float(atoms.get_temperature()),
                        float(atoms.get_potential_energy()),
                        atoms.get_stress(include_ideal_gas=True),
                    )
                    with (output_dir / "frame_metrics.jsonl").open("a", encoding="utf-8") as stream:
                        stream.write(json.dumps(criteria.monitor_record_to_json(frame), allow_nan=False) + "\n")
                    if state["stage"] == "D_CANDIDATE_1250K":
                        spacing = int(state["candidate_snapshot_spacing_steps"])
                        if (current_step - int(state["stage_d_start_step"])) % spacing == 0:
                            write(output_dir / "candidate_snapshots.extxyz", atoms, format="extxyz", append=True)
                            write(output_dir / "candidate_structure.extxyz", atoms, format="extxyz")
                            state["candidate_snapshot_steps"] = state.get("candidate_snapshot_steps", []) + [current_step]
                free = shutil.disk_usage(output_dir).free
                rss = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024
                if free < int(config["limits"]["minimum_free_bytes"]) + int(config["limits"]["checkpoint_reserve_bytes"]):
                    stop_reason = "DISK_RESOURCE_LIMIT"
                    break
                if rss > int(config["limits"]["maximum_rss_bytes"]):
                    stop_reason = "RSS_RESOURCE_LIMIT"
                    break
                if stop_reason:
                    break
            if stop_reason is None and time.monotonic() - chunk_start_wall >= 3600.0:
                # A step already in progress is allowed to finish, but this
                # chunk may not begin a new segment after crossing its wall cap.
                stop_reason = "CHUNK_WALLTIME_STOP"
            trajectory.close()
            trajectory = None
            state["global_step"] = current_step
            last_checkpoint_sha = _save_checkpoint(output_dir / "checkpoint.json", atoms, rngs, state, config, start_wall)
            if stop_reason:
                status = stop_reason
            else:
                status = "CHUNK_COMPLETE"
            if status == "CHUNK_COMPLETE":
                advance_completed_phase(state, current_step)
            monitor = None
            candidate_spacing = None
            if status == "CHUNK_COMPLETE" and current_step % 1000 == 0 and state["stage"] in {"A_HOT_1600K", "B_COOL_1600_TO_1250K", "C_HOLD_1250K", "D_CANDIDATE_1250K"}:
                state["last_monitor_step"] = current_step
                monitor = _monitor_at_boundary(output_dir, current_step, state, config, auth)
                state["last_monitor"] = monitor
                _atomic_json(output_dir / "QA_monitor_latest.json", monitor)
                if monitor.get("gate_status") == "PASS" and state["stage"] == "C_HOLD_1250K":
                    hold_start = max(int(state["stage_c_start_step"]), current_step - 8000)
                    metrics_for_tau = [m for m in _load_frame_metrics(output_dir) if hold_start < m.step <= current_step]
                    thermal_for_tau = [float(r["temperature_K"]) for r in _read_csv(output_dir / "thermal.csv") if hold_start < int(r["step"]) <= current_step and int(r["step"]) % 100 == 0]
                    energy_for_tau = [m.potential_energy_eV_per_atom for m in metrics_for_tau]
                    candidate_spacing = criteria.candidate_snapshot_spacing_steps(metrics_for_tau, thermal_for_tau, energy_for_tau)
            if status == "CHUNK_COMPLETE" and monitor is not None:
                status = apply_monitor_transition(state, current_step, monitor, candidate_spacing)
            status = apply_completed_chunk_limits(state, current_step, status, maximum_global)

            # Persist the post-QA phase decision and all monitor state in the
            # checkpoint that is about to be archived.
            state["global_step"] = current_step
            last_checkpoint_sha = _save_checkpoint(output_dir / "checkpoint.json", atoms, rngs, state, config, start_wall)
            _atomic_json(output_dir / "RUN_STATE.json", {**state, "status": status})
            _write_result_and_sums(output_dir, state, status, time.monotonic() - start_wall, last_checkpoint_sha)
            try:
                archive_result = _bundle_and_archive(repo, output_dir, current_step, config, auth)
            except BaseException as exc:
                status = "ARCHIVE_OR_RESTORE_FAILURE"
                state["stage"] = status
                state["archive_failure"] = f"{type(exc).__name__}: {exc}"
                _atomic_json(output_dir / "RUN_STATE.json", {**state, "status": status})
                _write_result_and_sums(output_dir, state, status, time.monotonic() - start_wall, last_checkpoint_sha)
                break
            if not archive_backend.segment_archive_passed(archive_result):
                status = "ARCHIVE_OR_RESTORE_FAILURE"
                state["stage"] = status
                state["archive_failure"] = "segment archive result did not pass every required gate"
                _atomic_json(output_dir / "RUN_STATE.json", {**state, "status": status})
                _write_result_and_sums(output_dir, state, status, time.monotonic() - start_wall, last_checkpoint_sha)
                break
            if status != "CHUNK_COMPLETE":
                break
            # Continue automatically only after upload, fresh fetch, restore,
            # and per-member SHA checks have all passed.
            trajectory = Trajectory(output_dir / "trajectory.traj", "a", properties=[])
    except BaseException as exc:
        # Preserve the last completed integrator boundary if analysis or I/O
        # fails after a successful step.  The failed bundle is retained locally.
        if trajectory is not None:
            trajectory.close()
            trajectory = None
        state["global_step"] = current_step
        state["failure"] = f"{type(exc).__name__}: {exc}"
        status = "POST_STEP_OR_CHECKPOINT_FAILURE"
        state["stage"] = status
        try:
            last_checkpoint_sha = _save_checkpoint(output_dir / "checkpoint.json", atoms, rngs, state, config, start_wall)
            _atomic_json(output_dir / "RUN_STATE.json", {**state, "status": status})
            _write_result_and_sums(output_dir, state, status, time.monotonic() - start_wall, last_checkpoint_sha)
        except BaseException as save_exc:
            state["checkpoint_save_failure"] = f"{type(save_exc).__name__}: {save_exc}"
            status = "POST_STEP_AND_CHECKPOINT_SAVE_FAILED"
    finally:
        if trajectory is not None:
            trajectory.close()
        for signum, previous in previous_handlers.items():
            signal.signal(signum, previous)
    return {
        "status": status,
        "actual_start_step": 1000,
        "actual_end_step": current_step,
        "time_ps": current_step * 0.0005,
        "checkpoint_sha256": last_checkpoint_sha,
        "archive": archive_result,
    }


def _monitor_at_boundary(out: Path, step: int, state: dict, config: dict, auth: dict) -> dict:
    metrics: list[criteria.FrameMetrics] = []
    metrics_path = out / "frame_metrics.jsonl"
    if metrics_path.is_file():
        for line in metrics_path.read_text(encoding="utf-8").splitlines():
            if line.strip():
                row = json.loads(line)
                for key in ("rdf_r_A", "rdf_g", "stress_eV_A3"):
                    row[key] = tuple(row[key])
                metrics.append(criteria.FrameMetrics(**row))
    recent = [m for m in metrics if step - 4000 < int(m.step) <= step]
    thermal = _read_csv(out / "thermal.csv")
    recent_rows = [r for r in thermal if step - 4000 < int(r["step"]) <= step]
    recent_t = [float(r["temperature_K"]) for r in recent_rows]
    if recent:
        mean_q6 = float(np.mean([m.global_q6 for m in recent]))
        crystal = float(np.mean([m.crystalline_fraction for m in recent]))
        coordination = float(np.mean([m.coordination_number for m in recent]))
    else:
        mean_q6 = crystal = coordination = math.nan
    monitor = {
        "step": step,
        "temperature_mean_recent_2ps_K": float(np.mean(recent_t)) if recent_t else None,
        "temperature_std_recent_2ps_K": float(np.std(recent_t, ddof=1)) if len(recent_t) > 1 else None,
        "global_Q6_mean_recent_2ps": mean_q6 if math.isfinite(mean_q6) else None,
        "crystalline_fraction_mean_recent_2ps": crystal if math.isfinite(crystal) else None,
        "MSD_A2": None,
        "coordination_mean_recent_2ps": coordination if math.isfinite(coordination) else None,
        "gate_status": "MONITOR_ONLY",
        "reason": None,
    }
    stage = state["stage"]
    if stage == "A_HOT_1600K":
        discard_end = int(state["stage_a_initial_transient_end_step"])
        target = 1600.0
    elif stage == "C_HOLD_1250K":
        discard_end = int(state["stage_c_start_step"]) + 4000
        target = 1250.0
    elif stage == "D_CANDIDATE_1250K":
        discard_end = int(state["stage_d_start_step"])
        target = 1250.0
    else:
        monitor["gate_status"] = "MONITOR_ONLY_STAGE_B_RAMP"
        return monitor
    if step < discard_end + 8000:
        monitor["gate_status"] = "WAITING_FOR_TWO_COMPLETE_2PS_WINDOWS"
        return monitor
    pair_end = discard_end + ((step - discard_end) // 4000) * 4000
    first_start, first_end = pair_end - 8000, pair_end - 4000
    second_start, second_end = first_end, pair_end
    thermal_windows = []
    energy_rows = _read_csv(out / "energy.csv")
    structure_windows = []
    for wstart, wend in ((first_start, first_end), (second_start, second_end)):
        tr = [r for r in thermal if wstart < int(r["step"]) <= wend]
        er = [r for r in energy_rows if wstart < int(r["step"]) <= wend]
        fr = [f for f in metrics if wstart < int(f.step) <= wend]
        if not tr or not er or not fr:
            monitor["gate_status"] = "WAITING_FOR_COMPLETE_WINDOW_DATA"
            return monitor
        thermal_windows.append({
            "temperature_K": [float(r["temperature_K"]) for r in tr],
            "time_ps": [int(r["step"]) * 0.0005 for r in tr],
            "potential_energy_eV_per_atom": [float(r["potential_energy_eV"]) / 64.0 for r in er],
        })
        structure_windows.append((fr, criteria.structure_window_summary(
            fr,
            config["criteria_thresholds"],
            rdf_policy=auth["rdf_shape_policy"],
        )))
    thermal_gate = criteria.paired_window_thermal_gate(
        thermal_windows[0], thermal_windows[1], config["thermal_thresholds"], target,
    )
    first_struct, second_struct = structure_windows[0][1], structure_windows[1][1]
    unwrapped_path = out / "unwrapped_positions.jsonl"
    unwrapped = []
    if unwrapped_path.is_file():
        for line in unwrapped_path.read_text(encoding="utf-8").splitlines():
            if line.strip():
                item = json.loads(line)
                if first_start <= int(item["step"]) <= second_end:
                    unwrapped.append((int(item["step"]), item["positions_A"]))
    unwrapped.sort(key=lambda x: x[0])
    if len(unwrapped) < 8:
        monitor["gate_status"] = "WAITING_FOR_MSD_DATA"
        return monitor
    msd = criteria.msd_summary([item[1] for item in unwrapped], frame_interval_ps=0.05)
    monitor["MSD_A2"] = msd["delta_msd_A2_over_full_interval"]
    all_pair_frames = [f for f in metrics if first_start < int(f.step) <= second_end]
    stress = criteria.stress_drift_gate(
        all_pair_frames,
        [f.step * 0.0005 for f in all_pair_frames],
    )
    monitor.update({
        "window_step_bounds": [[first_start, first_end], [second_start, second_end]],
        "thermal_gate": thermal_gate,
        "structure_gate_first": first_struct,
        "structure_gate_second": second_struct,
        "msd_gate": msd,
        "stress_drift_gate": stress,
    })
    gate_pass = (
        thermal_gate["temperature_gate_pass"]
        and thermal_gate["energy_gate_pass"]
        and first_struct["automatic_structure_gate_pass"]
        and second_struct["automatic_structure_gate_pass"]
        and msd["msd_gate_pass"]
        and stress["stress_drift_pass"]
    )
    monitor["gate_status"] = "PASS" if gate_pass else "FAIL_WINDOW"
    return monitor


def _write_result_and_sums(out: Path, state: dict, status: str, wall_s: float, checkpoint_sha: str) -> None:
    result = {
        "task_id": TASK_ID,
        "source_family_id": SOURCE_FAMILY_ID,
        "status": status,
        "global_step": int(state["global_step"]),
        "time_ps": float(state["global_step"] * 0.0005),
        "wall_seconds": float(wall_s),
        "checkpoint_sha256": checkpoint_sha,
        "liquid_onset_step": state.get("liquid_onset_step"),
        "execution_authorized": True,
    }
    _atomic_json(out / "RUN_RESULT.json", result)
    members = {
        p.name: sha256_file(p)
        for p in sorted(out.iterdir())
        if p.is_file() and p.name not in {"MANIFEST.json", "SHA256SUMS.txt"}
    }
    _atomic_json(out / "MANIFEST.json", {
        "owner": "D",
        "task_id": TASK_ID,
        "source_family_id": SOURCE_FAMILY_ID,
        "status": status,
        "global_step": int(state["global_step"]),
        "time_ps": float(state["global_step"] * 0.0005),
        "liquid_onset_step": state.get("liquid_onset_step"),
        "scientific_lock_sha256": SCIENTIFIC_LOCK_SHA256,
        "strategy_lock_sha256": STRATEGY_LOCK_SHA256,
        "source_step1000_checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
        "potential_sha256": EXPECTED_POTENTIAL_SHA256,
        "runner_sha256": state.get("runner_sha256"),
        "criteria_sha256": state.get("criteria_sha256"),
        "environment_sha256": EXPECTED_ENVIRONMENT_SHA256,
        "data_role": "SINGLE_SOURCE_FAMILY_SAMPLING_CANDIDATE_ONLY_NO_TRAIN_DEV_TEST",
        "member_sha256": members,
    })
    lines = []
    for path in sorted(out.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            lines.append(f"{sha256_file(path)}  {path.name}\n")
    (out / "SHA256SUMS.txt").write_text("".join(lines), encoding="ascii")


def run_authorized(config_path: Path, config_sha: str, authorization_path: Path, authorization_sha: str) -> dict:
    report = run_preflight(config_path, config_sha)
    if report["preflight_status"] != "PREPARED_HOLD":
        raise GateError("live preflight did not pass")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    auth = _authorization_check(authorization_path, authorization_sha, config_sha, config)
    auth["_release_sha256"] = authorization_sha
    repo = _git_root(PACKAGE_DIR)
    source = repo / config["source_checkpoint"]["archive_path"]
    checkpoint = json.loads((source / "checkpoint.json").read_text(encoding="utf-8"))
    atoms, rngs = _restore_atoms(checkpoint)
    state = {
        "stage": "A_HOT_1600K",
        "global_step": 1000,
        "stage_a_start_step": 1000,
        "stage_a_initial_transient_end_step": 5000,
        "liquid_onset_step": None,
        "stage_b_start_step": None,
        "stage_c_start_step": None,
        "stage_d_start_step": None,
        "unwrapped_positions_A": atoms.get_positions().tolist(),
        "qa_msd_origin_positions_A": None,
        "last_monitor_step": None,
        "last_monitor": None,
    }
    out = Path(config["output_directory"])
    if out.exists():
        raise GateError("output directory exists; no overwrite/reinitialization allowed")
    return _run_loop(config_path, config_sha, config, auth, atoms, rngs, state, out, repo)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=PACKAGE_DIR / CONFIG_NAME)
    parser.add_argument("--config-sha256", required=True)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--authorization", type=Path)
    parser.add_argument("--authorization-sha256")
    args = parser.parse_args()
    if args.preflight == args.run:
        parser.error("choose exactly one of --preflight or --run")
    try:
        if args.preflight:
            result = run_preflight(args.config, args.config_sha256)
            print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
            return 0 if result["preflight_status"] == "PREPARED_HOLD" else 2
        if not args.authorization or not args.authorization_sha256:
            raise GateError("--run requires a C-signed authorization file and its SHA256")
        result = run_authorized(args.config, args.config_sha256, args.authorization, args.authorization_sha256)
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
        return 0
    except Exception as exc:
        print(json.dumps({"status": "BLOCKED", "error": f"{type(exc).__name__}: {exc}"}, sort_keys=True), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
