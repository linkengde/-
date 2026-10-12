#!/usr/bin/env python3
"""Isolated Ag64 rho100 continuation from the immutable step-21000 archive.

This entry supports read-only preflight and a separately bounded user-approved
continuation. It never creates or reinitializes an Ag structure. This task runs
only tests and preflight, not MD.
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

import archive_backend_v2 as archive_backend
import liquid_criteria_v2 as criteria


PACKAGE_DIR = Path(__file__).resolve().parent
CONFIG_NAME = "FOLLOWUP_EXECUTION_CONFIG_STEP21000_60000.json"
RUNNER_NAME = "followup_runner_step21000_extension.py"
SCIENTIFIC_LOCK_SHA256 = "e1adaf12d542184c709d2eb2a5363cbfe9578c168ebc5e7ad260ee39a4faf09b"
STRATEGY_LOCK_SHA256 = "054f36dcda113e37b23b112e5919ba1c8201c669626db37b79bbd2281b020e1c"
QA_SCIENCE_LOCK_SHA256 = "848f8db72d83cea91dcf4cd42ee9fbbf49a80e7b5e3759cb4f535be8cc5873d2"
EXPECTED_CHECKPOINT_SHA256 = "c820437a2512c4b3d7dc9e60bb60ac3ced6588aa1f2a629d58cfdabb24f34046"
EXPECTED_POTENTIAL_SHA256 = "958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447"
SOURCE_CHECKPOINT_ENVIRONMENT_SHA256 = "f3f6d27f7e8b387bb5cc7b91048dba7fa26280f9ed54aa1572c73879946e867f"
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
    payload = json.dumps(_json_safe(value), sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
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


def can_start_step(current_step: int, segment_end: int, elapsed_wall_s: float,
                   maximum_global_step: int = 60000, maximum_wall_s: float = 3600.0) -> bool:
    """Fail closed before a new integration step at any signed boundary."""
    return bool(
        current_step < segment_end
        and current_step < maximum_global_step
        and math.isfinite(elapsed_wall_s)
        and 0.0 <= elapsed_wall_s < maximum_wall_s
    )


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
    elif gate == "QA_HOLD":
        status = "HOLD_QA_REQUIRED_INPUT_OR_METRIC"
        state["stage"] = status
    elif state["stage"] == "D_CANDIDATE_1250K" and gate == "FAIL_WINDOW":
        status = "HOLD_CANDIDATE_WINDOW_QA_FAILED"
        state["stage"] = status
    return status


def apply_completed_chunk_limits(state: dict, current_step: int, status: str, maximum_global: int,
                                 stage_a_maximum_global: int = 40000) -> str:
    """Apply the bounded Stage A extension and global stop gates."""
    if status != "CHUNK_COMPLETE":
        return status
    stage = state["stage"]
    if stage == "A_HOT_1600K" and current_step >= stage_a_maximum_global:
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


def _archive_allows_continuation(status: str, archive_result: dict | None) -> bool:
    """The production loop may begin another chunk only after every archive gate passes."""
    return status == "CHUNK_COMPLETE" and archive_backend.segment_archive_passed(archive_result)


def _should_continue_after_archive(status: str, archive_result: dict | None,
                                   stop_requested: bool) -> bool:
    """Do not revisit an immutable boundary after a signal arrived post-archive."""
    return not stop_requested and _archive_allows_continuation(status, archive_result)


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
    import importlib.metadata
    runtime = Path("/workspace/onboarding/runtime.sh")
    if not runtime.is_file():
        raise GateError("runtime.sh is missing")
    runtime_text = runtime.read_text(encoding="utf-8")
    if "export GPAW_MPI_BACKEND=cgpaw" not in runtime_text or "gpaw_data.datapath()" not in runtime_text:
        raise GateError("runtime.sh does not define the locked GPAW backend/setup path")
    mpi = subprocess.run(["mpiexec", "--version"], check=True, capture_output=True, text=True, timeout=10).stdout.splitlines()[0].strip()
    setup_path = str(gpaw_data.datapath())
    return {
        "python_version": ".".join(map(str, sys.version_info[:3])),
        "python_executable": sys.executable,
        "ase_version": ase.__version__,
        "numpy_version": np.__version__,
        "scipy_version": __import__("scipy").__version__,
        "ovito_version": importlib.metadata.version("ovito"),
        "runtime_sh": str(runtime),
        "runtime_sh_sha256": sha256_file(runtime),
        "GPAW_SETUP_PATH": setup_path,
        "GPAW_MPI_BACKEND": "cgpaw",
        "openmpi_version": mpi,
        "host_os": platform.platform(),
        "gpaw_version": gpaw.__version__,
        "gpaw_data_path": setup_path,
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


def _restore_atoms(checkpoint: dict, expected_global_step: int | None = None) -> tuple[Atoms, dict[str, np.random.Generator]]:
    if expected_global_step is not None and checkpoint.get("global_step") != expected_global_step:
        raise GateError("restored checkpoint global step differs from the expected segment boundary")
    if int(checkpoint.get("global_step", -1)) < 1000:
        raise GateError("continuation checkpoint predates the locked step-1000 source")
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


def _load_resume_checkpoint(repo: Path, config: dict) -> dict:
    """Verify and restore the immutable step-21000 archive without integration."""
    resume = config.get("resume_checkpoint")
    expected_step = 21000
    if not isinstance(resume, dict) or int(resume.get("global_step", -1)) != expected_step:
        raise GateError("configuration lacks the exact step-21000 resume binding")
    binding = _verify_json_sha(
        PACKAGE_DIR / resume["binding_file"], resume["binding_sha256"], "step-21000 checkpoint binding"
    )
    if binding.get("task_id") != TASK_ID or binding.get("owner") != "D":
        raise GateError("step-21000 checkpoint binding task/owner mismatch")
    bound_fields = (
        "archive_commit", "archive_path", "archive_sha256sums_sha256", "manifest_sha256",
        "segment_manifest_sha256", "checkpoint_sha256", "global_step", "producer_runner_sha256",
        "source_step1000_checkpoint_sha256", "potential_sha256", "environment_sha256",
        "parameter_sha256", "scientific_lock_sha256", "strategy_lock_sha256", "qa_science_lock_sha256",
    )
    for field in bound_fields:
        if binding.get(field) != resume.get(field):
            raise GateError(f"step-21000 binding/config mismatch: {field}")
    policy = resume.get("recovery_policy")
    if binding.get("recovery_policy") != policy or not isinstance(policy, dict):
        raise GateError("step-21000 continuation policy differs from its immutable binding")
    expected_policy = {
        "qa_status": "PENDING", "qa_stage": "A_HOT_1600K", "qa_target_K": 1600.0,
        "resume_stage": "A_HOT_1600K", "resume_target_K": 1600.0,
        "checkpoint_file_must_remain_immutable": True,
        "recorded_run_status": "CHUNK_COMPLETE",
        "requires_C_authorization": False,
    }
    if policy != expected_policy:
        raise GateError("step-21000 continuation policy does not match the approved Stage A continuation")
    if resume["source_step1000_checkpoint_sha256"] != config["source_checkpoint"]["checkpoint_sha256"]:
        raise GateError("step-21000 checkpoint is not bound to the frozen step-1000 source")
    if resume["potential_sha256"] != EXPECTED_POTENTIAL_SHA256:
        raise GateError("step-21000 binding potential SHA mismatch")
    if resume["environment_sha256"] != config["runtime"]["expected_environment_sha256"]:
        raise GateError("step-21000 binding environment SHA mismatch")
    archive_rel = PurePosixPath(resume["archive_path"])
    if archive_rel.is_absolute() or ".." in archive_rel.parts:
        raise GateError("step-21000 archive path is unsafe")
    archive_dir = repo / Path(*archive_rel.parts)
    archive_commit = resume["archive_commit"]
    if subprocess.run(["git", "-C", str(repo), "merge-base", "--is-ancestor", archive_commit, "HEAD"], check=False).returncode != 0:
        raise GateError("step-21000 archive commit is not an ancestor of the D branch")
    members = archive_backend.verify_sha256sums(archive_dir, filename="SHA256SUMS.txt")
    if archive_backend.sha256_file(archive_dir / "SHA256SUMS.txt") != resume["archive_sha256sums_sha256"]:
        raise GateError("step-21000 SHA256SUMS binding mismatch")
    for filename, member_key in (("checkpoint.json", "checkpoint_sha256"),
                                 ("MANIFEST.json", "manifest_sha256"),
                                 ("SEGMENT_MANIFEST.json", "segment_manifest_sha256")):
        expected = resume[member_key]
        if members.get(filename) != expected:
            raise GateError(f"step-21000 archive member SHA mismatch: {filename}")
        blob = subprocess.run(
            ["git", "-C", str(repo), "show", f"{archive_commit}:{archive_rel.as_posix()}/{filename}"],
            check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        if blob.returncode or hashlib.sha256(blob.stdout).hexdigest() != expected:
            raise GateError(f"step-21000 Git archive blob mismatch: {filename}")
    checkpoint = _verify_json_sha(archive_dir / "checkpoint.json", resume["checkpoint_sha256"], "step-21000 checkpoint")
    sidecar = (archive_dir / "checkpoint.json.sha256").read_text(encoding="ascii").split()
    if not sidecar or sidecar[0] != resume["checkpoint_sha256"]:
        raise GateError("step-21000 checkpoint sidecar SHA mismatch")
    if checkpoint.get("global_step") != expected_step or checkpoint.get("time_fs") != expected_step * 0.5:
        raise GateError("step-21000 checkpoint boundary/time mismatch")
    exact_checkpoint_bindings = {
        "source_step1000_checkpoint_sha256": config["source_checkpoint"]["checkpoint_sha256"],
        "lock_sha256": SCIENTIFIC_LOCK_SHA256,
        "potential_sha256": EXPECTED_POTENTIAL_SHA256,
        "environment_sha256": config["runtime"]["expected_environment_sha256"],
        "parameter_sha256": config["source_checkpoint"]["parameter_sha256"],
        "runner_sha256": resume["producer_runner_sha256"],
    }
    for key, expected in exact_checkpoint_bindings.items():
        if checkpoint.get(key) != expected:
            raise GateError(f"step-21000 checkpoint binding mismatch: {key}")
    atoms, rngs = _restore_atoms(checkpoint, expected_global_step=expected_step)
    manifest = json.loads((archive_dir / "MANIFEST.json").read_text(encoding="utf-8"))
    segment = json.loads((archive_dir / "SEGMENT_MANIFEST.json").read_text(encoding="utf-8"))
    qa = json.loads((archive_dir / "QA_monitor_latest.json").read_text(encoding="utf-8"))
    run_result = json.loads((archive_dir / "RUN_RESULT.json").read_text(encoding="utf-8"))
    if manifest.get("owner") != "D" or manifest.get("task_id") != TASK_ID:
        raise GateError("step-21000 manifest owner/task mismatch")
    if manifest.get("checkpoint_sha256") != resume["checkpoint_sha256"] or manifest.get("status") != expected_policy["recorded_run_status"]:
        raise GateError("step-21000 manifest status/checkpoint mismatch")
    if manifest.get("source_step1000_checkpoint_sha256") != config["source_checkpoint"]["checkpoint_sha256"]:
        raise GateError("step-21000 manifest source checkpoint mismatch")
    if segment.get("task_id") != TASK_ID or segment.get("step_start_exclusive") != 20000 or segment.get("step_end_inclusive") != expected_step:
        raise GateError("step-21000 segment boundary/task mismatch")
    if qa.get("gate_status") != "PENDING" or qa.get("stage") != "A_HOT_1600K" or qa.get("target_temperature_K") != 1600.0 or qa.get("step") != expected_step:
        raise GateError("step-21000 QA monitor does not support Stage A continuation")
    if run_result.get("global_step") != expected_step or run_result.get("checkpoint_sha256") != resume["checkpoint_sha256"] or run_result.get("status") != expected_policy["recorded_run_status"]:
        raise GateError("step-21000 run result/status/checkpoint mismatch")
    raw_state = checkpoint.get("stage_state")
    if not isinstance(raw_state, dict) or raw_state.get("global_step") != expected_step:
        raise GateError("step-21000 checkpoint stage state is missing or at a different boundary")
    if raw_state.get("stage") != "A_HOT_1600K" or raw_state.get("failure"):
        raise GateError("step-21000 checkpoint does not contain a clean Stage A continuation state")
    if raw_state.get("last_monitor") != qa or raw_state.get("last_monitor_step") != expected_step:
        raise GateError("step-21000 checkpoint QA monitor differs from its archived report")
    schedule = checkpoint.get("temperature_schedule_position", {})
    if schedule.get("stage") != "A_HOT_1600K" or schedule.get("next_integrator_step_target_K") != 1600.0:
        raise GateError("step-21000 schedule marker differs from the archived 1600 K stage")
    thermal_rows = _read_csv(archive_dir / "thermal.csv")
    stage_rows = [row for row in thermal_rows if 1000 < int(row["step"]) <= expected_step]
    if not stage_rows or int(stage_rows[-1]["step"]) != expected_step or any(float(row["target_temperature_K"]) != 1600.0 for row in stage_rows):
        raise GateError("step-21000 thermal log does not establish the locked 1600 K Stage A history")
    state = copy.deepcopy(raw_state)
    state["stage"] = "A_HOT_1600K"
    state["user_authorized_extension"] = True
    state["recovered_from_checkpoint_sha256"] = resume["checkpoint_sha256"]
    return {
        "archive_dir": archive_dir, "checkpoint": checkpoint,
        "checkpoint_sha256": resume["checkpoint_sha256"], "atoms": atoms,
        "rngs": rngs, "state": state, "qa_monitor": qa,
    }


def _restore_verified_segment(checkpoint: dict, current_atoms: Atoms, current_rngs: dict,
                              current_state: dict, calculator):
    """Reload the exact remotely restored boundary before permitting another chunk."""
    expected_step = int(current_state["global_step"])
    if checkpoint.get("global_step") != expected_step:
        raise GateError("fresh-fetch checkpoint is at the wrong global step")
    if canonical_sha(checkpoint.get("stage_state", {})) != canonical_sha(current_state):
        raise GateError("fresh-fetch checkpoint stage/RNG control state differs")
    if canonical_sha(checkpoint.get("rng_states", {})) != canonical_sha({name: current_rngs[name].bit_generator.state for name in RNG_NAMES}):
        raise GateError("fresh-fetch checkpoint RNG streams differ from the live integrator")
    if not np.array_equal(np.asarray(checkpoint["positions_A"]), current_atoms.get_positions()):
        raise GateError("fresh-fetch checkpoint positions differ from the live boundary")
    live_momenta = current_atoms.get_momenta() * fs
    if not np.allclose(np.asarray(checkpoint["momenta_amu_A_per_fs"]), live_momenta, rtol=1e-14, atol=1e-15):
        raise GateError("fresh-fetch checkpoint velocities/momenta differ from the live boundary")
    atoms, rngs = _restore_atoms(checkpoint, expected_global_step=expected_step)
    atoms.calc = calculator
    restored_constraints = atoms.constraints
    live_constraints = current_atoms.constraints
    if len(restored_constraints) != 1 or len(live_constraints) != 1 or type(restored_constraints[0]) is not FixCom or type(live_constraints[0]) is not FixCom:
        raise GateError("fresh-fetch checkpoint constraints differ")
    if current_state["stage"] == "A_HOT_1600K":
        target = 1600.0
    elif current_state["stage"] == "B_COOL_1600_TO_1250K":
        target = cooling_target_K(expected_step - int(current_state["stage_b_start_step"]))
    else:
        target = 1250.0
    dynamics = Langevin(atoms, timestep=0.5 * fs, temperature_K=target,
                        friction=0.01 / fs, rng=rngs["thermostat"], fixcm=False)
    return atoms, rngs, dynamics


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


def run_preflight(config_path: Path, expected_config_sha: str, *, allow_active_output: bool = False) -> dict:
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
    phase_a = config.get("schedule", {}).get("phase_1_hot_liquidification", {})
    step_policy = config.get("step_policy", {})
    resume_step = int(config.get("resume_checkpoint", {}).get("global_step", -1))
    maximum_global = int(step_policy.get("maximum_global_step", -1))
    maximum_added = int(step_policy.get("maximum_added_steps", -1))
    if (
        resume_step != 21000
        or phase_a.get("maximum_global_step") != maximum_global
        or phase_a.get("maximum_additional_steps_from_current_checkpoint") != maximum_added
        or phase_a.get("maximum_total_duration_from_step0_ps") != maximum_global * 0.0005
        or step_policy.get("resume_global_step") != resume_step
        or maximum_added != maximum_global - resume_step
        or maximum_global != 60000
        or step_policy.get("maximum_chunk_steps") != 1000
        or step_policy.get("maximum_chunk_wall_seconds") != 3600
        or phase_a.get("target_K") != 1600.0
    ):
        raise GateError("extension config does not match the step21000-to-60000 condition-bounded continuation")
    strategy = _verify_json_sha(PACKAGE_DIR / "C_STRATEGY_LOCK.json", STRATEGY_LOCK_SHA256, "C strategy lock")
    scientific = _verify_json_sha(PACKAGE_DIR / "C_SCIENTIFIC_LOCK.json", SCIENTIFIC_LOCK_SHA256, "C scientific lock")
    qa_lock = _verify_json_sha(PACKAGE_DIR / "LIQUID_QA_SCIENCE_LOCK_V2_1.json", QA_SCIENCE_LOCK_SHA256, "C QA v2.1 lock")
    if criteria.SCIENCE_LOCK_SHA256 != QA_SCIENCE_LOCK_SHA256 or qa_lock != criteria.LOCK:
        raise GateError("QA v2.1 criteria are not bound to the adjudicated lock")
    if strategy.get("execution_authorized") is not False:
        raise GateError("C frozen strategy unexpectedly grants execution")
    if scientific.get("potential", {}).get("sha256") != EXPECTED_POTENTIAL_SHA256:
        raise GateError("scientific lock potential identity mismatch")
    repo = _git_root(PACKAGE_DIR)
    branch = subprocess.run(["git", "-C", str(repo), "branch", "--show-current"], check=True, capture_output=True, text=True).stdout.strip()
    if branch != config["archive"]["branch"]:
        raise GateError("current checkout is not the D-owned archive branch")
    branch_head_check = archive_backend.verify_branch_head(repo, config["archive"]["remote"], branch)
    if not branch_head_check["pass"]:
        raise GateError("D local branch head differs from the real GitHub branch head")
    source_commit = config["source_checkpoint"]["D_source_commit"]
    if subprocess.run(["git", "-C", str(repo), "merge-base", "--is-ancestor", source_commit, "HEAD"], check=False).returncode != 0:
        raise GateError("D source archive commit is not an ancestor of the current D branch")
    for commit in (config["strategy"]["C_commit"], config["qa_science_lock"]["C_commit"]):
        if subprocess.run(["git", "-C", str(repo), "cat-file", "-e", commit], check=False).returncode != 0:
            raise GateError(f"fixed C commit is unavailable in this checkout: {commit}")
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
    if checkpoint.get("environment_sha256") != SOURCE_CHECKPOINT_ENVIRONMENT_SHA256:
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
    resume_context = _load_resume_checkpoint(repo, config) if "resume_checkpoint" in config else None
    actual_code_sha = {name: sha256_file(PACKAGE_DIR / name) for name in config["code_sha256"]}
    code_pass = actual_code_sha == config["code_sha256"]
    if not code_pass:
        raise GateError("runner, criteria, archive backend, or test code SHA differs from config")
    history_path = repo / config["q6_regression"]["trajectory_path"]
    history_sha = sha256_file(history_path) if history_path.is_file() else None
    history_q6 = {"status": "BLOCKED", "trajectory_sha256": history_sha}
    if history_sha == config["q6_regression"]["trajectory_sha256"]:
        history_frames = list(Trajectory(str(history_path)))
        historical_global = [criteria.global_q6(frame) for frame in history_frames]
        historical_local = [criteria.local_q6(frame)[0] for frame in history_frames]
        history_q6 = criteria.q6_regression_report(historical_global, historical_local)
        history_q6["trajectory_sha256"] = history_sha
    q6_regression_pass = history_q6.get("regression_pass") is True
    env = _environment_snapshot()
    env_sha = canonical_sha(env)
    expected_env = config["runtime"]["expected_environment"]
    env_compare = {k: {"expected": expected_env[k], "actual": env.get(k), "pass": env.get(k) == expected_env[k]} for k in expected_env}
    expected_env_sha = config["runtime"]["expected_environment_sha256"]
    env_pass = all(v["pass"] for v in env_compare.values()) and env_sha == expected_env_sha
    potential_path = Path(config["potential"]["path"])
    potential_sha = sha256_file(potential_path) if potential_path.is_file() else None
    potential_pass = potential_sha == EXPECTED_POTENTIAL_SHA256 and potential_path.stat().st_size == config["potential"]["bytes"] if potential_sha else False
    eam_parser_pass = False
    if potential_pass:
        try:
            EAM(potential=str(potential_path), form="alloy")
            eam_parser_pass = True
        except Exception as exc:
            eam_parser_error = f"{type(exc).__name__}: {exc}"
        else:
            eam_parser_error = None
    else:
        eam_parser_error = "potential hash or file-size check failed"
    resources = _resource_snapshot(Path("/workspace"))
    resource_pass = _resource_gate(resources, config["limits"])
    output_path = Path(config["output_directory"])
    if output_path.exists() and allow_active_output:
        # This option is internal to the already-authorized loop. The external
        # preflight and run entry points always require a fresh output path.
        output_boundary_pass = output_path.is_dir() and (output_path / "checkpoint.json").is_file()
    else:
        output_boundary_pass = not output_path.exists() and os.access(output_path.parent, os.W_OK)
    resume_pass = resume_context is not None
    preflight_pass = bool(potential_pass and eam_parser_pass and env_pass and resource_pass and q6_regression_pass and output_boundary_pass and resume_pass)
    return {
        "schema_version": 1,
        "owner": "D",
        "task_id": TASK_ID,
        "config_sha256": config_sha,
        "d_branch_head": branch_head_check["local_head"],
        "remote_branch_head": branch_head_check["remote_head"],
        "code_sha256": actual_code_sha,
        "code_sha_pass": code_pass,
        "strategy_lock_sha256": STRATEGY_LOCK_SHA256,
        "scientific_lock_sha256": SCIENTIFIC_LOCK_SHA256,
        "qa_science_lock_sha256": QA_SCIENCE_LOCK_SHA256,
        "checkpoint_step": checkpoint["global_step"],
        "checkpoint_sha256": sha256_file(checkpoint_path),
        "checkpoint_position_momentum_rng_restore_pass": True,
        "checkpoint_member_sha256_pass": True,
        "resume_checkpoint_step": resume_context["checkpoint"]["global_step"] if resume_context else None,
        "resume_checkpoint_sha256": resume_context["checkpoint_sha256"] if resume_context else None,
        "resume_checkpoint_position_momentum_rng_restore_pass": resume_pass,
        "resume_checkpoint_archive_manifest_pass": resume_pass,
        "potential_path": str(potential_path),
        "potential_sha256": potential_sha,
        "potential_sha_pass": potential_pass,
        "ase_eam_alloy_parse_pass": eam_parser_pass,
        "ase_eam_alloy_parse_error": eam_parser_error,
        "q6_historical_regression": history_q6,
        "q6_historical_regression_pass": q6_regression_pass,
        "output_directory": str(output_path),
        "output_directory_absent_and_parent_writable": output_boundary_pass,
        "resources": resources,
        "resource_preflight_pass": resource_pass,
        "environment": env,
        "environment_sha256": env_sha,
        "environment_expected_sha256": expected_env_sha,
        "environment_sha_pass": env_pass,
        "environment_field_comparison": env_compare,
        "execution_authorized": False,
        "preflight_status": "PREPARED_HOLD" if preflight_pass else "BLOCKED",
        "formal_md_started": False,
        "formal_steps": 0,
        "restored_rng_streams": sorted(restored_rngs),
        "current_task_only": True,
        "notes": [
            "Read-only check only; ASE parsed the EAM file, but no energy, force, dynamics, output directory, or checkpoint write was performed.",
            "Preparation-only report; execution remains separate and no MD was started.",
            "RDF peak width and higher-shell peaks are diagnostics only under QA v2.1; all locked numeric gates remain fail-closed.",
        ],
    }


def _authorization_check(auth_path: Path, auth_sha: str, config_sha: str, config: dict, repo: Path, release_commit: str, release_repo_path: str) -> dict:
    _verify_release_blob(repo, release_commit, release_repo_path, auth_sha, auth_path)
    auth = _verify_json_sha(auth_path, auth_sha, "C execution release")
    resume = config["resume_checkpoint"]
    required = {
        "execution_authorized": True,
        "owner": "D",
        "task_id": TASK_ID,
        "strategy_lock_sha256": STRATEGY_LOCK_SHA256,
        "scientific_lock_sha256": SCIENTIFIC_LOCK_SHA256,
        "qa_science_lock_sha256": QA_SCIENCE_LOCK_SHA256,
        "config_sha256": config_sha,
        "runner_sha256": config["code_sha256"][RUNNER_NAME],
        "criteria_sha256": config["code_sha256"]["liquid_criteria_v2.py"],
        "archive_backend_sha256": config["code_sha256"]["archive_backend_v2.py"],
        "test_suite_sha256": config["code_sha256"]["test_production_runner_v2_1.py"],
        "github_archive_test_harness_sha256": config["code_sha256"]["github_archive_integration_test.py"],
        "checkpoint_sha256": resume["checkpoint_sha256"],
        "checkpoint_binding_sha256": resume["binding_sha256"],
        "resume_global_step": resume["global_step"],
        "potential_sha256": EXPECTED_POTENTIAL_SHA256,
        "environment_sha256": config["runtime"]["expected_environment_sha256"],
        "maximum_global_step": 60000,
        "maximum_added_steps": config["step_policy"]["maximum_added_steps"],
        "maximum_chunk_steps": 1000,
        "maximum_chunk_wall_seconds": 3600,
        "d_branch_head": subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip(),
    }
    for key, value in required.items():
        if auth.get(key) != value:
            raise GateError(f"C execution release binding mismatch: {key}")
    if auth.get("execution_mode") not in (None, "bounded_automatic_continuation"):
        raise GateError("C execution release has an unrecognized execution mode")
    return auth


def _json_safe(value):
    """Convert NumPy diagnostics to JSON-native containers without hiding NaNs."""
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _atomic_json(path: Path, value: dict) -> None:
    temp = path.with_name(path.name + f".tmp-{os.getpid()}")
    with temp.open("x", encoding="utf-8") as stream:
        # RDF diagnostics include NumPy arrays.  Preserve their values as JSON
        # lists while keeping allow_nan=False so invalid statistics still fail closed.
        json.dump(_json_safe(value), stream, indent=2, sort_keys=True, allow_nan=False)
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


STRESS_CSV_COLUMNS = (
    "step", "time_fs", "xx_eV_A3", "yy_eV_A3", "zz_eV_A3",
    "yz_eV_A3", "xz_eV_A3", "xy_eV_A3",
)
THERMAL_CSV_COLUMNS = (
    "step", "time_fs", "target_temperature_K", "temperature_K",
    "potential_energy_eV", "kinetic_energy_eV", "total_energy_eV",
)


def _eligible_pressure_stress_samples(out: Path, target_K: float, discard_end_step: int,
                                      through_step: int) -> dict:
    """Join fixed-target stress.csv and thermal.csv rows after the phase discard."""
    stress_path = out / "stress.csv"
    thermal_path = out / "thermal.csv"
    if not stress_path.is_file() or not thermal_path.is_file():
        return {"input_error": "stress.csv or thermal.csv is missing"}
    try:
        with stress_path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != STRESS_CSV_COLUMNS:
                return {"input_error": "stress.csv does not contain the required six named Voigt components"}
            stress_rows = list(reader)
        with thermal_path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != THERMAL_CSV_COLUMNS:
                return {"input_error": "thermal.csv header is incomplete or inconsistent"}
            thermal_rows = list(reader)
        stress_by_step: dict[int, dict] = {}
        thermal_by_step: dict[int, dict] = {}
        for row in stress_rows:
            step = int(row["step"])
            if step in stress_by_step:
                return {"input_error": f"duplicate stress row at step {step}"}
            stress_by_step[step] = row
        for row in thermal_rows:
            step = int(row["step"])
            if step in thermal_by_step:
                return {"input_error": f"duplicate thermal row at step {step}"}
            thermal_by_step[step] = row

        eligible_thermal: dict[int, dict] = {}
        for step, row in thermal_by_step.items():
            if not (discard_end_step < step <= through_step):
                continue
            target = float(row["target_temperature_K"])
            if not math.isfinite(target):
                return {"input_error": f"non-finite target temperature at step {step}"}
            if math.isclose(target, target_K, rel_tol=0.0, abs_tol=1e-9):
                eligible_thermal[step] = row
        eligible_stress = {
            step: row for step, row in stress_by_step.items()
            if discard_end_step < step <= through_step
        }
        if set(eligible_thermal) != set(eligible_stress):
            missing_stress = sorted(set(eligible_thermal) - set(eligible_stress))[:8]
            missing_thermal = sorted(set(eligible_stress) - set(eligible_thermal))[:8]
            return {"input_error": "stress/thermal sample mismatch; "
                    f"missing_stress_steps={missing_stress}, missing_thermal_steps={missing_thermal}"}
        ordered_steps = sorted(eligible_thermal)
        times_ps, components = [], []
        for step in ordered_steps:
            thermal_row = eligible_thermal[step]
            stress_row = eligible_stress[step]
            time_fs_thermal = float(thermal_row["time_fs"])
            time_fs_stress = float(stress_row["time_fs"])
            expected_time_fs = step * 0.5
            if not all(math.isfinite(value) for value in (time_fs_thermal, time_fs_stress)):
                return {"input_error": f"non-finite sample time at step {step}"}
            if not (math.isclose(time_fs_thermal, expected_time_fs, rel_tol=0.0, abs_tol=1e-9)
                    and math.isclose(time_fs_stress, expected_time_fs, rel_tol=0.0, abs_tol=1e-9)):
                return {"input_error": f"stress/thermal time does not match step {step}"}
            try:
                row_values = [float(stress_row[name]) for name in STRESS_CSV_COLUMNS[2:]]
            except (TypeError, ValueError) as exc:
                return {"input_error": f"missing or invalid stress component at step {step}: {exc}"}
            if not np.isfinite(row_values).all():
                return {"input_error": f"non-finite stress component at step {step}"}
            times_ps.append(expected_time_fs * 0.001)
            components.append(row_values)
        return {"times_ps": times_ps, "stress_eV_A3": components,
                "sample_count": len(ordered_steps), "target_temperature_K": target_K,
                "discard_end_step_exclusive": discard_end_step, "through_step_inclusive": through_step}
    except (OSError, csv.Error, KeyError, TypeError, ValueError) as exc:
        return {"input_error": f"unable to parse stress/thermal samples: {type(exc).__name__}: {exc}"}


def _load_frame_metrics(directory: Path) -> list[criteria.FrameMetrics]:
    result = []
    path = directory / "frame_metrics.jsonl"
    if not path.is_file():
        return result
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        item = json.loads(line)
        result.append(criteria.frame_metrics_from_json(item))
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
        "runner_sha256": config["code_sha256"][RUNNER_NAME],
        "criteria_sha256": config["code_sha256"]["liquid_criteria_v2.py"],
        "environment_sha256": config["runtime"]["expected_environment_sha256"],
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
            "runner_sha256": config["code_sha256"][RUNNER_NAME],
            "criteria_sha256": config["code_sha256"]["liquid_criteria_v2.py"],
            "qa_science_lock_sha256": QA_SCIENCE_LOCK_SHA256,
        }, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        archive_relpath = f"{config['archive']['root']}/chunks/step_{global_step:05d}"
        restore_dir = Path(tempfile.gettempdir()) / f"ag64-rho100-restored-{global_step:05d}-{os.getpid()}"
        head = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip()
        try:
            result = archive_backend.archive_push_fresh_restore(
                repo_root=repo,
                source_dir=bundle,
                archive_relpath=archive_relpath,
                branch=config["archive"]["branch"],
                remote=config["archive"]["remote"],
                expected_remote_head=head,
                restore_dir=restore_dir,
                commit_message=f"Archive Ag64 rho100 verified segment through step {global_step}",
            )
            result["_verified_checkpoint_payload"] = json.loads((restore_dir / "checkpoint.json").read_text(encoding="utf-8"))
            return result
        finally:
            if restore_dir.exists():
                shutil.rmtree(restore_dir)


def _run_loop(config_path: Path, config_sha: str, config: dict, auth: dict, atoms: Atoms, rngs: dict, state: dict, output_dir: Path, repo: Path, initial_archive_dir: Path | None = None) -> dict:
    global _STOP_REQUESTED
    _STOP_REQUESTED = False
    potential = Path(config["potential"]["path"])
    if sha256_file(potential) != EXPECTED_POTENTIAL_SHA256:
        raise GateError("potential changed after preflight")
    if output_dir.exists():
        raise GateError("run output directory already exists; use only a separately verified resume archive")
    if initial_archive_dir is not None:
        # Resume from the exact archived boundary, preserving its trajectory,
        # thermodynamic/stress logs, unwrapped coordinates, metrics, and RNG checkpoint.
        archive_backend.verify_sha256sums(initial_archive_dir)
        shutil.copytree(initial_archive_dir, output_dir)
    else:
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
    execution_start_step = current_step
    maximum_global = min(int(config["step_policy"]["maximum_global_step"]), current_step + int(config["step_policy"]["maximum_added_steps"]))
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
                live_preflight = run_preflight(config_path, config_sha, allow_active_output=True)
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
                phase_boundary = min(int(config["schedule"]["phase_1_hot_liquidification"]["maximum_global_step"]), current_step + chunk_steps)
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
                if not can_start_step(current_step, segment_end, elapsed, maximum_global, 3600.0):
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
            status = apply_completed_chunk_limits(state, current_step, status, maximum_global, int(config["schedule"]["phase_1_hot_liquidification"]["maximum_global_step"]))

            # Persist the post-QA phase decision and all monitor state in the
            # checkpoint that is about to be archived.
            state["global_step"] = current_step
            last_checkpoint_sha = _save_checkpoint(output_dir / "checkpoint.json", atoms, rngs, state, config, start_wall)
            _atomic_json(output_dir / "RUN_STATE.json", {**state, "status": status})
            _write_result_and_sums(output_dir, state, status, time.monotonic() - start_wall, last_checkpoint_sha, config)
            try:
                archive_result = _bundle_and_archive(repo, output_dir, current_step, config, auth)
            except BaseException as exc:
                status = "ARCHIVE_OR_RESTORE_FAILURE"
                state["stage"] = status
                state["archive_failure"] = f"{type(exc).__name__}: {exc}"
                _atomic_json(output_dir / "RUN_STATE.json", {**state, "status": status})
                _write_result_and_sums(output_dir, state, status, time.monotonic() - start_wall, last_checkpoint_sha, config)
                break
            if not archive_backend.segment_archive_passed(archive_result):
                status = "ARCHIVE_OR_RESTORE_FAILURE"
                state["stage"] = status
                state["archive_failure"] = "segment archive result did not pass every required gate"
                _atomic_json(output_dir / "RUN_STATE.json", {**state, "status": status})
                _write_result_and_sums(output_dir, state, status, time.monotonic() - start_wall, last_checkpoint_sha, config)
                break
            restored_payload = archive_result.pop("_verified_checkpoint_payload", None)
            try:
                restored_atoms, restored_rngs, restored_dynamics = _restore_verified_segment(
                    restored_payload, atoms, rngs, state, calculator
                )
            except BaseException as exc:
                status = "ARCHIVE_OR_RESTORE_FAILURE"
                state["stage"] = status
                state["archive_failure"] = f"fresh-fetch restore mismatch: {type(exc).__name__}: {exc}"
                _atomic_json(output_dir / "RUN_STATE.json", {**state, "status": status})
                _write_result_and_sums(output_dir, state, status, time.monotonic() - start_wall, last_checkpoint_sha, config)
                break
            if not _should_continue_after_archive(status, archive_result, _STOP_REQUESTED):
                # A signal can arrive exactly after a chunk's last step. That
                # chunk is already archived and restored; do not write a second
                # immutable archive for the same step boundary.
                if _STOP_REQUESTED:
                    status = "GRACEFUL_SIGNAL_STOP"
                break
            # Continue automatically only after upload, fresh fetch, restore,
            # and per-member SHA checks have all passed.
            atoms, rngs, dynamics = restored_atoms, restored_rngs, restored_dynamics
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
            _write_result_and_sums(output_dir, state, status, time.monotonic() - start_wall, last_checkpoint_sha, config)
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
        "actual_start_step": execution_start_step,
        "actual_end_step": current_step,
        "time_ps": current_step * 0.0005,
        "checkpoint_sha256": last_checkpoint_sha,
        "archive": archive_result,
    }


def _monitor_at_boundary(out: Path, step: int, state: dict, config: dict, auth: dict) -> dict:
    metrics = _load_frame_metrics(out)
    stage = state["stage"]
    if stage == "A_HOT_1600K":
        discard_end, target = int(state["stage_a_initial_transient_end_step"]), 1600.0
    elif stage == "C_HOLD_1250K":
        discard_end, target = int(state["stage_c_start_step"]) + 4000, 1250.0
    elif stage == "D_CANDIDATE_1250K":
        discard_end, target = int(state["stage_d_start_step"]), 1250.0
    else:
        return {"step": step, "gate_status": "MONITOR_ONLY_STAGE_B_RAMP"}

    end_ps = step * 0.0005
    windows = [[end_ps - 4.0, end_ps - 2.0], [end_ps - 2.0, end_ps]]
    stationary_frames = [f for f in metrics if f.step > discard_end and f.step <= step]
    unwrapped_path = out / "unwrapped_positions.jsonl"
    lower_step = int(round(windows[0][0] / 0.0005))
    unwrapped = []
    if unwrapped_path.is_file():
        for line in unwrapped_path.read_text(encoding="utf-8").splitlines():
            if line.strip():
                item = json.loads(line)
                frame_step = int(item["step"])
                if lower_step <= frame_step <= step:
                    unwrapped.append(item["positions_A"])
    historical = _historical_q6_regression(_git_root(PACKAGE_DIR), config)
    assessment = criteria.assess_phase(
        stationary_frames,
        windows,
        unwrapped,
        target,
        historical,
        config["thermal_thresholds"],
        pressure_stress=_eligible_pressure_stress_samples(out, target, discard_end, step),
    )
    last_window = [f for f in stationary_frames if windows[1][0] <= f.time_ps < windows[1][1]]
    monitor = {
        "step": step,
        "stage": stage,
        "target_temperature_K": target,
        "stationary_discard_end_step": discard_end,
        "window_step_bounds": [[int(round(a / 0.0005)), int(round(b / 0.0005))] for a, b in windows],
        "stationary_frame_count": len(stationary_frames),
        "final_window_global_q6_mean": float(np.mean([f.global_q6 for f in last_window])) if last_window else None,
        "final_window_mean_local_q6_report_only": float(np.mean([f.mean_local_q6 for f in last_window])) if last_window else None,
        "final_window_ptm_crystalline_fraction_mean": float(np.mean([f.ptm_crystalline_fraction for f in last_window])) if last_window and all(f.ptm_crystalline_fraction is not None for f in last_window) else None,
        "final_window_coordination_number_mean": float(np.mean([f.coordination_number for f in last_window])) if last_window and all(f.coordination_number is not None for f in last_window) else None,
        "local_q6_role": "REPORT_ONLY",
        "q6_historical_regression": historical,
        "assessment": assessment,
        "gate_status": assessment["gate_status"],
        "reason": assessment.get("reason"),
    }
    if unwrapped:
        monitor["MSD_A2"] = assessment.get("msd", {}).get("delta_msd_A2")
    return monitor


def _historical_q6_regression(repo: Path, config: dict) -> dict:
    source = repo / config["q6_regression"]["trajectory_path"]
    actual_sha = sha256_file(source) if source.is_file() else None
    if actual_sha != config["q6_regression"]["trajectory_sha256"]:
        return {"status": "QA_HOLD", "regression_pass": False, "reason": "historical trajectory SHA mismatch"}
    frames = list(Trajectory(str(source)))
    global_values = [criteria.global_q6(frame) for frame in frames]
    local_values = [criteria.local_q6(frame)[0] for frame in frames]
    report = criteria.q6_regression_report(global_values, local_values)
    report["trajectory_sha256"] = actual_sha
    return report


def _write_result_and_sums(out: Path, state: dict, status: str, wall_s: float, checkpoint_sha: str, config: dict) -> None:
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
        "checkpoint_sha256": checkpoint_sha,
        "checkpoint_binding_sha256": config.get("resume_checkpoint", {}).get("binding_sha256"),
        "potential_sha256": EXPECTED_POTENTIAL_SHA256,
        "runner_sha256": state.get("runner_sha256"),
        "criteria_sha256": state.get("criteria_sha256"),
        "environment_sha256": config["runtime"]["expected_environment_sha256"],
        "data_role": "SINGLE_SOURCE_FAMILY_SAMPLING_CANDIDATE_ONLY_NO_TRAIN_DEV_TEST",
        "member_sha256": members,
    })
    archive_backend.write_sha256sums(out)
    archive_backend.verify_sha256sums(out)


def _user_extension_authorization_check(auth_path: Path, auth_sha: str, config_sha: str,
                                        config: dict, repo: Path) -> dict:
    auth = _verify_json_sha(auth_path, auth_sha, "direct user Stage A extension authorization")
    resume = config["resume_checkpoint"]
    expected = {
        "execution_authorized": True,
        "authorization_source": "direct user instruction in current session",
        "owner": "D", "task_id": TASK_ID, "config_sha256": config_sha,
        "runner_sha256": config["code_sha256"][RUNNER_NAME],
        "criteria_sha256": config["code_sha256"]["liquid_criteria_v2.py"],
        "archive_backend_sha256": config["code_sha256"]["archive_backend_v2.py"],
        "test_suite_sha256": config["code_sha256"]["test_step21000_extension_runner.py"],
        "checkpoint_sha256": resume["checkpoint_sha256"],
        "checkpoint_binding_sha256": resume["binding_sha256"],
        "resume_global_step": int(resume["global_step"]),
        "stage_a_maximum_global_step": int(config["schedule"]["phase_1_hot_liquidification"]["maximum_global_step"]),
        "maximum_global_step": int(config["step_policy"]["maximum_global_step"]),
        "maximum_added_steps": int(config["step_policy"]["maximum_added_steps"]),
        "maximum_chunk_steps": 1000, "maximum_chunk_wall_seconds": 3600,
        "potential_sha256": EXPECTED_POTENTIAL_SHA256,
        "environment_sha256": config["runtime"]["expected_environment_sha256"],
        "checkpoint_archive_commit": resume["archive_commit"],
    }
    for key, value in expected.items():
        if auth.get(key) != value:
            raise GateError(f"direct user extension authorization mismatch: {key}")
    if subprocess.run(["git", "-C", str(repo), "merge-base", "--is-ancestor", resume["archive_commit"], "HEAD"], check=False).returncode != 0:
        raise GateError("resume checkpoint archive commit is not an ancestor of the D branch")
    return auth


def run_authorized(config_path: Path, config_sha: str, authorization_path: Path,
                   authorization_sha: str) -> dict:
    report = run_preflight(config_path, config_sha)
    if report["preflight_status"] != "PREPARED_HOLD":
        raise GateError("live preflight did not pass")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    repo = _git_root(PACKAGE_DIR)
    auth = _user_extension_authorization_check(authorization_path, authorization_sha,
                                               config_sha, config, repo)
    auth["_release_sha256"] = authorization_sha
    resume_context = _load_resume_checkpoint(repo, config)
    atoms, rngs, state = resume_context["atoms"], resume_context["rngs"], resume_context["state"]
    out = Path(config["output_directory"])
    if out.exists():
        raise GateError("output directory exists; no overwrite/reinitialization allowed")
    return _run_loop(config_path, config_sha, config, auth, atoms, rngs, state, out, repo,
                     initial_archive_dir=resume_context["archive_dir"])


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
            raise GateError("--run requires the separately hash-bound direct user authorization file")
        result = run_authorized(args.config, args.config_sha256, args.authorization, args.authorization_sha256)
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
        return 0
    except Exception as exc:
        print(json.dumps({"status": "BLOCKED", "error": f"{type(exc).__name__}: {exc}"}, sort_keys=True), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
