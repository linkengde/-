#!/usr/bin/env python3
"""Owner-only, separately gated liquid-Ag sampling with reproducible restart.

C may use --help and static checks. Scientific execution belongs to D/E/F.
No default potential, invented hash, coordinate, or release is supplied.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import shutil
import signal
import socket
import sys
import time
import uuid

from inspect_eam import inspect_potential, sha256_file


THERMAL_FIELDS = ["step", "time_fs", "phase", "target_temperature_K", "T_K",
                  "sampler_energy_eV", "pressure_eV_A3", "COMspeed_A_fs"]
REQUIRED_COUNTS = {"warmup_steps": 1000, "hot_steps": 20000, "cool_steps": 20000,
                   "equil_steps": 40000, "production_steps": 40000}
PLACEHOLDERS = ("UNASSIGNED", "UNVERIFIED", "TO_BE_FROZEN", "NOT_COMPUTED", "PENDING")
RESOURCE_PARAMETER_KEYS = {"max_wall_s", "warmup_max_wall_s", "max_rss_bytes",
                           "min_free_bytes", "trajectory_archive_max_bytes",
                           "checkpoint_reserve_bytes"}


class Hold(RuntimeError):
    """A release, identity, resource or numerical stop; never retry implicitly."""


def canonical_sha(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                      allow_nan=False).encode()).hexdigest()


def reject_symlinks(path: Path) -> None:
    for component in (path, *path.parents):
        if component.is_symlink():
            raise Hold("SYMLINK_PATH_FORBIDDEN:" + str(component))


def owned_path(path: Path, owner_root: Path) -> Path:
    """Check lexical and resolved containment before every owner-side write."""
    owner_root = Path(os.path.abspath(owner_root))
    path = Path(os.path.abspath(path))
    reject_symlinks(owner_root)
    if path != owner_root and owner_root not in path.parents:
        raise Hold("PATH_ESCAPES_APPROVED_OWNER_ROOT")
    reject_symlinks(path)
    resolved = path.resolve(strict=False)
    if resolved != owner_root and owner_root not in resolved.parents:
        raise Hold("RESOLVED_PATH_ESCAPES_APPROVED_OWNER_ROOT")
    return resolved


def fsync_directory(path: Path) -> None:
    reject_symlinks(path)
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fsync_owned_file(path: Path, owner_root: Path) -> None:
    path = owned_path(path, owner_root)
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def copy_owned(source: Path, destination: Path, owner_root: Path) -> None:
    source = owned_path(source, owner_root)
    destination = owned_path(destination, owner_root)
    shutil.copyfile(source, destination)
    fsync_owned_file(destination, owner_root)
    fsync_directory(destination.parent)


def atomic_json(path: Path, value, owner_root: Path) -> None:
    path = owned_path(path, owner_root)
    temporary = path.with_name(path.name + ".tmp-" + uuid.uuid4().hex)
    owned_path(temporary, owner_root)
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    owned_path(path, owner_root)
    os.replace(temporary, path)
    fsync_directory(path.parent)


def require_true(value, message: str) -> None:
    if value is not True:
        raise Hold(message)


def fixed_source(value) -> bool:
    return isinstance(value, str) and bool(value.strip()) and not any(
        token in value.upper() for token in PLACEHOLDERS)


def overlaps(first: Path, second: Path) -> bool:
    return first == second or first in second.parents or second in first.parents


def validate_release(args):
    """All authorization checks occur before ASE import, coordinates or mkdir."""
    if not args.config.is_file():
        raise Hold("C_RELEASE_CONFIG_MISSING")
    expected_config_sha = args.config_sha256.lower()
    if len(expected_config_sha) != 64 or any(c not in "0123456789abcdef" for c in expected_config_sha):
        raise Hold("C_FROZEN_CONFIG_SHA_NOT_VALID")
    raw_config = args.config.read_bytes()
    actual_config_sha = hashlib.sha256(raw_config).hexdigest()
    if actual_config_sha != expected_config_sha:
        raise Hold("C_FROZEN_CONFIG_SHA_MISMATCH")
    args.verified_config_sha256 = actual_config_sha
    config = json.loads(raw_config.decode("utf-8"))
    for name in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_SIZE"):
        if int(os.environ.get(name, "1")) > 1:
            raise Hold("SAMPLING_IS_SINGLE_PROCESS; do not launch this script with mpiexec")
    if config.get("schema_version") != 1:
        raise Hold("UNSUPPORTED_RELEASE_SCHEMA")
    for gate in ("source_release", "storage_release", "warmup_release"):
        require_true(config.get(gate), gate.upper() + "_HOLD")
    require_true(config.get("density_reference_pass"), "DENSITY_REFERENCE_HOLD")
    if args.stage == "continue":
        require_true(config.get("continuation_release"), "CONTINUATION_RELEASE_HOLD")
        if args.task_id not in config.get("continuation_task_ids", []):
            raise Hold("TASK_WARMUP_NOT_SIGNED_OFF_BY_C")
    hosts = config.get("hosts", {})
    if args.owner not in hosts:
        raise Hold("OWNER_HOST_NOT_REGISTERED")
    host = hosts[args.owner]
    require_true(host.get("persistence_verified"), "PERSISTENT_WORK_ROOT_UNVERIFIED")
    require_true(host.get("isolation_verified"), "WORK_ROOT_ISOLATION_UNVERIFIED")
    if args.task_id not in host.get("allowed_task_ids", []):
        raise Hold("TASK_NOT_OWNED_BY_THIS_HOST")
    root_value = host.get("sampler_work_root")
    if not fixed_source(root_value) or not Path(root_value).is_absolute():
        raise Hold("SAMPLER_WORK_ROOT_NOT_FROZEN")
    reject_symlinks(Path(root_value))
    root = Path(root_value).resolve()
    repository = Path(__file__).resolve()
    for parent in repository.parents:
        if (parent / ".git").exists():
            if overlaps(root, parent):
                raise Hold("SAMPLER_ROOT_MUST_BE_OUTSIDE_GIT_CHECKOUT")
            break
    for owner, other in hosts.items():
        other_root = other.get("sampler_work_root")
        if owner != args.owner and fixed_source(other_root):
            if overlaps(root, Path(other_root).resolve()):
                raise Hold("SAMPLER_ROOT_OVERLAPS_ANOTHER_OWNER")
        dft_root = other.get("dft_work_root")
        if fixed_source(dft_root) and overlaps(root, Path(dft_root).resolve()):
            raise Hold("SAMPLER_ROOT_OVERLAPS_DFT_WORKSPACE")
    for env_name in ("DFT_WORK_ROOT", "SAMPLER_WORK_ROOT"):
        value = os.environ.get(env_name)
        if value and env_name == "DFT_WORK_ROOT" and overlaps(root, Path(value).resolve()):
            raise Hold("SAMPLER_ROOT_OVERLAPS_DFT_WORK_ROOT_ENV")
        if value and env_name == "SAMPLER_WORK_ROOT" and root != Path(value).resolve():
            raise Hold("SAMPLER_WORK_ROOT_ENV_DIFFERS_FROM_C_CONFIG")
    tasks_path = args.config.parent / config.get("tasks_file", "AUTHORIZED_SIX_TASKS.csv")
    if not tasks_path.is_file():
        raise Hold("AUTHORIZED_TASK_MATRIX_MISSING")
    task_sha = sha256_file(tasks_path)
    if task_sha != config.get("tasks_sha256"):
        raise Hold("AUTHORIZED_TASK_MATRIX_SHA_MISMATCH")
    with tasks_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 6 or len({row["task_id"] for row in rows}) != 6:
        raise Hold("TASK_MATRIX_MUST_HAVE_EXACTLY_SIX_UNIQUE_TASKS")
    own_rows = [row for row in rows if row["owner"] == args.owner]
    if len(own_rows) != 2:
        raise Hold("OWNER_MUST_HAVE_EXACTLY_TWO_TASKS")
    selected = [row for row in own_rows if row["task_id"] == args.task_id]
    if len(selected) != 1:
        raise Hold("TASK_IS_NOT_IN_OWNED_MATRIX_ROWS")
    row = selected[0]
    if int(row["atom_count"]) != 64 or row["PBC"].lower() != "true,true,true":
        raise Hold("ATOM_COUNT_OR_PBC_CHANGED")
    if int(row["RNG_seed"]) <= 0 or not fixed_source(row["trajectory_parent"]):
        raise Hold("SEED_OR_PARENT_MISSING")
    source_root = config.get("lineage", {}).get("source_root_by_task", {}).get(
        args.task_id, row.get("source_root"))
    if not fixed_source(source_root):
        raise Hold("REAL_SOURCE_ROOT_NOT_FROZEN")
    density, edge = float(row["density_g_cm3"]), float(row["box_edge_A"])
    if not 8.0 < density < 11.0 or not 10.0 < edge < 12.0:
        raise Hold("DENSITY_OR_BOX_OUT_OF_FROZEN_PILOT_RANGE")
    # Do not change either density or edge to resolve an inconsistency.
    expected_density = 64 * 107.8682 / 6.02214076e23 / (edge**3 * 1e-24)
    if abs(expected_density - density) > 1e-9:
        raise Hold("DENSITY_BOX_MASS_IDENTITY_MISMATCH")
    sampling = config.get("sampling", {})
    for key, value in REQUIRED_COUNTS.items():
        if sampling.get(key) != value:
            raise Hold("FROZEN_STEP_SCHEDULE_MISMATCH:" + key)
    required_scalars = {"dt_fs": 0.5, "T_hot_K": 1600, "T_target_K": 1250,
                        "min_initial_distance_A": 2.2, "write_every": 100,
                        "thermal_every": 20}
    for key, value in required_scalars.items():
        if sampling.get(key) != value:
            raise Hold("FROZEN_SAMPLING_PARAMETER_MISMATCH:" + key)
    for key in ("friction_per_fs", "checkpoint_every", "packing_max_attempts", "min_runtime_distance_A",
                "max_rss_bytes", "min_free_bytes", "trajectory_archive_max_bytes"):
        if not math.isfinite(float(sampling.get(key, 0))) or float(sampling.get(key, 0)) <= 0:
            raise Hold("RESOURCE_OR_SAMPLING_PARAMETER_NOT_FROZEN:" + key)
    potential = config.get("potential", {})
    require_true(potential.get("source_verified"), "POTENTIAL_SOURCE_UNVERIFIED")
    require_true(potential.get("license_verified"), "POTENTIAL_LICENSE_UNVERIFIED")
    potential_sha = potential.get("sha256", "")
    if len(potential_sha) != 64 or any(c not in "0123456789abcdef" for c in potential_sha):
        raise Hold("POTENTIAL_SHA_NOT_FROZEN")
    path_value = potential.get("absolute_path") or potential.get("relative_path")
    if not fixed_source(path_value):
        raise Hold("POTENTIAL_PATH_NOT_FROZEN")
    potential_path = Path(path_value)
    if not potential_path.is_absolute():
        potential_path = args.config.parent / potential_path
    if not potential_path.is_file():
        raise Hold("POTENTIAL_NOT_AVAILABLE; no alternative is loaded")
    if sha256_file(potential_path) != potential_sha:
        raise Hold("POTENTIAL_SHA_MISMATCH")
    return config, host, root, row, source_root, task_sha, potential_path


def tree_bytes(root: Path) -> int:
    total = 0
    owned_path(root, root)
    for path in root.rglob("*"):
        owned_path(path, root)
        if path.is_file():
            total += path.stat().st_size
    return total


def require_checkpoint_resume(config: dict, task_id: str, state_path: Path, owner_root: Path) -> None:
    state_path = owned_path(state_path, owner_root)
    if (task_id not in config.get("resume_task_ids", []) or not state_path.is_file()
            or config.get("resume_checkpoint_sha256_by_task", {}).get(task_id) != sha256_file(state_path)):
        raise Hold("PAUSED_CHECKPOINT_REQUIRES_C_RESUME_RELEASE_BOUND_TO_CURRENT_STATE_SHA")


def current_rss_bytes() -> int:
    with Path("/proc/self/status").open() as handle:
        for line in handle:
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) * 1024
    raise Hold("RSS_ACCOUNTING_UNAVAILABLE")


def nearest_distance(positions, edge, np) -> float:
    difference = positions[:, None, :] - positions[None, :, :]
    difference -= np.rint(difference / edge) * edge
    squared = (difference**2).sum(axis=2)
    np.fill_diagonal(squared, np.inf)
    return float(np.sqrt(squared.min()))


def pack_random(seed: int, edge: float, limit: int, np, budget_check=None):
    """Independent 64-particle periodic rejection packing, not a sliced solid."""
    packing_sequence, velocity_sequence, thermostat_sequence = np.random.SeedSequence(seed).spawn(3)
    packing_rng = np.random.Generator(np.random.PCG64(packing_sequence))
    velocity_rng = np.random.Generator(np.random.PCG64(velocity_sequence))
    thermostat_rng = np.random.Generator(np.random.PCG64(thermostat_sequence))
    placed = []
    attempts = 0
    while len(placed) < 64 and attempts < limit:
        if budget_check is not None and attempts % 1000 == 0:
            budget_check(check_archive=False)
        trial = packing_rng.uniform(0.0, edge, size=3)
        attempts += 1
        if placed:
            delta = np.asarray(placed) - trial
            delta -= np.rint(delta / edge) * edge
            if np.any((delta**2).sum(axis=1) < 2.2**2):
                continue
        placed.append(trial)
    if len(placed) != 64:
        raise Hold("RANDOM_PACKING_ATTEMPTS_EXHAUSTED; no distance or cell change allowed")
    return np.asarray(placed), packing_rng, velocity_rng, thermostat_rng, attempts


def schedule(step: int, sampling: dict):
    hot, cooling = sampling["hot_steps"], sampling["cool_steps"]
    equilibrating = sampling["equil_steps"]
    if step < hot:
        return "hot_melt_candidate", float(sampling["T_hot_K"])
    if step < hot + cooling:
        fraction = (step - hot) / cooling
        return "cooling", sampling["T_hot_K"] + fraction * (
            sampling["T_target_K"] - sampling["T_hot_K"])
    if step < hot + cooling + equilibrating:
        return "equilibration_candidate", float(sampling["T_target_K"])
    return "production_candidate", float(sampling["T_target_K"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--config-sha256", required=True,
                        help="C-published SHA-256 of the exact release-config bytes")
    parser.add_argument("--owner", choices=["D", "E", "F"], required=True)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--stage", choices=["warmup", "continue"], required=True)
    args = parser.parse_args()
    task_dir = None
    lock_path = None
    host_lock_path = None
    completed_step = 0
    checkpoint = None
    trajectory = None
    thermal_handle = None
    lock_acquired = False
    host_lock_acquired = False
    stop_requested = [None]
    previous_handlers = {}
    budget_started = False
    persist_budget = None
    try:
        config, host, root, row, source_root, task_sha, potential_path = validate_release(args)
        sampling = config["sampling"]
        # Only after all source, ownership and storage gates pass do imports occur.
        import numpy as np
        from ase import Atoms, units
        from ase.constraints import FixCom
        from ase.io.trajectory import Trajectory
        from ase.md.langevin import Langevin
        from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
        packages = {name: version(name) for name in ("ase", "numpy", "scipy")}
        if packages["ase"] != "3.29.0":
            raise Hold("ASE_VERSION_MISMATCH")
        for name, expected in config.get("expected_versions", {}).items():
            if version(name) != expected:
                raise Hold("PACKAGE_VERSION_MISMATCH:" + name)
        edge = float(row["box_edge_A"])
        inspection, calc = inspect_potential(potential_path, config["potential"]["format"], [edge])
        if inspection["format"] != config["potential"]["format"]:
            raise Hold("POTENTIAL_FORMAT_DIFFERS_FROM_FROZEN_CONFIG")
        if abs(inspection["cutoff_A"] - float(config["potential"].get("cutoff_A", -1))) > 1e-9:
            raise Hold("POTENTIAL_CUTOFF_DIFFERS_FROM_FROZEN_CONFIG")
        tail = inspection["radial_tail"]
        if tail["requires_explicit_C_policy"]:
            if config["potential"].get("radial_tail_policy") != "C_APPROVED_PUBLISHED_SETFL_ONE_BIN_ZERO_TAIL":
                raise Hold("PUBLISHED_SETFL_RADIAL_TAIL_C_REVIEW_REQUIRED")
            if (inspection["format"] != "alloy" or not tail["within_one_dr"]
                    or tail["max_abs_density"] > 1e-10 or tail["max_abs_phi_eV"] > 1e-10):
                raise Hold("PUBLISHED_SETFL_RADIAL_TAIL_EXCEEDS_APPROVED_ONE_BIN_ZERO_TAIL")
        code_sha = {path.name: sha256_file(path) for path in (
            Path(__file__), Path(__file__).with_name("inspect_eam.py"))}
        fixed_identity = {"potential": config["potential"],
                          "sampling": {key: value for key, value in sampling.items()
                                       if key not in RESOURCE_PARAMETER_KEYS},
                          "row": row, "source_root": source_root, "task_matrix_sha256": task_sha,
                          "packages": packages, "code_sha256": code_sha,
                          "sampler_work_root": str(root), "langevin_version": Langevin._lgv_version,
                          "hostname": socket.gethostname(), "python_version": sys.version,
                          "thread_environment": {key: os.environ.get(key) for key in (
                              "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
                          "thermostat": "ASE Langevin fixcm=False with FixCom; PCG64 Generator"}
        scientific_sha = canonical_sha(fixed_identity)
        raw_config_sha = args.verified_config_sha256
        owned_path(root, root)
        if not root.is_dir():
            raise Hold("APPROVED_PERSISTENT_OWNER_ROOT_MUST_BE_PRECREATED_AND_CHECKED")
        host_lock_path = owned_path(root / "HOST_SAMPLER_RUNNING.lock", root)
        try:
            descriptor = os.open(host_lock_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
        except FileExistsError as exc:
            raise Hold("ANOTHER_SAMPLER_TASK_OR_STALE_HOST_LOCK; no concurrent same-host trajectories") from exc
        with os.fdopen(descriptor, "w") as handle:
            json.dump({"owner": args.owner, "pid": os.getpid(), "hostname": socket.gethostname(),
                       "task_id": args.task_id}, handle)
            handle.flush()
            os.fsync(handle.fileno())
        fsync_directory(root)
        host_lock_acquired = True
        task_dir = owned_path(root / args.task_id, root)
        task_dir.mkdir(exist_ok=True)
        owned_path(task_dir, root)
        fsync_directory(root)
        lock_path = owned_path(task_dir / "RUNNING.lock", root)
        try:
            descriptor = os.open(lock_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
        except FileExistsError as exc:
            raise Hold("TASK_LOCK_EXISTS; report and review stale ownership, never auto-remove") from exc
        with os.fdopen(descriptor, "w") as handle:
            json.dump({"owner": args.owner, "pid": os.getpid(), "hostname": socket.gethostname(),
                       "stage": args.stage, "task_id": args.task_id}, handle)
            handle.flush()
            os.fsync(handle.fileno())
        fsync_directory(task_dir)
        lock_acquired = True
        for name in ("restart", "logs", "resume_journal"):
            directory = owned_path(task_dir / name, root)
            directory.mkdir(exist_ok=True)
            owned_path(directory, root)
        fsync_directory(task_dir)
        trajectory_path = owned_path(task_dir / "trajectory.traj", root)
        thermal_path = owned_path(task_dir / "thermal.csv", root)
        state_path = owned_path(task_dir / "sampler_state.json", root)
        source_manifest_path = owned_path(task_dir / "source_manifest.json", root)
        budget_path = owned_path(task_dir / "stage_wall_budget.json", root)
        start_time = time.monotonic()
        wall_value = sampling.get("max_wall_s", {})
        if isinstance(wall_value, dict):
            wall_limit = float(wall_value.get(args.stage, 0))
        else:
            wall_limit = float(sampling.get("warmup_max_wall_s", 0) if args.stage == "warmup" else wall_value)
        if not math.isfinite(wall_limit) or wall_limit <= 0:
            raise Hold("STAGE_WALL_BUDGET_NOT_FROZEN")
        archive_cap = int(sampling["trajectory_archive_max_bytes"])
        reserve = int(sampling.get("checkpoint_reserve_bytes", 32 * 1024**2))
        if reserve <= 0 or reserve >= archive_cap:
            raise Hold("CHECKPOINT_RESERVE_NOT_VALID")

        stage_base = {"warmup": 0.0, "continue": 0.0}
        unclean_charge_s = 0.0
        if budget_path.exists():
            owned_path(budget_path, root)
            previous_budget = json.loads(budget_path.read_text())
            if (previous_budget.get("task_id") != args.task_id or previous_budget.get("owner") != args.owner
                    or previous_budget.get("scientific_parameters_sha256") != scientific_sha):
                raise Hold("STAGE_WALL_LEDGER_IDENTITY_MISMATCH")
            for stage_name in stage_base:
                used = float(previous_budget["cumulative_wall_s_by_stage"][stage_name])
                if not math.isfinite(used) or used < 0:
                    raise Hold("INVALID_STAGE_WALL_ACCOUNTING")
                stage_base[stage_name] = used
            if previous_budget.get("active_invocation") is True:
                require_checkpoint_resume(config, args.task_id, state_path, root)
                unclean_charge_s = time.time() - float(previous_budget["updated_epoch_s"])
                if not math.isfinite(unclean_charge_s) or unclean_charge_s < 0:
                    raise Hold("UNCLEAN_WALL_CLOCK_ACCOUNTING_UNVERIFIED")
                active_stage = previous_budget["active_stage"]
                if active_stage not in stage_base:
                    raise Hold("UNKNOWN_ACTIVE_STAGE_IN_WALL_LEDGER")
                # An unclean interruption has no trustworthy stop timestamp.
                # Conservatively charge elapsed offline time until C reviews it;
                # a signed budget extension cannot reset already accrued usage.
                stage_base[active_stage] += unclean_charge_s
        elif state_path.exists():
            raise Hold("EXISTING_CHECKPOINT_HAS_NO_DURABLE_STAGE_WALL_LEDGER")

        def budget_snapshot(active: bool):
            cumulative = dict(stage_base)
            cumulative[args.stage] += time.monotonic() - start_time
            return {"task_id": args.task_id, "owner": args.owner,
                    "scientific_parameters_sha256": scientific_sha,
                    "cumulative_wall_s_by_stage": cumulative, "active_stage": args.stage,
                    "active_invocation": active, "updated_epoch_s": time.time(),
                    "unclean_elapsed_charge_s_this_resume": unclean_charge_s,
                    "C_release_config_sha256": raw_config_sha}

        def save_budget(active: bool):
            snapshot = budget_snapshot(active)
            atomic_json(budget_path, snapshot, root)
            return snapshot

        persist_budget = save_budget

        def used_stage_wall():
            return stage_base[args.stage] + time.monotonic() - start_time

        def resource_check(forecast=0, check_archive=True):
            if stop_requested[0]:
                raise Hold(stop_requested[0])
            if used_stage_wall() >= wall_limit:
                raise Hold("APPROVED_STAGE_WALL_BUDGET_EXCEEDED")
            if current_rss_bytes() > int(sampling["max_rss_bytes"]):
                raise Hold("APPROVED_RSS_BUDGET_EXCEEDED")
            if shutil.disk_usage(root).free < int(sampling["min_free_bytes"]) + forecast:
                raise Hold("SAMPLER_DISK_FREE_STOP")
            if check_archive and tree_bytes(task_dir) + forecast > archive_cap - reserve:
                raise Hold("TRAJECTORY_ARCHIVE_BUDGET_STOP; existing checkpoints are preserved")
            if budget_started:
                save_budget(True)

        resource_check()
        manifest = {
            "task_id": args.task_id, "owner": args.owner, "RNG_seed": int(row["RNG_seed"]),
            "trajectory_parent": row["trajectory_parent"], "source_root": source_root,
            "source_role": "liquid_pilot_candidate_not_training_or_independent_test",
            "sampler": "ASE EAM; exact frozen parameter file", "potential": inspection,
            "density_g_cm3": float(row["density_g_cm3"]), "cell_A": np.eye(3).tolist(),
            "box_edge_A": edge, "pbc": [True, True, True], "atom_ids": list(range(1, 65)),
            "initialization": "independent PCG64 periodic random packing; min distance 2.2 A",
            "scientific_parameters_sha256": scientific_sha, "task_matrix_sha256": task_sha,
            "initial_C_release_config_sha256": raw_config_sha, "fixed_identity": fixed_identity,
            "STRUCTURE_PASS": False, "liquid_QA": "NOT_PERFORMED", "DFT_LABELS": "ABSENT",
        }
        manifest["cell_A"] = (np.eye(3) * edge).tolist()
        if state_path.exists():
            owned_path(state_path, root)
            state = json.loads(state_path.read_text())
            if state["scientific_parameters_sha256"] != scientific_sha:
                raise Hold("CHECKPOINT_FIXED_PARAMETERS_OR_CODE_VERSION_MISMATCH")
            if state["task_id"] != args.task_id or state["owner"] != args.owner:
                raise Hold("CHECKPOINT_TASK_IDENTITY_MISMATCH")
            if state.get("status") not in {"WARMUP_COMPLETE_C_REVIEW_REQUIRED", "TRAJECTORY_COMPLETE_QA_REQUIRED"}:
                require_checkpoint_resume(config, args.task_id, state_path, root)
            if sha256_file(owned_path(source_manifest_path, root)) != state["source_manifest_sha256"]:
                raise Hold("SOURCE_MANIFEST_SHA_MISMATCH")
            bundle = owned_path(task_dir / state["restart_bundle"], task_dir)
            if set(state["bundle_sha256"]) != {"state_atoms.traj", "state_arrays.npz", "trajectory.traj", "thermal.csv"}:
                raise Hold("CHECKPOINT_REQUIRED_FILE_SET_MISMATCH")
            for name, expected in state["bundle_sha256"].items():
                if sha256_file(owned_path(bundle / name, bundle)) != expected:
                    raise Hold("CHECKPOINT_SHA_MISMATCH:" + name)
            # Exact restart uses raw binary arrays, avoiding reapplication of
            # momentum constraints by an Atoms/trajectory reader constructor.
            with np.load(bundle / "state_arrays.npz", allow_pickle=False) as arrays:
                atoms = Atoms(numbers=arrays["numbers"], positions=arrays["positions"],
                              cell=arrays["cell"], pbc=arrays["pbc"])
                atoms.set_masses(arrays["masses"])
                atoms.set_momenta(arrays["momenta"], apply_constraint=False)
                atoms.new_array("atom_ids", arrays["atom_ids"])
            completed_step = int(state["step"])
            if args.stage == "continue" and completed_step < sampling["warmup_steps"]:
                raise Hold("WARMUP_NOT_COMPLETE")
            if args.stage == "warmup" and completed_step >= sampling["warmup_steps"]:
                print(json.dumps({"status": "WARMUP_ALREADY_COMPLETE", "step": completed_step}))
                return 0
            if args.stage == "continue" and completed_step >= 120000:
                print(json.dumps({"status": "TRAJECTORY_ALREADY_COMPLETE_QA_REQUIRED", "step": completed_step}))
                return 0
            thermostat_rng = np.random.Generator(np.random.PCG64())
            thermostat_rng.bit_generator.state = state["rng_state"]["thermostat"]
            initial_rng_states = state["rng_state"]
            # Restore exact committed output snapshots; preserve uncommitted tails.
            for live, name in ((trajectory_path, "trajectory.traj"), (thermal_path, "thermal.csv")):
                owned_path(live, root)
                if live.exists() and sha256_file(live) != state["bundle_sha256"][name]:
                    archive = owned_path(task_dir / "resume_journal" / (live.name + ".uncommitted-" + uuid.uuid4().hex), root)
                    copy_owned(live, archive, root)
                temporary = owned_path(live.with_name(live.name + ".restore-" + uuid.uuid4().hex), root)
                copy_owned(bundle / name, temporary, root)
                owned_path(live, root)
                os.replace(temporary, live)
                fsync_directory(live.parent)
        else:
            if args.stage != "warmup":
                raise Hold("CONTINUE_REQUIRES_EXISTING_WARMUP_CHECKPOINT")
            if any(path.exists() for path in (trajectory_path, thermal_path, task_dir / "source_manifest.json")):
                raise Hold("UNREGISTERED_EXISTING_SAMPLER_OUTPUT; no restart or overwrite allowed")
            budget_started = True
            save_budget(True)
            positions, packing_rng, velocity_rng, thermostat_rng, attempts = pack_random(
                int(row["RNG_seed"]), edge, int(sampling["packing_max_attempts"]), np, resource_check)
            atoms = Atoms("Ag64", positions=positions, cell=np.eye(3) * edge, pbc=True)
            atoms.new_array("atom_ids", np.arange(1, 65, dtype=np.int64))
            atoms.set_masses(np.full(64, 107.8682))
            atoms.set_constraint(FixCom())
            MaxwellBoltzmannDistribution(atoms, temperature_K=sampling["T_hot_K"],
                                          force_temp=True, rng=velocity_rng)
            initial_rng_states = {"packing": packing_rng.bit_generator.state,
                                  "velocity": velocity_rng.bit_generator.state,
                                  "thermostat": thermostat_rng.bit_generator.state}
            manifest["packing_attempts"] = attempts
            atomic_json(source_manifest_path, manifest, root)
            owned_path(thermal_path, root)
            with thermal_path.open("x", newline="") as handle:
                csv.DictWriter(handle, fieldnames=THERMAL_FIELDS, lineterminator="\n").writeheader()
                handle.flush()
                os.fsync(handle.fileno())
            fsync_directory(task_dir)
        if (len(atoms) != 64 or not (atoms.numbers == 47).all() or not atoms.pbc.all()
                or not np.array_equal(atoms.get_array("atom_ids"), np.arange(1, 65))
                or not np.array_equal(atoms.cell.array, np.eye(3) * edge)):
            raise Hold("RESTART_COMPOSITION_CELL_OR_PBC_MISMATCH")
        momenta = atoms.get_momenta().copy()
        atoms.set_constraint(FixCom())
        atoms.set_momenta(momenta, apply_constraint=False)
        atoms.calc = calc
        _, temperature = schedule(completed_step, sampling)
        dyn = Langevin(atoms, timestep=sampling["dt_fs"] * units.fs,
                       temperature_K=temperature,
                       friction=sampling["friction_per_fs"] / units.fs,
                       fixcm=False, rng=thermostat_rng)
        dyn.nsteps = completed_step
        budget_started = True
        save_budget(True)
        owned_path(trajectory_path, root)
        trajectory = Trajectory(str(trajectory_path), "a", properties=[])
        owned_path(thermal_path, root)
        thermal_handle = thermal_path.open("a", newline="")
        thermal_writer = csv.DictWriter(thermal_handle, fieldnames=THERMAL_FIELDS, lineterminator="\n")

        def save_frame():
            owned_path(trajectory_path, root)
            phase, target = schedule(completed_step, sampling)
            frame = atoms.copy()
            frame.calc = None
            frame.info = {"step": completed_step, "time_fs": completed_step * sampling["dt_fs"],
                          "phase": phase, "target_temperature_K": target,
                          "task_id": args.task_id, "trajectory_parent": row["trajectory_parent"],
                          "source_root": source_root, "atom_ids": list(range(1, 65)),
                          "positions_are_unwrapped": True, "DFT_LABELS": "ABSENT"}
            trajectory.write(frame)

        def save_thermal():
            owned_path(thermal_path, root)
            energy = float(atoms.get_potential_energy())
            kinetic = float(atoms.get_kinetic_energy())
            stress = atoms.get_stress(voigt=False)
            pressure = -float(np.trace(stress)) / 3 + 2 * kinetic / (3 * atoms.get_volume())
            velocity = atoms.get_velocities()
            com_velocity = (atoms.get_masses()[:, None] * velocity).sum(axis=0) / atoms.get_masses().sum()
            phase, target = schedule(completed_step, sampling)
            data = {"step": completed_step, "time_fs": completed_step * sampling["dt_fs"],
                    "phase": phase, "target_temperature_K": target,
                    "T_K": float(atoms.get_temperature()), "sampler_energy_eV": energy,
                    "pressure_eV_A3": pressure,
                    "COMspeed_A_fs": float(np.linalg.norm(com_velocity)) * units.fs}
            if not all(math.isfinite(float(value)) for key, value in data.items() if key != "phase"):
                raise Hold("SAMPLER_THERMAL_NAN_OR_INF")
            thermal_writer.writerow(data)
            thermal_handle.flush()

        def save_checkpoint(status: str, reason: str | None = None):
            # A completed Langevin step is the only restart boundary; transient
            # random position/velocity increments are recomputed from saved RNG.
            thermal_handle.flush()
            os.fsync(thermal_handle.fileno())
            if trajectory.backend.fd is not None:
                trajectory.backend.fd.flush()
                os.fsync(trajectory.backend.fd.fileno())
            owned_path(trajectory_path, root)
            owned_path(thermal_path, root)
            suffix = f"step-{completed_step:09d}-" + uuid.uuid4().hex
            bundle = owned_path(task_dir / "restart" / suffix, root)
            bundle.mkdir()
            owned_path(bundle, root)
            fsync_directory(bundle.parent)
            state_atoms = atoms.copy()
            state_atoms.calc = None
            state_atoms.info = {"step": completed_step, "task_id": args.task_id,
                                "atom_ids": list(range(1, 65)), "positions_are_unwrapped": True}
            state_atoms_path = owned_path(bundle / "state_atoms.traj", root)
            with Trajectory(str(state_atoms_path), "w", properties=[]) as writer:
                writer.write(state_atoms)
            arrays_path = owned_path(bundle / "state_arrays.npz", root)
            with arrays_path.open("xb") as handle:
                np.savez_compressed(handle, numbers=atoms.numbers, positions=atoms.positions,
                                    momenta=atoms.get_momenta(), masses=atoms.get_masses(),
                                    cell=atoms.cell.array, pbc=atoms.pbc,
                                    atom_ids=atoms.get_array("atom_ids"))
                handle.flush()
                os.fsync(handle.fileno())
            copy_owned(trajectory_path, bundle / "trajectory.traj", root)
            copy_owned(thermal_path, bundle / "thermal.csv", root)
            rng_states = dict(initial_rng_states)
            rng_states["thermostat"] = thermostat_rng.bit_generator.state
            files = ["state_atoms.traj", "state_arrays.npz", "trajectory.traj", "thermal.csv"]
            for name in files:
                fsync_owned_file(bundle / name, root)
            active = status in {"INITIALIZED", "SAMPLING_IN_PROGRESS"}
            wall_snapshot = save_budget(active)
            state = {"task_id": args.task_id, "owner": args.owner, "step": completed_step,
                     "time_fs": completed_step * sampling["dt_fs"], "status": status,
                     "reason": reason, "restart_bundle": str(bundle.relative_to(task_dir)),
                     "rng_state": rng_states, "rng_class": "numpy.Generator(PCG64)",
                     "langevin_version": Langevin._lgv_version, "packages": packages,
                     "scientific_parameters_sha256": scientific_sha,
                     "C_release_config_sha256_at_checkpoint": raw_config_sha,
                     "task_matrix_sha256": task_sha,
                     "source_manifest_sha256": sha256_file(owned_path(source_manifest_path, root)),
                     "bundle_sha256": {name: sha256_file(bundle / name) for name in files},
                     "wall_s_this_invocation": time.monotonic() - start_time,
                     "cumulative_wall_s_by_stage": wall_snapshot["cumulative_wall_s_by_stage"],
                     "STRUCTURE_PASS": False, "QA_status": "NOT_PERFORMED_BY_SAMPLER"}
            atomic_json(bundle / "state.json", state, root)
            # Commit durability order: every bundle file, bundle directory,
            # parent restart directory, then the atomic current-state pointer.
            fsync_directory(bundle)
            fsync_directory(bundle.parent)
            atomic_json(state_path, state, root)
            fsync_directory(task_dir)
            lines = [f"{sha256_file(owned_path(path, root))}  {path.relative_to(task_dir)}" for path in (
                bundle / "state.json", *(bundle / name for name in files),
                task_dir / "source_manifest.json", state_path)]
            atomic_json(task_dir / "logs" / f"integrity-{suffix}.json", {"sha256_lines": lines}, root)

        checkpoint = save_checkpoint
        if completed_step == 0:
            save_frame()
            save_thermal()
            save_checkpoint("INITIALIZED")
        for signum in (signal.SIGTERM, signal.SIGINT):
            previous_handlers[signum] = signal.getsignal(signum)
            signal.signal(signum, lambda number, _frame: stop_requested.__setitem__(0, "STOP_SIGNAL_" + str(number)))
        end_step = sampling["warmup_steps"] if args.stage == "warmup" else sum(
            sampling[name] for name in ("hot_steps", "cool_steps", "equil_steps", "production_steps"))
        invocation_start_step = completed_step
        while completed_step < end_step:
            # Save a valid boundary in memory so an interrupted/failed force
            # evaluation cannot turn a half-step into a purported restart.
            if completed_step % sampling["thermal_every"] == 0:
                resource_check(check_archive=False)
            elif stop_requested[0] or used_stage_wall() >= wall_limit:
                raise Hold(stop_requested[0] or "APPROVED_STAGE_WALL_BUDGET_EXCEEDED")
            old_positions = atoms.positions.copy()
            old_momenta = atoms.get_momenta().copy()
            old_rng = json.loads(json.dumps(thermostat_rng.bit_generator.state))
            _, target = schedule(completed_step, sampling)
            dyn.set_temperature(temperature_K=target)
            try:
                forces = dyn.step()
                if not (np.isfinite(atoms.positions).all() and np.isfinite(atoms.get_momenta()).all()
                        and np.isfinite(forces).all()):
                    raise Hold("SAMPLER_STATE_NAN_OR_INF")
                density_values = np.asarray(calc.total_density)
                if (not np.isfinite(density_values).all() or density_values.min() < -1e-8
                        or density_values.max() > inspection["density_grid_max"] + 1e-8):
                    raise Hold("EAM_EMBEDDING_DENSITY_OUTSIDE_PARAMETER_GRID")
            except BaseException:
                atoms.set_positions(old_positions, apply_constraint=False)
                atoms.set_momenta(old_momenta, apply_constraint=False)
                thermostat_rng.bit_generator.state = old_rng
                calc.reset()
                raise
            completed_step += 1
            dyn.nsteps = completed_step
            if completed_step % sampling["thermal_every"] == 0:
                minimum = nearest_distance(atoms.positions, edge, np)
                threshold = float(sampling.get("min_runtime_distance_A", 0))
                if threshold <= 0:
                    raise Hold("MIN_RUNTIME_DISTANCE_STOP_NOT_FROZEN")
                if minimum < threshold:
                    raise Hold("SAMPLER_SHORT_DISTANCE_STOP")
                save_thermal()
            if completed_step % sampling["write_every"] == 0:
                save_frame()
            if completed_step % int(sampling["checkpoint_every"]) == 0 and completed_step < end_step:
                forecast = trajectory_path.stat().st_size + thermal_path.stat().st_size + 16384
                resource_check(forecast=forecast)
                save_checkpoint("SAMPLING_IN_PROGRESS")
        status = "WARMUP_COMPLETE_C_REVIEW_REQUIRED" if args.stage == "warmup" else "TRAJECTORY_COMPLETE_QA_REQUIRED"
        save_checkpoint(status)
        elapsed = time.monotonic() - start_time
        advanced = completed_step - invocation_start_step
        summary = {"task_id": args.task_id, "owner": args.owner, "stage": args.stage,
                   "status": status, "step": completed_step, "advanced_steps": advanced,
                   "wall_s": elapsed, "steps_per_s": advanced / elapsed if elapsed else None,
                   "cumulative_wall_s_by_stage": budget_snapshot(False)["cumulative_wall_s_by_stage"],
                   "projected_remaining_wall_s": (120000 - completed_step) * elapsed / advanced if advanced else None,
                   "rss_bytes": current_rss_bytes(), "archive_bytes": tree_bytes(task_dir),
                   "disk_free_bytes": shutil.disk_usage(root).free,
                   "STRUCTURE_PASS": False, "DFT_STARTED": False}
        atomic_json(task_dir / "logs" / (args.stage + "-performance-" + uuid.uuid4().hex + ".json"), summary, root)
        print(json.dumps(summary, indent=2))
        return 0
    except Exception as exc:
        reason = str(exc)
        if checkpoint is not None:
            try:
                checkpoint("STOPPED_REVIEW_REQUIRED", reason)
            except Exception as save_error:
                reason += "; checkpoint_save_failure=" + str(save_error)
        result = {"status": "SAMPLING_HOLD", "task_id": args.task_id, "owner": args.owner,
                  "stage": args.stage, "last_completed_step": completed_step, "reason": reason,
                  "STRUCTURE_PASS": False, "DFT_STARTED": False}
        if task_dir is not None and lock_acquired:
            try:
                atomic_json(task_dir / "logs" / ("failure-" + uuid.uuid4().hex + ".json"), result, root)
            except Exception as log_error:
                result["safe_log_write_failure"] = str(log_error)
        print(json.dumps(result, ensure_ascii=False, indent=2), file=sys.stderr)
        return 2
    finally:
        if budget_started and persist_budget is not None:
            try:
                persist_budget(False)
            except Exception as budget_error:
                print(json.dumps({"status": "WALL_LEDGER_CLOSE_HOLD", "reason": str(budget_error)}), file=sys.stderr)
        if trajectory is not None:
            trajectory.close()
        if thermal_handle is not None:
            thermal_handle.close()
        for signum, handler in previous_handlers.items():
            signal.signal(signum, handler)
        if lock_path is not None and lock_acquired:
            # Only the successfully acquired lock belongs to this invocation.
            try:
                owned_path(lock_path, root)
                identity = json.loads(lock_path.read_text())
                if identity.get("pid") == os.getpid() and identity.get("hostname") == socket.gethostname():
                    lock_path.unlink()
                    fsync_directory(lock_path.parent)
            except (FileNotFoundError, json.JSONDecodeError, Hold):
                pass
        if host_lock_path is not None and host_lock_acquired:
            try:
                owned_path(host_lock_path, root)
                identity = json.loads(host_lock_path.read_text())
                if identity.get("pid") == os.getpid() and identity.get("hostname") == socket.gethostname():
                    host_lock_path.unlink()
                    fsync_directory(host_lock_path.parent)
            except (FileNotFoundError, json.JSONDecodeError, Hold):
                pass


if __name__ == "__main__":
    raise SystemExit(main())
