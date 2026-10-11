"""Offline-only bounded-segment protocol candidate for the D Ag64 task.

This module deliberately has no ASE, EAM, MACE, GPAW, LAMMPS, or MD imports.
Its functional stepper interface is exercised only with synthetic state in
test_segment_guard.py. It is not connected to the production sampler and the
unsigned release files in this directory must fail the release gate.
"""

from __future__ import annotations

import copy
import hashlib
import io
import json
import math
import os
import subprocess
import tarfile
import tempfile
import time
import uuid
import shutil
from pathlib import Path, PurePosixPath
from typing import Any, Callable, Mapping


class GateError(RuntimeError):
    """A fail-closed identity, integrity, execution, or archive gate."""


EXPECTED_OWNER = "D"
EXPECTED_TASK = "C-D-AGL-1250K-rho097-seed01"
START_STEP = 1000
STOP_STEP = 2000
MAX_ADDED_STEPS = 1000
MAX_ARCHIVE_MEMBER_BYTES = 50 * 1024 * 1024
EXPECTED_ARCHIVE_SHA256 = "b62a818e58a2cf423e4cb1b9ed0468ae73ad578a75c4b2cb9fb3f28baf75e7b1"
EXPECTED_ARCHIVE_MANIFEST_SHA256 = "ab9f0dffa8ef145eac7792dfe1e14e07c1354df1c2152dafec97bad1dc3a57b1"
EXPECTED_CHECKPOINT_STATE_SHA256 = "e08574678b3f8231c0a52ab25febb13c9d06e815b28149eb306d4b451c83ef6f"
EXPECTED_SOURCE_MANIFEST_SHA256 = "9a694ea8f93655efa7a1c8695669ab2f17520c2d1b90f6df28115b2b1f9c2ad6"
EXPECTED_SCIENCE_SHA256 = "8f792d6206c1b7bae7b900e7ea81445d6cd3e837de2b5979c88bb977650ab7b2"
EXPECTED_POTENTIAL_SHA256 = "958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447"


def sha256_file(path: os.PathLike[str] | str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _json_bytes(path: os.PathLike[str] | str) -> tuple[dict[str, Any], bytes]:
    raw = Path(path).read_bytes()
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise GateError("JSON_ROOT_MUST_BE_OBJECT")
    return value, raw


def _safe_member(name: str) -> bool:
    p = PurePosixPath(name)
    return bool(name) and not p.is_absolute() and all(part not in {"", ".", ".."} for part in p.parts)


def _tar_member_bytes(tf: tarfile.TarFile, member: tarfile.TarInfo) -> bytes:
    stream = tf.extractfile(member)
    if stream is None:
        raise GateError(f"ARCHIVE_MEMBER_NOT_REGULAR_FILE:{member.name}")
    return stream.read()


def verify_original_archive(
    archive_path: os.PathLike[str] | str,
    manifest_path: os.PathLike[str] | str,
    *,
    expected_archive_sha256: str = EXPECTED_ARCHIVE_SHA256,
    expected_manifest_sha256: str = EXPECTED_ARCHIVE_MANIFEST_SHA256,
) -> dict[str, Any]:
    """Recompute the original D archive and every manifest member from bytes."""
    archive_path, manifest_path = Path(archive_path), Path(manifest_path)
    actual_archive_sha = sha256_file(archive_path)
    actual_manifest_sha = sha256_file(manifest_path)
    if actual_archive_sha != expected_archive_sha256:
        raise GateError("ORIGINAL_ARCHIVE_SHA256_MISMATCH")
    if actual_manifest_sha != expected_manifest_sha256:
        raise GateError("ORIGINAL_ARCHIVE_MANIFEST_SHA256_MISMATCH")
    manifest, _ = _json_bytes(manifest_path)
    if manifest.get("archive_sha256") != actual_archive_sha:
        raise GateError("MANIFEST_ARCHIVE_SHA256_MISMATCH")
    if manifest.get("archive_bytes") != archive_path.stat().st_size:
        raise GateError("MANIFEST_ARCHIVE_SIZE_MISMATCH")
    if (manifest.get("owner"), manifest.get("task_id"), manifest.get("step")) != (
        EXPECTED_OWNER, EXPECTED_TASK, START_STEP
    ):
        raise GateError("ORIGINAL_ARCHIVE_IDENTITY_MISMATCH")
    if manifest.get("potential_file_included") is not False:
        raise GateError("POTENTIAL_FILE_MUST_NOT_BE_IN_ARCHIVE")

    expected = {item["path"]: item for item in manifest.get("members", [])}
    if len(expected) != len(manifest.get("members", [])):
        raise GateError("DUPLICATE_MANIFEST_MEMBER")
    actual_names: set[str] = set()
    with tarfile.open(archive_path, "r:gz") as tf:
        for item in tf.getmembers():
            if not _safe_member(item.name):
                raise GateError(f"UNSAFE_ARCHIVE_PATH:{item.name}")
            if not item.isfile() or item.issym() or item.islnk():
                raise GateError(f"NONREGULAR_ARCHIVE_MEMBER:{item.name}")
            if item.name in actual_names:
                raise GateError(f"DUPLICATE_ARCHIVE_MEMBER:{item.name}")
            actual_names.add(item.name)
            if item.size > MAX_ARCHIVE_MEMBER_BYTES:
                raise GateError(f"ARCHIVE_MEMBER_OVER_50MIB:{item.name}")
            if item.name not in expected:
                raise GateError(f"UNMANIFESTED_ARCHIVE_MEMBER:{item.name}")
            data = _tar_member_bytes(tf, item)
            entry = expected[item.name]
            member_sha = sha256_bytes(data)
            if len(data) != entry.get("size_bytes") or member_sha != entry.get("sha256"):
                raise GateError(f"ARCHIVE_MEMBER_INTEGRITY_FAILURE:{item.name}")
        if actual_names != set(expected):
            raise GateError("ARCHIVE_MEMBER_SET_MISMATCH")
        by_name = {member.name: member for member in tf.getmembers()}
        sampler = json.loads(_tar_member_bytes(tf, by_name["sampler_state.json"]))
        source = _tar_member_bytes(tf, by_name["source_manifest.json"])
        restart_root = sampler.get("restart_bundle", "")
        arrays_bytes = _tar_member_bytes(tf, by_name[restart_root + "/state_arrays.npz"])
        state_atoms_bytes = _tar_member_bytes(tf, by_name[restart_root + "/state_atoms.traj"])
        checkpoint_state_path = sampler.get("restart_bundle", "") + "/state.json"
        checkpoint_state_bytes = _tar_member_bytes(tf, by_name[checkpoint_state_path])
        checkpoint_state = json.loads(checkpoint_state_bytes)

    source_json = json.loads(source)
    if sha256_bytes(source) != EXPECTED_SOURCE_MANIFEST_SHA256:
        raise GateError("ORIGINAL_SOURCE_MANIFEST_SHA256_MISMATCH")
    if sampler.get("source_manifest_sha256") != EXPECTED_SOURCE_MANIFEST_SHA256:
        raise GateError("SAMPLER_SOURCE_MANIFEST_IDENTITY_MISMATCH")
    if sampler.get("step") != START_STEP or sampler.get("task_id") != EXPECTED_TASK:
        raise GateError("SAMPLER_CHECKPOINT_STEP_OR_TASK_MISMATCH")
    if sampler.get("owner") != EXPECTED_OWNER or sampler.get("status") != "WARMUP_COMPLETE_C_REVIEW_REQUIRED":
        raise GateError("SAMPLER_CHECKPOINT_OWNER_OR_STATUS_MISMATCH")
    if sampler.get("rng_class") != "numpy.Generator(PCG64)":
        raise GateError("RNG_CLASS_MISMATCH")
    rng = sampler.get("rng_state")
    if not isinstance(rng, dict) or set(rng) != {"packing", "velocity", "thermostat"}:
        raise GateError("RNG_STATE_SET_INCOMPLETE")
    if any(not isinstance(v, dict) or v.get("bit_generator") != "PCG64" for v in rng.values()):
        raise GateError("RNG_STATE_NOT_PCG64")
    if checkpoint_state.get("step") != START_STEP or checkpoint_state.get("task_id") != EXPECTED_TASK:
        raise GateError("CHECKPOINT_STATE_IDENTITY_MISMATCH")
    if checkpoint_state.get("scientific_parameters_sha256") != EXPECTED_SCIENCE_SHA256:
        raise GateError("CHECKPOINT_SCIENTIFIC_PARAMETERS_CHANGED")
    if checkpoint_state.get("source_manifest_sha256") != EXPECTED_SOURCE_MANIFEST_SHA256:
        raise GateError("CHECKPOINT_SOURCE_MANIFEST_CHANGED")
    if checkpoint_state.get("rng_state") != rng:
        raise GateError("SAMPLER_AND_CHECKPOINT_RNG_STATE_DIFFER")
    if sha256_bytes(checkpoint_state_bytes) != EXPECTED_CHECKPOINT_STATE_SHA256:
        raise GateError("ORIGINAL_CHECKPOINT_STATE_SHA256_MISMATCH")
    expected_bundle = checkpoint_state.get("bundle_sha256", {})
    for basename, data in (("state_arrays.npz", arrays_bytes), ("state_atoms.traj", state_atoms_bytes)):
        if sha256_bytes(data) != expected_bundle.get(basename):
            raise GateError(f"CHECKPOINT_BUNDLE_SHA_MISMATCH:{basename}")
    try:
        import numpy as np
        with np.load(io.BytesIO(arrays_bytes), allow_pickle=False) as arrays:
            required = {"numbers", "positions", "momenta", "velocities_A_per_fs", "masses", "cell", "pbc", "atom_ids"}
            if set(arrays.files) != required:
                raise GateError("CHECKPOINT_ARRAY_MEMBERS_MISMATCH")
            finite_keys = ("positions", "momenta", "velocities_A_per_fs", "masses", "cell")
            if any(not np.isfinite(arrays[key]).all() for key in finite_keys):
                raise GateError("CHECKPOINT_ARRAYS_NONFINITE")
            if any(arrays[key].shape != (64, 3) for key in ("positions", "momenta", "velocities_A_per_fs")):
                raise GateError("CHECKPOINT_VECTOR_ARRAY_SHAPE_MISMATCH")
            if arrays["numbers"].shape != (64,) or not np.all(arrays["numbers"] == 47):
                raise GateError("CHECKPOINT_NOT_64_AG_ATOMS")
            if not np.array_equal(arrays["atom_ids"], np.asarray(source_json["atom_ids"])):
                raise GateError("CHECKPOINT_ATOM_ID_ORDER_CHANGED")
            if not np.allclose(arrays["cell"], np.asarray(source_json["cell_A"]), rtol=0, atol=1e-12):
                raise GateError("CHECKPOINT_CELL_CHANGED")
            if not np.array_equal(arrays["pbc"], np.asarray(source_json["pbc"])):
                raise GateError("CHECKPOINT_PBC_CHANGED")
    except GateError:
        raise
    except Exception as exc:
        raise GateError(f"CHECKPOINT_ARRAY_READ_FAILURE:{type(exc).__name__}") from exc
    if manifest.get("potential_sha256") != EXPECTED_POTENTIAL_SHA256:
        raise GateError("ORIGINAL_POTENTIAL_IDENTITY_MISMATCH")
    return {
        "status": "ORIGINAL_ARCHIVE_BYTES_VERIFIED_READ_ONLY",
        "archive_sha256": actual_archive_sha,
        "archive_bytes": archive_path.stat().st_size,
        "manifest_sha256": actual_manifest_sha,
        "member_count": len(expected),
        "checkpoint_step": START_STEP,
        "task_id": EXPECTED_TASK,
        "source_manifest_sha256": EXPECTED_SOURCE_MANIFEST_SHA256,
        "checkpoint_state_sha256": sha256_bytes(checkpoint_state_bytes),
        "state_arrays_sha256": sha256_bytes(arrays_bytes),
        "state_atoms_traj_sha256": sha256_bytes(state_atoms_bytes),
        "atom_count": 64,
        "elements": ["Ag"] * 64,
        "rng_streams": sorted(rng),
    }


def verify_release_inputs(
    *,
    proposal_path: os.PathLike[str] | str,
    release_lock_path: os.PathLike[str] | str,
    archive_path: os.PathLike[str] | str,
    archive_manifest_path: os.PathLike[str] | str,
    code_paths: Mapping[str, os.PathLike[str] | str],
    owner: str,
    task_id: str,
    requested_stop_step: int | None,
    expected_release_lock_sha256: str,
    potential_path: os.PathLike[str] | str,
) -> dict[str, Any]:
    """Fail-closed preflight; all file hashes are recomputed from local bytes.

    A real run must obtain release_lock from the specifically C-approved immutable
    commit. The proposal in this package is unsigned, so it must always HOLD.
    """
    config, config_raw = _json_bytes(proposal_path)
    lock, lock_raw = _json_bytes(release_lock_path)
    actual_lock_sha = sha256_bytes(lock_raw)
    if actual_lock_sha != expected_release_lock_sha256:
        raise GateError("RELEASE_LOCK_SHA256_MISMATCH")
    if config.get("release_status") != "SIGNED_BY_C" or lock.get("approval_status") != "SIGNED_BY_C":
        raise GateError("RELEASE_NOT_SIGNED_BY_C")
    if not lock.get("signed_commit"):
        raise GateError("SIGNED_RELEASE_COMMIT_MISSING")
    actual_config_sha = sha256_bytes(config_raw)
    if actual_config_sha != lock.get("proposal_sha256"):
        raise GateError("ACTUAL_CONFIG_SHA256_MISMATCH")
    if owner != EXPECTED_OWNER or task_id != EXPECTED_TASK:
        raise GateError("CALLER_OWNER_OR_TASK_NOT_ALLOWED")
    if (config.get("owner"), config.get("task_id")) != (owner, task_id):
        raise GateError("CONFIG_OWNER_OR_TASK_MISMATCH")
    if requested_stop_step != STOP_STEP or requested_stop_step != config.get("stop_at_step_required"):
        raise GateError("EXPLICIT_STOP_STEP_2000_REQUIRED")
    if (config.get("start_step"), config.get("stop_at_step_required"), config.get("max_added_steps")) != (
        START_STEP, STOP_STEP, MAX_ADDED_STEPS
    ):
        raise GateError("SEGMENT_BOUNDARY_MISMATCH")
    if config.get("stage_continue_default_endpoint_warning") != 120000:
        raise GateError("DEFAULT_120000_ENDPOINT_NOT_BOUND")
    actual_code = {name: sha256_file(path) for name, path in code_paths.items()}
    if actual_code != lock.get("code_sha256"):
        raise GateError("ACTUAL_EXECUTION_CODE_SHA256_MISMATCH")
    if actual_code != config.get("execution_code_sha256"):
        raise GateError("CONFIG_CODE_SHA256_MISMATCH")
    archive_result = verify_original_archive(archive_path, archive_manifest_path)
    if archive_result["archive_sha256"] != config.get("checkpoint_archive_sha256"):
        raise GateError("CONFIG_CHECKPOINT_ARCHIVE_SHA256_MISMATCH")
    if archive_result["manifest_sha256"] != config.get("checkpoint_archive_manifest_sha256"):
        raise GateError("CONFIG_CHECKPOINT_MANIFEST_SHA256_MISMATCH")
    actual_potential = sha256_file(potential_path)
    if actual_potential != EXPECTED_POTENTIAL_SHA256 or actual_potential != config.get("potential_sha256_must_remain"):
        raise GateError("LOCAL_POTENTIAL_SHA256_MISMATCH")
    if config.get("scientific_parameters_sha256_must_remain") != EXPECTED_SCIENCE_SHA256:
        raise GateError("SCIENTIFIC_PARAMETERS_SHA256_MISMATCH")
    return {
        "status": "PREFLIGHT_PASS_NOT_EXECUTION_AUTHORIZATION",
        "release_lock_sha256": actual_lock_sha,
        "proposal_sha256": actual_config_sha,
        "code_sha256": actual_code,
        "archive": archive_result,
        "potential_sha256": actual_potential,
        "owner": owner,
        "task_id": task_id,
        "start_step": START_STEP,
        "stop_step": STOP_STEP,
    }


def assert_step_allowed(completed_step: int, *, start_step: int, stop_at_step: int, max_added_steps: int) -> int:
    if start_step != START_STEP or stop_at_step != STOP_STEP or max_added_steps != MAX_ADDED_STEPS:
        raise GateError("ONLY_EXACT_1000_TO_2000_SEGMENT_SUPPORTED")
    if completed_step < start_step:
        raise GateError("CHECKPOINT_START_STEP_NOT_REACHED")
    if completed_step >= stop_at_step:
        raise GateError("STOP_AT_STEP_REACHED_NO_NEXT_INTEGRATION")
    return completed_step + 1


def validate_resume_state(state: Mapping[str, Any]) -> None:
    """Require the complete immutable step-1000 identity before stepping."""
    if (state.get("owner"), state.get("task_id"), state.get("step")) != (EXPECTED_OWNER, EXPECTED_TASK, START_STEP):
        raise GateError("RESUME_STATE_OWNER_TASK_OR_STEP_MISMATCH")
    if state.get("scientific_parameters_sha256") != EXPECTED_SCIENCE_SHA256:
        raise GateError("RESUME_STATE_SCIENCE_SHA_MISMATCH")
    if state.get("source_manifest_sha256") != EXPECTED_SOURCE_MANIFEST_SHA256:
        raise GateError("RESUME_STATE_SOURCE_SHA_MISMATCH")
    required_arrays = ("positions_A", "momenta", "velocities_A_per_fs")
    for name in required_arrays:
        value = state.get(name)
        if not isinstance(value, list) or len(value) != 64 or any(not isinstance(row, list) or len(row) != 3 for row in value):
            raise GateError(f"RESUME_STATE_{name.upper()}_MISSING_OR_BAD_SHAPE")
    if state.get("atom_ids") != list(range(1, 65)):
        raise GateError("RESUME_STATE_ATOM_IDS_OR_ORDER_CHANGED")
    if state.get("pbc") != [True, True, True]:
        raise GateError("RESUME_STATE_PBC_CHANGED")
    expected_cell = 10.826824658074
    cell = state.get("cell_A")
    if not isinstance(cell, list) or len(cell) != 3 or any(not isinstance(row, list) or len(row) != 3 for row in cell):
        raise GateError("RESUME_STATE_CELL_MISSING_OR_BAD_SHAPE")
    for i in range(3):
        for j in range(3):
            expected = expected_cell if i == j else 0.0
            if not math.isclose(float(cell[i][j]), expected, rel_tol=0, abs_tol=1e-12):
                raise GateError("RESUME_STATE_CELL_CHANGED")
    rng = state.get("rng_state")
    if not isinstance(rng, dict) or set(rng) != {"packing", "velocity", "thermostat"}:
        raise GateError("RESUME_STATE_RNG_STREAMS_INCOMPLETE")
    if any(not isinstance(item, dict) or item.get("bit_generator") != "PCG64" for item in rng.values()):
        raise GateError("RESUME_STATE_RNG_NOT_PCG64")
    if not _all_finite(state):
        raise GateError("RESUME_STATE_NONFINITE")


def _all_finite(value: Any) -> bool:
    if isinstance(value, bool) or value is None or isinstance(value, (str, bytes)):
        return True
    if isinstance(value, (int, float)):
        return math.isfinite(float(value))
    if isinstance(value, Mapping):
        return all(_all_finite(v) for v in value.values())
    if isinstance(value, (list, tuple)):
        return all(_all_finite(v) for v in value)
    # Numpy arrays and similar test/adapter state use .tolist() without importing them here.
    tolist = getattr(value, "tolist", None)
    return _all_finite(tolist()) if callable(tolist) else True


def run_bounded_segment(
    initial_state: Mapping[str, Any],
    advance_one: Callable[[dict[str, Any]], None],
    save_checkpoint: Callable[[dict[str, Any]], None],
    *,
    stop_at_step: int,
    max_added_steps: int,
    max_wall_seconds: float,
    clock: Callable[[], float] = time.monotonic,
    cancel_requested: Callable[[], bool] = lambda: False,
    save_every: int = 100,
) -> dict[str, Any]:
    """Run a functional adapter step-by-step; this is only exercised by fakes.

    Each trial step is applied to a deep copy, so failed/NaN mutations cannot
    corrupt the last verified checkpoint state. The production sampler is not
    wired to this function.
    """
    state = copy.deepcopy(dict(initial_state))
    validate_resume_state(state)
    if stop_at_step != STOP_STEP or max_added_steps != MAX_ADDED_STEPS:
        raise GateError("SEGMENT_LIMIT_MUST_BE_1000_STEPS")
    if stop_at_step - state["step"] != max_added_steps:
        raise GateError("SEGMENT_STEP_BUDGET_MISMATCH")
    if save_every <= 0 or (stop_at_step - state["step"]) % save_every != 0:
        raise GateError("CHECKPOINT_CADENCE_MUST_DIVIDE_SEGMENT")
    started = clock()
    last_saved_step = state["step"]
    while state["step"] < stop_at_step:
        if cancel_requested():
            raise GateError("SEGMENT_CANCELLED; LAST_CHECKPOINT_RETAINED")
        if clock() - started >= max_wall_seconds:
            raise GateError("SEGMENT_WALL_LIMIT; LAST_CHECKPOINT_RETAINED")
        next_step = assert_step_allowed(
            state["step"], start_step=START_STEP,
            stop_at_step=stop_at_step, max_added_steps=max_added_steps,
        )
        trial = copy.deepcopy(state)
        advance_one(trial)
        # If a completed fake/adapter step consumes the remaining budget, do
        # not commit its candidate state or start the next integration step.
        if clock() - started >= max_wall_seconds:
            raise GateError("SEGMENT_WALL_LIMIT_AFTER_STEP; LAST_CHECKPOINT_RETAINED")
        trial["step"] = next_step
        if not _all_finite(trial):
            raise GateError("NONFINITE_STATE; LAST_CHECKPOINT_RETAINED")
        state = trial
        if state["step"] % save_every == 0 or state["step"] == stop_at_step:
            save_checkpoint(copy.deepcopy(state))
            last_saved_step = state["step"]
    if state["step"] != STOP_STEP or last_saved_step != STOP_STEP:
        raise GateError("EXACT_FINAL_CHECKPOINT_NOT_SAVED")
    return {"status": "SYNTHETIC_SEGMENT_BOUNDARY_PASS_ONLY", "start_step": START_STEP,
            "end_step": state["step"], "added_steps": state["step"] - START_STEP}


def atomic_write_bytes(path: os.PathLike[str] | str, data: bytes) -> None:
    """Write a fresh file atomically; existing outputs are immutable and refused."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise GateError("REFUSE_TO_OVERWRITE_EXISTING_ARCHIVE_OR_CHECKPOINT")
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temp_name, path)
        os.unlink(temp_name)
        dir_fd = os.open(path.parent, os.O_DIRECTORY)
        try:
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)
    except Exception:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def write_checkpoint_callback(output_root: os.PathLike[str] | str) -> Callable[[dict[str, Any]], None]:
    root = Path(output_root)
    root.mkdir(parents=True, exist_ok=True)

    def save(state: dict[str, Any]) -> None:
        payload = json.dumps(state, sort_keys=True, separators=(",", ":"), allow_nan=False).encode() + b"\n"
        atomic_write_bytes(root / f"state-step-{state['step']:07d}.json", payload)
    return save


def atomic_publish_checkpoint_bundle(
    checkpoint_root: os.PathLike[str] | str,
    final_name: str,
    members: Mapping[str, bytes],
) -> dict[str, Any]:
    """Publish a complete checkpoint directory atomically, never overwrite.

    Expected production adapter members include sampler_state.json, state.json,
    state_arrays.npz, state_atoms.traj, thermal.csv, and trajectory.traj. Tests
    use only synthetic bytes. A stale lock is intentionally a hard stop.
    """
    root = Path(checkpoint_root)
    root.mkdir(parents=True, exist_ok=True)
    if not final_name or Path(final_name).name != final_name or final_name in {".", ".."}:
        raise GateError("INVALID_CHECKPOINT_BUNDLE_NAME")
    final = root / final_name
    lock = root / f".{final_name}.publish-lock"
    try:
        lock.mkdir()
    except FileExistsError as exc:
        raise GateError("CHECKPOINT_PUBLISH_LOCK_EXISTS; STOP") from exc
    if final.exists():
        lock.rmdir()
        raise GateError("REFUSE_TO_OVERWRITE_CHECKPOINT_BUNDLE")
    stage = root / f".{final_name}.{uuid.uuid4().hex}.staging"
    try:
        stage.mkdir()
        payload: dict[str, bytes] = {}
        for name, data in members.items():
            if not _safe_member(name) or name in {"MANIFEST.json", "SHA256SUMS.txt"}:
                raise GateError(f"INVALID_CHECKPOINT_MEMBER_PATH:{name}")
            if len(data) > MAX_ARCHIVE_MEMBER_BYTES:
                raise GateError(f"CHECKPOINT_MEMBER_OVER_50MIB:{name}")
            payload[name] = data
        if not payload:
            raise GateError("EMPTY_CHECKPOINT_BUNDLE")
        manifest = {"schema": 1, "members": [
            {"path": name, "size_bytes": len(data), "sha256": sha256_bytes(data)}
            for name, data in sorted(payload.items())
        ]}
        manifest_bytes = (json.dumps(manifest, sort_keys=True, indent=2) + "\n").encode()
        checksum_payload = {**payload, "MANIFEST.json": manifest_bytes}
        sums_bytes = "".join(f"{sha256_bytes(data)}  {name}\n" for name, data in sorted(checksum_payload.items())).encode()
        checksum_payload["SHA256SUMS.txt"] = sums_bytes
        for name, data in checksum_payload.items():
            target = stage.joinpath(*PurePosixPath(name).parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            with open(target, "xb") as stream:
                stream.write(data)
                stream.flush()
                os.fsync(stream.fileno())
        for current, dirs, _files in os.walk(stage, topdown=False):
            for dirname in dirs:
                fd = os.open(Path(current) / dirname, os.O_DIRECTORY)
                try:
                    os.fsync(fd)
                finally:
                    os.close(fd)
            fd = os.open(current, os.O_DIRECTORY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
        if final.exists():
            raise GateError("CHECKPOINT_BUNDLE_APPEARED_DURING_PUBLISH")
        os.rename(stage, final)
        parent_fd = os.open(root, os.O_DIRECTORY)
        try:
            os.fsync(parent_fd)
        finally:
            os.close(parent_fd)
        lock.rmdir()
        return {"status": "ATOMIC_CHECKPOINT_BUNDLE_PUBLISHED", "path": final_name,
                "member_count": len(checksum_payload), "sha256": {
                    name: sha256_file(final.joinpath(*PurePosixPath(name).parts)) for name in checksum_payload
                }}
    except Exception:
        if stage.exists():
            shutil.rmtree(stage)
        # A normal exception releases the lock; a process crash leaves it as a fail-closed marker.
        try:
            lock.rmdir()
        except OSError:
            pass
        raise


def make_archive_bytes(payload: Mapping[str, bytes], *, task_id: str, owner: str, final_step: int) -> bytes:
    """Build a deterministic test/archive tar with README, MANIFEST and sums."""
    if task_id != EXPECTED_TASK or owner != EXPECTED_OWNER or final_step != STOP_STEP:
        raise GateError("SUCCESSOR_ARCHIVE_IDENTITY_OR_STEP_MISMATCH")
    normalized: dict[str, bytes] = {}
    for name, data in payload.items():
        if not _safe_member(name) or name in {"README.md", "MANIFEST.json", "SHA256SUMS.txt"}:
            raise GateError(f"INVALID_OR_RESERVED_ARCHIVE_PATH:{name}")
        if "Ag.eam.alloy" in PurePosixPath(name).name:
            raise GateError("POTENTIAL_FILE_MUST_NOT_BE_ARCHIVED")
        if len(data) > MAX_ARCHIVE_MEMBER_BYTES:
            raise GateError(f"ARCHIVE_MEMBER_OVER_50MIB:{name}")
        normalized[name] = data
    generated: dict[str, bytes] = {
        "README.md": (
            f"# {task_id} segment archive\n\nOwner: {owner}. Global step: {final_step}.\n"
            "This is a synthetic test archive unless its provenance manifest says otherwise.\n"
            "Potential parameter files are excluded.\n"
        ).encode()
    }
    generated.update(normalized)
    manifest = {
        "schema": 1, "owner": owner, "task_id": task_id, "final_step": final_step,
        "members": [{"path": name, "size_bytes": len(data), "sha256": sha256_bytes(data)}
                    for name, data in sorted(generated.items())],
        "potential_file_included": False,
    }
    generated["MANIFEST.json"] = (json.dumps(manifest, sort_keys=True, indent=2) + "\n").encode()
    sums = "".join(f"{sha256_bytes(data)}  {name}\n" for name, data in sorted(generated.items()))
    generated["SHA256SUMS.txt"] = sums.encode()
    result = bytearray()
    with tempfile.TemporaryFile() as fileobj:
        with tarfile.open(fileobj=fileobj, mode="w:gz", format=tarfile.PAX_FORMAT) as tf:
            for name, data in sorted(generated.items()):
                info = tarfile.TarInfo(name)
                info.size = len(data)
                info.mtime = 0
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                import io
                tf.addfile(info, io.BytesIO(data))
        fileobj.seek(0)
        result.extend(fileobj.read())
    return bytes(result)


def verify_successor_archive(archive_bytes: bytes) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="ag64-archive-qa-") as tmp:
        archive = Path(tmp) / "candidate.tar.gz"
        archive.write_bytes(archive_bytes)
        with tarfile.open(archive, "r:gz") as tf:
            members = tf.getmembers()
            by_name: dict[str, bytes] = {}
            for member in members:
                if not _safe_member(member.name) or not member.isfile() or member.size > MAX_ARCHIVE_MEMBER_BYTES:
                    raise GateError("SUCCESSOR_ARCHIVE_MEMBER_POLICY_FAILURE")
                if member.name in by_name:
                    raise GateError("SUCCESSOR_ARCHIVE_DUPLICATE_MEMBER")
                by_name[member.name] = _tar_member_bytes(tf, member)
        for required in ("README.md", "MANIFEST.json", "SHA256SUMS.txt"):
            if required not in by_name:
                raise GateError(f"SUCCESSOR_ARCHIVE_MISSING_{required}")
        if any(PurePosixPath(name).name == "Ag.eam.alloy" for name in by_name):
            raise GateError("POTENTIAL_FILE_MUST_NOT_BE_ARCHIVED")
        manifest = json.loads(by_name["MANIFEST.json"])
        expected = {item["path"]: item for item in manifest["members"]}
        if set(expected) != set(by_name) - {"MANIFEST.json", "SHA256SUMS.txt"}:
            raise GateError("SUCCESSOR_ARCHIVE_MANIFEST_MEMBER_SET_MISMATCH")
        for name, entry in expected.items():
            data = by_name[name]
            if len(data) != entry["size_bytes"] or sha256_bytes(data) != entry["sha256"]:
                raise GateError(f"SUCCESSOR_ARCHIVE_MEMBER_SHA_FAILURE:{name}")
        expected_sums = "".join(
            f"{sha256_bytes(by_name[name])}  {name}\n"
            for name in sorted(set(by_name) - {"SHA256SUMS.txt"})
        ).encode()
        if by_name["SHA256SUMS.txt"] != expected_sums:
            raise GateError("SUCCESSOR_SHA256SUMS_CONTENT_MISMATCH")
        if (manifest.get("task_id"), manifest.get("owner"), manifest.get("final_step")) != (
            EXPECTED_TASK, EXPECTED_OWNER, STOP_STEP
        ):
            raise GateError("SUCCESSOR_ARCHIVE_IDENTITY_MISMATCH")
        return {"archive_sha256": sha256_bytes(archive_bytes), "archive_bytes": len(archive_bytes),
                "member_count": len(by_name), "final_step": STOP_STEP}


def local_git_roundtrip(archive_bytes: bytes, work_parent: os.PathLike[str] | str) -> dict[str, Any]:
    """Non-scientific local bare-Git push/fresh-clone/readback simulation."""
    parent = Path(work_parent)
    parent.mkdir(parents=True, exist_ok=True)
    bare = parent / "remote.git"
    work = parent / "writer"
    fresh = parent / "fresh-reader"
    for path in (bare, work, fresh):
        if path.exists():
            raise GateError("ROUNDTRIP_PATH_MUST_BE_NEW")
    subprocess.run(["git", "init", "--bare", str(bare)], check=True, capture_output=True)
    subprocess.run(["git", "clone", str(bare), str(work)], check=True, capture_output=True)
    subprocess.run(["git", "-C", str(work), "config", "user.name", "offline-test"], check=True)
    subprocess.run(["git", "-C", str(work), "config", "user.email", "offline-test@example.invalid"], check=True)
    bundle = work / "archive" / "segment.tar.gz"
    bundle.parent.mkdir(parents=True)
    bundle.write_bytes(archive_bytes)
    subprocess.run(["git", "-C", str(work), "add", "archive/segment.tar.gz"], check=True)
    subprocess.run(["git", "-C", str(work), "commit", "-m", "synthetic offline archive"], check=True,
                   capture_output=True)
    subprocess.run(["git", "-C", str(work), "branch", "-M", "test-segment-archive"], check=True)
    subprocess.run(["git", "-C", str(work), "push", "-u", "origin", "test-segment-archive"], check=True,
                   capture_output=True)
    subprocess.run(["git", "clone", "--branch", "test-segment-archive", str(bare), str(fresh)], check=True,
                   capture_output=True)
    fetched = (fresh / "archive" / "segment.tar.gz").read_bytes()
    if sha256_bytes(fetched) != sha256_bytes(archive_bytes):
        raise GateError("FRESH_GIT_READBACK_ARCHIVE_SHA_FAILURE")
    internal = verify_successor_archive(fetched)
    return {"status": "LOCAL_GIT_FRESH_CLONE_READBACK_PASS_NOT_GITHUB", "commit": subprocess.check_output(
        ["git", "-C", str(fresh), "rev-parse", "HEAD"], text=True).strip(),
        "archive_sha256": sha256_bytes(fetched), "internal": internal}
