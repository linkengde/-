"""Bounded ASE Langevin restart adapter for the D Ag64 pilot.

The production path is locked to global steps 1000 -> 2000.  Tests inject a
virtual calculator; this module never chooses or downloads a potential.
"""
from __future__ import annotations

import argparse
import copy
import ctypes
import csv
import errno
import hashlib
import json
import os
import tarfile
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from ase import Atoms, units
from ase.constraints import FixCom
from ase.md.langevin import Langevin


TASK_ID = "C-D-AGL-1250K-rho097-seed01"
OWNER = "D"
START_STEP = 1000
STOP_AT_STEP = 2000
MAX_ADDED_STEPS = STOP_AT_STEP - START_STEP
EXPECTED_ARCHIVE_SHA256 = "b62a818e58a2cf423e4cb1b9ed0468ae73ad578a75c4b2cb9fb3f28baf75e7b1"
EXPECTED_ARCHIVE_MANIFEST_SHA256 = "ab9f0dffa8ef145eac7792dfe1e14e07c1354df1c2152dafec97bad1dc3a57b1"
EXPECTED_STATE_SHA256 = "e08574678b3f8231c0a52ab25febb13c9d06e815b28149eb306d4b451c83ef6f"
EXPECTED_ARRAYS_SHA256 = "39da7d05b186ce61134e14dcc53c19311cfbda183a2d1403e1451d38d0205466"
EXPECTED_SOURCE_MANIFEST_SHA256 = "9a694ea8f93655efa7a1c8695669ab2f17520c2d1b90f6df28115b2b1f9c2ad6"
EXPECTED_SCIENTIFIC_PARAMETERS_SHA256 = "8f792d6206c1b7bae7b900e7ea81445d6cd3e837de2b5979c88bb977650ab7b2"
EXPECTED_OLD_RUNNER_SHA256 = "6b54d094c26bb9d4a5e4eb5997a31955eb8829a2975042e23fe10527cb943583"
EXPECTED_POTENTIAL_SHA256 = "958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447"
EXPECTED_CELL_EDGE_A = 10.826824658074
DT_FS = 0.5
FRICTION_PER_FS = 0.01
T_HOT_K = 1600.0
T_TARGET_K = 1250.0
HOT_STEPS = 20000
COOL_STEPS = 20000
EQUIL_STEPS = 40000
MAX_ARCHIVE_MEMBER_BYTES = 50 * 1024 * 1024
SIGNED_SAMPLING_PARAMETERS = {
    "atom_count": 64,
    "density_g_cm3": 9.032733418760,
    "box_edge_A": EXPECTED_CELL_EDGE_A,
    "seed": 125009701,
    "T_hot_K": T_HOT_K,
    "T_target_K": T_TARGET_K,
    "dt_fs": DT_FS,
    "friction_per_fs": FRICTION_PER_FS,
    "hot_steps": HOT_STEPS,
    "cool_steps": COOL_STEPS,
    "equil_steps": EQUIL_STEPS,
    "production_steps": 40000,
    "checkpoint_every": 100,
    "thermal_every": 20,
    "write_every": 100,
    "max_wall_s": 600,
}


class SegmentRefused(RuntimeError):
    """Identity, signed-boundary, resource, or numerical gate refused work."""


def _rename_noreplace(source: Path, destination: Path) -> None:
    """Atomically publish a directory without replacing any existing path."""
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    _require(renameat2 is not None, "ATOMIC_NOREPLACE_RENAME_UNAVAILABLE")
    renameat2.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    renameat2.restype = ctypes.c_int
    result = renameat2(-100, os.fsencode(source), -100, os.fsencode(destination), 1)
    if result != 0:
        error = ctypes.get_errno()
        if error == errno.EEXIST:
            raise SegmentRefused("OUTPUT_BUNDLE_ALREADY_EXISTS_NO_OVERWRITE")
        raise OSError(error, os.strerror(error), str(destination))


def _fsync_dir(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")


def schedule_target_K(global_step: int) -> tuple[str, float]:
    """Exact global-step schedule from the original sampler."""
    if global_step < HOT_STEPS:
        return "hot_melt_candidate", T_HOT_K
    if global_step < HOT_STEPS + COOL_STEPS:
        fraction = (global_step - HOT_STEPS) / COOL_STEPS
        return "cooling", T_HOT_K + fraction * (T_TARGET_K - T_HOT_K)
    if global_step < HOT_STEPS + COOL_STEPS + EQUIL_STEPS:
        return "equilibration_candidate", T_TARGET_K
    return "production_candidate", T_TARGET_K


def _require(condition: bool, code: str) -> None:
    if not condition:
        raise SegmentRefused(code)


@dataclass
class RestoredCheckpoint:
    atoms: Atoms
    thermostat_rng: np.random.Generator
    rng_state_all: dict
    global_step: int
    sampler_state: dict
    arrays_sha256: str
    state_sha256: str
    thermal_rows: list[dict] | None = None
    trajectory_frames: list[dict] | None = None


def restore_checkpoint(arrays_path: str | Path, state_path: str | Path) -> RestoredCheckpoint:
    """Restore raw ASE arrays and RNG without drawing or projecting velocities."""
    arrays_path, state_path = Path(arrays_path), Path(state_path)
    arrays_sha = sha256_file(arrays_path)
    state_sha = sha256_file(state_path)
    _require(arrays_sha == EXPECTED_ARRAYS_SHA256, "STEP1000_ARRAYS_SHA_MISMATCH")
    _require(state_sha == EXPECTED_STATE_SHA256, "STEP1000_STATE_SHA_MISMATCH")
    state = json.loads(state_path.read_text(encoding="utf-8"))
    _require(state.get("task_id") == TASK_ID and state.get("owner") == OWNER,
             "STEP1000_TASK_OR_OWNER_MISMATCH")
    _require(state.get("step") == START_STEP, "STEP1000_GLOBAL_STEP_MISMATCH")
    _require(state.get("status") == "WARMUP_COMPLETE_C_REVIEW_REQUIRED",
             "STEP1000_STATE_STATUS_MISMATCH")
    _require(state.get("scientific_parameters_sha256") == EXPECTED_SCIENTIFIC_PARAMETERS_SHA256,
             "SOURCE_SCIENTIFIC_IDENTITY_MISMATCH")
    _require(state.get("source_manifest_sha256") == EXPECTED_SOURCE_MANIFEST_SHA256,
             "SOURCE_MANIFEST_SHA_MISMATCH")
    _require(state.get("bundle_sha256", {}).get("state_arrays.npz") == arrays_sha,
             "STEP1000_STATE_ARRAYS_INNER_SHA_MISMATCH")
    rng_states = state.get("rng_state")
    _require(isinstance(rng_states, dict) and set(rng_states) == {"packing", "velocity", "thermostat"},
             "PCG64_STREAM_SET_MISMATCH")
    _require(all(rng_states[name].get("bit_generator") == "PCG64" for name in rng_states),
             "PCG64_GENERATOR_TYPE_MISMATCH")

    required = {"numbers", "positions", "momenta", "velocities_A_per_fs", "masses",
                "cell", "pbc", "atom_ids"}
    with np.load(arrays_path, allow_pickle=False) as stored:
        _require(set(stored.files) == required, "STEP1000_ARRAY_KEY_SET_MISMATCH")
        arrays = {name: stored[name].copy() for name in required}
    numbers = arrays["numbers"]
    positions = arrays["positions"]
    momenta = arrays["momenta"]
    masses = arrays["masses"]
    cell = arrays["cell"]
    pbc = arrays["pbc"]
    atom_ids = arrays["atom_ids"]
    velocities = arrays["velocities_A_per_fs"]
    _require(numbers.shape == (64,) and np.all(numbers == 47), "STEP1000_COMPOSITION_MISMATCH")
    _require(positions.shape == momenta.shape == velocities.shape == (64, 3),
             "STEP1000_VECTOR_ARRAY_SHAPE_MISMATCH")
    _require(masses.shape == (64,) and np.isfinite(masses).all() and np.all(masses > 0),
             "STEP1000_MASS_ARRAY_INVALID")
    _require(cell.shape == (3, 3) and pbc.shape == (3,), "STEP1000_CELL_PBC_SHAPE_MISMATCH")
    _require(np.array_equal(pbc, np.ones(3, dtype=bool)), "STEP1000_PBC_MISMATCH")
    _require(np.array_equal(cell, np.eye(3) * EXPECTED_CELL_EDGE_A), "STEP1000_CELL_MISMATCH")
    _require(np.array_equal(atom_ids, np.arange(1, 65, dtype=np.int64)), "STEP1000_ATOM_ID_ORDER_MISMATCH")
    _require(np.isfinite(positions).all() and np.isfinite(momenta).all() and np.isfinite(velocities).all(),
             "STEP1000_STATE_NONFINITE")
    _require(np.allclose(momenta, masses[:, None] * velocities, rtol=2e-15, atol=2e-15),
             "STEP1000_MOMENTA_VELOCITY_INCONSISTENCY")

    atoms = Atoms(numbers=numbers, positions=positions, cell=cell, pbc=pbc)
    atoms.set_masses(masses)
    atoms.set_momenta(momenta, apply_constraint=False)
    atoms.new_array("atom_ids", atom_ids)
    # This reproduces the old restart path: FixCom with Langevin fixcm=False.
    # Existing momenta are deliberately not reprojected on restoration.
    atoms.set_constraint(FixCom())
    atoms.set_momenta(momenta, apply_constraint=False)
    _require(np.array_equal(atoms.positions, positions), "RESTORE_CHANGED_POSITIONS")
    _require(np.array_equal(atoms.get_momenta(), momenta), "RESTORE_CHANGED_MOMENTA")
    _require(np.array_equal(atoms.get_masses(), masses), "RESTORE_CHANGED_MASSES")
    _require(np.array_equal(atoms.get_velocities(), velocities), "RESTORE_CHANGED_VELOCITIES")

    thermostat_rng = np.random.Generator(np.random.PCG64())
    thermostat_rng.bit_generator.state = copy.deepcopy(rng_states["thermostat"])
    # Validate, then preserve all three source streams. Only thermostat is used
    # by ASE Langevin; packing and velocity streams remain provenance state.
    for stream in rng_states.values():
        probe = np.random.Generator(np.random.PCG64())
        probe.bit_generator.state = copy.deepcopy(stream)
    return RestoredCheckpoint(atoms, thermostat_rng, copy.deepcopy(rng_states), START_STEP,
                              state, arrays_sha, state_sha)


def restore_segment_bundle(bundle_dir: str | Path) -> RestoredCheckpoint:
    """Restore a later immutable segment bundle after checking every digest."""
    root = Path(bundle_dir)
    manifest = json.loads((root / "MANIFEST.json").read_text(encoding="utf-8"))
    _require(manifest.get("task_id") == TASK_ID and manifest.get("owner") == OWNER,
             "SEGMENT_BUNDLE_TASK_OR_OWNER_MISMATCH")
    sums = {}
    for line in (root / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        sums[name] = digest
    expected_sum_names = set(manifest.get("members_sha256", {})) | {"MANIFEST.json"}
    _require(set(sums) == expected_sum_names, "SEGMENT_SUMS_MEMBER_SET_MISMATCH")
    _require(sums.get("MANIFEST.json") == sha256_file(root / "MANIFEST.json"),
             "SEGMENT_MANIFEST_SHA_MISMATCH")
    for name, expected in manifest.get("members_sha256", {}).items():
        path = root / name
        _require(path.is_file() and sha256_file(path) == expected == sums.get(name),
                 "SEGMENT_MEMBER_SHA_MISMATCH:" + name)
    state = json.loads((root / "state.json").read_text(encoding="utf-8"))
    step = state.get("step")
    _require(state.get("task_id") == TASK_ID and state.get("owner") == OWNER,
             "SEGMENT_STATE_TASK_OR_OWNER_MISMATCH")
    _require(isinstance(step, int) and START_STEP < step <= STOP_AT_STEP,
             "SEGMENT_RESTART_STEP_OUTSIDE_1000_2000")
    _require(manifest.get("global_step") == step, "SEGMENT_MANIFEST_STEP_MISMATCH")
    _require(state.get("execution_code_sha256") == sha256_file(Path(__file__)),
             "SEGMENT_RESTART_EXECUTOR_SHA_MISMATCH")
    with np.load(root / "state_arrays.npz", allow_pickle=False) as stored:
        arrays = {key: stored[key].copy() for key in stored.files}
    required = {"numbers", "positions", "momenta", "velocities_ASE_internal_units", "masses",
                "cell", "pbc", "atom_ids"}
    _require(set(arrays) == required, "SEGMENT_RESTART_ARRAY_SET_MISMATCH")
    _require(arrays["numbers"].shape == (64,) and np.all(arrays["numbers"] == 47),
             "SEGMENT_RESTART_COMPOSITION_MISMATCH")
    _require(arrays["positions"].shape == arrays["momenta"].shape == (64, 3)
             and arrays["masses"].shape == (64,), "SEGMENT_RESTART_ARRAY_SHAPE_MISMATCH")
    _require(arrays["cell"].shape == (3, 3) and arrays["pbc"].shape == (3,),
             "SEGMENT_RESTART_CELL_PBC_SHAPE_MISMATCH")
    _require(np.array_equal(arrays["cell"], np.eye(3) * EXPECTED_CELL_EDGE_A)
             and np.array_equal(arrays["pbc"], np.ones(3, dtype=bool)),
             "SEGMENT_RESTART_CELL_OR_PBC_MISMATCH")
    _require(np.array_equal(arrays["atom_ids"], np.arange(1, 65, dtype=np.int64)),
             "SEGMENT_RESTART_ATOM_ID_ORDER_MISMATCH")
    _require(np.isfinite(arrays["positions"]).all() and np.isfinite(arrays["momenta"]).all()
             and np.isfinite(arrays["masses"]).all() and np.all(arrays["masses"] > 0),
             "SEGMENT_RESTART_NONFINITE_OR_INVALID_ARRAY")
    rng_state = json.loads((root / "rng_state.json").read_text(encoding="utf-8"))
    _require(set(rng_state) == {"packing", "velocity", "thermostat"},
             "SEGMENT_RESTART_RNG_STREAM_SET_MISMATCH")
    atoms = Atoms(numbers=arrays["numbers"], positions=arrays["positions"],
                  cell=arrays["cell"], pbc=arrays["pbc"])
    atoms.set_masses(arrays["masses"])
    atoms.set_momenta(arrays["momenta"], apply_constraint=False)
    atoms.new_array("atom_ids", arrays["atom_ids"])
    _require(np.array_equal(atoms.get_velocities(), arrays["velocities_ASE_internal_units"]),
             "SEGMENT_RESTART_VELOCITY_MISMATCH")
    atoms.set_constraint(FixCom())
    atoms.set_momenta(arrays["momenta"], apply_constraint=False)
    thermostat_rng = np.random.Generator(np.random.PCG64())
    thermostat_rng.bit_generator.state = copy.deepcopy(rng_state["thermostat"])
    for stream in rng_state.values():
        probe = np.random.Generator(np.random.PCG64())
        probe.bit_generator.state = copy.deepcopy(stream)
    with (root / "thermal.csv").open(newline="", encoding="utf-8") as handle:
        thermal_rows = []
        for row in csv.DictReader(handle):
            for name in ("step",):
                row[name] = int(row[name])
            for name in ("time_fs", "target_temperature_K", "T_K", "sampler_energy_eV",
                         "kinetic_energy_eV", "pressure_eV_A3"):
                row[name] = float(row[name])
            thermal_rows.append(row)
    with np.load(root / "trajectory.npz", allow_pickle=False) as trajectory:
        _require(trajectory["positions"].shape[1:] == (64, 3)
                 and trajectory["steps"].shape == (trajectory["positions"].shape[0],),
                 "SEGMENT_TRAJECTORY_HISTORY_SHAPE_MISMATCH")
        trajectory_frames = [{"step": int(step), "positions": positions.copy()}
                             for step, positions in zip(trajectory["steps"], trajectory["positions"])]
    return RestoredCheckpoint(atoms, thermostat_rng, rng_state, step, state,
                              sha256_file(root / "state_arrays.npz"),
                              sha256_file(root / "state.json"), thermal_rows,
                              trajectory_frames)


def snapshot_state(atoms: Atoms, global_step: int,
                   rng_state_all: dict, thermostat_rng: np.random.Generator,
                   source_provenance: dict) -> dict:
    all_states = copy.deepcopy(rng_state_all)
    all_states["thermostat"] = copy.deepcopy(thermostat_rng.bit_generator.state)
    return {
        "task_id": TASK_ID,
        "owner": OWNER,
        "step": int(global_step),
        "positions": atoms.get_positions().copy(),
        "momenta": atoms.get_momenta().copy(),
        "velocities_A_per_fs": atoms.get_velocities().copy(),
        "masses": atoms.get_masses().copy(),
        "cell": atoms.cell.array.copy(),
        "pbc": atoms.pbc.copy(),
        "numbers": atoms.numbers.copy(),
        "atom_ids": atoms.get_array("atom_ids").copy(),
        "rng_state": all_states,
        "source_provenance": copy.deepcopy(source_provenance),
    }


class BoundedLangevin(Langevin):
    """ASE Langevin whose actual step method owns and enforces the hard stop."""

    def __init__(self, atoms: Atoms, thermostat_rng: np.random.Generator,
                 start_step: int = START_STEP, stop_at_step: int = STOP_AT_STEP,
                 initial_step: int = START_STEP, on_step: Callable | None = None,
                 synthetic_test_mode: bool = False,
                 release_lock_path: str | Path | None = None,
                 release_lock_sha256: str | None = None,
                 execution_code_sha256: str | None = None):
        _require(start_step == START_STEP and stop_at_step == STOP_AT_STEP,
                 "ONLY_SIGNED_1000_TO_2000_SEGMENT_SUPPORTED")
        _require(START_STEP <= initial_step <= STOP_AT_STEP,
                 "SEGMENT_MUST_START_WITHIN_GLOBAL_STEP_1000_2000")
        _require(ase_version_is_329(), "ASE_VERSION_MUST_BE_3_29_0")
        if synthetic_test_mode:
            _require(bool(getattr(atoms.calc, "codex_synthetic_test_calculator", False)),
                     "SYNTHETIC_MODE_REQUIRES_MARKED_VIRTUAL_CALCULATOR")
            self.release_lock = None
        else:
            _require(release_lock_path is not None and release_lock_sha256 is not None
                     and execution_code_sha256 is not None,
                     "SIGNED_RELEASE_LOCK_REQUIRED_BEFORE_INTEGRATOR_CONSTRUCTION")
            self.release_lock = validate_release_lock(release_lock_path, release_lock_sha256,
                                                      execution_code_sha256)
        self.synthetic_test_mode = synthetic_test_mode
        super().__init__(atoms, timestep=DT_FS * units.fs, temperature_K=T_HOT_K,
                         friction=FRICTION_PER_FS / units.fs, fixcm=False, rng=thermostat_rng)
        self.nsteps = int(initial_step)
        self.segment_start = int(start_step)
        self.stop_at_step = int(stop_at_step)
        self.on_step = on_step

    def step(self, forces=None):
        current = int(self.nsteps)
        _require(self.segment_start <= current < self.stop_at_step,
                 "SIGNED_STOP_AT_STEP_REACHED_OR_STEP_OUT_OF_RANGE")
        phase, target = schedule_target_K(current)
        _require(phase == "hot_melt_candidate" and target == T_HOT_K,
                 "GLOBAL_TEMPERATURE_SCHEDULE_MISMATCH")
        self.set_temperature(temperature_K=target)
        old_positions = self.atoms.get_positions().copy()
        old_momenta = self.atoms.get_momenta().copy()
        old_rng_state = copy.deepcopy(self.rng.bit_generator.state)
        try:
            new_forces = super().step(forces)
            _require(np.isfinite(new_forces).all() and np.isfinite(self.atoms.positions).all()
                     and np.isfinite(self.atoms.get_momenta()).all(), "INTEGRATOR_NONFINITE_STATE")
        except BaseException:
            self.atoms.set_positions(old_positions, apply_constraint=False)
            self.atoms.set_momenta(old_momenta, apply_constraint=False)
            self.rng.bit_generator.state = old_rng_state
            if self.atoms.calc is not None:
                self.atoms.calc.reset()
            raise
        self.nsteps = current + 1
        if self.on_step is not None:
            self.on_step(self, phase, target, new_forces)
        return new_forces

    def run(self, *args, **kwargs):
        raise SegmentRefused("ASE_RUN_API_DISABLED_USE_BOUNDED_STEP_LOOP")

    def irun(self, *args, **kwargs):
        raise SegmentRefused("ASE_IRUN_API_DISABLED_USE_BOUNDED_STEP_LOOP")


def ase_version_is_329() -> bool:
    import ase
    return ase.__version__ == "3.29.0"


def advance_to_stop(dyn: BoundedLangevin, *, stop_at_step: int = STOP_AT_STEP,
                    budget_check: Callable[[int], None] | None = None,
                    max_wall_s: float | None = None,
                    clock: Callable[[], float] = time.monotonic) -> int:
    _require(stop_at_step == STOP_AT_STEP, "STOP_BOUNDARY_NOT_2000")
    _require(START_STEP <= dyn.nsteps <= stop_at_step, "CONTINUATION_GLOBAL_STEP_MISMATCH")
    start = clock()
    while dyn.nsteps < stop_at_step:
        if budget_check is not None:
            budget_check(dyn.nsteps)
        if max_wall_s is not None and clock() - start >= max_wall_s:
            raise SegmentRefused("SEGMENT_WALLTIME_BUDGET_EXCEEDED")
        dyn.step()
    _require(dyn.nsteps == stop_at_step, "SEGMENT_DID_NOT_STOP_EXACTLY_AT_2000")
    return dyn.nsteps


def publish_and_verify(bundle_dir: str | Path, archive_path: str | Path,
                       remote_roundtrip: Callable[[Path, str], None]) -> str:
    """Publish archive, then require caller's independent remote readback.

    Any upload/readback exception stops the caller. The verified local bundle,
    archive, and original step-1000 checkpoint are left untouched.
    """
    archive_sha = create_archive(bundle_dir, archive_path)
    remote_roundtrip(Path(archive_path), archive_sha)
    return archive_sha


def write_checkpoint_bundle(output_root: str | Path, label: str, snapshot: dict,
                           thermal_rows: list[dict] | None = None,
                           trajectory_frames: list[dict] | None = None) -> Path:
    """Write one immutable, internally hashed directory; never overwrite."""
    root = Path(output_root)
    root.mkdir(parents=True, exist_ok=True)
    step = int(snapshot["step"])
    _require(step > START_STEP and step <= STOP_AT_STEP, "OUTPUT_STEP_OUTSIDE_SIGNED_SEGMENT")
    destination = root / f"{label}-step-{step:09d}"
    _require(not destination.exists(), "OUTPUT_BUNDLE_ALREADY_EXISTS_NO_OVERWRITE")
    staging = Path(tempfile.mkdtemp(prefix=".segment-staging-", dir=root))
    try:
        arrays = staging / "state_arrays.npz"
        with arrays.open("xb") as handle:
            payload = {key: snapshot[key] for key in (
                "numbers", "positions", "momenta", "masses", "cell", "pbc", "atom_ids")}
            payload["velocities_ASE_internal_units"] = snapshot["velocities_A_per_fs"]
            np.savez_compressed(handle, **payload)
            handle.flush()
            os.fsync(handle.fileno())
        rng_path = staging / "rng_state.json"
        rng_path.write_bytes(canonical_json_bytes(snapshot["rng_state"]))
        thermal_rows = thermal_rows or []
        thermal_path = staging / "thermal.csv"
        fields = ["step", "time_fs", "phase", "target_temperature_K", "T_K",
                  "sampler_energy_eV", "kinetic_energy_eV", "pressure_eV_A3"]
        with thermal_path.open("x", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(thermal_rows)
            handle.flush()
            os.fsync(handle.fileno())
        trajectory_frames = trajectory_frames or []
        trajectory_path = staging / "trajectory.npz"
        frame_positions = np.stack([frame["positions"] for frame in trajectory_frames]) if trajectory_frames else np.empty((0, 64, 3))
        frame_steps = np.asarray([frame["step"] for frame in trajectory_frames], dtype=np.int64)
        with trajectory_path.open("xb") as handle:
            np.savez_compressed(handle, steps=frame_steps, positions=frame_positions,
                                cell=np.asarray(snapshot["cell"]), pbc=np.asarray(snapshot["pbc"]),
                                atom_ids=np.asarray(snapshot["atom_ids"]))
            handle.flush()
            os.fsync(handle.fileno())
        state = {key: snapshot[key] for key in ("task_id", "owner", "step")}
        state["source_provenance"] = snapshot["source_provenance"]
        state.update({"global_time_fs": step * DT_FS,
                      "scientific_parameters_sha256_source": EXPECTED_SCIENTIFIC_PARAMETERS_SHA256,
                      "execution_code_sha256": sha256_file(Path(__file__)),
                      "velocity_unit_semantics": "raw ASE Atoms.get_velocities() internal units; no conversion",
                      "execution_status": "SYNTHETIC_OR_UNSIGNED_CANDIDATE_ONLY",
                      "rng_state_file": "rng_state.json"})
        state_path = staging / "state.json"
        state_path.write_bytes(canonical_json_bytes(state))
        readme = ("# Ag64 segment checkpoint\n\n"
                  "Immutable per-step bundle for global step 1000->2000.\n"
                  "The scientific source identity is preserved as provenance; this bundle's "
                  "execution code hash is recorded separately.\n"
                  "No EAM potential bytes are included.\n").encode()
        (staging / "README.md").write_bytes(readme)
        file_hashes = {p.name: sha256_file(p) for p in sorted(staging.iterdir()) if p.is_file()}
        manifest = {"schema": 1, "task_id": TASK_ID, "owner": OWNER,
                    "global_step": step, "members_sha256": file_hashes,
                    "potential_file_included": False}
        manifest_path = staging / "MANIFEST.json"
        manifest_path.write_bytes(canonical_json_bytes(manifest))
        sums_lines = [f"{file_hashes[name]}  {name}" for name in sorted(file_hashes)]
        sums_lines.append(f"{sha256_file(manifest_path)}  MANIFEST.json")
        sums_path = staging / "SHA256SUMS.txt"
        sums_path.write_text("\n".join(sums_lines) + "\n", encoding="utf-8")
        for path in staging.iterdir():
            _require(path.stat().st_size <= MAX_ARCHIVE_MEMBER_BYTES,
                     "OUTPUT_MEMBER_EXCEEDS_50_MIB_ARCHIVE_LIMIT")
            with path.open("rb") as handle:
                os.fsync(handle.fileno())
        _fsync_dir(staging)
        _require(not destination.exists(), "OUTPUT_BUNDLE_ALREADY_EXISTS_NO_OVERWRITE")
        _rename_noreplace(staging, destination)
        _fsync_dir(root)
        return destination
    except BaseException:
        # Failed staging is not a checkpoint; original source is never touched.
        if staging.exists():
            import shutil
            shutil.rmtree(staging)
        raise


def create_archive(bundle_dir: str | Path, archive_path: str | Path) -> str:
    bundle_dir, archive_path = Path(bundle_dir), Path(archive_path)
    _require(not archive_path.exists(), "ARCHIVE_ALREADY_EXISTS_NO_OVERWRITE")
    temp_path = archive_path.with_name("." + archive_path.name + ".tmp")
    _require(not temp_path.exists(), "ARCHIVE_TEMP_PATH_ALREADY_EXISTS")
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(temp_path, "x:gz") as archive:
            archive.add(bundle_dir, arcname=bundle_dir.name, recursive=True)
        _require(temp_path.stat().st_size <= MAX_ARCHIVE_MEMBER_BYTES,
                 "OUTPUT_ARCHIVE_EXCEEDS_50_MIB_LIMIT")
        with temp_path.open("rb") as handle:
            os.fsync(handle.fileno())
        # Hard-link publication is atomic and cannot overwrite an existing name.
        os.link(temp_path, archive_path)
        temp_path.unlink()
        _fsync_dir(archive_path.parent)
    except BaseException:
        if temp_path.exists():
            temp_path.unlink()
        raise
    return sha256_file(archive_path)


def advance_segment_with_outputs(dyn: BoundedLangevin, *, rng_state_all: dict,
                                 source_provenance: dict, output_root: str | Path,
                                 remote_roundtrip: Callable[[Path, str], None],
                                 existing_thermal_rows: list[dict] | None = None,
                                 existing_trajectory_frames: list[dict] | None = None,
                                 checkpoint_every: int = 100, max_wall_s: float = 600.0,
                                 clock: Callable[[], float] = time.monotonic) -> list[dict]:
    """Execute only through the integrator-owned stop and archive each boundary.

    The callback is mandatory: a failed upload/fresh-fetch raises before another
    ASE step. This helper was exercised only with a virtual calculator here.
    """
    _require(START_STEP <= dyn.nsteps < STOP_AT_STEP, "OUTPUT_SEGMENT_START_OUTSIDE_1000_2000")
    _require(checkpoint_every == 100, "CHECKPOINT_INTERVAL_MUST_REMAIN_100_STEPS")
    _require(remote_roundtrip is not None, "REMOTE_ARCHIVE_CALLBACK_REQUIRED")
    started = clock()
    if dyn.nsteps == START_STEP:
        _require(existing_thermal_rows is None and existing_trajectory_frames is None,
                 "FRESH_SEGMENT_MUST_NOT_IMPORT_OTHER_OUTPUT_HISTORY")
        thermal_rows: list[dict] = []
        trajectory_frames: list[dict] = [{"step": START_STEP, "positions": dyn.atoms.positions.copy()}]
    else:
        _require(existing_thermal_rows is not None and existing_trajectory_frames is not None,
                 "RESUME_REQUIRES_VERIFIED_OUTPUT_HISTORY")
        thermal_rows = copy.deepcopy(existing_thermal_rows)
        trajectory_frames = copy.deepcopy(existing_trajectory_frames)
        _require(trajectory_frames and trajectory_frames[-1]["step"] == dyn.nsteps,
                 "RESUME_TRAJECTORY_ENDPOINT_MISMATCH")
    outputs: list[dict] = []
    original_callback = dyn.on_step

    def on_step(current_dyn, phase, target, _forces):
        step = current_dyn.nsteps
        if original_callback is not None:
            original_callback(current_dyn, phase, target, _forces)
        # Match source sampler order: successful integration; thermal row every
        # 20; trajectory frame every 100; immutable checkpoint/archive every 100.
        if step % 20 == 0:
            atoms = current_dyn.atoms
            potential_e = float(atoms.get_potential_energy())
            kinetic_e = float(atoms.get_kinetic_energy())
            stress = atoms.get_stress(voigt=False)
            pressure = -float(np.trace(stress)) / 3.0 + 2.0 * kinetic_e / (3.0 * atoms.get_volume())
            row = {"step": step, "time_fs": step * DT_FS, "phase": phase,
                   "target_temperature_K": target, "T_K": float(atoms.get_temperature()),
                   "sampler_energy_eV": potential_e, "kinetic_energy_eV": kinetic_e,
                   "pressure_eV_A3": pressure}
            _require(all(np.isfinite(float(row[key])) for key in fields_for_finite(row)),
                     "SEGMENT_THERMAL_LOG_NONFINITE")
            thermal_rows.append(row)
        if step % 100 == 0:
            trajectory_frames.append({"step": step, "positions": current_dyn.atoms.positions.copy()})
        if step % checkpoint_every == 0 or step == STOP_AT_STEP:
            snapshot = snapshot_state(current_dyn.atoms, step, rng_state_all,
                                      current_dyn.rng, source_provenance)
            bundle = write_checkpoint_bundle(output_root, "segment", snapshot,
                                             thermal_rows, trajectory_frames)
            archive_path = Path(output_root) / "archives" / f"step-{step:09d}.tar.gz"
            archive_sha = publish_and_verify(bundle, archive_path, remote_roundtrip)
            outputs.append({"step": step, "bundle": str(bundle),
                            "archive": str(archive_path), "archive_sha256": archive_sha})

    dyn.on_step = on_step
    while dyn.nsteps < STOP_AT_STEP:
        if clock() - started >= max_wall_s:
            raise SegmentRefused("SEGMENT_WALLTIME_BUDGET_EXCEEDED")
        dyn.step()
    _require(dyn.nsteps == STOP_AT_STEP and outputs and outputs[-1]["step"] == STOP_AT_STEP,
             "FINAL_STEP2000_ARCHIVE_NOT_VERIFIED")
    return outputs


def fields_for_finite(row: dict) -> list[str]:
    return [key for key in row if key != "phase"]


def validate_release_lock(lock_path: str | Path, expected_lock_sha256: str,
                          expected_code_sha256: str) -> dict:
    """Real execution is possible only with a separately signed future lock."""
    lock_path = Path(lock_path)
    _require(lock_path.is_file(), "SIGNED_RELEASE_LOCK_MISSING")
    actual_lock_sha = sha256_file(lock_path)
    _require(actual_lock_sha == expected_lock_sha256, "SIGNED_RELEASE_LOCK_SHA_MISMATCH")
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    _require(expected_code_sha256 == sha256_file(Path(__file__)), "EXECUTION_CODE_SHA_MISMATCH")
    _require(lock.get("execution_code_sha256") == expected_code_sha256,
             "RELEASE_LOCK_DOES_NOT_BIND_THIS_EXECUTOR")
    _require(lock.get("release_status") == "C_SIGNED_FOR_SEGMENT_1000_2000"
             and lock.get("execution_authorized") is True,
             "RELEASE_LOCK_UNSIGNED_OR_NOT_AUTHORIZED")
    _require(lock.get("task_id") == TASK_ID and lock.get("owner") == OWNER,
             "RELEASE_LOCK_TASK_OR_OWNER_MISMATCH")
    _require(lock.get("start_step") == START_STEP and lock.get("stop_at_step") == STOP_AT_STEP
             and lock.get("max_added_steps") == MAX_ADDED_STEPS,
             "RELEASE_LOCK_BOUNDARY_MISMATCH")
    _require(lock.get("source_scientific_parameters_sha256") == EXPECTED_SCIENTIFIC_PARAMETERS_SHA256,
             "RELEASE_LOCK_SOURCE_SCIENTIFIC_IDENTITY_MISMATCH")
    required_pins = {
        "source_archive_sha256": EXPECTED_ARCHIVE_SHA256,
        "source_archive_manifest_sha256": EXPECTED_ARCHIVE_MANIFEST_SHA256,
        "source_checkpoint_state_sha256": EXPECTED_STATE_SHA256,
        "source_state_arrays_sha256": EXPECTED_ARRAYS_SHA256,
        "source_manifest_sha256": EXPECTED_SOURCE_MANIFEST_SHA256,
        "source_executor_code_sha256": EXPECTED_OLD_RUNNER_SHA256,
        "potential_sha256": EXPECTED_POTENTIAL_SHA256,
    }
    for key, expected in required_pins.items():
        _require(lock.get(key) == expected, "RELEASE_LOCK_PIN_MISMATCH:" + key)
    _require(lock.get("sampling") == SIGNED_SAMPLING_PARAMETERS,
             "RELEASE_LOCK_SCIENTIFIC_OR_RESOURCE_PARAMETERS_MISMATCH")
    _require(lock.get("remote_archive_verified") is True
             and lock.get("remote_archive_method") == "GITHUB_PUSH_FRESH_FETCH_SHA256",
             "REMOTE_GITHUB_FRESH_FETCH_RELEASE_GATE_NOT_PASSED")
    return lock


def verify_source_archive(archive_path: str | Path, archive_manifest_path: str | Path) -> dict:
    """Read-only verification of the immutable D step-1000 archive bytes."""
    archive_path = Path(archive_path)
    archive_manifest_path = Path(archive_manifest_path)
    _require(archive_path.is_file(), "SOURCE_ARCHIVE_NOT_FOUND")
    _require(archive_manifest_path.is_file(), "SOURCE_ARCHIVE_MANIFEST_NOT_FOUND")
    archive_sha = sha256_file(archive_path)
    _require(archive_sha == EXPECTED_ARCHIVE_SHA256, "SOURCE_ARCHIVE_SHA_MISMATCH")
    outer_manifest_sha = sha256_file(archive_manifest_path)
    _require(outer_manifest_sha == EXPECTED_ARCHIVE_MANIFEST_SHA256,
             "SOURCE_ARCHIVE_OUTER_MANIFEST_SHA_MISMATCH")
    outer_manifest = json.loads(archive_manifest_path.read_text(encoding="utf-8"))
    _require(outer_manifest.get("archive_sha256") == archive_sha
             and outer_manifest.get("task_id") == TASK_ID
             and outer_manifest.get("owner") == OWNER
             and outer_manifest.get("step") == START_STEP,
             "SOURCE_ARCHIVE_OUTER_MANIFEST_IDENTITY_MISMATCH")
    _require(outer_manifest.get("potential_sha256") == EXPECTED_POTENTIAL_SHA256
             and outer_manifest.get("source_code_sha256", {}).get("run_sampling_codex.py")
             == EXPECTED_OLD_RUNNER_SHA256,
             "SOURCE_ARCHIVE_OLD_EXECUTOR_OR_POTENTIAL_IDENTITY_MISMATCH")
    with tarfile.open(archive_path, "r:gz") as archive:
        members = archive.getmembers()
        by_name = {member.name: member for member in members}
        _require(len(by_name) == len(members), "SOURCE_ARCHIVE_DUPLICATE_PATHS")
        for name in by_name:
            p = Path(name)
            _require(not p.is_absolute() and ".." not in p.parts, "SOURCE_ARCHIVE_UNSAFE_PATH")
            _require(by_name[name].isfile(), "SOURCE_ARCHIVE_NONREGULAR_MEMBER")
        listed_members = {item["path"]: item for item in outer_manifest.get("members", [])}
        _require(set(listed_members) == set(by_name), "SOURCE_ARCHIVE_OUTER_MEMBER_SET_MISMATCH")
        for name, member in by_name.items():
            content = archive.extractfile(member).read()
            item = listed_members[name]
            _require(len(content) == item.get("size_bytes")
                     and sha256_bytes(content) == item.get("sha256"),
                     "SOURCE_ARCHIVE_MEMBER_HASH_MISMATCH:" + name)
        arrays_names = [n for n in by_name if n.endswith(
            "restart/step-000001000-db6a96b6a11f49d4a97a7fa00db2b2a6/state_arrays.npz")]
        state_names = [n for n in by_name if n == "sampler_state.json"]
        source_names = [n for n in by_name if n == "source_manifest.json"]
        _require(len(arrays_names) == len(state_names) == len(source_names) == 1,
                 "SOURCE_ARCHIVE_REQUIRED_MEMBER_SET_MISMATCH")
        arrays_data = archive.extractfile(by_name[arrays_names[0]]).read()
        state_data = archive.extractfile(by_name[state_names[0]]).read()
        source_data = archive.extractfile(by_name[source_names[0]]).read()
    _require(sha256_bytes(arrays_data) == EXPECTED_ARRAYS_SHA256, "SOURCE_ARCHIVE_ARRAYS_SHA_MISMATCH")
    _require(sha256_bytes(state_data) == EXPECTED_STATE_SHA256, "SOURCE_ARCHIVE_STATE_SHA_MISMATCH")
    _require(sha256_bytes(source_data) == EXPECTED_SOURCE_MANIFEST_SHA256,
             "SOURCE_ARCHIVE_MANIFEST_SHA_MISMATCH")
    state = json.loads(state_data)
    _require(state.get("step") == START_STEP and state.get("task_id") == TASK_ID
             and state.get("owner") == OWNER, "SOURCE_ARCHIVE_STATE_IDENTITY_MISMATCH")
    return {"archive_sha256": archive_sha, "outer_manifest_sha256": outer_manifest_sha,
            "state_sha256": sha256_bytes(state_data),
            "arrays_sha256": sha256_bytes(arrays_data),
            "source_manifest_sha256": sha256_bytes(source_data),
            "step": START_STEP, "task_id": TASK_ID, "owner": OWNER}


def restore_source_from_archive(archive_path: str | Path) -> RestoredCheckpoint:
    """Verify and reconstruct source state in a disposable read-only temp dir."""
    with tarfile.open(archive_path, "r:gz") as archive:
        arrays_member = next(m for m in archive.getmembers() if m.name.endswith(
            "restart/step-000001000-db6a96b6a11f49d4a97a7fa00db2b2a6/state_arrays.npz"))
        state_member = next(m for m in archive.getmembers() if m.name == "sampler_state.json")
        arrays_bytes = archive.extractfile(arrays_member).read()
        state_bytes = archive.extractfile(state_member).read()
    with tempfile.TemporaryDirectory(prefix="codex-c-step1000-restore-") as tmp:
        root = Path(tmp)
        arrays_path = root / "state_arrays.npz"
        state_path = root / "sampler_state.json"
        arrays_path.write_bytes(arrays_bytes)
        state_path.write_bytes(state_bytes)
        return restore_checkpoint(arrays_path, state_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-lock", type=Path, required=True)
    parser.add_argument("--release-lock-sha256", required=True)
    parser.add_argument("--execution-code-sha256", required=True)
    parser.add_argument("--owner", required=True)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--start-step", type=int, required=True)
    parser.add_argument("--stop-at-step", type=int, required=True)
    parser.add_argument("--source-archive", type=Path, required=True)
    parser.add_argument("--source-archive-manifest", type=Path, required=True)
    parser.add_argument("--source-archive-sha256", required=True)
    parser.add_argument("--potential", type=Path, required=True)
    parser.add_argument("--potential-sha256", required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    try:
        _require(args.owner == OWNER and args.task_id == TASK_ID, "CLI_TASK_OR_OWNER_MISMATCH")
        _require(args.start_step == START_STEP and args.stop_at_step == STOP_AT_STEP,
                 "CLI_STEP_BOUNDARY_MISMATCH")
        _require(args.source_archive_sha256 == EXPECTED_ARCHIVE_SHA256,
                 "CLI_SOURCE_ARCHIVE_SHA_EXPECTATION_MISMATCH")
        source_evidence = verify_source_archive(args.source_archive, args.source_archive_manifest)
        _require(source_evidence["archive_sha256"] == args.source_archive_sha256,
                 "CLI_SOURCE_ARCHIVE_SHA_MISMATCH")
        _require(args.potential_sha256 == EXPECTED_POTENTIAL_SHA256,
                 "CLI_POTENTIAL_SHA_EXPECTATION_MISMATCH")
        lock = validate_release_lock(args.release_lock, args.release_lock_sha256,
                                     args.execution_code_sha256)
        restored = restore_source_from_archive(args.source_archive)
        _require(restored.global_step == START_STEP and restored.atoms.calc is None,
                 "SOURCE_RESTORE_PREFLIGHT_FAILED")
        if not args.execute:
            print(json.dumps({"status": "PREFLIGHT_PASS_NO_SCIENCE",
                              "lock_sha256": sha256_file(args.release_lock),
                              "restored_step": restored.global_step,
                              "atom_count": len(restored.atoms)}))
            return 0
        _require(args.potential.is_file(), "POTENTIAL_FILE_NOT_PRESENT")
        _require(sha256_file(args.potential) == EXPECTED_POTENTIAL_SHA256,
                 "POTENTIAL_FILE_SHA_MISMATCH")
        # The committed unsigned candidate does not reach this point. Remote
        # GitHub archive orchestration is intentionally a separate unresolved
        # production gate; no local-only execution can satisfy storage release.
        raise SegmentRefused("VERIFIED_GITHUB_ARCHIVE_ORCHESTRATION_NOT_IMPLEMENTED")
    except (SegmentRefused, OSError, ValueError, KeyError) as exc:
        print(json.dumps({"status": "HOLD", "reason": str(exc)}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
