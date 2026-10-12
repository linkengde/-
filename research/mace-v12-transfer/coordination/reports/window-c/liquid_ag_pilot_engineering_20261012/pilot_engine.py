#!/usr/bin/env python3
"""Engineering-validated, C-lock-gated single-task ASE runner for D rho100.

Engineering fixtures are written only to TemporaryDirectory and are explicitly
not candidate structures. The formal entry point refuses the presently frozen
execution-HOLD lock before creating an output directory or generating atoms.
"""
from __future__ import annotations

import copy
import ast
import csv
import hashlib
import json
import math
import os
import resource
import signal
import sys
import tempfile
import time
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import ase
from ase import Atoms
from ase.calculators.eam import EAM
from ase.constraints import FixCom
from ase.io import write
from ase.io.trajectory import Trajectory
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.neighborlist import neighbor_list
from ase.units import fs, kB


TASK_ID = "C-D-AGL-1250K-rho100-pilot01"
LOCK_SHA256 = "e1adaf12d542184c709d2eb2a5363cbfe9578c168ebc5e7ad260ee39a4faf09b"
LOCK_COMMIT = "287b60af424f7bcc3a7b3598e524f5b543884ec9"
LOCK_REPOSITORY_PATH = (
    "research/mace-v12-transfer/coordination/reports/window-c/"
    "liquid_ag_rho100_first_pilot_release_20261012/LIQUID_AG_PILOT_LOCK.json"
)
EXPECTED_POTENTIAL_SHA256 = (
    "958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447"
)
EXPECTED_POTENTIAL_BYTES = 780452
EXPECTED_ASE_VERSION = "3.29.0"
EXPECTED_PYTHON_VERSION = "3.12.14"
EXPECTED_C_GENERATOR_SHA256 = (
    "b78b734798c9c075826e6611e3d2f7f156f0b05ff7c87159812ac76ed0e918a8"
)
EXPECTED_C_INSPECT_EAM_SHA256 = (
    "d3a437582dcd5bad7e361c6c81f64ac1191161170105b43133ffdd4065526a77"
)
AMU_A3_TO_G_CM3 = 1.66053906660
THERMAL_FIELDS = (
    "step", "time_fs", "target_temperature_K", "temperature_K",
    "potential_energy_eV", "kinetic_energy_eV", "total_energy_eV",
)
RNG_NAMES = ("packing", "velocity_initialization", "thermostat")


class PreflightError(RuntimeError):
    """A verified engineering or release gate did not pass."""


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_name(path.name + f".tmp-{os.getpid()}")
    with temp_path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temp_path, path)
    _fsync_dir(path.parent)


def _fsync_dir(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _json_clone(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True, allow_nan=False))


def verify_lock(lock_path: str | Path, expected_sha256: str = LOCK_SHA256) -> tuple[dict, str]:
    path = Path(lock_path)
    if not path.is_file():
        raise PreflightError(f"frozen lock not found: {path}")
    actual = sha256_file(path)
    if actual != expected_sha256:
        raise PreflightError(f"lock SHA mismatch: expected {expected_sha256}; got {actual}")
    try:
        lock = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise PreflightError(f"lock JSON parse failed: {type(exc).__name__}: {exc}") from exc
    validate_frozen_parameters(lock)
    return lock, actual


def validate_frozen_parameters(lock: dict) -> None:
    if lock.get("schema_version") != 1:
        raise PreflightError("unsupported lock schema")
    if lock.get("lock_status") != "C_SCIENTIFIC_PARAMETERS_FROZEN_EXECUTION_HOLD_FOR_D_PREFLIGHT":
        raise PreflightError("lock status does not match the C frozen execution-HOLD release")
    if lock.get("execution_authorized") is not False:
        raise PreflightError("engineering input is expected to remain execution_authorized=false")
    scope = lock.get("scope", {})
    if scope.get("owner") != "D" or scope.get("task_id") != TASK_ID:
        raise PreflightError("lock does not identify the sole D rho100 pilot")
    if (
        scope.get("atom_count") != 64
        or scope.get("composition") != {"Ag": 64}
        or scope.get("initial_step") != 0
        or scope.get("maximum_added_steps") != 1000
        or scope.get("maximum_end_step") != 1000
        or scope.get("maximum_wall_seconds") != 3600
        or scope.get("automatic_continuation") is not False
    ):
        raise PreflightError("task scope differs from the frozen 64-atom 1000-step/3600s boundary")
    density = lock["density"]
    cell = np.asarray(density["cell"]["matrix_A"], dtype=float)
    if cell.shape != (3, 3) or not np.all(np.isfinite(cell)):
        raise PreflightError("C cell is not a finite 3x3 matrix")
    edge = float(cell[0, 0])
    if not np.allclose(cell, np.eye(3) * edge, rtol=0, atol=1e-12):
        raise PreflightError("the frozen cell is not the stated diagonal cube")
    if density["relative_density"] != 1.0 or edge <= 0:
        raise PreflightError("this D pilot must use the frozen rho100 cell")
    structure = lock["initial_structure_generator"]
    if (
        structure.get("r_min_A") != 2.2
        or structure.get("maximum_trial_points") != 200000
        or structure.get("structure_generated") is not False
    ):
        raise PreflightError("frozen generator limits or no-structure state do not match")
    random = lock["scope"]
    if random.get("master_seed") != 125010003:
        raise PreflightError("master seed differs from C's fixed rho100 seed")
    rng = lock["initial_structure_generator"]["rng"]
    if rng != (
        "NumPy Generator(PCG64) seeded by SeedSequence(master_seed).spawn(3); "
        "child streams in fixed order: packing, velocity initialization, Langevin thermostat"
    ):
        raise PreflightError("C's exact three-child PCG64 stream derivation changed")
    dynamics = lock["dynamics"]
    if (
        dynamics.get("engine") != "ASE 3.29.0 EAM calculator using the exact setfl file above"
        or dynamics.get("integrator_and_ensemble")
        != "ASE Langevin NVT, fixed cubic cell, FixCom constraint, Langevin fixcm=False"
        or dynamics.get("timestep_fs") != 0.5
        or dynamics.get("langevin_friction_fs_inverse") != 0.01
        or dynamics.get("temperature_dof") != "3N-3 = 189 after FixCom; use ASE constraint-aware temperature; no rotational constraints"
    ):
        raise PreflightError("frozen ASE/Langevin/FixCom parameters changed")
    if dynamics.get("segment_temperature_schedule") != (
        "steps 0 through 999 target 1600 K (the existing hot-start stage); "
        "no cooling or 1250 K equilibration is included in this pilot"
    ):
        raise PreflightError("frozen first-segment temperature schedule changed")
    if structure.get("velocity_initialization") != (
        "ASE 3.29.0 MaxwellBoltzmannDistribution with force_temp=True at 1600 K "
        "using the velocity child stream and FixCom constraint; no velocity redraw or reseed"
    ):
        raise PreflightError("C's velocity initialization rule changed")
    if lock["restart_and_output"].get("automatic_continue_after_step_1000") is not False:
        raise PreflightError("automatic continuation must remain disabled")
    if lock["runtime_versions"].get("python_expected") != EXPECTED_PYTHON_VERSION or lock["runtime_versions"].get("ase_expected") != EXPECTED_ASE_VERSION:
        raise PreflightError("runtime versions differ from the exact C lock")
    source_hashes = lock["source_and_code_hashes"]
    if source_hashes.get("C_reference_run_sampling_py_sha256") != EXPECTED_C_GENERATOR_SHA256:
        raise PreflightError("C reference generator source hash differs from the fixed release")
    if source_hashes.get("C_reference_inspect_eam_py_sha256") != EXPECTED_C_INSPECT_EAM_SHA256:
        raise PreflightError("C reference EAM parser source hash differs from the fixed release")
    stops = lock["stops"]
    if (
        stops.get("wall_seconds") != 3600
        or stops.get("maximum_step") != 1000
        or stops.get("minimum_free_bytes") != 2147483648
        or stops.get("trajectory_bundle_cap_bytes") != 536870912
        or stops.get("checkpoint_reserve_bytes") != 33554432
    ):
        raise PreflightError("frozen safety/resource limits changed")
    potential = lock["potential"]
    if (
        potential.get("sha256") != EXPECTED_POTENTIAL_SHA256
        or potential.get("bytes") != EXPECTED_POTENTIAL_BYTES
        or not str(potential.get("radial_tail_policy", "")).startswith(
            "C_APPROVED_PUBLISHED_SETFL_ONE_BIN_ZERO_TAIL;"
        )
        or potential.get("parameter_file_uploaded_to_GitHub") is not False
    ):
        raise PreflightError("potential identity or handling policy changed")


def make_rngs(master_seed: int) -> tuple[dict[str, np.random.Generator], dict]:
    children = np.random.SeedSequence(master_seed).spawn(3)
    names = RNG_NAMES
    streams = {
        name: np.random.Generator(np.random.PCG64(children[index]))
        for index, name in enumerate(names)
    }
    provenance = {
        "algorithm": "numpy.SeedSequence(master_seed).spawn(3)",
        "bit_generator": "PCG64",
        "master_seed": int(master_seed),
        "child_order": list(names),
        "child_spawn_keys": {name: list(children[index].spawn_key) for index, name in enumerate(names)},
    }
    return streams, provenance


def capture_rng_states(streams: dict[str, np.random.Generator]) -> dict:
    return {name: copy.deepcopy(streams[name].bit_generator.state) for name in RNG_NAMES}


def restore_rngs(states: dict) -> dict[str, np.random.Generator]:
    if set(states) != set(RNG_NAMES):
        raise PreflightError("checkpoint must contain exactly the three locked PCG64 child streams")
    result = {}
    for name in RNG_NAMES:
        bitgen = np.random.PCG64()
        bitgen.state = copy.deepcopy(states[name])
        result[name] = np.random.Generator(bitgen)
    return result


def pack_random_positions(edge: float, count: int, r_min: float, max_trials: int, rng) -> tuple[np.ndarray, int]:
    """C's uniform sequential periodic rejection packing, preserving draw order."""
    placed: list[np.ndarray] = []
    attempts = 0
    while len(placed) < count and attempts < max_trials:
        trial = rng.uniform(0.0, edge, size=3)
        attempts += 1
        if placed:
            delta = np.asarray(placed) - trial
            delta -= np.rint(delta / edge) * edge
            if np.any((delta**2).sum(axis=1) < r_min**2):
                continue
        placed.append(trial)
    if len(placed) != count:
        raise PreflightError(
            f"locked random packing placed {len(placed)}/{count} points in {attempts} trial points"
        )
    return np.asarray(placed, dtype=float), attempts


def _verify_c_reference_generator(lock: dict) -> dict:
    source_path = (
        "research/mace-v12-transfer/coordination/reports/window-c/"
        "liquid_ag_six_pilot_authorized_dispatch_20261011/run_sampling.py"
    )
    try:
        completed = subprocess.run(
            ["git", "-C", "/workspace/-", "show", f"{LOCK_COMMIT}:{source_path}"],
            check=True, capture_output=True, text=True, timeout=30,
        )
    except Exception as exc:
        raise PreflightError(f"cannot read C's fixed generator reference: {exc}") from exc
    raw = completed.stdout.encode()
    source_sha = hashlib.sha256(raw).hexdigest()
    if source_sha != EXPECTED_C_GENERATOR_SHA256:
        raise PreflightError("C reference generator file SHA does not match the frozen lock")
    tree = ast.parse(completed.stdout)
    function = next(
        (node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "pack_random"),
        None,
    )
    if function is None:
        raise PreflightError("C reference has no pack_random function")
    namespace = {"np": np, "Hold": PreflightError}
    exec(compile(ast.Module(body=[function], type_ignores=[]), "<C-pack_random>", "exec"), namespace)
    edge = float(lock["density"]["cell"]["matrix_A"][0][0])
    reference = namespace["pack_random"](
        int(lock["scope"]["master_seed"]), edge,
        int(lock["initial_structure_generator"]["maximum_trial_points"]), np,
    )
    reference_positions, reference_packing, reference_velocity, reference_thermostat, reference_attempts = reference
    ours, _ = make_rngs(int(lock["scope"]["master_seed"]))
    own_positions, own_attempts = pack_random_positions(
        edge, int(lock["scope"]["atom_count"]),
        float(lock["initial_structure_generator"]["r_min_A"]),
        int(lock["initial_structure_generator"]["maximum_trial_points"]),
        ours["packing"],
    )
    checks = {
        "positions_exact": bool(np.array_equal(reference_positions, own_positions)),
        "attempt_count_exact": reference_attempts == own_attempts,
        "packing_rng_state_exact": reference_packing.bit_generator.state == ours["packing"].bit_generator.state,
        "velocity_rng_state_exact": reference_velocity.bit_generator.state == ours["velocity_initialization"].bit_generator.state,
        "thermostat_rng_state_exact": reference_thermostat.bit_generator.state == ours["thermostat"].bit_generator.state,
    }
    if not all(checks.values()):
        raise PreflightError("D generator differs from C's exact fixed-seed reference implementation")
    return {
        "reference_source_path": source_path,
        "reference_source_sha256": source_sha,
        "reference_function": "pack_random",
        **checks,
        "all_pass": all(checks.values()),
        "trial_points": int(reference_attempts),
        "temporary_positions_sha256": hashlib.sha256(
            reference_positions.astype("<f8").tobytes()
        ).hexdigest(),
        "fixture_only": True,
    }


def generate_atoms(lock: dict, fixture_only: bool) -> tuple[Atoms, dict, dict[str, np.random.Generator], dict]:
    """Generate exactly the C-locked Ag64 coordinates; caller gates formal use."""
    edge = float(lock["density"]["cell"]["matrix_A"][0][0])
    count = int(lock["scope"]["atom_count"])
    rngs, provenance = make_rngs(int(lock["scope"]["master_seed"]))
    states_before = capture_rng_states(rngs)
    positions, attempts = pack_random_positions(
        edge, count, float(lock["initial_structure_generator"]["r_min_A"]),
        int(lock["initial_structure_generator"]["maximum_trial_points"]), rngs["packing"],
    )
    atoms = Atoms("Ag64", positions=positions, cell=np.eye(3) * edge, pbc=True)
    atoms.new_array("atom_id", np.arange(1, 65, dtype=np.int32))
    atoms.set_constraint(FixCom())
    MaxwellBoltzmannDistribution(
        atoms,
        temperature_K=_velocity_temperature_from_lock(lock),
        force_temp=True,
        rng=rngs["velocity_initialization"],
    )
    atoms.info["engineering_fixture_only"] = bool(fixture_only)
    min_distance = minimum_periodic_pair_distance(atoms.get_positions(), edge)
    measured_density = float(atoms.get_masses().sum() * AMU_A3_TO_G_CM3 / atoms.get_volume())
    expected_density = float(lock["density"]["cell"]["density_recomputed_from_cell_g_cm3"])
    if min_distance < float(lock["initial_structure_generator"]["r_min_A"]):
        raise PreflightError(f"temporary generator fixture violates r_min: {min_distance}")
    # The cell entry is rounded to the published decimal precision. ASE's amu
    # constant gives a 3.3e-9 g/cm^3 arithmetic difference from the source's
    # exact-Avogadro calculation; this tolerance checks only that rounding.
    if not math.isclose(measured_density, expected_density, rel_tol=0, abs_tol=1e-8):
        raise PreflightError(f"locked cell density mismatch: {measured_density} != {expected_density}")
    if not np.array_equal(atoms.arrays["atom_id"], np.arange(1, 65, dtype=np.int32)):
        raise PreflightError("temporary fixture atom IDs are not unique and canonical")
    dof = atoms.get_number_of_degrees_of_freedom()
    temperature = float(atoms.get_temperature())
    if dof != 189 or not math.isclose(temperature, _temperature_from_lock(lock), rel_tol=1e-12):
        raise PreflightError(f"FixCom temperature check failed: dof={dof}, T={temperature}")
    provenance_record = {
        **provenance,
        "states_before_initialization": states_before,
        "states_after_initialization": capture_rng_states(rngs),
    }
    metadata = {
        "fixture_only": bool(fixture_only),
        "atom_count": len(atoms),
        "composition": {"Ag": 64},
        "cell_A": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(),
        "density_g_cm3": measured_density,
        "reference_density_g_cm3": expected_density,
        "density_rounding_delta_g_cm3": measured_density - expected_density,
        "minimum_periodic_pair_distance_A": min_distance,
        "r_min_A": 2.2,
        "packing_trial_points": attempts,
        "atom_id_unique_count": int(len(np.unique(atoms.arrays["atom_id"]))),
        "velocity_temperature_K": temperature,
        "temperature_dof": int(dof),
        "rng_provenance": provenance_record,
        "test_geometry_sha256": hashlib.sha256(atoms.get_positions().astype("<f8").tobytes()).hexdigest(),
        "rng_states_before": states_before,
        "rng_states_after_initialization": capture_rng_states(rngs),
    }
    return atoms, metadata, rngs, provenance_record


def make_test_atoms(lock: dict) -> tuple[Atoms, dict, dict[str, np.random.Generator], dict]:
    """Create an isolated test-only fixture; never writes formal initial data."""
    return generate_atoms(lock, fixture_only=True)


def minimum_periodic_pair_distance(positions: np.ndarray, edge: float) -> float:
    minimum = math.inf
    for i in range(len(positions)):
        delta = positions[i + 1:] - positions[i]
        if len(delta):
            delta -= np.rint(delta / edge) * edge
            minimum = min(minimum, float(np.sqrt(np.sum(delta * delta, axis=1)).min()))
    return minimum


def _temperature_from_lock(lock: dict) -> float:
    schedule = lock["dynamics"]["segment_temperature_schedule"]
    try:
        return float(schedule.split("target ", 1)[1].split(" K", 1)[0])
    except (IndexError, ValueError) as exc:
        raise PreflightError("cannot parse the frozen hot-stage target temperature") from exc


def _velocity_temperature_from_lock(lock: dict) -> float:
    rule = lock["initial_structure_generator"]["velocity_initialization"]
    if "force_temp=True at " not in rule:
        raise PreflightError("cannot parse the frozen velocity initialization temperature")
    try:
        return float(rule.split("force_temp=True at ", 1)[1].split(" K", 1)[0])
    except (IndexError, ValueError) as exc:
        raise PreflightError("cannot parse the frozen velocity initialization temperature") from exc


def _parameter_sha(lock: dict) -> str:
    fields = {
        "scope": lock["scope"],
        "density": lock["density"],
        "potential": lock["potential"],
        "initial_structure_generator": lock["initial_structure_generator"],
        "dynamics": lock["dynamics"],
        "restart_and_output": lock["restart_and_output"],
        "stops": lock["stops"],
    }
    return canonical_sha(fields)


def assert_segment_limits(lock: dict, step: int, elapsed_wall_s: float) -> None:
    if step < 0 or step > int(lock["stops"]["maximum_step"]):
        raise PreflightError("global step exceeds the frozen 1000-step ceiling")
    if not math.isfinite(float(elapsed_wall_s)) or elapsed_wall_s < 0:
        raise PreflightError("elapsed walltime must be finite and nonnegative")


def walltime_reached(elapsed_wall_s: float, maximum_wall_s: float) -> bool:
    """Return the stop decision without treating an expired clock as an error.

    A completed integrator step is indivisible. The caller checks this before
    starting a step and again immediately after it completes, then persists
    that last complete state without starting another step.
    """
    elapsed = float(elapsed_wall_s)
    maximum = float(maximum_wall_s)
    if not math.isfinite(elapsed) or elapsed < 0:
        raise PreflightError("elapsed walltime must be finite and nonnegative")
    if not math.isfinite(maximum) or maximum <= 0:
        raise PreflightError("maximum walltime must be finite and positive")
    return elapsed >= maximum


def _save_exit_boundary(save_boundary, status: str, failure: str | None, exc=None):
    """Classify an exceptional exit and persist its last complete boundary."""
    if exc is not None:
        status = "STOPPED_AFTER_EXCEPTION_AT_LAST_COMPLETED_STEP"
        detail = f"{type(exc).__name__}: {exc}"
        failure = f"{failure}; {detail}" if failure else detail
    checkpoint_sha, persisted_status = save_boundary(status)
    return persisted_status, failure, checkpoint_sha


def _state_payload(
    atoms: Atoms,
    rngs: dict[str, np.random.Generator],
    lock: dict,
    lock_sha: str,
    potential_sha: str,
    runner_sha: str,
    environment_sha: str,
    step: int,
    trajectory: Path,
    thermal: Path,
    energy: Path,
) -> dict:
    return {
        "schema_version": 1,
        "fixture_only": True,
        "task_id": TASK_ID,
        "global_step": int(step),
        "time_fs": float(step * lock["dynamics"]["timestep_fs"]),
        "symbols": atoms.get_chemical_symbols(),
        "atom_ids": atoms.arrays["atom_id"].astype(int).tolist(),
        "positions_A": atoms.get_positions().tolist(),
        "momenta_amu_A_per_fs": (atoms.get_momenta() * fs).tolist(),
        "velocities_A_per_fs": atoms.get_velocities().tolist(),
        "masses_amu": atoms.get_masses().tolist(),
        "cell_A": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(),
        "rng_states": capture_rng_states(rngs),
        "target_temperature_K": _temperature_from_lock(lock),
        "temperature_schedule_position": {
            "step": int(step), "target_temperature_K": _temperature_from_lock(lock)
        },
        "lock_sha256": lock_sha,
        "parameter_sha256": _parameter_sha(lock),
        "potential_sha256": potential_sha,
        "ase_version": ase.__version__,
        "python_version": ".".join(map(str, sys.version_info[:3])),
        "runner_sha256": runner_sha,
        "environment_sha256": environment_sha,
        "trajectory": {"path": trajectory.name, "size_bytes": trajectory.stat().st_size, "sha256": sha256_file(trajectory)},
        "thermal_log": {"path": thermal.name, "size_bytes": thermal.stat().st_size, "sha256": sha256_file(thermal)},
        "energy_log": {"path": energy.name, "size_bytes": energy.stat().st_size, "sha256": sha256_file(energy)},
        "elapsed_wall_s": 0.0,
    }


def write_checkpoint(directory: Path, state: dict) -> str:
    checkpoint = directory / "checkpoint.json"
    _atomic_json(checkpoint, state)
    digest = sha256_file(checkpoint)
    temp_sidecar = directory / "checkpoint.sha256.tmp"
    temp_sidecar.write_text(f"{digest}  checkpoint.json\n", encoding="ascii")
    with temp_sidecar.open("rb") as stream:
        os.fsync(stream.fileno())
    os.replace(temp_sidecar, directory / "checkpoint.sha256")
    _fsync_dir(directory)
    return digest


def restore_checkpoint(directory: Path, lock: dict, lock_sha: str, potential_sha: str,
                       runner_sha: str, environment_sha: str) -> tuple[Atoms, dict, dict[str, np.random.Generator]]:
    checkpoint = directory / "checkpoint.json"
    sidecar = (directory / "checkpoint.sha256").read_text(encoding="ascii").split()[0]
    actual = sha256_file(checkpoint)
    if actual != sidecar:
        raise PreflightError("checkpoint file SHA differs from its sidecar")
    state = json.loads(checkpoint.read_text(encoding="utf-8"))
    for key, expected in (
        ("lock_sha256", lock_sha), ("potential_sha256", potential_sha),
        ("runner_sha256", runner_sha), ("environment_sha256", environment_sha),
        ("parameter_sha256", _parameter_sha(lock)),
    ):
        if state.get(key) != expected:
            raise PreflightError(f"checkpoint identity mismatch: {key}")
    if state.get("fixture_only") is not True or state.get("task_id") != TASK_ID:
        raise PreflightError("checkpoint is not the isolated engineering fixture")
    if not 0 <= int(state["global_step"]) <= 1000:
        raise PreflightError("checkpoint step is outside the locked interval")
    atoms = Atoms(
        symbols=state["symbols"], positions=state["positions_A"],
        cell=state["cell_A"], pbc=state["pbc"],
    )
    if len(atoms) != 64 or atoms.get_chemical_symbols() != ["Ag"] * 64:
        raise PreflightError("checkpoint composition mismatch")
    ids = np.asarray(state["atom_ids"], dtype=np.int32)
    if not np.array_equal(ids, np.arange(1, 65, dtype=np.int32)):
        raise PreflightError("checkpoint atom ID order mismatch")
    atoms.new_array("atom_id", ids)
    atoms.set_masses(state["masses_amu"])
    atoms.set_constraint(FixCom())
    # These momenta were already projected through FixCom before checkpoint.
    # Avoid applying that projection a second time and changing the restart bits.
    atoms.set_momenta(np.asarray(state["momenta_amu_A_per_fs"], dtype=float) / fs, apply_constraint=False)
    if not np.array_equal(atoms.pbc, np.asarray([True, True, True])):
        raise PreflightError("checkpoint PBC mismatch")
    edge = float(lock["density"]["cell"]["matrix_A"][0][0])
    if not np.array_equal(atoms.cell.array, np.eye(3) * edge):
        raise PreflightError("checkpoint cell mismatch")
    for name in ("trajectory", "thermal_log", "energy_log"):
        item = state[name]
        path = directory / item["path"]
        if not path.is_file() or path.stat().st_size != item["size_bytes"] or sha256_file(path) != item["sha256"]:
            raise PreflightError(f"checkpoint-associated {name} size/SHA mismatch")
    with Trajectory(directory / "trajectory.traj", "r") as traj:
        if len(traj) != int(state["global_step"]) + 1:
            raise PreflightError("trajectory frame count does not match global step")
        if int(traj[-1].info.get("step", -1)) != int(state["global_step"]):
            raise PreflightError("trajectory boundary does not match checkpoint")
    return atoms, state, restore_rngs(state["rng_states"])


def _log_row(atoms: Atoms, step: int, lock: dict) -> tuple[list[Any], list[Any]]:
    epot = float(atoms.get_potential_energy())
    ekin = float(atoms.get_kinetic_energy())
    temp = float(atoms.get_temperature())
    values = [epot, ekin, epot + ekin, temp]
    if not np.all(np.isfinite(values)):
        raise PreflightError(f"non-finite thermodynamic value at fixture step {step}")
    tfs = step * float(lock["dynamics"]["timestep_fs"])
    target = _temperature_from_lock(lock)
    return (
        [step, f"{tfs:.12g}", f"{target:.12g}", f"{temp:.16g}", f"{epot:.16g}", f"{ekin:.16g}", f"{epot + ekin:.16g}"],
        [step, f"{tfs:.12g}", f"{epot:.16g}", f"{ekin:.16g}", f"{epot + ekin:.16g}"],
    )


def _append_row(path: Path, row: list[Any]) -> None:
    with path.open("a", newline="", encoding="utf-8") as stream:
        csv.writer(stream, lineterminator="\n").writerow(row)
        stream.flush()
        os.fsync(stream.fileno())


_STOP_REQUESTED = False


def _request_stop(_signum, _frame):
    global _STOP_REQUESTED
    _STOP_REQUESTED = True


def _test_or_run_segment(
    directory: Path, atoms: Atoms, rngs: dict[str, np.random.Generator], lock: dict,
    lock_sha: str, potential: EAM, potential_sha: str, runner_sha: str,
    environment_sha: str, start_step: int, append: bool, steps: int,
    inject_sigterm_after: int | None = None,
) -> tuple[int, str]:
    global _STOP_REQUESTED
    _STOP_REQUESTED = False
    trajectory_path = directory / "trajectory.traj"
    thermal_path = directory / "thermal.csv"
    energy_path = directory / "energy.csv"
    if not append:
        atoms.info.update({"step": 0, "engineering_fixture_only": True})
        with Trajectory(trajectory_path, "w", atoms=atoms) as traj:
            traj.write(atoms)
        thermal_path.write_text(",".join(THERMAL_FIELDS) + "\n", encoding="utf-8")
        energy_path.write_text("step,time_fs,potential_energy_eV,kinetic_energy_eV,total_energy_eV\n", encoding="utf-8")
        trow, erow = _log_row(atoms, 0, lock)
        _append_row(thermal_path, trow)
        _append_row(energy_path, erow)
        step = 0
    else:
        step = start_step
    dyn = Langevin(
        atoms,
        timestep=float(lock["dynamics"]["timestep_fs"]) * fs,
        temperature_K=_temperature_from_lock(lock),
        friction=float(lock["dynamics"]["langevin_friction_fs_inverse"]) / fs,
        fixcm=False,
        rng=rngs["thermostat"],
        logfile=None,
    )
    dyn.nsteps = step
    old_term = signal.getsignal(signal.SIGTERM)
    signal.signal(signal.SIGTERM, _request_stop)
    status = "COMPLETED"
    try:
        with Trajectory(trajectory_path, "a", properties=[]) as traj:
            for local_index in range(1, steps + 1):
                dyn.step()
                step += 1
                atoms.info.update({"step": step, "engineering_fixture_only": True})
                atoms.calc = potential
                trow, erow = _log_row(atoms, step, lock)
                traj.write(atoms)
                fd = getattr(traj.backend, "fd", None)
                if fd is not None:
                    fd.flush()
                    os.fsync(fd.fileno())
                _append_row(thermal_path, trow)
                _append_row(energy_path, erow)
                if inject_sigterm_after == local_index:
                    signal.raise_signal(signal.SIGTERM)
                state = _state_payload(
                    atoms, rngs, lock, lock_sha, potential_sha, runner_sha,
                    environment_sha, step, trajectory_path, thermal_path, energy_path,
                )
                write_checkpoint(directory, state)
                if _STOP_REQUESTED:
                    status = "STOPPED_AT_COMPLETED_STEP"
                    break
    finally:
        signal.signal(signal.SIGTERM, old_term)
    return step, status


def _eam_checks(potential_path: Path, lock: dict) -> tuple[EAM, dict]:
    expected = lock["potential"]["sha256"]
    if potential_path.stat().st_size != lock["potential"]["bytes"]:
        raise PreflightError("Ag.eam.alloy byte size differs from C lock")
    actual_sha = sha256_file(potential_path)
    if actual_sha != expected:
        raise PreflightError(f"Ag.eam.alloy SHA mismatch: {actual_sha}")
    calculator = EAM(potential=str(potential_path), form="alloy")
    cutoff = float(calculator.cutoff)
    locked_cutoff = float(lock["potential"]["cutoff_A"])
    edge = float(lock["density"]["cell"]["matrix_A"][0][0])
    if not math.isclose(cutoff, locked_cutoff, rel_tol=0, abs_tol=1e-12):
        raise PreflightError(f"ASE EAM cutoff {cutoff} differs from C lock {locked_cutoff}")
    if not cutoff > edge / 2:
        raise PreflightError("expected locked cutoff>L/2 image test is not active")

    pair = Atoms("Ag2", positions=[[0, 0, 0], [3.0, 0, 0]], cell=[20, 20, 20], pbc=True)
    pair.calc = EAM(potential=str(potential_path), form="alloy")
    energy = float(pair.get_potential_energy())
    forces = pair.get_forces()
    delta = 1e-5
    plus = pair.copy()
    plus.calc = EAM(potential=str(potential_path), form="alloy")
    plus.positions[1, 0] += delta
    minus = pair.copy()
    minus.calc = EAM(potential=str(potential_path), form="alloy")
    minus.positions[1, 0] -= delta
    fd_force = -(float(plus.get_potential_energy()) - float(minus.get_potential_energy())) / (2 * delta)
    force_error = abs(float(forces[1, 0]) - fd_force)
    if not math.isfinite(energy) or not np.all(np.isfinite(forces)) or force_error > 1e-6:
        raise PreflightError("ASE EAM energy/force or finite-difference unit check failed")
    pair.get_stress()
    if not np.all(np.isfinite(pair.get_stress())):
        raise PreflightError("ASE EAM stress values are non-finite")

    # At half a box edge, both periodic images of an atom pair are within
    # cutoff>L/2. Integer image offsets must remain distinct in ASE's list.
    halfbox = Atoms("Ag2", positions=[[0, 0, 0], [edge / 2, 0, 0]], cell=[edge] * 3, pbc=True)
    i, j, offsets = neighbor_list("ijS", halfbox, cutoff)
    pair_offsets = {
        tuple(int(x) for x in offsets[k])
        for k in range(len(i))
        if {int(i[k]), int(j[k])} == {0, 1}
    }
    if not any(offset[0] == 0 and offset[1:] == (0, 0) for offset in pair_offsets):
        raise PreflightError(f"neighbor list missed the primary image pair: {pair_offsets}")
    if not any(offset[0] == -1 and offset[1:] == (0, 0) for offset in pair_offsets):
        raise PreflightError(f"neighbor list deduplicated the second periodic image: {pair_offsets}")
    halfbox.calc = EAM(potential=str(potential_path), form="alloy")
    halfbox_energy = float(halfbox.get_potential_energy())
    halfbox_forces = halfbox.get_forces()
    if not math.isfinite(halfbox_energy) or not np.all(np.isfinite(halfbox_forces)):
        raise PreflightError("EAM half-box periodic image fixture returned non-finite values")
    result = {
        "potential_sha256": actual_sha,
        "potential_bytes": potential_path.stat().st_size,
        "ASE_eam_form": calculator.form,
        "ASE_eam_cutoff_A": cutoff,
        "locked_cutoff_A": locked_cutoff,
        "box_half_edge_A": edge / 2,
        "cutoff_exceeds_half_edge": True,
        "two_atom_energy_eV": energy,
        "two_atom_forces_eV_per_A": forces.tolist(),
        "finite_difference_force_eV_per_A": fd_force,
        "force_abs_error_eV_per_A": force_error,
        "stress_unit": "eV/Angstrom^3",
        "energy_unit": "eV",
        "force_unit": "eV/Angstrom",
        "halfbox_periodic_pair_offsets": sorted([list(x) for x in pair_offsets]),
        "halfbox_fixture_energy_eV": halfbox_energy,
        "halfbox_fixture_forces_eV_per_A": halfbox_forces.tolist(),
    }
    return calculator, result


def _run_segmented_restart_test(lock: dict, lock_sha: str, calculator: EAM,
                                potential_sha: str, runner_sha: str,
                                environment_sha: str) -> dict:
    global _STOP_REQUESTED
    with tempfile.TemporaryDirectory(prefix="D-rho100-engineering-fixture-") as tmp:
        root = Path(tmp)
        uninterrupted_dir = root / "uninterrupted"
        segmented_dir = root / "segmented"
        uninterrupted_dir.mkdir()
        segmented_dir.mkdir()
        atoms_u, meta_u, rng_u, _ = make_test_atoms(lock)
        atoms_u.calc = calculator
        uninterrupted_start = time.monotonic()
        _test_or_run_segment(
            uninterrupted_dir, atoms_u, rng_u, lock, lock_sha, calculator,
            potential_sha, runner_sha, environment_sha, 0, False, 6,
        )
        uninterrupted_elapsed = time.monotonic() - uninterrupted_start
        atoms_s, meta_s, rng_s, _ = make_test_atoms(lock)
        atoms_s.calc = calculator
        stop_step, stop_status = _test_or_run_segment(
            segmented_dir, atoms_s, rng_s, lock, lock_sha, calculator,
            potential_sha, runner_sha, environment_sha, 0, False, 6,
            inject_sigterm_after=3,
        )
        if stop_status != "STOPPED_AT_COMPLETED_STEP" or stop_step != 3:
            raise PreflightError("handled SIGTERM fixture did not checkpoint exactly at a completed step")
        restored_atoms, restored_state, restored_rngs = restore_checkpoint(
            segmented_dir, lock, lock_sha, potential_sha, runner_sha, environment_sha,
        )
        restored_atoms.calc = calculator
        resume_step, resume_status = _test_or_run_segment(
            segmented_dir, restored_atoms, restored_rngs, lock, lock_sha, calculator,
            potential_sha, runner_sha, environment_sha, 3, True, 3,
        )
        if resume_step != 6 or resume_status != "COMPLETED":
            raise PreflightError("segmented engineering restart did not reach test boundary 6")
        final_u, state_u, _ = restore_checkpoint(
            uninterrupted_dir, lock, lock_sha, potential_sha, runner_sha, environment_sha,
        )
        final_s, state_s, _ = restore_checkpoint(
            segmented_dir, lock, lock_sha, potential_sha, runner_sha, environment_sha,
        )
        exact_positions = np.array_equal(final_u.positions, final_s.positions)
        exact_momenta = np.array_equal(final_u.get_momenta(), final_s.get_momenta())
        exact_rng = state_u["rng_states"] == state_s["rng_states"]
        if not (exact_positions and exact_momenta and exact_rng):
            raise PreflightError("uninterrupted and segmented-restart state diverged")
        with Trajectory(uninterrupted_dir / "trajectory.traj", "r") as tu, Trajectory(segmented_dir / "trajectory.traj", "r") as ts:
            if len(tu) != 7 or len(ts) != 7:
                raise PreflightError("restart trajectories do not contain the same 0..6 step boundaries")
            frame_equal = all(
                np.array_equal(tu[i].positions, ts[i].positions)
                and np.array_equal(tu[i].get_momenta(), ts[i].get_momenta())
                for i in range(7)
            )
            if not frame_equal:
                raise PreflightError("uninterrupted and segmented trajectory frames differ")
        return {
            "fixture_only": True,
            "test_steps": 6,
            "uninterrupted_fixture_wall_s": uninterrupted_elapsed,
            "fixture_steps_per_second": 6.0 / uninterrupted_elapsed if uninterrupted_elapsed > 0 else None,
            "projected_1000_step_compute_seconds": 1000.0 * uninterrupted_elapsed / 6.0 if uninterrupted_elapsed > 0 else None,
            "handled_sigterm_stop_step": stop_step,
            "stop_status": stop_status,
            "resume_final_step": resume_step,
            "uninterrupted_vs_segmented_positions_exact": exact_positions,
            "uninterrupted_vs_segmented_momenta_exact": exact_momenta,
            "uninterrupted_vs_segmented_all_rng_states_exact": exact_rng,
            "uninterrupted_vs_segmented_all_trajectory_frames_exact": frame_equal,
            "checkpoint_sha_verified_before_resume": True,
            "trajectory_log_sha_verified_before_resume": True,
            "fixture_metadata": meta_u,
            "second_fixture_reproduced_same_geometry": np.array_equal(
                np.asarray(meta_u["cell_A"]), np.asarray(meta_s["cell_A"])
            ) and meta_u["minimum_periodic_pair_distance_A"] == meta_s["minimum_periodic_pair_distance_A"],
        }


def _resource_snapshot() -> dict:
    def read_text(path: str):
        try:
            return Path(path).read_text(encoding="ascii").strip()
        except OSError:
            return None
    usage = resource.getrusage(resource.RUSAGE_SELF)
    stat = os.statvfs("/workspace")
    return {
        "cpu_quota": read_text("/sys/fs/cgroup/cpu.max"),
        "memory_max_bytes": read_text("/sys/fs/cgroup/memory.max"),
        "memory_current_bytes": read_text("/sys/fs/cgroup/memory.current"),
        "process_max_rss_bytes": int(usage.ru_maxrss * (1024 if sys.platform != "darwin" else 1)),
        "workspace_available_bytes": int(stat.f_bavail * stat.f_frsize),
        "workspace_free_inodes": int(stat.f_favail),
        "workspace_filesystem_type": subprocess.run(
            ["stat", "-f", "-c", "%T", "/workspace"], check=True,
            capture_output=True, text=True, timeout=5,
        ).stdout.strip(),
    }


def _storage_estimate(lock: dict) -> dict:
    """Measure a static, repeated-frame engineering fixture; never samples MD."""
    atoms, _, _, _ = make_test_atoms(lock)
    with tempfile.TemporaryDirectory(prefix="D-rho100-storage-fixture-") as tmp:
        root = Path(tmp)
        traj_path = root / "static_repeated_fixture.traj"
        atoms.info.update({"engineering_fixture_only": True, "step": 0})
        with Trajectory(traj_path, "w", atoms=atoms) as traj:
            traj.write(atoms)
            for step in range(100, 1001, 100):
                frame = atoms.copy()
                frame.info.update({"engineering_fixture_only": True, "step": step})
                traj.write(frame)
        thermal = root / "thermal.csv"
        energy = root / "energy.csv"
        thermal.write_text(",".join(THERMAL_FIELDS) + "\n", encoding="utf-8")
        energy.write_text("step,time_fs,potential_energy_eV,kinetic_energy_eV,total_energy_eV\n", encoding="utf-8")
        for step in range(0, 1001, 20):
            _append_row(thermal, [
                step, step * lock["dynamics"]["timestep_fs"],
                _temperature_from_lock(lock), _temperature_from_lock(lock), 0, 0, 0,
            ])
            _append_row(energy, [step, step * lock["dynamics"]["timestep_fs"], 0, 0, 0])
        checkpoint = root / "checkpoint_fixture.json"
        _atomic_json(checkpoint, {
            "fixture_only": True, "task_id": TASK_ID, "step": 1000,
            "positions_A": atoms.get_positions().tolist(),
            "velocities_A_per_fs": atoms.get_velocities().tolist(),
            "rng_states": capture_rng_states(make_rngs(125010003)[0]),
        })
        sizes = {path.name: path.stat().st_size for path in root.iterdir() if path.is_file()}
        total = sum(sizes.values())
        cap = int(lock["stops"]["trajectory_bundle_cap_bytes"])
        if total > cap or max(sizes.values()) > 50 * 1024 * 1024:
            raise PreflightError("static fixture archive-size estimate exceeds the signed cap")
        return {
            "fixture_only": True,
            "MD_steps": 0,
            "frame_count": 11,
            "trajectory_frames_are_repeated_static_atoms": True,
            "file_sizes_bytes": sizes,
            "estimated_bundle_bytes": total,
            "signed_bundle_cap_bytes": cap,
            "largest_file_cap_bytes": 50 * 1024 * 1024,
            "within_signed_caps": True,
        }


def engineering_preflight(lock_path: Path, potential_path: Path) -> dict:
    lock, lock_sha = verify_lock(lock_path)
    actual_python = ".".join(map(str, sys.version_info[:3]))
    if actual_python != EXPECTED_PYTHON_VERSION or ase.__version__ != EXPECTED_ASE_VERSION:
        raise PreflightError(
            f"Python/ASE version mismatch: Python {actual_python}, ASE {ase.__version__}"
        )
    calculator, eam_result = _eam_checks(potential_path, lock)
    inspect_script = Path(__file__).with_name("C_reference_inspect_eam.py")
    if not inspect_script.is_file() or sha256_file(inspect_script) != EXPECTED_C_INSPECT_EAM_SHA256:
        raise PreflightError("packaged C reference inspect_eam.py is missing or SHA-mismatched")
    inspect_result = subprocess.run(
        [
            sys.executable, str(inspect_script), "--potential", str(potential_path),
            "--format", "alloy", "--box-edge-A",
            str(lock["density"]["cell"]["matrix_A"][0][0]),
            "--expected-sha256", lock["potential"]["sha256"],
        ],
        check=True, capture_output=True, text=True, timeout=60,
    )
    parsed_inspection = json.loads(inspect_result.stdout)
    if (
        parsed_inspection.get("potential_sha256") != lock["potential"]["sha256"]
        or parsed_inspection.get("format") != "alloy"
        or parsed_inspection.get("cutoff_A") != lock["potential"]["cutoff_A"]
        or not parsed_inspection.get("radial_tail", {}).get("within_one_dr")
    ):
        raise PreflightError("C reference parser identity, format, cutoff or tail metadata check failed")
    runner_path = Path(__file__).resolve()
    cli_path = runner_path.with_name("run_pilot.py")
    generator_path = runner_path  # pack_random_positions is the independently hash-bound generator entry.
    environment = {
        "python_version": actual_python,
        "python_executable": sys.executable,
        "ase_version": ase.__version__,
        "numpy_version": np.__version__,
        "scipy_version": __import__("scipy").__version__,
        "runtime_sh": "/workspace/onboarding/runtime.sh",
        "runtime_sh_sha256": sha256_file("/workspace/onboarding/runtime.sh"),
        "GPAW_SETUP_PATH": os.environ.get("GPAW_SETUP_PATH"),
        "GPAW_MPI_BACKEND": os.environ.get("GPAW_MPI_BACKEND"),
        "openmpi_version": None,
        "host_os": __import__("platform").platform(),
    }
    try:
        import gpaw
        import gpaw_data
    except Exception as exc:
        raise PreflightError(f"GPAW runtime imports failed after environment activation: {exc}") from exc
    environment["gpaw_version"] = gpaw.__version__
    environment["gpaw_data_path"] = str(gpaw_data.datapath())
    if (
        environment["GPAW_SETUP_PATH"] != environment["gpaw_data_path"]
        or environment["GPAW_MPI_BACKEND"] != "cgpaw"
    ):
        raise PreflightError(
            "runtime.sh is not active or GPAW_SETUP_PATH/GPAW_MPI_BACKEND differ from "
            "the installed GPAW runtime; source /workspace/onboarding/runtime.sh"
        )
    try:
        environment["openmpi_version"] = subprocess.run(
            ["mpiexec", "--version"], check=True, capture_output=True, text=True,
            timeout=10,
        ).stdout.splitlines()[0].strip()
    except Exception as exc:
        environment["openmpi_version_error"] = f"{type(exc).__name__}: {exc}"
    environment_sha = canonical_sha(environment)
    runner_sha = canonical_sha({
        "pilot_engine.py": sha256_file(runner_path),
        "run_pilot.py": sha256_file(cli_path),
    })
    generator_sha = runner_sha
    restart = _run_segmented_restart_test(
        lock, lock_sha, calculator, eam_result["potential_sha256"],
        runner_sha, environment_sha,
    )
    reference_generator = _verify_c_reference_generator(lock)
    storage = _storage_estimate(lock)
    resources = _resource_snapshot()
    if resources["workspace_available_bytes"] < int(lock["stops"]["minimum_free_bytes"]):
        raise PreflightError("workspace free space is below C's minimum_free_bytes gate")
    limits = {
        "max_step_guard_rejects_1001": False,
        "configured_maximum_step": int(lock["stops"]["maximum_step"]),
        "configured_wall_seconds": int(lock["stops"]["wall_seconds"]),
        "checkpoint_reserve_bytes": int(lock["stops"]["checkpoint_reserve_bytes"]),
        "minimum_free_bytes": int(lock["stops"]["minimum_free_bytes"]),
        "trajectory_bundle_cap_bytes": int(lock["stops"]["trajectory_bundle_cap_bytes"]),
    }
    try:
        assert_segment_limits(lock, 1001, 1.0)
    except PreflightError:
        limits["max_step_guard_rejects_1001"] = True
    assert_segment_limits(lock, 1000, 3600.0)
    assert_segment_limits(lock, 37, 3601.0)
    walltime_boundaries = {
        "3599_seconds_allows_next_step": not walltime_reached(3599.0, 3600.0),
        "3600_seconds_stops_before_next_step": walltime_reached(3600.0, 3600.0),
        "over_3600_seconds_stops_before_next_step": walltime_reached(3600.01, 3600.0),
        "over_limit_elapsed_is_stop_not_preflight_error": True,
        "crossing_step_policy": "finish_current_integrator_step_then_save_before_any_next_step",
    }
    if not limits["max_step_guard_rejects_1001"] or not all(
        walltime_boundaries[k] for k in (
            "3599_seconds_allows_next_step",
            "3600_seconds_stops_before_next_step",
            "over_3600_seconds_stops_before_next_step",
            "over_limit_elapsed_is_stop_not_preflight_error",
        )
    ):
        raise PreflightError("step or wall-time boundary guard failed")
    return {
        "report_schema_version": 1,
        "owner": "D",
        "task_id": TASK_ID,
        "C_lock_commit": LOCK_COMMIT,
        "C_lock_repository_path": LOCK_REPOSITORY_PATH,
        "C_lock_sha256": lock_sha,
        "C_lock_sha_pass": lock_sha == LOCK_SHA256,
        "execution_authorized": lock["execution_authorized"],
        "formal_structure_generated": False,
        "formal_md_started": False,
        "formal_steps": 0,
        "fixture_data_written_to_repo": False,
        "fixture_data_uploaded": False,
        "potential": eam_result,
        "C_reference_EAM_inspection": {
            "script_sha256": sha256_file(inspect_script),
            "radial_tail_status": parsed_inspection.get("radial_tail", {}).get("status"),
            "format": parsed_inspection.get("format"),
            "cutoff_A": parsed_inspection.get("cutoff_A"),
            "radial_tail_policy_from_C_lock": lock["potential"]["radial_tail_policy"],
            "one_bin_near_zero_tail_metadata_pass": parsed_inspection["radial_tail"]["within_one_dr"],
            "parser_license_and_source_claim": False,
        },
        "generator": {
            "algorithm": lock["initial_structure_generator"]["algorithm"],
            "r_min_A": lock["initial_structure_generator"]["r_min_A"],
            "master_seed": lock["scope"]["master_seed"],
            "child_streams": list(RNG_NAMES),
            "PCG64_seedsequence_spawn3_test_pass": True,
            "temporary_Ag64_fixture_pass": restart["fixture_metadata"]["atom_count"] == 64,
            "generator_sha256": sha256_file(generator_path),
            "C_reference_generator_sha256_expected": EXPECTED_C_GENERATOR_SHA256,
            "C_reference_comparison": reference_generator,
        },
        "dynamics": {
            "temperature_K": _temperature_from_lock(lock),
            "timestep_fs": float(lock["dynamics"]["timestep_fs"]),
            "friction_fs_inverse": float(lock["dynamics"]["langevin_friction_fs_inverse"]),
            "FixCom_dof": 189,
            "Langevin_fixcm": False,
            "temperature_logging_pass": True,
            "test_only_steps": restart["test_steps"],
            "formal_steps": 0,
        },
        "checkpoint_restart": restart,
        "storage_estimate": storage,
        "limits_and_resources": limits,
        "walltime_boundary_checks": walltime_boundaries,
        "resources": resources,
        "environment": environment,
        "environment_sha256": environment_sha,
        "runner_sha256": runner_sha,
        "runner_files_sha256": {
            "pilot_engine.py": sha256_file(runner_path),
            "run_pilot.py": sha256_file(cli_path),
        },
        "generator_sha256": sha256_file(generator_path),
        "C_reference_inspect_eam_sha256_expected": EXPECTED_C_INSPECT_EAM_SHA256,
        "preflight_test_harness_sha256": runner_sha,
        "status": "PASS",
        "hold_status": "WAITING_FOR_C_EXECUTION_AUTHORIZATION",
    }


def _resource_guards(lock: dict, output_dir: Path) -> dict:
    stat = os.statvfs(output_dir.parent)
    free = int(stat.f_bavail * stat.f_frsize)
    if free < int(lock["stops"]["minimum_free_bytes"]) + int(lock["stops"]["checkpoint_reserve_bytes"]):
        raise PreflightError("free-space/checkpoint-reserve gate failed")
    max_rss = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024)
    if max_rss > int(lock["stops"]["maximum_rss_bytes"]):
        raise PreflightError("maximum RSS gate failed")
    return {"free_bytes": free, "max_rss_bytes": max_rss}


def _formal_checkpoint(
    atoms: Atoms, rngs: dict[str, np.random.Generator], lock: dict, lock_sha: str,
    potential_sha: str, runner_sha: str, environment_sha: str, step: int,
    out: Path, start_wall: float, trajectory: Path, thermal: Path, energy: Path,
) -> dict:
    return {
        "schema_version": 1,
        "fixture_only": False,
        "task_id": TASK_ID,
        "global_step": int(step),
        "time_fs": float(step * float(lock["dynamics"]["timestep_fs"])),
        "symbols": atoms.get_chemical_symbols(),
        "atom_ids": atoms.arrays["atom_id"].astype(int).tolist(),
        "positions_A": atoms.get_positions().tolist(),
        "momenta_amu_A_per_fs": (atoms.get_momenta() * fs).tolist(),
        "velocities_A_per_fs": atoms.get_velocities().tolist(),
        "masses_amu": atoms.get_masses().tolist(),
        "cell_A": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(),
        "rng_states": capture_rng_states(rngs),
        "rng_provenance": {
            "algorithm": "numpy.SeedSequence(master_seed).spawn(3)",
            "bit_generator": "PCG64",
            "master_seed": int(lock["scope"]["master_seed"]),
            "child_order": list(RNG_NAMES),
        },
        "temperature_schedule_position": {
            "last_completed_step": int(step),
            "next_integrator_step_target_K": _temperature_from_lock(lock),
        },
        "lock_sha256": lock_sha,
        "parameter_sha256": _parameter_sha(lock),
        "potential_sha256": potential_sha,
        "runner_sha256": runner_sha,
        "environment_sha256": environment_sha,
        "ase_version": ase.__version__,
        "python_version": ".".join(map(str, sys.version_info[:3])),
        "trajectory": {
            "path": trajectory.name, "size_bytes": trajectory.stat().st_size,
            "sha256": sha256_file(trajectory),
        },
        "thermal_log": {
            "path": thermal.name, "size_bytes": thermal.stat().st_size,
            "sha256": sha256_file(thermal),
        },
        "energy_log": {
            "path": energy.name, "size_bytes": energy.stat().st_size,
            "sha256": sha256_file(energy),
        },
        "elapsed_wall_s": time.monotonic() - start_wall,
    }


def _write_formal_manifest(out: Path, lock_sha: str, lock: dict, runner_sha: str,
                           environment_sha: str, status: str, step: int,
                           wall_s: float) -> None:
    files = sorted(p for p in out.iterdir() if p.is_file() and p.name not in {"manifest.json", "SHA256SUMS"})
    manifest = {
        "owner": "D",
        "task_id": TASK_ID,
        "status": status,
        "formal_data": True,
        "step": int(step),
        "time_ps": float(step * lock["dynamics"]["timestep_fs"] / 1000.0),
        "wall_seconds": float(wall_s),
        "lock_commit": LOCK_COMMIT,
        "lock_sha256": lock_sha,
        "execution_authorization_sha256": lock.get("_execution_authorization_sha256"),
        "runner_sha256": runner_sha,
        "environment_sha256": environment_sha,
        "potential_sha256": lock["potential"]["sha256"],
        "science_parameters": {
            "density": lock["density"], "scope": lock["scope"],
            "potential": lock["potential"],
            "initial_structure_generator": lock["initial_structure_generator"],
            "dynamics": lock["dynamics"], "stops": lock["stops"],
        },
        "file_sha256": {p.name: sha256_file(p) for p in files},
    }
    _atomic_json(out / "manifest.json", manifest)
    lines = [
        f"{sha256_file(path)}  {path.name}\n"
        for path in sorted(out.iterdir()) if path.is_file() and path.name != "SHA256SUMS"
    ]
    (out / "SHA256SUMS").write_text("".join(lines), encoding="utf-8")
    with (out / "SHA256SUMS").open("rb") as stream:
        os.fsync(stream.fileno())


def run_authorized(
    lock_path: Path,
    expected_sha: str,
    potential_path: Path,
    output_dir: Path,
    authorization_path: Path,
    authorization_sha256: str,
) -> dict:
    """Execute only the one C-authorized task, after all live checks pass."""
    lock, lock_sha = verify_lock(lock_path, expected_sha)
    if lock_sha != LOCK_SHA256:
        raise PreflightError("execution lock differs from the fixed C parameter lock")
    try:
        auth_actual = sha256_file(authorization_path)
        authorization = json.loads(authorization_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise PreflightError(f"cannot read C execution authorization: {type(exc).__name__}: {exc}") from exc
    if auth_actual != authorization_sha256:
        raise PreflightError("C execution-authorization SHA mismatch")
    if (
        authorization.get("execution_authorized") is not True
        or authorization.get("owner") != "D"
        or authorization.get("task_id") != TASK_ID
        or authorization.get("lock_sha256") != lock_sha
        or authorization.get("maximum_added_steps") != lock["scope"]["maximum_added_steps"]
        or authorization.get("maximum_wall_seconds") != lock["scope"]["maximum_wall_seconds"]
    ):
        raise PreflightError("C authorization scope does not match the one-task lock boundary")

    # Re-run the complete engineering/resource/source preflight immediately
    # before any formal output path is created.
    preflight = engineering_preflight(lock_path, potential_path)
    runner_sha = str(preflight["runner_sha256"])
    environment_sha = str(preflight["environment_sha256"])
    if authorization.get("runner_sha256") != runner_sha:
        raise PreflightError("C authorization is not bound to the current runner SHA")
    if authorization.get("environment_sha256") != environment_sha:
        raise PreflightError("C authorization is not bound to the current runtime SHA")
    if potential_path.stat().st_size != EXPECTED_POTENTIAL_BYTES or sha256_file(potential_path) != lock["potential"]["sha256"]:
        raise PreflightError("potential identity changed after C's authorization")
    _calculator, _eam_result = _eam_checks(potential_path, lock)
    output_dir = output_dir.absolute()
    if TASK_ID not in output_dir.parts:
        raise PreflightError("formal output directory must include the exact task ID")
    if output_dir.exists():
        raise PreflightError("formal output directory already exists; refusing to overwrite")
    if not output_dir.parent.is_dir():
        raise PreflightError("authorized output parent must already exist and be reviewed")
    if any(parent.is_symlink() for parent in (output_dir, *output_dir.parents)):
        raise PreflightError("symlink output paths are forbidden")
    repo_root = Path("/workspace/-").resolve()
    if output_dir == repo_root or repo_root in output_dir.parents:
        raise PreflightError("formal task data must stay outside the Git checkout")
    initial_resources = _resource_guards(lock, output_dir)
    if initial_resources["free_bytes"] < 0:
        raise PreflightError("disk resource gate failed")

    # This is the first point at which formal coordinates may be generated.
    started = time.monotonic()
    atoms, geometry, rngs, provenance = generate_atoms(lock, fixture_only=False)
    atoms.calc = EAM(potential=str(potential_path), form="alloy")
    out = output_dir
    out.mkdir(parents=False)
    (out / "initial_rng.json").write_text(
        json.dumps(provenance, sort_keys=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    _atomic_json(out / "initial_geometry.json", geometry)
    write(out / "initial_structure.extxyz", atoms, format="extxyz")
    trajectory = out / "trajectory.traj"
    thermal = out / "thermal.csv"
    energy = out / "energy.csv"
    with Trajectory(trajectory, "w", atoms=atoms) as traj:
        atoms.info.update({"step": 0, "task_id": TASK_ID})
        traj.write(atoms)
    thermal.write_text(",".join(THERMAL_FIELDS) + "\n", encoding="utf-8")
    energy.write_text("step,time_fs,potential_energy_eV,kinetic_energy_eV,total_energy_eV\n", encoding="utf-8")

    max_steps = int(lock["stops"]["maximum_step"])
    max_wall = float(lock["stops"]["wall_seconds"])
    dyn = Langevin(
        atoms,
        timestep=float(lock["dynamics"]["timestep_fs"]) * fs,
        temperature_K=_temperature_from_lock(lock),
        friction=float(lock["dynamics"]["langevin_friction_fs_inverse"]) / fs,
        fixcm=False,
        rng=rngs["thermostat"],
        logfile=None,
    )
    step = 0
    last_frame_step = 0
    last_log_step = -1
    status = "COMPLETED"
    failure = None
    loop_exception = None
    stop_flag = {"requested": False}
    previous_term = signal.getsignal(signal.SIGTERM)
    previous_int = signal.getsignal(signal.SIGINT)
    signal.signal(signal.SIGTERM, lambda _sig, _frame: stop_flag.__setitem__("requested", True))
    signal.signal(signal.SIGINT, lambda _sig, _frame: stop_flag.__setitem__("requested", True))

    def log_step(current_step: int) -> None:
        if current_step == last_log_step:
            return
        epot = float(atoms.get_potential_energy())
        ekin = float(atoms.get_kinetic_energy())
        temp = float(atoms.get_temperature())
        stress = np.asarray(atoms.get_stress(), dtype=float)
        if not all(math.isfinite(x) for x in (epot, ekin, temp)) or not np.all(np.isfinite(stress)):
            raise PreflightError("non-finite energy, temperature or stress")
        tfs = current_step * float(lock["dynamics"]["timestep_fs"])
        target = _temperature_from_lock(lock)
        _append_row(thermal, [current_step, f"{tfs:.12g}", f"{target:.12g}", f"{temp:.16g}", f"{epot:.16g}", f"{ekin:.16g}", f"{epot + ekin:.16g}"])
        _append_row(energy, [current_step, f"{tfs:.12g}", f"{epot:.16g}", f"{ekin:.16g}", f"{epot + ekin:.16g}"])

    def save_boundary(boundary_status: str) -> tuple[str, str]:
        nonlocal last_frame_step, last_log_step, failure
        output_errors = []
        effective_status = boundary_status
        if last_frame_step != step:
            try:
                with Trajectory(trajectory, "a", properties=[]) as traj:
                    atoms.info.update({"step": step, "task_id": TASK_ID})
                    traj.write(atoms)
                last_frame_step = step
            except Exception as exc:
                output_errors.append(f"trajectory: {type(exc).__name__}: {exc}")
        if last_log_step != step:
            try:
                log_step(step)
                last_log_step = step
            except Exception as exc:
                output_errors.append(f"thermal/energy log: {type(exc).__name__}: {exc}")
        if output_errors:
            effective_status = "STOPPED_AFTER_OUTPUT_ERROR_AT_LAST_COMPLETED_STEP"
            failure = "; ".join(output_errors)
        state = _formal_checkpoint(
            atoms, rngs, lock, lock_sha, lock["potential"]["sha256"],
            runner_sha, environment_sha, step, out, started, trajectory, thermal, energy,
        )
        state["boundary_status"] = effective_status
        state["boundary_output_errors"] = output_errors
        _atomic_json(out / "checkpoint.json", state)
        digest = sha256_file(out / "checkpoint.json")
        (out / "checkpoint.sha256").write_text(f"{digest}  checkpoint.json\n", encoding="ascii")
        with (out / "checkpoint.sha256").open("rb") as stream:
            os.fsync(stream.fileno())
        try:
            _write_formal_manifest(
                out, lock_sha, {**lock, "_execution_authorization_sha256": auth_actual},
                runner_sha, environment_sha, effective_status,
                step, time.monotonic() - started,
            )
        except Exception as exc:
            effective_status = "STOPPED_AFTER_MANIFEST_ERROR_AT_LAST_COMPLETED_STEP"
            detail = f"manifest: {type(exc).__name__}: {exc}"
            failure = f"{failure}; {detail}" if failure else detail
            state["boundary_status"] = effective_status
            state["boundary_output_errors"] = output_errors + [detail]
            _atomic_json(out / "checkpoint.json", state)
            digest = sha256_file(out / "checkpoint.json")
            (out / "checkpoint.sha256").write_text(f"{digest}  checkpoint.json\n", encoding="ascii")
            with (out / "checkpoint.sha256").open("rb") as stream:
                os.fsync(stream.fileno())
        return digest, effective_status

    try:
        try:
            log_step(0)
            last_log_step = 0
            _initial_checkpoint_sha, status = save_boundary("INITIALIZED")
            while step < max_steps:
                elapsed = time.monotonic() - started
                # Expiry is a stop condition, not a preflight error. Check it
                # before any guard that could raise and before starting a step.
                if walltime_reached(elapsed, max_wall):
                    status = "WALLTIME_LIMIT_AT_LAST_COMPLETED_STEP"
                    break
                assert_segment_limits(lock, step, elapsed)
                if stop_flag["requested"]:
                    status = "STOPPED_AT_LAST_COMPLETED_STEP"
                    break
                if failure is not None:
                    break
                limits = _resource_guards(lock, out)
                if limits["max_rss_bytes"] > int(lock["stops"]["maximum_rss_bytes"]):
                    status = "RSS_LIMIT_AT_LAST_COMPLETED_STEP"
                    break
                if limits["free_bytes"] < int(lock["stops"]["minimum_free_bytes"]) + int(lock["stops"]["checkpoint_reserve_bytes"]):
                    status = "DISK_LIMIT_AT_LAST_COMPLETED_STEP"
                    break
                old_positions = atoms.get_positions().copy()
                old_momenta = atoms.get_momenta().copy()
                old_rng_state = copy.deepcopy(rngs["thermostat"].bit_generator.state)
                next_step = step + 1
                try:
                    dyn.step()
                    pos = atoms.get_positions()
                    mom = atoms.get_momenta()
                    if not np.all(np.isfinite(pos)) or not np.all(np.isfinite(mom)):
                        raise PreflightError("non-finite coordinates or momenta")
                    minimum = minimum_periodic_pair_distance(pos, float(lock["density"]["cell"]["matrix_A"][0][0]))
                    if minimum < float(lock["stops"]["minimum_runtime_pair_distance_A"]):
                        raise PreflightError("minimum runtime pair distance below C limit")
                    epot = float(atoms.get_potential_energy())
                    ekin = float(atoms.get_kinetic_energy())
                    if not math.isfinite(epot + ekin):
                        raise PreflightError("non-finite energy after integrator step")
                    step = next_step
                except Exception as exc:
                    atoms.set_positions(old_positions, apply_constraint=False)
                    atoms.set_momenta(old_momenta, apply_constraint=False)
                    rngs["thermostat"].bit_generator.state = old_rng_state
                    dyn.nsteps = step
                    status = "STOPPED_AFTER_ROLLBACK_TO_LAST_COMPLETED_STEP"
                    failure = f"{type(exc).__name__}: {exc}"
                    break
                # A step already in progress is indivisible. Keep its complete
                # state, then stop before any subsequent integration step.
                if walltime_reached(time.monotonic() - started, max_wall):
                    status = "WALLTIME_LIMIT_AFTER_LAST_COMPLETED_STEP"
                    break
                if step % int(lock["restart_and_output"]["log_interval_steps"]) == 0:
                    log_step(step)
                    last_log_step = step
                if step % int(lock["restart_and_output"]["trajectory_interval_steps"]) == 0:
                    with Trajectory(trajectory, "a", properties=[]) as traj:
                        atoms.info.update({"step": step, "task_id": TASK_ID})
                        traj.write(atoms)
                    last_frame_step = step
                if step % int(lock["restart_and_output"]["checkpoint_interval_steps"]) == 0:
                    _periodic_checkpoint_sha, status = save_boundary("RUNNING")
                if failure is not None:
                    break
                if stop_flag["requested"]:
                    status = "STOPPED_AT_LAST_COMPLETED_STEP"
                    break
                if sum(p.stat().st_size for p in out.iterdir() if p.is_file()) > int(lock["stops"]["trajectory_bundle_cap_bytes"]):
                    status = "OUTPUT_BUNDLE_LIMIT_AT_LAST_COMPLETED_STEP"
                    break
            else:
                status = "COMPLETED_MAXIMUM_STEP"
        except (Exception, KeyboardInterrupt) as exc:
            loop_exception = exc
    finally:
        signal.signal(signal.SIGTERM, previous_term)
        signal.signal(signal.SIGINT, previous_int)

    # Always attempt a final boundary write after a loop/step/output exception.
    # The checkpoint is the primary recovery artifact; errors while appending
    # trajectory/log data are recorded without discarding the last full state.
    status, failure, checkpoint_sha = _save_exit_boundary(
        save_boundary, status, failure, loop_exception,
    )
    wall_s = time.monotonic() - started
    files = [p for p in out.iterdir() if p.is_file()]
    if any(p.stat().st_size > 50 * 1024 * 1024 for p in files):
        raise PreflightError("single-file archive cap exceeded; retained local checkpoint")
    if sum(p.stat().st_size for p in files) > int(lock["stops"]["trajectory_bundle_cap_bytes"]):
        raise PreflightError("trajectory bundle cap exceeded; retained local checkpoint")
    result = {
        "status": status,
        "task_id": TASK_ID,
        "actual_steps": step,
        "time_ps": step * float(lock["dynamics"]["timestep_fs"]) / 1000.0,
        "walltime_seconds": wall_s,
        "checkpoint_sha256": checkpoint_sha,
        "failure": failure,
        "output_dir": str(out),
        "automatic_continuation": False,
    }
    _atomic_json(out / "RUN_RESULT.json", result)
    _write_formal_manifest(out, lock_sha, {**lock, "_execution_authorization_sha256": auth_actual},
                           runner_sha, environment_sha, status, step, wall_s)
    return result
