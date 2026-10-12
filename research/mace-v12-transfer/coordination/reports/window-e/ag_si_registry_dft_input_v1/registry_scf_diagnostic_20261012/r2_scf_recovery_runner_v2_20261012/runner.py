#!/usr/bin/env python3
"""Hash-bound, R2-only SCF diagnostic runner. Default mode is read-only preflight."""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import sys
import threading
import time
import traceback
from typing import Any, Callable


HERE = Path(__file__).resolve().parent
CONFIG_PATH = HERE / "R2_RECOVERY_CONFIG.json"
EXPECTED_LOCK_PATH = (
    "research/mace-v12-transfer/coordination/reports/window-e/"
    "ag_si_registry_dft_input_v1/registry_phase2_result_collector_20261012/"
    "PHASE2_PARAMETER_LOCK.json"
)
EXPECTED_LOCK_SHA = "d97de3fe5a26ab6529d5e6e8ceef15adc1f3a2f548f03b59f9f092361c2c0c78"
WFS_BYTES = 1024913408
PROJECTOR_BYTES = 384342528


class HoldError(RuntimeError):
    pass


class CheckpointFailure(RuntimeError):
    pass


class SafeStop(RuntimeError):
    def __init__(self, reason: str, step: int):
        super().__init__(reason)
        self.reason = reason
        self.step = step


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def find_repo_root(start: Path = HERE) -> Path:
    for parent in (start, *start.parents):
        if (parent / ".git").exists():
            return parent
    raise HoldError("cannot locate the Git repository root")


REPO_ROOT = find_repo_root()


def load_config() -> tuple[dict[str, Any], str]:
    raw = CONFIG_PATH.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise HoldError(message)


def verify_parameter_lock(config: dict[str, Any]) -> tuple[dict[str, Any], Path]:
    source = config["source_bindings"]
    lock_rel = source["parameter_lock_path"]
    require(lock_rel == EXPECTED_LOCK_PATH,
            "parameter lock path is not the Phase2 collector lock")
    lock_path = REPO_ROOT / lock_rel
    require(lock_path.is_file(), f"Phase2 parameter lock missing: {lock_rel}")
    observed_lock_sha = sha256_file(lock_path)
    require(observed_lock_sha == EXPECTED_LOCK_SHA == source["parameter_lock_sha256"],
            "Phase2 parameter lock SHA mismatch")
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    method = config["method"]
    base = lock["method"]
    require(lock["structures"]["R2"]["input_sha256"] == source["input_sha256"],
            "R2 input SHA in lock and recovery config differ")
    require(lock["structures"]["R2"]["structure_sha256"] == source["structure_sha256"],
            "R2 structure SHA in lock and recovery config differ")
    require(lock["kpoint_mesh"] == method["kpoint_mesh"] == [8, 8, 1],
            "R2 k-point lock mismatch")
    require(lock["role"] == "REGISTRY_SCREENING_ONLY",
            "unexpected base task role")
    require(base["xc"] == method["xc"] == "PBE", "XC mismatch")
    require(base["cutoff_eV"] == method["cutoff_eV"] == 500, "cutoff mismatch")
    require(base["occupations"] == method["occupations"], "smearing mismatch")
    require(base["convergence"] == method["convergence"],
            "frozen SCF convergence thresholds changed")
    require(base["mixer"] == method["mixer"], "mixer changed from Phase2 lock")
    require(base["max_scf_iterations"] == method["baseline_max_scf_iterations"] == 160,
            "wrong baseline maxiter; older top-level lock may have been used")
    require(method["max_scf_iterations"] == 200, "recovery maxiter must be 200")
    require(method["eigensolver"] == {"name": "PPCG", "niter": 4},
            "PPCG niter=4 is not locked")
    require(lock["resources"] == {"mpi_ranks": 4, "omp_threads": 1},
            "MPI/OMP values differ from Phase2 lock")
    require(lock["paw_sha256_by_element"] == config["paw_sha256_by_element"],
            "PAW hashes differ from the Phase2 parameter lock")
    require(method["mpi_ranks"] == 4 and method["omp_threads"] == 1,
            "recovery MPI/OMP configuration mismatch")
    require(method["fixed_geometry"] is True and method["relax_atoms"] is False
            and method["relax_cell"] is False
            and method["coordinates_cell_pbc_changes_allowed"] is False,
            "geometry must remain fixed")
    return lock, lock_path


def verify_sources(config: dict[str, Any]) -> dict[str, Any]:
    from ase.io import read

    source = config["source_bindings"]
    input_path = REPO_ROOT / source["input_path"]
    structure_path = REPO_ROOT / source["structure_path"]
    require(input_path.is_file(), "R2 input.py is missing")
    require(structure_path.is_file(), "frozen R2 structure is missing")
    observed_input = sha256_file(input_path)
    observed_structure = sha256_file(structure_path)
    require(observed_input == source["input_sha256"], "R2 input SHA mismatch")
    require(observed_structure == source["structure_sha256"],
            "R2 structure SHA mismatch")
    atoms = read(structure_path, format="extxyz")
    require(len(atoms) == 20, "frozen R2 structure must contain 20 atoms")
    require(atoms.pbc.tolist() == [True, True, False], "frozen R2 PBC mismatch")
    return {
        "input_path": str(input_path),
        "input_sha256": observed_input,
        "structure_path": str(structure_path),
        "structure_sha256": observed_structure,
        "atom_count": len(atoms),
        "cell_A": atoms.cell.array.tolist(),
        "pbc": atoms.pbc.tolist(),
    }


def read_cgroup_int(path: Path) -> int | None:
    try:
        value = path.read_text().strip()
        return None if value == "max" else int(value)
    except (OSError, ValueError):
        return None


def collect_resources(work_parent: Path) -> dict[str, Any]:
    cpu_quota = None
    try:
        quota_text, period_text = Path("/sys/fs/cgroup/cpu.max").read_text().split()
        if quota_text != "max":
            cpu_quota = int(quota_text) / int(period_text)
    except (OSError, ValueError, ZeroDivisionError):
        pass
    memory_max = read_cgroup_int(Path("/sys/fs/cgroup/memory.max"))
    memory_current = read_cgroup_int(Path("/sys/fs/cgroup/memory.current"))
    mem_available_kib = _mem_available_kib()
    usage = shutil.disk_usage(work_parent)
    stat = os.statvfs(work_parent)
    method = load_config()[0]["method"]
    budget = load_config()[0]["checkpoint_space_budget"]
    raw_cp = WFS_BYTES + PROJECTOR_BYTES
    exact_per = int(raw_cp * budget["per_checkpoint_estimate_multiplier"])
    estimated_working_set = 2 * raw_cp
    memory_headroom = (memory_max - memory_current
                       if memory_max is not None and memory_current is not None else None)
    memory_pass = (memory_headroom is not None
                   and memory_headroom >= estimated_working_set
                   and mem_available_kib is not None
                   and mem_available_kib * 1024 >= estimated_working_set)
    required = (
        (budget["maximum_periodic_checkpoints"]
         + budget["simultaneous_write_temporary_copies"]) * exact_per
        + budget["log_allowance_bytes"]
        + budget["free_space_reserve_bytes"]
    )
    return {
        "cpu": {
            "os_cpu_count": os.cpu_count(),
            "cgroup_quota_cores": cpu_quota,
            "required_mpi_ranks": method["mpi_ranks"],
            "omp_threads": method["omp_threads"],
            "pass": cpu_quota is not None and cpu_quota >= method["mpi_ranks"],
        },
        "memory": {
            "cgroup_limit_bytes": memory_max,
            "cgroup_current_bytes": memory_current,
            "cgroup_headroom_bytes": memory_headroom,
            "estimated_r2_working_set_bytes": estimated_working_set,
            "mem_available_kib": mem_available_kib,
            "pass": memory_pass,
        },
        "storage": {
            "checked_path": str(work_parent),
            "filesystem_total_bytes": usage.total,
            "filesystem_free_bytes": usage.free,
            "inode_free": stat.f_favail,
            "checkpoint_bytes_per_estimate": exact_per,
            "maximum_checkpoint_count": budget["maximum_periodic_checkpoints"],
            "temporary_checkpoint_count": budget["simultaneous_write_temporary_copies"],
            "log_allowance_bytes": budget["log_allowance_bytes"],
            "free_space_reserve_bytes": budget["free_space_reserve_bytes"],
            "required_free_bytes": required,
            "pass": usage.free >= required and stat.f_favail > 1000,
        },
    }


def _mem_available_kib() -> int | None:
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1])
    except (OSError, ValueError):
        pass
    return None


def verify_runtime(check_mpi: bool) -> dict[str, Any]:
    import ase
    import gpaw
    import gpaw_data
    from importlib.metadata import version

    config, _ = load_config()
    method = config["method"]
    require(ase.__version__ == config["engine_versions"]["ase"], "ASE version mismatch")
    require(gpaw.__version__ == config["engine_versions"]["gpaw"], "GPAW version mismatch")
    require(version("gpaw-data") == config["engine_versions"]["gpaw_data"],
            "gpaw-data version mismatch")
    expected_setup_path = Path(gpaw_data.datapath()).resolve()
    observed_setup_path = os.environ.get("GPAW_SETUP_PATH")
    require(observed_setup_path is not None, "GPAW_SETUP_PATH is unset")
    require(Path(observed_setup_path).resolve() == expected_setup_path,
            "GPAW_SETUP_PATH does not resolve to gpaw_data.datapath()")
    require(os.environ.get("OMP_NUM_THREADS") == "1", "OMP_NUM_THREADS must equal 1")
    paw_observed = {}
    from gpaw.setup_data import search_for_file
    for element, expected in config["paw_sha256_by_element"].items():
        filename, _ = search_for_file(f"{element}.PBE")
        actual = sha256_file(Path(filename))
        require(actual == expected, f"PAW hash mismatch for {element}")
        paw_observed[element] = {"path": str(filename), "sha256": actual}
    mpi_size = None
    if check_mpi:
        from gpaw.mpi import world
        mpi_size = world.size
        require(mpi_size == method["mpi_ranks"], "MPI world size must equal 4")
    return {
        "python": sys.version.split()[0],
        "ase": ase.__version__,
        "gpaw": gpaw.__version__,
        "gpaw_data": version("gpaw-data"),
        "mpi_world_size": mpi_size,
        "omp_threads": os.environ.get("OMP_NUM_THREADS"),
        "gpaw_setup_path": str(expected_setup_path),
        "paw_files": paw_observed,
    }


def run_preflight(work_parent: Path) -> dict[str, Any]:
    config, config_sha = load_config()
    _, lock_path = verify_parameter_lock(config)
    sources = verify_sources(config)
    runtime = verify_runtime(check_mpi=False)
    resources = collect_resources(work_parent)
    return {
        "status": "PREFLIGHT_PASS" if resources["cpu"]["pass"]
        and resources["memory"]["pass"] and resources["storage"]["pass"] else "HOLD",
        "execution_authorized": False,
        "r2_dft_started": False,
        "config_sha256": config_sha,
        "parameter_lock_path": str(lock_path),
        "parameter_lock_sha256": sha256_file(lock_path),
        "sources": sources,
        "runtime": runtime,
        "resources": resources,
    }


def validate_authorization(auth_path: Path, auth_sha_expected: str,
                           config: dict[str, Any], config_sha: str) -> dict[str, Any]:
    require(auth_path.is_file(), "C execution authorization file is missing")
    actual_auth_sha = sha256_file(auth_path)
    require(actual_auth_sha == auth_sha_expected,
            "C execution authorization SHA-256 mismatch")
    auth = json.loads(auth_path.read_text(encoding="utf-8"))
    approval = auth.get("C_approval", {})
    require(auth.get("execution_authorized") is True,
            "C authorization does not set execution_authorized=true")
    require(approval.get("signed_by") == "C" and approval.get("approval_id"),
            "C approval signature fields missing")
    require(auth.get("scope") == "R2_ONLY_SINGLE_DIAGNOSTIC",
            "authorization scope must cover R2 only")
    bindings = auth.get("bindings", {})
    expected = {
        "config_sha256": config_sha,
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "collector_sha256": sha256_file(HERE / "collector.py"),
        "parameter_lock_sha256": config["source_bindings"]["parameter_lock_sha256"],
        "input_sha256": config["source_bindings"]["input_sha256"],
        "structure_sha256": config["source_bindings"]["structure_sha256"],
    }
    for key, value in expected.items():
        require(bindings.get(key) == value, f"authorization binding mismatch: {key}")
    require(auth.get("limits") == {
        "registry_id": "R2",
        "max_scf_iterations": 200,
        "max_wall_clock_seconds": 14400,
        "mpi_ranks": 4,
        "omp_threads": 1,
    }, "C authorization limits do not exactly match R2 recovery config")
    return auth


class CheckpointManager:
    """Collective full-GPW checkpoint writer with verified immutable snapshots."""

    def __init__(self, calculator: Any, atoms: Any, work_root: Path,
                 world: Any, config: dict[str, Any], config_sha: str,
                 lock_sha: str, clock: Callable[[], float] = time.monotonic,
                 started_at: float | None = None,
                 stop_event: threading.Event | None = None,
                 gpw_reader: Callable[..., Any] | None = None):
        self.calculator = calculator
        self.atoms = atoms
        self.work_root = work_root
        self.directory = work_root / "checkpoint"
        self.world = world
        self.config = config
        self.config_sha = config_sha
        self.lock_sha = lock_sha
        self.clock = clock
        self.started_at = clock() if started_at is None else started_at
        self.stop_event = stop_event or threading.Event()
        self.gpw_reader = gpw_reader
        self.saved_steps: set[int] = set()
        self.latest_checkpoint: Path | None = None

    def _reader(self):
        if self.gpw_reader is not None:
            return self.gpw_reader
        from gpaw.new.gpw import read_gpw
        return read_gpw

    def _verify_gp_w(self, path: Path) -> None:
        if self.gpw_reader is not None:
            restored = self.gpw_reader(path)
            require(bool(restored.get("wavefunctions")),
                    "synthetic GPW readback did not contain wavefunctions")
            return
        from gpaw.mpi import serial_comm
        restored_atoms, restored_dft, _ = self._reader()(path, comm=serial_comm)
        require((restored_atoms.numbers == self.atoms.numbers).all(),
                "checkpoint atom numbers differ")
        require((restored_atoms.pbc == self.atoms.pbc).all(), "checkpoint PBC differs")
        require((restored_atoms.positions == self.atoms.positions).all(),
                "checkpoint coordinates differ")
        require((restored_atoms.cell.array == self.atoms.cell.array).all(),
                "checkpoint cell differs")
        require(restored_dft.ibzwfs.has_wave_functions(),
                "GPW did not restore wavefunctions")
        require(restored_dft.scf_loop.eigensolver.niter == 4,
                "GPW readback PPCG niter is not 4")
        del restored_atoms, restored_dft
        gc.collect()

    def _checkpoint_budget_ok(self) -> None:
        estimate = self.config["checkpoint_space_budget"]["per_checkpoint_estimate_bytes"]
        reserve = self.config["checkpoint_space_budget"]["free_space_reserve_bytes"]
        current = shutil.disk_usage(self.work_root).free
        require(current >= estimate + reserve,
                f"remaining disk below one checkpoint plus reserve: {current} bytes")

    def save(self, step: int) -> Path:
        require(step > 0, "checkpoint iteration must be positive")
        require(step not in self.saved_steps, f"duplicate checkpoint iteration {step}")
        self.directory.mkdir(parents=True, exist_ok=True)
        final_path = self.directory / f"checkpoint_iter_{step:06d}.gpw"
        temporary_path = self.directory / f".checkpoint_iter_{step:06d}.gpw.partial"
        local_error = None
        try:
            self._checkpoint_budget_ok()
            require(not final_path.exists() and not temporary_path.exists(),
                    "checkpoint output collision; refusing overwrite")
            # GPAW's new ASECalculator write API is collective. On a rank-0
            # filesystem error, rendezvous once here so peers leave its final
            # write barrier, then rendezvous again below before broadcasting.
            self.calculator.write(str(temporary_path), mode="all")
        except Exception as exc:
            local_error = f"{type(exc).__name__}: {exc}"
            self.world.barrier()
        self.world.barrier()
        errors = self.world.sum_scalar(int(local_error is not None))
        payload = None
        if self.world.rank == 0:
            if errors:
                payload = {"ok": False, "error": local_error or "MPI rank GPW write failed"}
            else:
                try:
                    require(temporary_path.is_file(), "GPAW returned without a GPW file")
                    size = temporary_path.stat().st_size
                    require(size > 0, "GPW file is empty")
                    estimate = int(self.config["checkpoint_space_budget"]["per_checkpoint_estimate_bytes"])
                    require(size <= estimate,
                            f"GPW size {size} exceeds the checkpoint budget estimate {estimate}")
                    actual_sha = sha256_file(temporary_path)
                    self._verify_gp_w(temporary_path)
                    # Test-fixture readers can validate metadata as well as bytes.
                    if self.gpw_reader is not None:
                        fake_readback = self.gpw_reader(temporary_path)
                        require(fake_readback.get("iteration") == step,
                                "synthetic checkpoint iteration mismatch")
                    os.replace(temporary_path, final_path)
                    metadata = {
                        "status": "CHECKPOINT_VALID",
                        "registry_id": "R2",
                        "completed_scf_iteration": step,
                        "checkpoint_sha256": actual_sha,
                        "checkpoint_size_bytes": size,
                        "checkpoint_path": final_path.name,
                        "wavefunctions_included": True,
                        "readback_verified": True,
                        "config_sha256": self.config_sha,
                        "parameter_lock_sha256": self.lock_sha,
                        "input_sha256": self.config["source_bindings"]["input_sha256"],
                        "structure_sha256": self.config["source_bindings"]["structure_sha256"],
                        "scf_iteration_binding": "completed iteration stored in this sidecar and bound to this GPW SHA-256",
                    }
                    sidecar = final_path.with_suffix(".json")
                    write_json_atomic(sidecar, metadata)
                    payload = {"ok": True, "path": str(final_path), "sha256": actual_sha}
                except Exception as exc:
                    payload = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
        from gpaw.mpi import broadcast
        payload = broadcast(payload, 0, comm=self.world)
        if not payload["ok"]:
            raise CheckpointFailure(payload["error"])
        self.saved_steps.add(step)
        self.latest_checkpoint = Path(payload["path"])
        return self.latest_checkpoint

    def on_scf_step(self, ctx: Any) -> None:
        step = int(ctx.niter)
        interval = int(self.config["checkpoint_policy"]["interval_completed_scf_iterations"])
        if step % interval == 0:
            self.save(step)
        elapsed = self.clock() - self.started_at
        wall_limit = int(self.config["method"]["max_wall_clock_seconds"])
        exit_margin = int(self.config["method"]["wall_clock_exit_margin_seconds"])
        safe_stop_after = wall_limit - exit_margin
        signal_requested = self.stop_event.is_set()
        if signal_requested or elapsed >= safe_stop_after:
            if step not in self.saved_steps:
                self.save(step)
            reason = "SIGNAL_SAFE_STOP" if signal_requested else "WALL_CLOCK_SAFE_STOP"
            raise SafeStop(reason, step)


def build_calculator_kwargs(config: dict[str, Any], txt: Path) -> dict[str, Any]:
    from gpaw import FermiDirac, PW
    from gpaw.mixer import Mixer
    method = config["method"]
    return {
        "mode": PW(method["cutoff_eV"]),
        "xc": method["xc"],
        "kpts": tuple(method["kpoint_mesh"]),
        "occupations": FermiDirac(method["occupations"]["width_eV"]),
        "convergence": method["convergence"],
        "maxiter": method["max_scf_iterations"],
        "mixer": Mixer(method["mixer"]["beta"],
                       method["mixer"]["nmaxold"], method["mixer"]["weight"]),
        # The new GPAW ASE interface accepts the parameter mapping; the
        # legacy gpaw.dft.PPCG object is intentionally not used here.
        "eigensolver": {"name": "ppcg", "niter": method["eigensolver"]["niter"]},
        "symmetry": method["symmetry"],
        "poissonsolver": method["poissonsolver"],
        "parallel": method["parallel"],
        "txt": str(txt),
    }


def current_scf_steps(calculator: Any) -> int | None:
    try:
        value = calculator.dft.scf_loop.niter
    except (AttributeError, RuntimeError):
        return None
    return int(value) if value is not None else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true",
                        help="run R2 only; requires a hash-bound C authorization")
    parser.add_argument("--authorization", type=Path)
    parser.add_argument("--authorization-sha256")
    parser.add_argument("--work-root", type=Path,
                        help="new run directory outside the Git repository")
    parser.add_argument("--preflight-path", type=Path, default=Path("/workspace"),
                        help="existing path used for read-only disk/inode checks")
    args = parser.parse_args()
    try:
        config, config_sha = load_config()
        lock, lock_path = verify_parameter_lock(config)
        sources = verify_sources(config)
        if not args.run:
            runtime = verify_runtime(check_mpi=False)
            resources = collect_resources(args.preflight_path.resolve())
            print(json.dumps({
                "status": "PREFLIGHT_PASS" if resources["cpu"]["pass"]
                and resources["memory"]["pass"] and resources["storage"]["pass"] else "HOLD",
                "execution_authorized": False,
                "r2_dft_started": False,
                "config_sha256": config_sha,
                "parameter_lock_sha256": sha256_file(lock_path),
                "sources": sources,
                "runtime": runtime,
                "resources": resources,
            }, indent=2))
            return 0 if resources["cpu"]["pass"] and resources["memory"]["pass"] \
                and resources["storage"]["pass"] else 2

        require(args.authorization is not None and args.authorization_sha256,
                "--run requires a C authorization file and its SHA-256")
        require(args.work_root is not None, "--run requires --work-root")
        auth = validate_authorization(args.authorization, args.authorization_sha256,
                                      config, config_sha)
        runtime = verify_runtime(check_mpi=True)
        parent = args.work_root.expanduser().resolve().parent
        require(parent.is_dir(), "work-root parent must already exist")
        require(not args.work_root.expanduser().resolve().is_relative_to(REPO_ROOT),
                "runtime output must be outside the Git repository")
        resources = collect_resources(parent)
        require(resources["cpu"]["pass"], "CPU quota does not cover four MPI ranks")
        require(resources["memory"]["pass"], "R2 memory preflight failed")
        require(resources["storage"]["pass"], "checkpoint disk budget preflight failed")
        from gpaw.mpi import world
        if world.rank == 0:
            require(not args.work_root.exists(),
                    "work-root already exists; refusing to overwrite or resume")
            args.work_root.mkdir(parents=False, exist_ok=False)
        world.barrier()

        from ase.io import read
        # Use GPAW's new ASE interface: the legacy top-level calculator does
        # not expose the per-SCF-step hook required for verified checkpoints.
        from gpaw.new.ase_interface import GPAW
        atoms = read(sources["structure_path"], format="extxyz")
        positions0 = atoms.positions.copy()
        cell0 = atoms.cell.array.copy()
        pbc0 = atoms.pbc.copy()
        status_path = args.work_root / "RUN_STATUS.json"
        gpaw_log = args.work_root / "gpaw.log"
        stop_event = threading.Event()
        signal.signal(signal.SIGTERM, lambda *_: stop_event.set())
        signal.signal(signal.SIGINT, lambda *_: stop_event.set())
        calc = GPAW(**build_calculator_kwargs(config, gpaw_log))
        manager = CheckpointManager(
            calc, atoms, args.work_root, world, config, config_sha,
            sha256_file(lock_path), stop_event=stop_event)
        calc.hooks["scf_step"] = manager.on_scf_step
        atoms.calc = calc
        start = time.monotonic()
        failure = None
        result = None
        try:
            total_energy = float(atoms.get_potential_energy(force_consistent=False))
            free_energy = float(atoms.get_potential_energy(force_consistent=True))
            forces = atoms.get_forces(apply_constraint=False)
            stress = atoms.get_stress(voigt=False, apply_constraint=False)
            require((atoms.positions == positions0).all(), "atomic coordinates changed")
            require((atoms.cell.array == cell0).all(), "cell changed")
            require((atoms.pbc == pbc0).all(), "PBC changed")
            require(math.isfinite(total_energy) and math.isfinite(free_energy),
                    "energy output is not finite")
            require(all(math.isfinite(float(x)) for x in forces.ravel()),
                    "force output is not finite")
            require(all(math.isfinite(float(x)) for x in stress.ravel()),
                    "stress output is not finite")
            scf_steps = int(calc.dft.scf_loop.niter)
            result = {
                "status": "SUCCESS",
                "registry_id": "R2",
                "role": config["scope"]["role"],
                "total_energy_eV": total_energy,
                "force_consistent_free_energy_eV": free_energy,
                "forces_eV_A": forces.tolist(),
                "stress_eV_A3": stress.tolist(),
                "scf_converged": True,
                "scf_steps": scf_steps,
                "input_sha256": sources["input_sha256"],
                "structure_sha256": sources["structure_sha256"],
                "parameter_lock_sha256": sha256_file(lock_path),
                "config_sha256": config_sha,
                "authorization_sha256": args.authorization_sha256,
                "authorization_id": auth["C_approval"]["approval_id"],
                "geometry_unchanged": True,
                "mpi_ranks": runtime["mpi_world_size"],
                "omp_threads": runtime["omp_threads"],
                "elapsed_seconds": time.monotonic() - start,
            }
        except SafeStop as exc:
            failure = {"run_status": exc.reason, "failure_reason": str(exc),
                       "scf_steps": exc.step, "scf_converged": False}
        except BaseException as exc:
            failure = {"run_status": "FAILED_SCF" if type(exc).__name__ == "KohnShamConvergenceError"
                       else "FAILED_RUNTIME",
                       "failure_reason": f"{type(exc).__name__}: {exc}",
                       "scf_steps": current_scf_steps(calc),
                       "scf_converged": False,
                       "traceback": traceback.format_exc()}
        if world.rank == 0:
            if result is not None:
                write_json_atomic(args.work_root / "RESULTS_RAW.json", result)
                status = {
                    "run_status": "SUCCESS",
                    "scf_converged": True,
                    "scf_steps": result["scf_steps"],
                    "elapsed_seconds": result["elapsed_seconds"],
                    "latest_verified_checkpoint": str(manager.latest_checkpoint) if manager.latest_checkpoint else None,
                    "authorization_sha256": args.authorization_sha256,
                    "authorization_id": auth["C_approval"]["approval_id"],
                }
            else:
                status = {
                    **(failure or {"run_status": "FAILED_RUNTIME",
                                   "failure_reason": "no result or failure receipt"}),
                    "latest_verified_checkpoint": str(manager.latest_checkpoint) if manager.latest_checkpoint else None,
                    "execution_authorized": True,
                    "authorization_sha256": args.authorization_sha256,
                }
            write_json_atomic(status_path, status)
        world.barrier()
        if result is None:
            return 5
        return 0
    except Exception as exc:
        print(f"HOLD: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
