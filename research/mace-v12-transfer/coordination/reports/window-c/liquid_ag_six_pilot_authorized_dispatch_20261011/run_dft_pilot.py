#!/usr/bin/env python3
"""Fail-closed, owner-specific six-liquid-Ag GPAW runner; no implicit launch.

--help and --show-template use only the standard library.  A frozen C-approved
config SHA, real QA-approved input, verified host root and explicit --execute
are required before GPAW is imported.  No A/B task registration is imported.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import sys
import time

TASK_OWNERS = {
    "C-D-AGL-1250K-rho097-seed01": "D",
    "C-D-AGL-1250K-rho103-seed02": "D",
    "C-E-AGL-1250K-rho097-seed02": "E",
    "C-E-AGL-1250K-rho100-seed01": "E",
    "C-F-AGL-1250K-rho100-seed02": "F",
    "C-F-AGL-1250K-rho103-seed01": "F",
}
CONTROL_TASK = "C-D-AGL-1250K-rho097-seed01"
METHOD = {
    "xc": "PBE", "cutoff_eV": 500, "smearing_eV": 0.10,
    "convergence": {"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
    "maxiter": 160, "mixer": {"beta": 0.05, "nmaxold": 8, "weight": 100},
    "parallel": {"sl_auto": True}, "poissonsolver": {},
    "charge": 0, "spinpol": False,
    "software": {"gpaw": "26.7.0", "ase": "3.29.0", "gpaw-data": "1.2.1"},
    "energy_key": "GPAW energy (force_consistent=False)",
    "force_consistent_energy_key": "GPAW free_energy (force_consistent=True)",
    "kpoint_convention": "standard Monkhorst-Pack; tuple size; no Gamma centering offset",
}
MIN_MEMORY = 16 * 1024**3
MIN_DISK_GUIDANCE = 30_000_000_000
MIN_RESERVE = 5_000_000_000


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def helper_geometry_sha(atoms):
    import numpy as np
    if "atom_ids" not in atoms.arrays:
        raise ValueError("Input requires immutable atom_ids")
    return json_sha({"numbers": atoms.numbers.tolist(),
                     "positions": np.round(atoms.positions, 12).tolist(),
                     "cell": np.round(atoms.cell.array, 12).tolist(),
                     "pbc": atoms.pbc.tolist(),
                     "atom_ids": atoms.arrays["atom_ids"].astype(int).tolist()})


def template():
    pending = {"path": "NOT_READY", "sha256": "NOT_COMPUTED"}
    return {
        "schema_version": 1, "task_id": CONTROL_TASK, "owner": "D",
        "run_kind": "kpoint_control", "kpts": [1, 1, 1],
        "input": dict(pending, geometry_sha256="NOT_COMPUTED", box_edge_A=None),
        "source_manifest": dict(pending), "qa_manifest": dict(pending),
        "owner_preflight": dict(pending), "work_root": "NOT_VERIFIED",
        "release": {"pilot_authorized": True, "source_pass": False,
                    "structure_pass": False, "provenance_pass": False,
                    "host_storage_pass": False, "C_run_approved": False,
                    "kpoint_pass": False, "disk_exception_approved": False},
        "method": dict(METHOD, paw_setup_sha256={}),
        "checkpoint_budget_bytes": None,
        "memory_budget_bytes": None, "memory_reserve_bytes": 2 * 1024**3,
        "miscellaneous_budget_bytes": 1_000_000_000,
        "reserve_bytes": MIN_RESERVE, "max_elapsed_s": 43_200,
        "pause_after_s": 28_800, "checkpoint_interval_scf": 20,
        "checkpoint_rotation": "retain_current_and_previous_atomic_replace",
        "resume": "NOT_AUTHORIZED; existing run directories are never reused",
    }


def require(condition, message):
    if not condition:
        raise ValueError(message)


def bound_json(entry, base):
    p = Path(entry["path"])
    if not p.is_absolute():
        p = base / p
    p = p.resolve(strict=True)
    require(len(entry["sha256"]) == 64 and sha(p) == entry["sha256"],
            f"Bound manifest SHA mismatch: {p.name}")
    return json.loads(p.read_text()), p


def machine_fingerprint_sha(work_root):
    return json_sha({"hostname": socket.gethostname(),
                     "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
                     "cgroup": Path("/proc/self/cgroup").read_text().strip(),
                     "work_root": str(Path(work_root).resolve()),
                     "work_root_device": os.stat(work_root).st_dev})


def cgroup_resource_snapshot(work_root):
    """Visible cgroup ancestors constrain the host; bound ranks share cpuset."""
    cg = Path("/sys/fs/cgroup")
    relative = next((line[3:].lstrip("/") for line in Path("/proc/self/cgroup")
                     .read_text().splitlines() if line.startswith("0::")), "")
    current = (cg / relative).resolve()
    if not current.exists() or not current.is_relative_to(cg): current = cg
    ancestors = []
    while current == cg or cg in current.parents:
        ancestors.append(current)
        if current == cg: break
        current = current.parent
    quotas, cpus, memory_heads, memory_currents = [], [], [], []
    for directory in ancestors:
        quota, period = (directory / "cpu.max").read_text().split()
        if quota != "max": quotas.append(int(quota) / int(period))
        cpuset = (directory / "cpuset.cpus.effective").read_text().strip()
        if cpuset:
            cpus.append(sum(int(t.split("-")[1]) - int(t.split("-")[0]) + 1
                            if "-" in t else 1 for t in cpuset.split(",")))
        max_memory = (directory / "memory.max").read_text().strip()
        memory_current = int((directory / "memory.current").read_text())
        memory_currents.append(memory_current)
        if max_memory != "max": memory_heads.append(max(0, int(max_memory) - memory_current))
    require(cpus and quotas and memory_heads, "CPU/memory quota metadata is UNVERIFIED")
    cpu = min(cpus + quotas)
    host_available = next(int(x.split()[1]) * 1024 for x in Path("/proc/meminfo")
                          .read_text().splitlines() if x.startswith("MemAvailable:"))
    effective_memory = min([host_available] + memory_heads)
    st = os.statvfs(work_root)
    return {"cpu_quota_effective": cpu, "cgroup_cpuset_count": min(cpus),
            "cpu_quota_visible_ancestors": quotas,
            "memory_available_bytes": effective_memory,
            "memory_current_visible_ancestors": memory_currents,
            "disk_available_bytes": st.f_bavail * st.f_frsize,
            "disk_available_inodes": st.f_favail,
            "filesystem_device": os.stat(work_root).st_dev}


def common_validate(config, args, base):
    require(config.get("schema_version") == 1, "Unsupported run config schema")
    task = config.get("task_id")
    require(task in TASK_OWNERS and TASK_OWNERS[task] == args.owner == config["owner"],
            "This owner may execute only its two authorized pilot IDs")
    rel = config["release"]
    for key in ("pilot_authorized", "source_pass", "structure_pass", "provenance_pass",
                "host_storage_pass", "C_run_approved"):
        require(rel.get(key) is True, f"HOLD: C release {key} has not passed")
    require(config["kpts"] in ([1, 1, 1], [2, 2, 2]),
            "3x3x3/denser meshes require separate approval; prohibited by this package")
    require(config["run_kind"] in ("kpoint_control", "final_label"), "Unknown run kind")
    if config["run_kind"] == "kpoint_control":
        require(args.owner == "D" and task == CONTROL_TASK,
                "Only D's designated frozen snapshot may run the approved mesh controls")
    else:
        require(rel.get("kpoint_pass") is True, "HOLD: final label needs C k-point release")
    method = config["method"]
    require({k: method.get(k) for k in METHOD} == METHOD,
            "Frozen scientific parameters/version/energy conventions changed")
    require(set(method["paw_setup_sha256"]) == {"Ag.PBE.gz", "Ti.PBE.gz", "Si.PBE.gz", "C.PBE.gz"},
            "All four project PAW setup hashes must be frozen")
    require(config["max_elapsed_s"] == 43_200 and config["pause_after_s"] == 28_800,
            "12-hour stop and 8-hour subsequent-job pause are immutable")
    require(config["checkpoint_interval_scf"] == 20, "Checkpoint interval must be frozen at 20 SCF")
    require(config.get("checkpoint_rotation") == "retain_current_and_previous_atomic_replace",
            "C must approve the new-pilot-only two-generation checkpoint rotation")
    require(type(config["checkpoint_budget_bytes"]) is int and config["checkpoint_budget_bytes"] > 0,
            "Actual geometry/mesh checkpoint budget is missing")
    require(type(config["memory_budget_bytes"]) is int
            and config["memory_budget_bytes"] >= MIN_MEMORY, "Memory budget below pilot screen")
    require(config["reserve_bytes"] >= MIN_RESERVE, "At least 5 GB disk reserve required")
    require(config["memory_reserve_bytes"] >= 2 * 1024**3, "At least 2GiB runtime memory headroom required")
    require(config["miscellaneous_budget_bytes"] >= 1_000_000_000,
            "At least 1 GB input/log/compact-output allowance required")
    source, _ = bound_json(config["source_manifest"], base)
    qa, _ = bound_json(config["qa_manifest"], base)
    host, _ = bound_json(config["owner_preflight"], base)
    require(source.get("SOURCE_PASS") is True, "Source manifest has no C-approved SOURCE_PASS")
    require(qa.get("STRUCTURE_PASS") is True and qa.get("PROVENANCE_PASS") is True,
            "Real trajectory QA/provenance has not passed")
    require(qa.get("task_id") == task and qa.get("owner") == args.owner,
            "QA report belongs to another task or owner")
    require(qa.get("input_sha256") == config["input"]["sha256"]
            and qa.get("geometry_sha256") == config["input"]["geometry_sha256"],
            "C-approved QA does not bind this exact input and geometry")
    root = Path(config["work_root"]).resolve(strict=True)
    require(root.is_dir() and not root.is_relative_to(Path("/tmp"))
            and str(root) != "/workspace", "Verified independent persistent task root required")
    require(host.get("owner") == args.owner and Path(host["work_root"]).resolve() == root,
            "Owner preflight/work-root mismatch")
    require(host.get("persistence_verified") is True
            and host.get("independent_resources_verified") is True,
            "Persistence/resource independence remains UNVERIFIED")
    require(host.get("machine_fingerprint_sha256") == machine_fingerprint_sha(root),
            "Current machine differs from C-bound owner preflight; --owner alone is insufficient")
    runtime = Path(host["runtime_path"]).resolve(strict=True)
    require(sha(runtime) == host["runtime_sha256"], "Verified runtime script changed")
    for name, version in METHOD["software"].items():
        require(importlib.metadata.version(name) == version, f"Wrong installed {name} version")
    require(os.environ.get("GPAW_MPI_BACKEND") == "cgpaw", "Source verified runtime first; MPI backend must be cgpaw")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        require(os.environ.get(name) == "1", f"{name} must be 1; no oversubscription")
    setups = Path(os.environ.get("GPAW_SETUP_PATH", "NOT_SET")).resolve(strict=True)
    for name, digest in method["paw_setup_sha256"].items():
        require(len(digest) == 64 and sha(setups / name) == digest, f"PAW SHA mismatch: {name}")
    inp = Path(config["input"]["path"])
    inp = inp.resolve(strict=True) if inp.is_absolute() else (base / inp).resolve(strict=True)
    require(sha(inp) == config["input"]["sha256"], "Input file SHA differs from C approval")
    require(not (root / "PAUSE_FURTHER.json").exists(), "Later jobs paused by previous 8-hour/failed job")
    return root, inp


def atomic_json(path, value):
    path = Path(path)
    temp = path.with_name(path.name + f".tmp-{os.getpid()}")
    with temp.open("x") as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write("\n"); f.flush(); os.fsync(f.fileno())
    temp.replace(path)
    fd = os.open(path.parent, os.O_RDONLY)
    try: os.fsync(fd)
    finally: os.close(fd)


def mark_supervisor_failure(args, config, root, exc, started, child=None):
    """Every unsuccessful approved owner attempt pauses subsequent launches."""
    record = {"task_id": config["task_id"], "owner": args.owner,
              "C_config_sha256": args.config_sha256, "status": "PAUSE_FURTHER",
              "elapsed_s": time.monotonic() - started,
              "reason": f"{type(exc).__name__}: {exc}",
              "stage": "supervisor_or_child_failure", "automatic_restart": False,
              "existing_run_logs_and_checkpoints": "PRESERVED",
              "child_exit_code": child.poll() if child is not None else None}
    try:
        atomic_json(root / "PAUSE_FURTHER.json", record)
    except Exception as marker_exc:
        # ENOSPC/permission errors cannot be repaired by deleting checkpoints.
        # Expose failure loudly; caller still exits nonzero and must keep HOLD.
        print(f"CRITICAL: could not persist PAUSE_FURTHER; keep all later jobs HOLD: {marker_exc}",
              file=sys.stderr)


def supervisor(args, config, config_path, root):
    """Bounded watchdog for only the process group created for this exact job."""
    child = None
    started = time.monotonic()
    try:
        stats = cgroup_resource_snapshot(root)
        require(stats["cpu_quota_effective"] >= 4
                and stats["memory_available_bytes"] >= config["memory_budget_bytes"],
                "Prelaunch CPU/memory budget failed")
        needed = 3 * config["checkpoint_budget_bytes"] + config["miscellaneous_budget_bytes"] + config["reserve_bytes"]
        require(stats["disk_available_bytes"] >= needed, "Prelaunch checkpoint triple/reserve budget failed")
        if stats["disk_available_bytes"] < MIN_DISK_GUIDANCE:
            require(config["release"].get("disk_exception_approved") is True,
                    "Below unchanged 30GB screen without C-approved limited-pilot exception")
        launcher = shutil.which("mpiexec")
        require(launcher is not None, "Verified OpenMPI launcher missing")
        info = subprocess.run([launcher, "--version"], check=True, capture_output=True, text=True).stdout
        require("5.0.7" in info and ("Open MPI" in info or "OpenRTE" in info),
                "Expected reviewed OpenMPI5.0.7 launcher")
        env = os.environ.copy()
        env["C_PILOT_WATCHDOG_CONFIG_SHA"] = args.config_sha256
        command = [launcher, "--host", "localhost", "--bind-to", "core", "--map-by", "core", "-n", "4",
                   sys.executable, str(Path(__file__).resolve()), "--config", str(config_path),
                   "--config-sha256", args.config_sha256, "--owner", args.owner,
                   "--execute", "--mpi-child"]
        child = subprocess.Popen(command, env=env, start_new_session=True)
        # child.pid is this supervisor's new session/group, never A/B's group.
        paused = False
        while child.poll() is None:
            elapsed_s = time.monotonic() - started
            if elapsed_s >= config["pause_after_s"] and not paused:
                atomic_json(root / "PAUSE_FURTHER.json", {
                    "task_id": config["task_id"], "owner": args.owner,
                    "C_config_sha256": args.config_sha256, "status": "PAUSE_FURTHER",
                    "elapsed_s": elapsed_s, "reason": "8h: hold subsequent jobs, current job retains 12h limit"})
                paused = True
            if elapsed_s >= config["max_elapsed_s"]:
                raise ValueError("Approved 12h science budget exhausted; last promoted checkpoint retained")
            try: child.wait(timeout=min(20, max(0.01, config["max_elapsed_s"] - elapsed_s)))
            except subprocess.TimeoutExpired: pass
        require(child.returncode == 0, f"This pilot MPI job failed (exit {child.returncode}); no automatic restart")
    except BaseException as exc:
        # Includes Popen/version/resource errors, immediate rank/geometry/env
        # failures and nonzero child exits, not only failures inside SCF.
        mark_supervisor_failure(args, config, root, exc, started, child)
        if child is not None and child.poll() is None:
            try: os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError: pass
            try: child.wait(timeout=120)
            except subprocess.TimeoutExpired:
                try: os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError: pass
                child.wait()
        raise


def mark_early_validated_owner_failure(args, config, base, exc):
    """Pause only a C-bound real owner root; false release templates never write."""
    needed = ("pilot_authorized", "source_pass", "structure_pass", "provenance_pass",
              "host_storage_pass", "C_run_approved")
    if not (args.launch and args.execute
            and TASK_OWNERS.get(config.get("task_id")) == config.get("owner") == args.owner
            and all(config.get("release", {}).get(k) is True for k in needed)):
        return
    try:
        host, _ = bound_json(config["owner_preflight"], base)
        root = Path(config["work_root"]).resolve(strict=True)
        require(root.is_dir() and str(root) != "/workspace" and not root.is_relative_to(Path("/tmp")),
                "No verified isolated owner root")
        require(host.get("owner") == args.owner and Path(host["work_root"]).resolve() == root
                and host.get("persistence_verified") is True
                and host.get("independent_resources_verified") is True
                and host.get("machine_fingerprint_sha256") == machine_fingerprint_sha(root),
                "Early failure root ownership is unverified")
        mark_supervisor_failure(args, config, root, exc, time.monotonic())
    except Exception as root_exc:
        print(f"No pause marker written outside a verified owner root; all jobs stay HOLD: {root_exc}",
              file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--show-template", action="store_true")
    ap.add_argument("--config", type=Path)
    ap.add_argument("--config-sha256")
    ap.add_argument("--owner", choices=("D", "E", "F"))
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--launch", action="store_true", help="Start this owner job in a new supervised process group")
    ap.add_argument("--mpi-child", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.show_template:
        print(json.dumps(template(), indent=2)); return 0
    require(args.config is not None and args.config_sha256 and args.owner,
            "--config, --config-sha256 and --owner are required")
    config_path = args.config.resolve(strict=True)
    require(sha(config_path) == args.config_sha256, "Frozen C run-config SHA mismatch")
    config = json.loads(config_path.read_text())
    try:
        root, inp = common_validate(config, args, config_path.parent)
    except BaseException as exc:
        mark_early_validated_owner_failure(args, config, config_path.parent, exc)
        raise
    require(args.execute, "Checks only; explicit --execute is required for DFT")
    if args.launch:
        require(not args.mpi_child and "OMPI_COMM_WORLD_SIZE" not in os.environ,
                "Supervisor must start outside MPI")
        supervisor(args, config, config_path, root)
        return 0
    require(args.mpi_child and os.environ.get("C_PILOT_WATCHDOG_CONFIG_SHA") == args.config_sha256,
            "Use --launch --execute; unsupervised direct MPI launch forbidden")
    require(int(os.environ.get("OMPI_COMM_WORLD_SIZE", "1")) == 4,
            "Launch exactly four MPI ranks from the verified owner runtime")
    # All scientific imports are below C-bound source/structure/storage gates.
    from gpaw.mpi import world, broadcast
    require(world.size == 4, "GPAW MPI backend did not initialize four ranks")
    import numpy as np
    from ase.io import read
    from gpaw import GPAW, PW
    from gpaw.mixer import Mixer
    from gpaw.occupations import FermiDirac
    import fcntl

    def root_call(callback):
        status = None
        if world.rank == 0:
            try: status = {"ok": True, "value": callback()}
            except Exception as exc: status = {"ok": False, "reason": f"{type(exc).__name__}: {exc}"}
        status = broadcast(status, comm=world)
        require(status["ok"], status.get("reason", "Rank-zero operation failed"))
        return status.get("value")

    atoms = read(inp)
    require(len(atoms) == 64 and set(atoms.get_chemical_symbols()) == {"Ag"}, "Input must be exactly 64 Ag")
    require(atoms.pbc.tolist() == [True, True, True], "Liquid bulk PBC must be TTT")
    require(np.isfinite(atoms.positions).all() and np.isfinite(atoms.cell.array).all(), "Nonfinite geometry")
    require(atoms.arrays["atom_ids"].dtype.kind in "iu" and len(set(atoms.arrays["atom_ids"].tolist())) == 64,
            "Immutable atom_ids must be unique integers")
    edge = config["input"]["box_edge_A"]
    require(isinstance(edge, (int, float)) and edge > 0
            and np.allclose(atoms.cell.array, np.eye(3) * edge, rtol=0, atol=1e-9),
            "Frozen density-derived cubic box mismatch")
    geom_sha = helper_geometry_sha(atoms)
    require(geom_sha == config["input"]["geometry_sha256"], "Geometry/order/IDs SHA mismatch")
    forbidden_info = {"energy", "free_energy", "REF_energy", "PW_PBE_energy_eV"}
    require(not forbidden_info.intersection(atoms.info)
            and not {"forces", "REF_forces", "PW_PBE_forces"}.intersection(atoms.arrays)
            and atoms.calc is None, "Candidate input contains a preexisting calculator or force/energy label")
    frozen_geometry = atoms.copy()
    mesh = "x".join(str(k) for k in config["kpts"])
    run_id = f"{config['run_kind']}-k{mesh}-{args.config_sha256[:12]}"
    task_dir = root / config["task_id"]
    out = task_dir / run_id
    lock_handle = None

    def launch_directory():
        nonlocal lock_handle
        lock_handle = (root / ".pilot.lock").open("a")
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        require(not out.exists(), "Existing run preserved; automatic or unapproved resume forbidden")
        if task_dir.exists():
            for previous_result in task_dir.glob("*/output/result.json"):
                require(previous_result.resolve().is_relative_to(task_dir.resolve()),
                        "Result path leaves this owner's approved task directory")
                old = json.loads(previous_result.read_text())
                duplicate = (old.get("scf_converged") is True
                             and old.get("geometry_sha256") == geom_sha
                             and old.get("kpts") == config["kpts"]
                             and old.get("method_without_kpoints_sha256") == json_sha(config["method"]))
                require(not duplicate,
                        "Completed same-geometry/method/mesh result exists: preserve and ask C to reuse; no duplicate SCF")
        task_dir.mkdir(exist_ok=True)
        out.mkdir(exist_ok=False)
        for name in ("input", "checkpoint", "logs", "output"):
            (out / name).mkdir()
        shutil.copy2(inp, out / "input" / inp.name)
        shutil.copy2(config_path, out / "input" / "C_run_config.json")
        return str(out)

    def disk_and_memory_gate(additional_checkpoint=True, initial=False):
        stats = cgroup_resource_snapshot(root)
        require(stats["cpu_quota_effective"] >= 4, "Effective CPU quota below four ranks")
        memory_target = config["memory_budget_bytes"] if initial else config["memory_reserve_bytes"]
        require(stats["memory_available_bytes"] >= memory_target,
                "Effective available cgroup memory below approved headroom")
        require(stats["disk_available_inodes"] >= 100, "Insufficient filesystem inodes")
        if stats["disk_available_bytes"] < MIN_DISK_GUIDANCE:
            require(config["release"].get("disk_exception_approved") is True,
                    "Below unchanged 30GB screen without C-approved limited-pilot exception")
        # free bytes already exclude all retained files on this host; account for
        # up to three additional generations rather than double-counting retention.
        active_retained = sum(p.stat().st_size for p in
                              (out / "checkpoint" / "state.gpw", out / "checkpoint" / "state.previous.gpw")
                              if p.exists())
        needed = (max(0, 3 * config["checkpoint_budget_bytes"] - active_retained)
                  if additional_checkpoint else 0)
        needed += config["miscellaneous_budget_bytes"] + config["reserve_bytes"]
        require(stats["disk_available_bytes"] >= needed, "Checkpoint triple plus reserve does not fit")
        return stats

    initial_stats = root_call(lambda: disk_and_memory_gate(initial=True))
    root_call(launch_directory)
    started = time.monotonic()
    progress_path = out / "output" / "progress.json"
    pause_path = root / "PAUSE_FURTHER.json"
    checkpoint_dir = out / "checkpoint"
    converged = {"value": False}
    dims = {}
    elapsed = lambda: time.monotonic() - started

    def update_progress(iteration, status="running", reason=None):
        record = {"schema_version": 1, "task_id": config["task_id"], "owner": args.owner,
                  "run_id": run_id, "status": status, "scf_iteration": iteration,
                  "elapsed_s": elapsed(), "input_sha256": config["input"]["sha256"],
                  "geometry_sha256": geom_sha, "C_config_sha256": args.config_sha256,
                  "mpi_ranks": world.size}
        if reason: record["reason"] = reason
        atomic_json(progress_path, record)
        if elapsed() >= config["pause_after_s"] or status in ("failed", "stopped"):
            atomic_json(pause_path, dict(record, status="PAUSE_FURTHER"))
        return record

    calc = GPAW(mode=PW(500), xc="PBE", kpts=tuple(config["kpts"]),
                occupations=FermiDirac(0.10), mixer=Mixer(0.05, 8, 100),
                convergence=dict(METHOD["convergence"]), maxiter=160,
                parallel={"sl_auto": True}, poissonsolver={}, charge=0,
                spinpol=False, txt=str(out / "logs" / "gpaw.log"))
    require(type(calc).__module__ == "gpaw.new.ase_interface", "New GPAW API required")
    atoms.calc = calc

    def atomic_checkpoint(iteration, final=False):
        root_call(disk_and_memory_gate)
        temp = checkpoint_dir / f".state-iter{iteration}-{'final' if final else 'scf'}.gpw.tmp"
        root_call(lambda: require(not temp.exists(), "Existing temporary checkpoint preserved"))
        # This must execute on every rank: write_gpw gathers arrays and barriers.
        calc.write(temp, mode="all", precision="double")
        def promote():
            size = temp.stat().st_size
            require(size <= config["checkpoint_budget_bytes"], "Checkpoint exceeds C-approved estimate; preserve temp")
            with temp.open("rb") as f: os.fsync(f.fileno())
            from ase.io.ulm import Reader
            with Reader(temp) as reader:
                require(reader.get_tag() == "gpaw" and reader.version == 7,
                        "Incomplete or unexpected GPAW ULM checkpoint")
            digest = sha(temp)
            current = checkpoint_dir / "state.gpw"
            previous = checkpoint_dir / "state.previous.gpw"
            # A third checkpoint generation is never silently deleted; changing
            # retention requires an explicit C-approved policy.  A bounded
            # two-generation rotation is explicitly frozen in each run config.
            require(config.get("checkpoint_rotation") == "retain_current_and_previous_atomic_replace",
                    "No C-approved bounded checkpoint retention policy")
            if current.exists():
                current.replace(previous)
                old_meta = checkpoint_dir / "state.gpw.json"
                if old_meta.exists(): old_meta.replace(checkpoint_dir / "state.previous.gpw.json")
            temp.replace(current)
            atomic_json(checkpoint_dir / "state.gpw.json", {
                "task_id": config["task_id"], "run_id": run_id, "iteration": iteration,
                "final": final, "sha256": digest, "bytes": size,
                "input_sha256": config["input"]["sha256"], "geometry_sha256": geom_sha,
                "C_config_sha256": args.config_sha256, "method": config["method"], "kpts": config["kpts"],
                "policy": "new pilot only; retain current+previous; no other task checkpoint touched"})
            return {"bytes": size, "sha256": digest}
        return root_call(promote)

    def collect_dimensions():
        ibzwfs = calc.dft.ibzwfs
        shape = ibzwfs.get_max_shape(global_shape=True)  # collective
        return {"nbands": int(ibzwfs.nbands), "ncomponents": int(ibzwfs.ncomponents),
                "nspins": int(ibzwfs.nspins), "NIBZ": len(ibzwfs.ibz),
                "wavefunction_shape_max": [int(x) for x in shape],
                "wavefunction_dtype": str(np.dtype(ibzwfs.dtype)),
                "wavefunction_storage_dtype": str(ibzwfs._wfs_u[0].psit_nX.data.dtype),
                "density_FFT_shape": [int(x) for x in calc.dft.density.nt_sR.desc.global_shape()],
                "projector_count_per_band": sum(int(setup.ni) for setup in calc.dft.setups),
                "BZ_kpoints": ibzwfs.ibz.bz.kpt_Kc.tolist(),
                "IBZ_kpoints": ibzwfs.ibz.kpt_kc.tolist(),
                "IBZ_weights": ibzwfs.ibz.weight_k.tolist(),
                "expected_default_bands_estimate": 657,
                "nbands_explicitly_overridden": False}

    def safe_hook(ctx):
        iteration = int(ctx.niter)
        if not dims:
            dims.update(collect_dimensions())
            wave_bytes = (math.prod(dims["wavefunction_shape_max"]) * dims["nbands"] * dims["NIBZ"]
                          * dims["nspins"] * np.dtype(dims["wavefunction_storage_dtype"]).itemsize)
            dims["estimated_wavefunction_component_bytes"] = int(wave_bytes)
            root_call(lambda: atomic_json(out / "output" / "actual_dimensions.json", dims))
            root_call(lambda: require(config["checkpoint_budget_bytes"] >= 2 * wave_bytes + 1_000_000_000,
                                     "Actual initialization exceeds conservative checkpoint budget; stop before first write"))
        def checks():
            disk_and_memory_gate(False)
            if (out / "STOP_REQUESTED").exists() or elapsed() >= config["max_elapsed_s"]:
                update_progress(iteration, "stopped", "Stop requested or approved 12h budget exhausted")
                raise ValueError("Safe SCF stop: preserve last complete checkpoint")
            if iteration % 10 == 0 or elapsed() >= config["pause_after_s"]:
                update_progress(iteration)
        root_call(checks)
        if iteration % config["checkpoint_interval_scf"] == 0:
            atomic_checkpoint(iteration)

    calc.hooks["scf_step"] = safe_hook
    calc.hooks["converged"] = lambda: converged.__setitem__("value", True)
    root_call(lambda: update_progress(0))
    try:
        native = float(atoms.get_potential_energy(force_consistent=False))
        free = float(atoms.get_potential_energy(force_consistent=True))
        forces = np.asarray(atoms.get_forces(), dtype=float)
        require(converged["value"], "SCF did not report convergence")
        require(math.isfinite(native) and math.isfinite(free)
                and forces.shape == (64, 3) and np.isfinite(forces).all(), "Nonfinite or incomplete DFT result")
        require(helper_geometry_sha(atoms) == geom_sha
                and np.array_equal(atoms.positions, frozen_geometry.positions)
                and np.array_equal(atoms.cell.array, frozen_geometry.cell.array), "Fixed geometry changed")
        iteration = int(calc.get_number_of_iterations())
        require(iteration <= 160 and elapsed() < config["max_elapsed_s"], "Approved runtime/SCF budget exceeded")
        checkpoint = atomic_checkpoint(iteration, final=True)
        result = {"schema_version": 1, "task_id": config["task_id"], "owner": args.owner,
                  "run_id": run_id, "run_kind": config["run_kind"], "kpts": config["kpts"],
                  "input_sha256": config["input"]["sha256"], "geometry_sha256": geom_sha,
                  "C_config_sha256": args.config_sha256, "method": config["method"],
                  "method_without_kpoints_sha256": json_sha(config["method"]),
                  "atom_count": 64, "numbers": atoms.numbers.tolist(),
                  "atom_ids": atoms.arrays["atom_ids"].tolist(), "positions_A": atoms.positions.tolist(),
                  "cell_A": atoms.cell.array.tolist(), "pbc": atoms.pbc.tolist(),
                  "native_energy_eV": native, "free_energy_eV": free, "forces_eV_A": forces.tolist(),
                  "energy_key": METHOD["energy_key"],
                  "force_consistent_energy_key": METHOD["force_consistent_energy_key"],
                  "force_vector_rms_eV_A": float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
                  "max_atom_force_eV_A": float(np.linalg.norm(forces, axis=1).max()),
                  "scf_converged": True, "scf_iterations": iteration, "mpi_ranks": 4,
                  "elapsed_s": elapsed(), "actual_dimensions": dims, "initial_resources": initial_stats,
                  "runtime_versions": {"python": sys.version, "numpy": importlib.metadata.version("numpy"),
                                       "scipy": importlib.metadata.version("scipy")},
                  "checkpoint": checkpoint, "dataset_role": "liquid_Ag_pilot_candidate_pending_C_acceptance",
                  "PILOT_ACCEPTED": False, "scientific_validation": "static numerical pilot only; not phase-change validation"}
        root_call(lambda: atomic_json(out / "output" / "result.json", result))
        root_call(lambda: update_progress(iteration, "complete"))
        root_call(lambda: print(json.dumps({"status": "COMPLETE_PENDING_C_REVIEW", "run_dir": str(out)}, indent=2)))
    except Exception as exc:
        # Root-only IO failures have already been broadcast.  Collective GPAW
        # failures may strand a rank; the owner launcher MUST enforce the outer
        # job-group timeout, preserving the last promoted checkpoint.
        if world.rank == 0:
            try: update_progress(int(getattr(calc.dft.scf_loop, "niter", 0)), "failed", f"{type(exc).__name__}: {exc}")
            except Exception: pass
        raise
    finally:
        if world.rank == 0 and lock_handle is not None:
            lock_handle.close()
    return 0


if __name__ == "__main__":
    try: raise SystemExit(main())
    except Exception as exc:
        print(f"HOLD/STOP: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(2)
