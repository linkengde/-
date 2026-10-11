"""Real backend-v2 upload/clone/verify/restore through a local bare Git remote."""
from __future__ import annotations

import hashlib
import io
import json
import subprocess
import tarfile
import tempfile
from pathlib import Path

import numpy as np

import sys
DEV_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DEV_ROOT))
import github_archive_backend as backend  # noqa: E402

TASK_ID = "C-D-AGL-1250K-rho097-seed01"
SOURCE_STATE_SHA = "5d03864f07a98f177c3c0cbdcb850158ed89a64fa2a6ab7f879e8c7ab003be2a"


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git(args: list[str], *, cwd: Path | None = None) -> str:
    result = subprocess.run(["git", *args], cwd=cwd, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            check=False)
    if result.returncode:
        raise AssertionError(f"git {args[0]} failed: {result.stderr.strip()}")
    return result.stdout.strip()


def _build_transport_archive(path: Path) -> dict:
    root = "synthetic-transport-step-000003000"
    rng_states = {}
    for name, seed in (("packing", 901), ("velocity", 902), ("thermostat", 903)):
        rng_states[name] = np.random.Generator(np.random.PCG64(seed)).bit_generator.state
    arrays_buffer = io.BytesIO()
    np.savez_compressed(
        arrays_buffer,
        positions=np.arange(192, dtype=float).reshape(64, 3),
        velocities=np.zeros((64, 3), dtype=float),
        atom_ids=np.arange(1, 65, dtype=np.int64),
    )
    trajectory_buffer = io.BytesIO()
    np.savez_compressed(trajectory_buffer,
                        steps=np.asarray([2000, 3000], dtype=np.int64),
                        positions=np.zeros((2, 64, 3), dtype=float))
    members = {
        "README.md": b"Synthetic transport fixture only; not an MD result.\n",
        "state.json": json.dumps(
            {"task_id": TASK_ID, "owner": "D", "step": 3000,
             "source_step": 2000, "synthetic_transport_fixture": True},
            sort_keys=True, separators=(",", ":")).encode(),
        "state_arrays.npz": arrays_buffer.getvalue(),
        "rng_state.json": json.dumps(rng_states, sort_keys=True,
                                      separators=(",", ":")).encode(),
        "thermal.csv": b"step,time_fs\n3000,1500.0\n",
        "trajectory.npz": trajectory_buffer.getvalue(),
    }
    manifest = {
        "schema": 1, "task_id": TASK_ID, "owner": "D", "global_step": 3000,
        "members_sha256": {name: _sha(payload) for name, payload in members.items()},
        "potential_file_included": False,
    }
    manifest_bytes = json.dumps(manifest, sort_keys=True,
                                separators=(",", ":")).encode()
    sums = [f"{manifest['members_sha256'][name]}  {name}"
            for name in sorted(members)]
    sums.append(f"{_sha(manifest_bytes)}  MANIFEST.json")
    members["MANIFEST.json"] = manifest_bytes
    members["SHA256SUMS.txt"] = ("\n".join(sums) + "\n").encode()
    with tarfile.open(path, "w:gz") as archive:
        directory = tarfile.TarInfo(root + "/")
        directory.type = tarfile.DIRTYPE
        directory.mode = 0o755
        directory.mtime = 0
        archive.addfile(directory)
        for name, payload in sorted(members.items()):
            info = tarfile.TarInfo(f"{root}/{name}")
            info.size = len(payload)
            info.mode = 0o644
            info.mtime = 0
            archive.addfile(info, io.BytesIO(payload))
    return {"manifest": manifest, "manifest_sha256": _sha(manifest_bytes),
            "rng_states": rng_states,
            "arrays": np.load(io.BytesIO(arrays_buffer.getvalue()))["positions"].copy()}


def run_backend_v2_roundtrip(temp_root: Path | None = None) -> dict:
    own_temp = tempfile.TemporaryDirectory(prefix="ag64-backend-v2-bare-") if temp_root is None else None
    root = Path(own_temp.name) if own_temp is not None else temp_root
    root.mkdir(parents=True, exist_ok=True)
    try:
        bare = root / "archive-remote.git"
        work = root / "seed-work"
        fresh = root / "explicit-fresh-clone"
        archive = root / "synthetic-step3000.tar.gz"
        restored = root / "restored-checkpoint"
        _git(["init", "--quiet", "--initial-branch=main", str(work)])
        (work / "README.md").write_text("isolated backend-v2 test remote\n", encoding="utf-8")
        _git(["-C", str(work), "config", "user.name", "Backend Test"])
        _git(["-C", str(work), "config", "user.email", "backend-test@example.invalid"])
        _git(["-C", str(work), "add", "README.md"])
        _git(["-C", str(work), "commit", "--quiet", "-m", "seed local bare remote"])
        _git(["init", "--bare", "--quiet", "--initial-branch=main", str(bare)])
        _git(["-C", str(work), "remote", "add", "origin", bare.as_uri()])
        _git(["-C", str(work), "push", "--quiet", "origin", "HEAD:refs/heads/main"])
        base_head = _git(["-C", str(work), "rev-parse", "HEAD"])
        fixture = _build_transport_archive(archive)
        archive_sha = backend.sha256_file(archive)

        client = backend.GitHubArchiveBackend(
            bare.as_uri(), "main", base_head, "test-archives",
            allow_local_test_remote=True)
        receipt = client.upload(
            archive, archive_sha, start_step=2000, end_step=3000,
            max_added_steps=1000, source_checkpoint_sha256=SOURCE_STATE_SHA)
        upload_status = "PASS"

        _git(["clone", "--quiet", "--depth", "1", "--single-branch", "--no-tags",
              "--branch", "main", bare.as_uri(), str(fresh)])
        fresh_head = _git(["-C", str(fresh), "rev-parse", "HEAD"])
        remote_receipt_path = fresh / receipt["receipt_path"]
        remote_archive_path = fresh / receipt["archive_path"]
        receipt_bytes = remote_receipt_path.read_bytes()
        remote_receipt = json.loads(receipt_bytes)
        remote_manifest, remote_files = backend._read_archive_members(
            remote_archive_path, archive_sha)
        verify_result = client.verify(receipt)
        verify_status = "PASS"

        restore_result = client.restore(receipt, restored)
        restore_status = "PASS"
        restored_state = json.loads((restored / "state.json").read_text(encoding="utf-8"))
        restored_arrays = np.load(restored / "state_arrays.npz", allow_pickle=False)["positions"]
        restored_rng = json.loads((restored / "rng_state.json").read_text(encoding="utf-8"))

        assert fresh_head == receipt["commit_sha"]
        assert remote_receipt["start_step"] == 2000
        assert remote_receipt["end_step"] == 3000
        assert remote_receipt["source_checkpoint_sha256"] == SOURCE_STATE_SHA
        assert remote_receipt["checkpoint_sha256"] == _sha(remote_files[
            next(name for name in remote_files if name.endswith("/state.json"))])
        assert restored_state["step"] == 3000
        np.testing.assert_array_equal(restored_arrays, fixture["arrays"])
        assert restored_rng == fixture["rng_states"]
        assert backend.sha256_file(restored / "state_arrays.npz") == remote_manifest[
            "members_sha256"]["state_arrays.npz"]
        assert backend.sha256_file(restored / "rng_state.json") == remote_manifest[
            "members_sha256"]["rng_state.json"]

        return {
            "UPLOAD_STATUS": upload_status,
            "VERIFY_STATUS": verify_status,
            "RESTORE_STATUS": restore_status,
            "backend_v2_remote_restore_pass": True,
            "test_remote_type": "local_bare_git_file_remote",
            "production_backend_modified": False,
            "real_backend_methods": ["upload", "verify", "restore"],
            "fresh_clone_status": "PASS",
            "base_head": base_head,
            "fresh_clone_head": fresh_head,
            "upload_commit_sha": receipt["commit_sha"],
            "archive_sha256": archive_sha,
            "manifest_sha256": fixture["manifest_sha256"],
            "receipt_sha256": _sha(receipt_bytes),
            "checkpoint_state_sha256": remote_receipt["checkpoint_sha256"],
            "source_checkpoint_sha256": SOURCE_STATE_SHA,
            "arrays_sha256": backend.sha256_file(restored / "state_arrays.npz"),
            "rng_state_sha256": backend.sha256_file(restored / "rng_state.json"),
            "manifest_step": remote_manifest["global_step"],
            "receipt_segment": {key: remote_receipt[key] for key in (
                "start_step", "end_step", "max_added_steps",
                "source_checkpoint_sha256", "checkpoint_sha256")},
            "verify_result": verify_result,
            "restore_result": {**restore_result,
                               "restore_path": "<temporary-local-test-directory>"},
            "arrays_equal_after_restore": True,
            "rng_equal_after_restore": True,
            "synthetic_transport_fixture_only": True,
            "potential_loaded": False,
            "md_run": False,
        }
    finally:
        if own_temp is not None:
            own_temp.cleanup()


if __name__ == "__main__":
    print(json.dumps(run_backend_v2_roundtrip(), sort_keys=True, indent=2))
