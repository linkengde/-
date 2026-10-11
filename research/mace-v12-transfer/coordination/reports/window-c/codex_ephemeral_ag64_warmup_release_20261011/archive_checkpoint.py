#!/usr/bin/env python3
"""Archive one completed 1000-step checkpoint to its owner's GitHub branch and read it back."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tarfile
import tempfile
from datetime import datetime, timezone

from codex_release_gate import GateRefusal, MAX_GITHUB_FILE_BYTES, sha256_file


class ArchiveHold(RuntimeError):
    pass


def _run(args: list[str], cwd: Path | None = None) -> str:
    result = subprocess.run(args, cwd=cwd, check=False, capture_output=True, text=True,
                            timeout=120, env={**os.environ, "GIT_TERMINAL_PROMPT": "0"})
    if result.returncode != 0:
        # Avoid emitting remotes or environment credentials in exceptions.
        raise ArchiveHold("GIT_OPERATION_FAILED:" + args[1])
    return result.stdout.strip()


def _regular_payload_files(task_dir: Path, cap: int) -> list[Path]:
    if not task_dir.is_dir():
        raise ArchiveHold("TASK_OUTPUT_DIRECTORY_MISSING")
    files = []
    for path in sorted(task_dir.rglob("*")):
        if path.is_symlink():
            raise ArchiveHold("ARCHIVE_SYMLINK_FORBIDDEN")
        if path.is_dir():
            continue
        if path.name.endswith(".lock") or path.name == "HOST_SAMPLER_RUNNING.lock":
            continue
        if path.name.lower() == "ag.eam.alloy" or path.suffix.lower() == ".pdf":
            raise ArchiveHold("FORBIDDEN_POTENTIAL_OR_DOCUMENT_IN_TASK_OUTPUT")
        if not path.is_file():
            raise ArchiveHold("NON_REGULAR_ARCHIVE_MEMBER")
        size = path.stat().st_size
        if size > cap:
            raise ArchiveHold("SINGLE_FILE_EXCEEDS_50MIB_GITHUB_CAP")
        files.append(path)
    required = {
        "source_manifest.json", "sampler_state.json", "trajectory.traj", "thermal.csv",
    }
    rels = {path.relative_to(task_dir).as_posix() for path in files}
    allowed_top = {"source_manifest.json", "sampler_state.json", "trajectory.traj", "thermal.csv",
                   "stage_wall_budget.json"}
    for rel in rels:
        top = PurePosixPath(rel).parts[0]
        if top not in allowed_top and top not in {"restart", "logs", "resume_journal"}:
            raise ArchiveHold("UNEXPECTED_TASK_FILE_NOT_APPROVED_FOR_ARCHIVE:" + top)
    if not required.issubset(rels):
        raise ArchiveHold("REQUIRED_SAMPLING_STATE_MISSING")
    state = json.loads((task_dir / "sampler_state.json").read_text())
    if state.get("step") != 1000 or state.get("status") != "WARMUP_COMPLETE_C_REVIEW_REQUIRED":
        raise ArchiveHold("ARCHIVE_REQUIRES_COMPLETED_1000_STEP_WARMUP")
    bundle = task_dir / state.get("restart_bundle", "")
    bundle_files = {p.relative_to(task_dir).as_posix() for p in files}
    for name in ("state.json", "state_atoms.traj", "state_arrays.npz", "trajectory.traj", "thermal.csv"):
        if (bundle / name).relative_to(task_dir).as_posix() not in bundle_files:
            raise ArchiveHold("RESTART_BUNDLE_MEMBER_MISSING:" + name)
    if not any(p.relative_to(task_dir).as_posix().startswith("logs/") and "performance" in p.name for p in files):
        raise ArchiveHold("WARMUP_PERFORMANCE_LOG_MISSING")
    if (not state.get("rng_state") or not state.get("C_release_config_sha256_at_checkpoint")
            or not state.get("task_matrix_sha256") or not state.get("source_manifest_sha256")):
        raise ArchiveHold("RNG_OR_CONFIG_OR_SOURCE_IDENTITY_MISSING")
    source_manifest = task_dir / "source_manifest.json"
    if sha256_file(source_manifest) != state["source_manifest_sha256"]:
        raise ArchiveHold("SOURCE_MANIFEST_CHECKPOINT_SHA_MISMATCH")
    expected_bundle = {"state_atoms.traj", "state_arrays.npz", "trajectory.traj", "thermal.csv"}
    if set(state.get("bundle_sha256", {})) != expected_bundle:
        raise ArchiveHold("CHECKPOINT_BUNDLE_SHA_SET_MISMATCH")
    for member_name, expected_sha in state["bundle_sha256"].items():
        if sha256_file(bundle / member_name) != expected_sha:
            raise ArchiveHold("CHECKPOINT_BUNDLE_SHA_MISMATCH:" + member_name)
    bundle_state = json.loads((bundle / "state.json").read_text())
    if bundle_state != state:
        raise ArchiveHold("CURRENT_AND_BUNDLE_CHECKPOINT_STATE_DIFFER")
    with __import__("numpy").load(bundle / "state_arrays.npz", allow_pickle=False) as arrays:
        required_arrays = {"numbers", "positions", "velocities_A_per_fs", "momenta", "masses", "cell", "pbc", "atom_ids"}
        if not required_arrays.issubset(set(arrays.files)):
            raise ArchiveHold("RESTART_ARRAYS_DO_NOT_CONTAIN_FULL_ATOMIC_STATE")
    return files


def _member_hashes(task_dir: Path, files: list[Path]) -> list[dict]:
    return [{"path": p.relative_to(task_dir).as_posix(), "size_bytes": p.stat().st_size,
             "sha256": sha256_file(p)} for p in files]


def _safe_extract(archive: Path, destination: Path) -> None:
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            rel = PurePosixPath(member.name)
            if rel.is_absolute() or ".." in rel.parts or not (member.isfile() or member.isdir()):
                raise ArchiveHold("UNSAFE_TAR_MEMBER")
        tar.extractall(destination, filter="data")


def _readback(repo: Path, remote: str, branch: str, commit: str, rel_archive: str,
              rel_manifest: str, work: Path, expected_task: str, expected_config_sha: str,
              rel_receipt: str | None = None, expected_receipt: dict | None = None) -> dict:
    checkout = work / ("fresh-readback-" + next(tempfile._get_candidate_names()))
    _run(["git", "clone", "--depth", "1", "--branch", branch, "--single-branch", remote, str(checkout)])
    head = _run(["git", "rev-parse", "HEAD"], cwd=checkout)
    if head != commit:
        raise ArchiveHold("REMOTE_READBACK_COMMIT_MISMATCH")
    archive = checkout / rel_archive
    manifest_path = checkout / rel_manifest
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("task_id") != expected_task or manifest.get("config_sha256") != expected_config_sha:
        raise ArchiveHold("REMOTE_ARCHIVE_TASK_OR_CONFIG_IDENTITY_MISMATCH")
    if sha256_file(archive) != manifest.get("archive_sha256"):
        raise ArchiveHold("REMOTE_ARCHIVE_OUTER_SHA_MISMATCH")
    restore = work / ("restored-" + next(tempfile._get_candidate_names()))
    restore.mkdir()
    _safe_extract(archive, restore)
    for entry in manifest["members"]:
        path = restore / entry["path"]
        if not path.is_file() or path.stat().st_size != entry["size_bytes"] or sha256_file(path) != entry["sha256"]:
            raise ArchiveHold("REMOTE_ARCHIVE_MEMBER_SHA_MISMATCH:" + entry["path"])
    actual_files = {p.relative_to(restore).as_posix() for p in restore.rglob("*") if p.is_file()}
    expected_files = {entry["path"] for entry in manifest["members"]}
    if actual_files != expected_files:
        raise ArchiveHold("RESTORED_ARCHIVE_FILE_SET_MISMATCH")
    state = json.loads((restore / "sampler_state.json").read_text())
    if state.get("task_id") != expected_task or state.get("step") != 1000:
        raise ArchiveHold("RESTORED_CHECKPOINT_IDENTITY_MISMATCH")
    if state.get("C_release_config_sha256_at_checkpoint") != expected_config_sha:
        raise ArchiveHold("RESTORED_CHECKPOINT_CONFIG_SHA_MISMATCH")
    receipt_sha = None
    if rel_receipt is not None:
        receipt_path = checkout / rel_receipt
        if not receipt_path.is_file():
            raise ArchiveHold("REMOTE_READBACK_RECEIPT_MISSING")
        actual_receipt = json.loads(receipt_path.read_text())
        if actual_receipt != expected_receipt:
            raise ArchiveHold("REMOTE_READBACK_RECEIPT_CONTENT_MISMATCH")
        receipt_sha = sha256_file(receipt_path)
    return {"readback_commit": head, "archive_sha256": manifest["archive_sha256"],
            "manifest_sha256": sha256_file(manifest_path), "receipt_sha256": receipt_sha,
            "restored_member_count": len(manifest["members"]),
            "restore_path": str(restore), "all_member_sha256_pass": True}


def archive_and_restore(task_dir: Path, staging_root: Path, config: dict, host: dict,
                        task_id: str, config_sha: str, source_code_sha: dict) -> dict:
    policy = config["codex_archive_policy"]
    cap = int(policy["max_single_file_bytes"])
    if cap != MAX_GITHUB_FILE_BYTES or int(policy["max_archive_bytes"]) != MAX_GITHUB_FILE_BYTES:
        raise ArchiveHold("CODEX_ARCHIVE_CAP_NOT_50MIB")
    archive_cfg = host["remote_archive"]
    if archive_cfg.get("synthetic_fixture_verified") is not True:
        raise ArchiveHold("SYNTHETIC_REMOTE_RESTORE_GATE_HOLD")
    staging_root.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="ag64-archive-", dir=staging_root))
    try:
        files = _regular_payload_files(task_dir, cap)
        estimated_bytes = sum(p.stat().st_size for p in files)
        anchor = staging_root
        while not anchor.exists():
            anchor = anchor.parent
        min_free = int(config["sampling"]["min_free_bytes"])
        reserve = int(config["sampling"]["checkpoint_reserve_bytes"])
        if shutil.disk_usage(anchor).free < min_free + reserve + 2 * estimated_bytes + 1024 * 1024:
            raise ArchiveHold("INSUFFICIENT_TEMP_SPACE_FOR_ARCHIVE_AND_RECOVERY_COPY")
        members = _member_hashes(task_dir, files)
        partial_archive = work / "warmup.tar.gz.partial"
        archive_path = work / "warmup.tar.gz"
        with tarfile.open(partial_archive, "w:gz", compresslevel=6) as tar:
            for path in files:
                tar.add(path, arcname=path.relative_to(task_dir).as_posix(), recursive=False)
        with partial_archive.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(partial_archive, archive_path)
        directory_fd = os.open(work, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        archive_bytes = archive_path.stat().st_size
        if archive_bytes > int(policy["max_archive_bytes"]):
            raise ArchiveHold("ARCHIVE_EXCEEDS_50MIB_GITHUB_CAP")
        if shutil.disk_usage(anchor).free < min_free + reserve:
            raise ArchiveHold("TEMP_SPACE_RESERVE_BREACHED_DURING_ARCHIVE")
        sha = sha256_file(archive_path)
        rel_dir = Path("codex_warmup_archive") / task_id / sha[:16]
        rel_archive = (rel_dir / ("warmup-" + sha[:16] + ".tar.gz")).as_posix()
        rel_manifest = (rel_dir / "archive-manifest.json").as_posix()
        manifest = {
            "schema": 1, "task_id": task_id, "owner": host["owner"],
            "step": 1000, "time_fs": 500.0, "config_sha256": config_sha,
            "task_matrix_sha256": config["tasks_sha256"], "archive_sha256": sha,
            "archive_bytes": archive_bytes, "members": members,
            "source_code_sha256": source_code_sha,
            "potential_sha256": config["potential"]["sha256"],
            "potential_file_included": False,
            "unlicensed_documents_included": False,
            "created_utc": datetime.now(timezone.utc).isoformat(),
        }
        manifest_path = work / "archive-manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        if manifest_path.stat().st_size > cap:
            raise ArchiveHold("MANIFEST_EXCEEDS_50MIB_GITHUB_CAP")

        repo = work / "push-repository"
        remote = config["repository"]["url"]
        branch = archive_cfg["branch"]
        _run(["git", "clone", "--depth", "1", "--branch", branch, "--single-branch", remote, str(repo)])
        original_head = _run(["git", "rev-parse", "HEAD"], cwd=repo)
        if original_head != archive_cfg["expected_remote_head"]:
            raise ArchiveHold("REMOTE_BRANCH_HEAD_MOVED_SINCE_SIGNED_PREFLIGHT")
        destination = repo / rel_dir
        destination.mkdir(parents=True, exist_ok=False)
        shutil.copy2(archive_path, destination / Path(rel_archive).name)
        shutil.copy2(manifest_path, destination / "archive-manifest.json")
        _run(["git", "add", rel_archive, rel_manifest], cwd=repo)
        _run(["git", "-c", "user.name=Codex Cloud C archive", "-c", "user.email=codex-cloud-c@users.noreply.github.com",
              "commit", "-m", "C: archive Ag64 1000-step warmup " + task_id], cwd=repo)
        archive_commit = _run(["git", "rev-parse", "HEAD"], cwd=repo)
        _run(["git", "push", "origin", "HEAD:refs/heads/" + branch], cwd=repo)

        readback = _readback(repo, remote, branch, archive_commit, rel_archive, rel_manifest,
                             work, task_id, config_sha)
        receipt = {
            "schema": 1, "task_id": task_id, "owner": host["owner"],
            "config_sha256": config_sha, "archive_sha256": sha,
            "manifest_sha256": readback["manifest_sha256"],
            "archive_commit": archive_commit, "receipt_written_after_fresh_fetch_restore": True,
            "outer_sha_pass": True, "inner_sha_pass": True,
            "restored_member_count": readback["restored_member_count"],
            "source_code_sha256": source_code_sha,
            "created_utc": datetime.now(timezone.utc).isoformat(),
        }
        receipt_rel = (rel_dir / "REMOTE_READBACK_RECEIPT.json").as_posix()
        receipt_path = repo / receipt_rel
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        _run(["git", "add", receipt_rel], cwd=repo)
        _run(["git", "-c", "user.name=Codex Cloud C archive", "-c", "user.email=codex-cloud-c@users.noreply.github.com",
              "commit", "-m", "C: verify fresh remote restore " + task_id], cwd=repo)
        receipt_commit = _run(["git", "rev-parse", "HEAD"], cwd=repo)
        _run(["git", "push", "origin", "HEAD:refs/heads/" + branch], cwd=repo)
        final_readback = _readback(repo, remote, branch, receipt_commit, rel_archive, rel_manifest,
                                   work, task_id, config_sha, receipt_rel, receipt)
        return {**receipt, "receipt_commit": receipt_commit,
                "remote_receipt_sha256": final_readback["receipt_sha256"],
                "remote_branch": branch, "remote_path": rel_dir.as_posix(),
                "final_fresh_fetch_commit": final_readback["readback_commit"],
                "final_restore_member_count": final_readback["restored_member_count"],
                "STORAGE_RELEASE_PROOF": "PASS_FOR_THIS_COMPLETED_WARMUP_ONLY"}
    finally:
        # Archive workspace is a staging copy; keep the signed science output only in owner task root and GitHub.
        shutil.rmtree(work, ignore_errors=True)
