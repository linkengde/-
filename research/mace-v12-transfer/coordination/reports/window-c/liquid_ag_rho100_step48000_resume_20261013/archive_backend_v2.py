#!/usr/bin/env python3
"""Immutable per-segment Git archive with fresh-fetch restore verification.

The backend never stages files outside the explicit archive path and refuses
to archive a potential file.  It is intentionally not called by preflight.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tempfile
from typing import Iterable


class ArchiveError(RuntimeError):
    pass


def segment_archive_passed(result: dict | None) -> bool:
    """Return whether a segment is safe to continue past.

    This predicate is shared by the real runner and synthetic gate tests so a
    missing remote upload, fresh fetch, restore, checkpoint restore, or SHA
    match always stops automatic continuation.
    """
    required = (
        "archive_upload_pass",
        "archive_fresh_fetch_pass",
        "archive_restore_pass",
        "checkpoint_restore_pass",
        "sha256_match_pass",
    )
    return isinstance(result, dict) and all(result.get(key) is True for key in required)


FORBIDDEN_NAMES = {"ag.eam.alloy", "ag.eam", "potential.eam", "potential.setfl"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _git(repo: Path, *args: str, check: bool = True) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo), *args],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if check and result.returncode:
        raise ArchiveError(f"git {' '.join(args[:2])} failed ({result.returncode})")
    return result.stdout.strip()


def write_sha256sums(directory: Path, filename: str = "SHA256SUMS.txt") -> Path:
    directory = directory.resolve()
    sum_path = directory / filename
    entries = []
    for path in sorted(directory.rglob("*")):
        if not path.is_file() or path == sum_path:
            continue
        if path.is_symlink():
            raise ArchiveError(f"symlinks are forbidden in an archive: {path.relative_to(directory)}")
        if path.name.casefold() in FORBIDDEN_NAMES:
            raise ArchiveError("potential files are forbidden in the public-safe archive")
        rel = path.relative_to(directory).as_posix()
        entries.append(f"{sha256_file(path)}  {rel}")
    sum_path.write_text("\n".join(entries) + "\n", encoding="ascii")
    return sum_path


def verify_sha256sums(directory: Path, filename: str = "SHA256SUMS.txt") -> dict[str, str]:
    directory = directory.resolve()
    manifest = directory / filename
    if not manifest.is_file():
        raise ArchiveError(f"missing SHA file: {filename}")
    expected: dict[str, str] = {}
    for line in manifest.read_text(encoding="ascii").splitlines():
        if not line.strip():
            continue
        try:
            digest, rel = line.split("  ", 1)
        except ValueError as exc:
            raise ArchiveError(f"malformed SHA line: {line!r}") from exc
        pure = PurePosixPath(rel)
        if pure.is_absolute() or ".." in pure.parts or rel in expected:
            raise ArchiveError(f"unsafe or duplicate archive path: {rel}")
        target = (directory / Path(*pure.parts)).resolve()
        if directory not in target.parents or not target.is_file():
            raise ArchiveError(f"missing or out-of-bundle file: {rel}")
        if target.name.casefold() in FORBIDDEN_NAMES:
            raise ArchiveError("potential files are forbidden in the public-safe archive")
        actual = sha256_file(target)
        if actual != digest:
            raise ArchiveError(f"SHA mismatch for {rel}")
        expected[rel] = digest
    actual_files = {
        p.relative_to(directory).as_posix()
        for p in directory.rglob("*")
        if p.is_file() and p != manifest
    }
    if actual_files != set(expected):
        raise ArchiveError("archive members do not match SHA256SUMS.txt")
    return expected


def _remote_head(repo: Path, remote: str, branch: str) -> str:
    url = _git(repo, "remote", "get-url", remote)
    result = subprocess.run(
        ["git", "ls-remote", url, f"refs/heads/{branch}"],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode:
        raise ArchiveError("could not read remote branch head")
    fields = result.stdout.strip().split()
    return fields[0] if fields else ""


def remote_branch_head(repo: Path, remote: str, branch: str) -> str:
    """Read the live remote branch head without fetching or changing the checkout."""
    return _remote_head(repo.resolve(), remote, branch)


def verify_branch_head(repo: Path, remote: str, branch: str) -> dict:
    local = _git(repo.resolve(), "rev-parse", "HEAD")
    remote_sha = remote_branch_head(repo, remote, branch)
    return {"local_head": local, "remote_head": remote_sha, "pass": bool(local and local == remote_sha)}


def archive_push_fresh_restore(
    *,
    repo_root: Path,
    source_dir: Path,
    archive_relpath: str,
    branch: str,
    remote: str,
    expected_remote_head: str,
    restore_dir: Path,
    commit_message: str,
) -> dict:
    """Commit one immutable archive, push, clone/fetch, restore, compare SHAs."""
    repo = repo_root.resolve()
    source = source_dir.resolve()
    destination_rel = PurePosixPath(archive_relpath)
    if destination_rel.is_absolute() or ".." in destination_rel.parts:
        raise ArchiveError("archive path must stay inside the repository")
    if _git(repo, "branch", "--show-current") != branch:
        raise ArchiveError(f"wrong owner branch; expected {branch}")
    local_head = _git(repo, "rev-parse", "HEAD")
    remote_head = _remote_head(repo, remote, branch)
    if remote_head != expected_remote_head or local_head != expected_remote_head:
        raise ArchiveError("local/remote branch head changed since the signed base")
    destination = repo / Path(*destination_rel.parts)
    if destination.exists():
        raise ArchiveError("immutable archive destination already exists")
    if not source.is_dir():
        raise ArchiveError("source bundle directory does not exist")
    write_sha256sums(source)
    expected = verify_sha256sums(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, destination, copy_function=shutil.copy2)
    try:
        _git(repo, "add", "-f", "--", destination_rel.as_posix())
        _git(repo, "commit", "--only", "-m", commit_message, "--", destination_rel.as_posix())
        commit = _git(repo, "rev-parse", "HEAD")
        pushed = subprocess.run(
            ["git", "-C", str(repo), "push", "--porcelain", remote, f"HEAD:refs/heads/{branch}"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if pushed.returncode:
            raise ArchiveError("push rejected; no force-push was attempted")
        with tempfile.TemporaryDirectory(prefix="ag64-fresh-fetch-") as temp_name:
            clone_dir = Path(temp_name) / "clone"
            clone_dir.mkdir()
            _git(clone_dir, "init", "--quiet")
            _git(clone_dir, "remote", "add", "origin", _git(repo, "remote", "get-url", remote))
            fetch = subprocess.run(
                ["git", "-C", str(clone_dir), "fetch", "--depth=1", "--no-tags", "--filter=blob:none", "origin", f"refs/heads/{branch}"],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if fetch.returncode:
                raise ArchiveError("fresh fetch failed")
            fetched_head = _git(clone_dir, "rev-parse", "FETCH_HEAD")
            if fetched_head != commit:
                raise ArchiveError("fresh fetch did not retrieve the pushed commit")
            _git(clone_dir, "checkout", "--detach", fetched_head)
            restored = clone_dir / Path(*destination_rel.parts)
            fetched_manifest = verify_sha256sums(restored)
            if fetched_manifest != expected:
                raise ArchiveError("fresh-fetch member hashes differ from the local bundle")
            restore = restore_dir.resolve()
            if restore.exists():
                raise ArchiveError("restore destination already exists")
            shutil.copytree(restored, restore, copy_function=shutil.copy2)
            restored_manifest = verify_sha256sums(restore)
            if restored_manifest != expected:
                raise ArchiveError("restored bundle hashes differ from the source bundle")
            checkpoint_restore_pass = True
            if "checkpoint.json" in expected:
                try:
                    source_checkpoint = json.loads((source / "checkpoint.json").read_text(encoding="utf-8"))
                    restored_checkpoint = json.loads((restore / "checkpoint.json").read_text(encoding="utf-8"))
                except Exception as exc:
                    raise ArchiveError(f"restored checkpoint JSON is unreadable: {type(exc).__name__}") from exc
                required_checkpoint_fields = {"global_step", "positions_A", "momenta_amu_A_per_fs", "rng_states"}
                if not required_checkpoint_fields.issubset(restored_checkpoint):
                    raise ArchiveError("restored checkpoint lacks position, momentum, step, or RNG state")
                if source_checkpoint != restored_checkpoint:
                    raise ArchiveError("restored checkpoint payload differs from the archived source payload")
                checkpoint_restore_pass = True
            else:
                checkpoint_restore_pass = None
        return {
            "archive_upload_pass": True,
            "archive_fresh_fetch_pass": True,
            "archive_restore_pass": True,
            "checkpoint_restore_pass": checkpoint_restore_pass,
            "sha256_match_pass": True,
            "commit": commit,
            "member_count": len(expected),
            "manifest_sha256": sha256_file(restore_dir / "SHA256SUMS.txt"),
        }
    except Exception:
        # Leave the immutable local directory/commit available for diagnosis;
        # never rewrite or delete a possibly valid archive.
        raise
