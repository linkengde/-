"""Append-only GitHub Git archive transport with fresh-clone verification.

This module uses the configured HTTPS Git remote and Git credential proxy; it
never handles or stores a token. It does not run scientific calculations.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path, PurePosixPath

MAX_ARCHIVE_BYTES = 50 * 1024 * 1024
_SAFE_COMPONENT = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


class ArchiveBackendError(RuntimeError):
    pass


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_component(value: str, label: str) -> str:
    if not isinstance(value, str) or not _SAFE_COMPONENT.fullmatch(value) or value in {".", ".."}:
        raise ArchiveBackendError(f"UNSAFE_{label.upper()}")
    return value


def _safe_relative(value: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or not path.parts or any(p in {"", ".", ".."} for p in path.parts):
        raise ArchiveBackendError("ARCHIVE_PATH_TRAVERSAL")
    return path


def _read_archive_members(archive_path: str | Path, expected_sha256: str) -> tuple[dict, dict[str, bytes]]:
    path = Path(archive_path)
    if not path.is_file():
        raise ArchiveBackendError("ARCHIVE_NOT_FOUND")
    size = path.stat().st_size
    if size > MAX_ARCHIVE_BYTES:
        raise ArchiveBackendError("ARCHIVE_EXCEEDS_50_MIB")
    actual_sha = sha256_file(path)
    if actual_sha != expected_sha256:
        raise ArchiveBackendError("ARCHIVE_SHA_MISMATCH")
    files: dict[str, bytes] = {}
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        if not members:
            raise ArchiveBackendError("EMPTY_ARCHIVE")
        for member in members:
            rel = _safe_relative(member.name)
            if member.isdir():
                continue
            if not member.isfile():
                raise ArchiveBackendError("NONREGULAR_ARCHIVE_MEMBER")
            if rel.as_posix() in files:
                raise ArchiveBackendError("DUPLICATE_ARCHIVE_PATH")
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ArchiveBackendError("UNREADABLE_ARCHIVE_MEMBER")
            blob = extracted.read()
            if len(blob) != member.size or len(blob) > MAX_ARCHIVE_BYTES:
                raise ArchiveBackendError("ARCHIVE_MEMBER_SIZE_MISMATCH_OR_LIMIT")
            files[rel.as_posix()] = blob
    manifests = [name for name in files if PurePosixPath(name).name == "MANIFEST.json"]
    sums_files = [name for name in files if PurePosixPath(name).name == "SHA256SUMS.txt"]
    if len(manifests) != 1 or len(sums_files) != 1:
        raise ArchiveBackendError("ARCHIVE_MANIFEST_SET_INVALID")
    manifest_name = PurePosixPath(manifests[0])
    root = manifest_name.parent
    manifest = json.loads(files[manifests[0]])
    member_hashes = manifest.get("members_sha256")
    if not isinstance(member_hashes, dict) or not member_hashes:
        raise ArchiveBackendError("MANIFEST_MEMBERS_INVALID")
    expected_files = {str(root / name) for name in member_hashes}
    expected_files |= {str(root / "MANIFEST.json"), str(root / "SHA256SUMS.txt")}
    actual_files = set(files)
    if actual_files != expected_files:
        raise ArchiveBackendError("ARCHIVE_MEMBER_SET_MISMATCH")
    for name, expected in member_hashes.items():
        rel = _safe_relative(name)
        full = str(root / rel)
        if hashlib.sha256(files[full]).hexdigest() != expected:
            raise ArchiveBackendError(f"INNER_MEMBER_SHA_MISMATCH:{name}")
    sums = {}
    for line in files[str(root / "SHA256SUMS.txt")].decode("utf-8").splitlines():
        digest, name = line.split("  ", 1)
        _safe_relative(name)
        if name in sums:
            raise ArchiveBackendError("DUPLICATE_SHA256SUMS_ENTRY")
        sums[name] = digest
    expected_sums = dict(member_hashes)
    expected_sums["MANIFEST.json"] = hashlib.sha256(files[str(root / "MANIFEST.json")]).hexdigest()
    if sums != expected_sums:
        raise ArchiveBackendError("SHA256SUMS_MANIFEST_MISMATCH")
    for key in ("task_id", "owner", "global_step"):
        if key not in manifest:
            raise ArchiveBackendError(f"MANIFEST_IDENTITY_MISSING:{key}")
    return manifest, files


def _verify_restored_directory(root: Path) -> dict:
    manifest_path = root / "MANIFEST.json"
    sums_path = root / "SHA256SUMS.txt"
    if not manifest_path.is_file() or not sums_path.is_file():
        raise ArchiveBackendError("RESTORED_MANIFEST_MISSING")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    member_hashes = manifest.get("members_sha256", {})
    expected_names = set(member_hashes) | {"MANIFEST.json", "SHA256SUMS.txt"}
    actual_names = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}
    if actual_names != expected_names:
        raise ArchiveBackendError("RESTORED_MEMBER_SET_MISMATCH")
    for name, expected in member_hashes.items():
        _safe_relative(name)
        if sha256_file(root / name) != expected:
            raise ArchiveBackendError(f"RESTORED_MEMBER_SHA_MISMATCH:{name}")
    expected_sums = dict(member_hashes)
    expected_sums["MANIFEST.json"] = sha256_file(manifest_path)
    actual_sums = {}
    for line in sums_path.read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        actual_sums[name] = digest
    if actual_sums != expected_sums:
        raise ArchiveBackendError("RESTORED_SHA256SUMS_MISMATCH")
    return manifest


class GitHubArchiveBackend:
    """GitHub branch archive API: upload(), verify(), restore().

    Every upload is a new hash-addressed commit on one configured branch.
    The expected branch head is pinned and updated only after a successful
    normal fast-forward push. Verification and restoration each use a new
    independent clone. Existing remote objects are never overwritten.
    """

    def __init__(self, remote_url: str, branch: str, expected_head: str,
                 archive_prefix: str, *, author_name: str = "Codex Cloud D",
                 author_email: str = "codex-cloud-d@users.noreply.github.com",
                 git_timeout_s: int = 120, allow_local_test_remote: bool = False):
        valid_remote = remote_url.startswith("https://github.com/")
        local_test_remote = allow_local_test_remote and remote_url.startswith("file://")
        if not (valid_remote or local_test_remote):
            raise ArchiveBackendError("REMOTE_MUST_BE_GITHUB_HTTPS")
        self.remote_url = remote_url
        self.branch = _safe_component(branch, "branch")
        if not re.fullmatch(r"[0-9a-f]{40}", expected_head):
            raise ArchiveBackendError("EXPECTED_REMOTE_HEAD_INVALID")
        self.expected_head = expected_head
        prefix = PurePosixPath(archive_prefix)
        if prefix.is_absolute() or any(p in {"", ".", ".."} for p in prefix.parts):
            raise ArchiveBackendError("ARCHIVE_PREFIX_UNSAFE")
        self.archive_prefix = prefix.as_posix()
        self.author_name = author_name
        self.author_email = author_email
        self.git_timeout_s = int(git_timeout_s)

    def _run(self, args: list[str], *, cwd: Path | None = None, check: bool = True) -> subprocess.CompletedProcess:
        env = os.environ.copy()
        env["GIT_TERMINAL_PROMPT"] = "0"
        env["GIT_OPTIONAL_LOCKS"] = "0"
        try:
            result = subprocess.run(args, cwd=cwd, env=env, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    timeout=self.git_timeout_s, check=False)
        except subprocess.TimeoutExpired as exc:
            raise ArchiveBackendError("GIT_OPERATION_TIMEOUT") from exc
        if check and result.returncode != 0:
            # Do not echo command arguments or environment; report bounded diagnostics.
            detail = result.stderr.strip().splitlines()[-1:] or ["git command failed"]
            raise ArchiveBackendError("GIT_OPERATION_FAILED:" + detail[0][:240])
        return result

    def _remote_head(self) -> str:
        result = self._run(["git", "ls-remote", self.remote_url, f"refs/heads/{self.branch}"])
        rows = [line.split() for line in result.stdout.splitlines() if line.strip()]
        if len(rows) != 1 or len(rows[0]) != 2:
            raise ArchiveBackendError("REMOTE_BRANCH_HEAD_UNAVAILABLE_OR_AMBIGUOUS")
        return rows[0][0]

    def _clone_branch(self, destination: Path) -> Path:
        if destination.exists():
            raise ArchiveBackendError("CLONE_DESTINATION_ALREADY_EXISTS")
        self._run(["git", "clone", "--quiet", "--depth", "1", "--single-branch",
                   "--no-tags", "--branch", self.branch, self.remote_url, str(destination)])
        return destination

    def upload(self, archive_path: str | Path, expected_sha256: str) -> dict:
        source = Path(archive_path)
        manifest, _ = _read_archive_members(source, expected_sha256)
        task_id = _safe_component(str(manifest["task_id"]), "task_id")
        owner = _safe_component(str(manifest["owner"]), "owner")
        step = int(manifest["global_step"])
        if step < 1000 or step > 2000 or step % 100 != 0:
            raise ArchiveBackendError("ARCHIVE_STEP_OUTSIDE_SIGNED_100_STEP_SEGMENT_BOUNDARIES")
        if self._remote_head() != self.expected_head:
            raise ArchiveBackendError("REMOTE_BRANCH_MOVED_SINCE_BOUND_HEAD")
        with tempfile.TemporaryDirectory(prefix="codex-github-archive-upload-") as tmp:
            clone = self._clone_branch(Path(tmp) / "repo")
            base = self._run(["git", "rev-parse", "HEAD"], cwd=clone).stdout.strip()
            if base != self.expected_head:
                raise ArchiveBackendError("CLONED_BRANCH_HEAD_MISMATCH")
            rel = PurePosixPath(self.archive_prefix) / owner / task_id / f"step-{step:09d}-{expected_sha256}.tar.gz"
            rel = _safe_relative(rel.as_posix())
            archive_dest = clone.joinpath(*rel.parts)
            archive_dest.parent.mkdir(parents=True, exist_ok=True)
            if archive_dest.exists():
                raise ArchiveBackendError("HASHED_ARCHIVE_PATH_ALREADY_EXISTS")
            shutil.copyfile(source, archive_dest)
            sidecar = archive_dest.with_suffix(archive_dest.suffix + ".sha256")
            sidecar.write_text(f"{expected_sha256}  {archive_dest.name}\n", encoding="utf-8")
            receipt_rel = rel.with_suffix(rel.suffix + ".receipt.json")
            receipt_path = clone.joinpath(*receipt_rel.parts)
            receipt_path.write_text(json.dumps({
                "schema_version": 1, "task_id": task_id, "owner": owner,
                "global_step": step, "archive_sha256": expected_sha256,
                "archive_path": rel.as_posix(), "base_head": base,
                "backend": "GitHubArchiveBackend.upload-v1"
            }, sort_keys=True, indent=2) + "\n", encoding="utf-8")
            self._run(["git", "add", "--", rel.as_posix(),
                       rel.as_posix() + ".sha256", receipt_rel.as_posix()], cwd=clone)
            self._run(["git", "-c", f"user.name={self.author_name}",
                       "-c", f"user.email={self.author_email}", "commit", "--quiet",
                       "-m", f"Archive {task_id} step {step} sha256 {expected_sha256[:16]}"], cwd=clone)
            commit = self._run(["git", "rev-parse", "HEAD"], cwd=clone).stdout.strip()
            if self._remote_head() != base:
                raise ArchiveBackendError("REMOTE_BRANCH_MOVED_BEFORE_PUSH")
            self._run(["git", "push", "origin", f"HEAD:refs/heads/{self.branch}"], cwd=clone)
        self.expected_head = commit
        return {"branch": self.branch, "commit_sha": commit, "base_head": base,
                "archive_path": rel.as_posix(), "receipt_path": receipt_rel.as_posix(),
                "archive_sha256": expected_sha256, "task_id": task_id,
                "owner": owner, "global_step": step}

    def verify(self, receipt: dict) -> dict:
        if receipt.get("branch") != self.branch or receipt.get("archive_sha256") is None:
            raise ArchiveBackendError("RECEIPT_BRANCH_OR_SHA_MISMATCH")
        if receipt.get("global_step", 2001) > 2000:
            raise ArchiveBackendError("RECEIPT_STEP_ABOVE_SIGNED_CEILING")
        with tempfile.TemporaryDirectory(prefix="codex-github-archive-verify-") as tmp:
            clone = self._clone_branch(Path(tmp) / "fresh-clone")
            head = self._run(["git", "rev-parse", "HEAD"], cwd=clone).stdout.strip()
            ancestor = self._run(["git", "merge-base", "--is-ancestor",
                                  receipt["commit_sha"], head], cwd=clone, check=False)
            if ancestor.returncode != 0:
                raise ArchiveBackendError("UPLOAD_COMMIT_NOT_IN_FRESH_REMOTE_BRANCH")
            archive_path = clone.joinpath(*_safe_relative(receipt["archive_path"]).parts)
            sidecar_path = clone.joinpath(*_safe_relative(receipt["archive_path"] + ".sha256").parts)
            remote_receipt_path = clone.joinpath(*_safe_relative(receipt["receipt_path"]).parts)
            if not archive_path.is_file() or not sidecar_path.is_file() or not remote_receipt_path.is_file():
                raise ArchiveBackendError("REMOTE_ARCHIVE_OR_RECEIPT_MISSING")
            sidecar_hash = sidecar_path.read_text(encoding="utf-8").split()[0]
            if sidecar_hash != receipt["archive_sha256"]:
                raise ArchiveBackendError("REMOTE_ARCHIVE_SIDECAR_SHA_MISMATCH")
            manifest, _ = _read_archive_members(archive_path, receipt["archive_sha256"])
            remote_receipt = json.loads(remote_receipt_path.read_text(encoding="utf-8"))
            for key in ("task_id", "owner", "global_step", "archive_sha256", "archive_path", "base_head"):
                if key == "global_step":
                    if remote_receipt.get(key) != receipt.get(key):
                        raise ArchiveBackendError("REMOTE_RECEIPT_IDENTITY_MISMATCH")
                elif remote_receipt.get(key) != receipt.get(key):
                    raise ArchiveBackendError("REMOTE_RECEIPT_IDENTITY_MISMATCH")
            if manifest.get("task_id") != receipt["task_id"] or manifest.get("owner") != receipt["owner"] or manifest.get("global_step") != receipt["global_step"]:
                raise ArchiveBackendError("REMOTE_ARCHIVE_CONTENT_IDENTITY_MISMATCH")
            return {"status": "REMOTE_FRESH_FETCH_VERIFIED", "fresh_head": head,
                    "upload_commit": receipt["commit_sha"],
                    "archive_sha256": receipt["archive_sha256"],
                    "member_count": len(manifest["members_sha256"])}

    def restore(self, receipt: dict, destination: str | Path) -> dict:
        self.verify(receipt)
        dest = Path(destination)
        if dest.exists():
            raise ArchiveBackendError("RESTORE_DESTINATION_ALREADY_EXISTS_NO_OVERWRITE")
        dest.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="codex-github-archive-restore-") as tmp:
            clone = self._clone_branch(Path(tmp) / "fresh-clone")
            archive_path = clone.joinpath(*_safe_relative(receipt["archive_path"]).parts)
            _read_archive_members(archive_path, receipt["archive_sha256"])
            try:
                with tarfile.open(archive_path, "r:gz") as archive:
                    members = archive.getmembers()
                    roots = {PurePosixPath(m.name).parts[0] for m in members}
                    if len(roots) != 1:
                        raise ArchiveBackendError("ARCHIVE_MUST_HAVE_ONE_ROOT_DIRECTORY")
                    root_name = next(iter(roots))
                    for member in members:
                        rel = _safe_relative(member.name)
                        if rel.parts[0] != root_name:
                            raise ArchiveBackendError("ARCHIVE_ROOT_MISMATCH")
                        suffix = rel.parts[1:]
                        if not suffix:
                            if not member.isdir():
                                raise ArchiveBackendError("ARCHIVE_ROOT_NOT_DIRECTORY")
                            continue
                        target = dest.joinpath(*suffix)
                        if member.isdir():
                            target.mkdir(parents=True, exist_ok=True)
                            continue
                        if not member.isfile():
                            raise ArchiveBackendError("NONREGULAR_ARCHIVE_MEMBER")
                        target.parent.mkdir(parents=True, exist_ok=True)
                        if target.exists():
                            raise ArchiveBackendError("RESTORE_MEMBER_ALREADY_EXISTS")
                        payload = archive.extractfile(member)
                        if payload is None:
                            raise ArchiveBackendError("RESTORE_MEMBER_UNREADABLE")
                        with target.open("xb") as stream:
                            shutil.copyfileobj(payload, stream)
                manifest = _verify_restored_directory(dest)
            except BaseException:
                if dest.exists():
                    shutil.rmtree(dest)
                raise
        return {"status": "REMOTE_RESTORE_VERIFIED", "restore_path": str(dest),
                "task_id": manifest["task_id"], "owner": manifest["owner"],
                "global_step": manifest["global_step"],
                "archive_sha256": receipt["archive_sha256"],
                "member_count": len(manifest["members_sha256"])}
