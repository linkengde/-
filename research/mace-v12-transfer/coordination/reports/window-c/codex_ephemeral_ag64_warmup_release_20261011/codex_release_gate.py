#!/usr/bin/env python3
"""Fail-closed storage and owner gate for the separate Codex warmup runner."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import socket
import subprocess

MODE_EPHEMERAL = "CODEX_EPHEMERAL_WITH_VERIFIED_REMOTE_ARCHIVE"
MODE_PERSISTENT = "PERSISTENT_WORK_ROOT"
MAX_GITHUB_FILE_BYTES = 50 * 1024 * 1024
EXPECTED_POTENTIAL_SHA256 = "958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447"
SHA_RE = re.compile(r"^[0-9a-f]{64}$")
GIT_OBJECT_RE = re.compile(r"^(?:[0-9a-f]{40}|[0-9a-f]{64})$")


class GateRefusal(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_signed_config(path: Path, expected_sha256: str) -> tuple[dict, str]:
    if not Path(path).is_file():
        raise GateRefusal("CODEX_CONFIG_MISSING")
    expected = (expected_sha256 or "").lower()
    if not SHA_RE.fullmatch(expected):
        raise GateRefusal("CODEX_CONFIG_SHA_INVALID")
    raw = Path(path).read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != expected:
        raise GateRefusal("CODEX_CONFIG_SHA_MISMATCH")
    try:
        value = json.loads(raw.decode("utf-8"))
    except Exception as exc:
        raise GateRefusal("CODEX_CONFIG_JSON_INVALID") from exc
    if not isinstance(value, dict):
        raise GateRefusal("CODEX_CONFIG_ROOT_NOT_OBJECT")
    return value, actual


def machine_fingerprint(workspace: Path = Path("/workspace")) -> str:
    """Task-session binding only; it does not establish exclusive physical host use."""
    root = Path(workspace).resolve(strict=True)
    try:
        identity = {
            "hostname": socket.gethostname(),
            "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
            "cgroup": Path("/proc/self/cgroup").read_text().strip(),
            "work_root": str(root),
            "work_root_device": os.stat(root).st_dev,
        }
    except OSError as exc:
        raise GateRefusal("NODE_IDENTITY_UNAVAILABLE") from exc
    raw = json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _no_symlink_components(path: Path) -> None:
    for item in (path, *path.parents):
        if item.is_symlink():
            raise GateRefusal("SYMLINK_PATH_FORBIDDEN:" + str(item))


def _overlap(a: Path, b: Path) -> bool:
    a, b = Path(os.path.abspath(a)), Path(os.path.abspath(b))
    return a == b or a in b.parents or b in a.parents


def _inside(path: Path, parent: Path) -> bool:
    path, parent = Path(os.path.abspath(path)), Path(os.path.abspath(parent))
    return path == parent or parent in path.parents


def _checkout_root(script_path: Path) -> Path:
    for parent in Path(script_path).resolve().parents:
        if (parent / ".git").exists():
            return parent
    raise GateRefusal("GIT_CHECKOUT_ROOT_NOT_FOUND")


def storage_mode_precheck(mode: str | None, host: dict) -> str:
    """Legacy persistent mode still requires verified persistence; ephemeral never flips it."""
    if mode in (None, MODE_PERSISTENT):
        if host.get("persistence_verified") is not True:
            raise GateRefusal("PERSISTENT_WORK_ROOT_UNVERIFIED")
        return MODE_PERSISTENT
    if mode != MODE_EPHEMERAL:
        raise GateRefusal("UNSUPPORTED_STORAGE_MODE")
    if host.get("persistence_verified") is not False:
        raise GateRefusal("EPHEMERAL_MODE_MUST_KEEP_PERSISTENCE_FALSE")
    return MODE_EPHEMERAL


def validate_ephemeral_release(
    config: dict,
    *, owner: str,
    task_id: str,
    stage: str,
    repository_root: Path,
    current_fingerprint: str,
    potential_path: Path | None,
    remote_branch_head: str | None,
) -> tuple[dict, Path, Path]:
    """Validate all release gates without making directories or generating atoms."""
    if config.get("storage_mode") != MODE_EPHEMERAL:
        raise GateRefusal("CODEX_EPHEMERAL_MODE_REQUIRED")
    hosts = config.get("hosts")
    if not isinstance(hosts, dict) or owner not in hosts:
        raise GateRefusal("OWNER_NOT_REGISTERED")
    host = hosts[owner]
    storage_mode_precheck(MODE_EPHEMERAL, host)
    if host.get("owner") != owner:
        raise GateRefusal("NODE_OWNER_MISMATCH")
    if stage != "warmup" or config.get("release_scope") != "WARMUP_1000_ONLY":
        raise GateRefusal("WARMUP_SCOPE_ONLY")
    if config.get("config_status") != "C_SIGNED_FOR_WARMUP":
        raise GateRefusal("CODEX_CONFIG_NOT_C_SIGNED")
    if config.get("release_decision") != "APPROVE_WARMUP_ONLY":
        raise GateRefusal("C_RELEASE_DECISION_NOT_APPROVE")
    for key in ("source_release", "density_reference_pass", "storage_release", "warmup_release"):
        if config.get(key) is not True:
            raise GateRefusal(key.upper() + "_HOLD")
    if config.get("continuation_release") is not False:
        raise GateRefusal("CONTINUATION_MUST_REMAIN_HOLD")
    if host.get("C_host_release") is not True:
        raise GateRefusal("OWNER_HOST_RELEASE_HOLD")
    if task_id not in config.get("warmup_task_ids", []):
        raise GateRefusal("TASK_NOT_IN_SIGNED_WARMUP_SCOPE")
    if task_id not in host.get("allowed_task_ids", []):
        raise GateRefusal("TASK_NOT_ALLOWED_FOR_OWNER")
    if host.get("machine_fingerprint_sha256") != current_fingerprint:
        raise GateRefusal("NODE_FINGERPRINT_MISMATCH")
    if not SHA_RE.fullmatch(current_fingerprint or ""):
        raise GateRefusal("NODE_FINGERPRINT_INVALID")

    workspace = Path("/workspace").resolve(strict=True)
    repo = Path(repository_root).resolve(strict=True)
    root_value, stage_value = host.get("sampler_work_root"), host.get("archive_staging_root")
    if not isinstance(root_value, str) or not isinstance(stage_value, str):
        raise GateRefusal("SAMPLER_OR_ARCHIVE_ROOT_NOT_FROZEN")
    root, staging = Path(root_value), Path(stage_value)
    if not root.is_absolute() or not staging.is_absolute():
        raise GateRefusal("SAMPLER_OR_ARCHIVE_ROOT_NOT_ABSOLUTE")
    _no_symlink_components(root)
    _no_symlink_components(staging)
    root, staging = root.resolve(strict=False), staging.resolve(strict=False)
    if not _inside(root, workspace) or _inside(root, repo):
        raise GateRefusal("SAMPLER_ROOT_MUST_BE_IN_WORKSPACE_OUTSIDE_GIT")
    if not _inside(staging, workspace) or _inside(staging, repo):
        raise GateRefusal("ARCHIVE_STAGING_ROOT_MUST_BE_IN_WORKSPACE_OUTSIDE_GIT")
    if _overlap(root, staging):
        raise GateRefusal("SAMPLER_AND_ARCHIVE_ROOT_OVERLAP")

    all_roots = []
    for other_owner, other in hosts.items():
        if not isinstance(other, dict):
            raise GateRefusal("HOST_CONFIG_INVALID")
        for field in ("sampler_work_root", "archive_staging_root"):
            value = other.get(field)
            if value is not None:
                path = Path(value)
                if not path.is_absolute() or not _inside(path, workspace) or _inside(path, repo):
                    raise GateRefusal("OWNER_ROOT_OUTSIDE_ALLOWED_WORKSPACE")
                all_roots.append((other_owner, field, path.resolve(strict=False)))
        dft = other.get("dft_work_root")
        if isinstance(dft, str) and dft:
            dft_path = Path(dft).resolve(strict=False)
            if _overlap(root, dft_path) or _overlap(staging, dft_path):
                raise GateRefusal("CODEX_ROOT_OVERLAPS_DFT_ROOT")
    for index, (owner_a, field_a, path_a) in enumerate(all_roots):
        for owner_b, field_b, path_b in all_roots[index + 1:]:
            if _overlap(path_a, path_b):
                raise GateRefusal(f"OWNER_ROOTS_OVERLAP:{owner_a}:{field_a}:{owner_b}:{field_b}")

    archive = host.get("remote_archive")
    if not isinstance(archive, dict):
        raise GateRefusal("REMOTE_ARCHIVE_CHANNEL_MISSING")
    if archive.get("roundtrip_status") != "PASS" or archive.get("synthetic_fixture_verified") is not True:
        raise GateRefusal("REMOTE_SYNTHETIC_ROUNDTRIP_NOT_VERIFIED")
    if archive.get("repository") != "linkengde/-" or archive.get("repository_visibility") != "PUBLIC":
        raise GateRefusal("REMOTE_REPOSITORY_VISIBILITY_UNVERIFIED")
    if archive.get("scientific_output_upload_authorized") is not True:
        raise GateRefusal("SCIENTIFIC_ARCHIVE_UPLOAD_NOT_AUTHORIZED")
    if archive.get("branch") != f"cloud-{owner.lower()}-b1-dft":
        raise GateRefusal("REMOTE_BRANCH_OWNER_MISMATCH")
    expected_head = archive.get("expected_remote_head")
    if not GIT_OBJECT_RE.fullmatch(expected_head or "") or remote_branch_head != expected_head:
        raise GateRefusal("REMOTE_ARCHIVE_HEAD_CHANGED_OR_UNVERIFIED")
    limits = archive.get("limits", {})
    if limits.get("max_single_file_bytes") != MAX_GITHUB_FILE_BYTES:
        raise GateRefusal("REMOTE_SINGLE_FILE_CAP_MISMATCH")
    if limits.get("max_archive_bytes") != MAX_GITHUB_FILE_BYTES:
        raise GateRefusal("REMOTE_ARCHIVE_CAP_MISMATCH")
    if archive.get("exclude_potential_and_unlicensed_documents") is not True:
        raise GateRefusal("ARCHIVE_EXCLUSION_POLICY_MISSING")

    potential = config.get("potential", {})
    if potential.get("sha256") != EXPECTED_POTENTIAL_SHA256:
        raise GateRefusal("FROZEN_POTENTIAL_SHA_MISMATCH")
    if potential_path is None or not potential_path.is_file():
        raise GateRefusal("POTENTIAL_NOT_AVAILABLE_NO_ALTERNATIVE")
    if potential_path.name != "Ag.eam.alloy":
        raise GateRefusal("POTENTIAL_FILENAME_MISMATCH_NO_ALTERNATIVE")
    if sha256_file(potential_path) != EXPECTED_POTENTIAL_SHA256:
        raise GateRefusal("POTENTIAL_SHA_MISMATCH_NO_ALTERNATIVE")
    if host.get("potential_local_sha_verified") is not True:
        raise GateRefusal("NODE_LOCAL_POTENTIAL_RECEIPT_MISSING")
    receipt = host.get("potential_local_sha_receipt")
    if not isinstance(receipt, dict) or receipt.get("sha256") != EXPECTED_POTENTIAL_SHA256:
        raise GateRefusal("NODE_LOCAL_POTENTIAL_RECEIPT_INVALID")
    return host, root, staging


def remote_branch_head(remote_url: str, branch: str) -> str:
    """Read a branch SHA without printing credentials or changing any ref."""
    result = subprocess.run(
        ["git", "ls-remote", "--heads", remote_url, "refs/heads/" + branch],
        check=False, capture_output=True, text=True, timeout=30,
        env={**os.environ, "GIT_TERMINAL_PROMPT": "0"},
    )
    if result.returncode != 0:
        raise GateRefusal("REMOTE_ARCHIVE_CHANNEL_UNREACHABLE")
    rows = [line.split() for line in result.stdout.splitlines() if line.strip()]
    if len(rows) != 1 or len(rows[0]) < 2:
        raise GateRefusal("REMOTE_ARCHIVE_BRANCH_NOT_FOUND")
    if rows[0][1] != "refs/heads/" + branch or not GIT_OBJECT_RE.fullmatch(rows[0][0]):
        raise GateRefusal("REMOTE_ARCHIVE_BRANCH_RESPONSE_INVALID")
    return rows[0][0]
