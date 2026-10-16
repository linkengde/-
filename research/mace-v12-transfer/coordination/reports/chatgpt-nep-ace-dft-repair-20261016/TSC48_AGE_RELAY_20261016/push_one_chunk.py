#!/usr/bin/env python3
"""Age-encrypt and push one bounded checkpoint chunk on a dedicated Git branch.

One invocation creates at most one <=32 MiB ciphertext file. It commits and
pushes that file, verifies the immutable remote commit/blob through GitHub API
and a raw HTTPS readback, records the commit in a local state file, then removes
the temporary worktree ciphertext. Git history remains append-only by design.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from urllib.request import Request, urlopen

MIB = 1024 * 1024
GIB = 1024 * MIB
REPO = "linkengde/-"
SCHEMA = "h2-age-gpw-relay-v1"


def run(args: list[str], cwd: Path | None = None, capture: bool = True) -> str:
    p = subprocess.run(args, cwd=cwd, check=True, text=True,
                       stdout=subprocess.PIPE if capture else None,
                       stderr=subprocess.PIPE if capture else None)
    return p.stdout.strip() if capture else ""


def hash_path(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(MIB), b""):
            digest.update(block)
            size += len(block)
    return digest.hexdigest(), size


def atomic_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".writing")
    with temp.open("x", encoding="utf-8") as out:
        json.dump(data, out, sort_keys=True, indent=2)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
    os.replace(temp, path)


def git_object_bytes(gitdir: Path) -> int:
    total = 0
    root = gitdir / "objects"
    for base, _, names in os.walk(root):
        for name in names:
            try:
                total += (Path(base) / name).stat().st_size
            except FileNotFoundError:
                pass
    return total


def verify_repo(repo_dir: Path, branch: str) -> None:
    current = run(["git", "branch", "--show-current"], cwd=repo_dir)
    if current != branch:
        raise RuntimeError("WRONG_TRANSFER_BRANCH: " + current)
    url = run(["git", "remote", "get-url", "origin"], cwd=repo_dir)
    if "github.com/linkengde/-.git" not in url and "github.com/linkengde/-" not in url:
        raise RuntimeError("WRONG_GITHUB_ORIGIN")
    if run(["git", "config", "--bool", "gc.auto"], cwd=repo_dir) != "false":
        raise RuntimeError("GIT_AUTO_GC_MUST_BE_DISABLED")


def init_state(source: Path, state_path: Path, recipient: str,
               expected_sha: str, chunk_bytes: int) -> dict:
    before = source.stat()
    actual_sha, actual_bytes = hash_path(source)
    after = source.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        raise RuntimeError("SOURCE_CHANGED_DURING_HASH")
    if actual_sha != expected_sha:
        raise RuntimeError("SOURCE_SHA_MISMATCH")
    state = {"schema": SCHEMA, "repo": REPO, "source_path": str(source),
             "source_bytes": actual_bytes, "source_sha256": actual_sha,
             "recipient": recipient, "chunk_bytes": chunk_bytes,
             "source_stat": {"dev": before.st_dev, "ino": before.st_ino,
                             "size": before.st_size, "mtime_ns": before.st_mtime_ns},
             "parts": []}
    atomic_json(state_path, state)
    return state


def remote_verify(commit: str, path: str, expected_git_blob: str,
                  expected_sha: str, expected_bytes: int) -> dict:
    api_commit = run(["gh", "api", f"repos/{REPO}/commits/{commit}", "--jq", ".sha"])
    if api_commit != commit:
        raise RuntimeError("REMOTE_COMMIT_MISMATCH")
    jq = f'.tree[] | select(.path == "{path}") | .sha'
    remote_blob = run(["gh", "api", f"repos/{REPO}/git/trees/{commit}?recursive=1",
                       "--jq", jq])
    if remote_blob != expected_git_blob:
        raise RuntimeError("REMOTE_TREE_BLOB_MISMATCH")
    url = f"https://raw.githubusercontent.com/{REPO}/{commit}/{path}"
    digest = hashlib.sha256()
    size = 0
    with urlopen(Request(url, headers={"User-Agent": "TSC48-age-relay/1"}), timeout=120) as response:
        while True:
            block = response.read(MIB)
            if not block:
                break
            size += len(block)
            if size > expected_bytes:
                raise RuntimeError("REMOTE_READBACK_OVERSIZE")
            digest.update(block)
    if size != expected_bytes or digest.hexdigest() != expected_sha:
        raise RuntimeError("REMOTE_READBACK_SHA_MISMATCH")
    return {"commit": commit, "path": path, "git_blob": remote_blob,
            "age_bytes": size, "age_sha256": digest.hexdigest(),
            "remote_sha256_readback": "PASS"}


def push_one(source: Path, repo_dir: Path, state_path: Path,
             branch: str, recipient: str, age_binary: str,
             max_objects: int, min_free: int) -> dict:
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if state.get("schema") != SCHEMA or state.get("repo") != REPO:
        raise RuntimeError("STATE_SCHEMA_MISMATCH")
    if state.get("recipient") != recipient or state.get("source_path") != str(source):
        raise RuntimeError("STATE_SOURCE_OR_RECIPIENT_MISMATCH")
    verify_repo(repo_dir, branch)
    source_stat = source.stat()
    saved_stat = state["source_stat"]
    if {"dev": source_stat.st_dev, "ino": source_stat.st_ino,
        "size": source_stat.st_size, "mtime_ns": source_stat.st_mtime_ns} != saved_stat:
        raise RuntimeError("SOURCE_STAT_CHANGED_SINCE_INIT")
    index = len(state["parts"])
    offset = index * state["chunk_bytes"]
    if offset >= state["source_bytes"]:
        return {"status": "ALL_PARTS_UPLOADED", "parts": index}
    objects = git_object_bytes(repo_dir / ".git")
    if objects + state["chunk_bytes"] > max_objects:
        raise RuntimeError("LOCAL_GIT_OBJECT_BUDGET_EXCEEDED")
    if shutil.disk_usage(repo_dir).free < min_free + state["chunk_bytes"]:
        raise RuntimeError("LOCAL_FREE_SPACE_RESERVE_WOULD_BE_BREACHED")
    with source.open("rb") as stream:
        stream.seek(offset)
        raw = stream.read(state["chunk_bytes"])
    if not raw:
        raise RuntimeError("UNEXPECTED_EMPTY_SOURCE_CHUNK")
    plain_sha = hashlib.sha256(raw).hexdigest()
    payload_dir = repo_dir / "payload"
    payload_dir.mkdir(exist_ok=True)
    path = f"payload/part_{index:06d}.age"
    staged_cipher = repo_dir / (path + ".writing")
    final_cipher = repo_dir / path
    staged_cipher.unlink(missing_ok=True)
    try:
        with staged_cipher.open("xb") as out:
            proc = subprocess.run([age_binary, "-r", recipient], input=raw,
                                  stdout=out, stderr=subprocess.PIPE, check=False)
            out.flush()
            os.fsync(out.fileno())
        if proc.returncode != 0:
            raise RuntimeError("AGE_ENCRYPT_FAILED: " + proc.stderr.decode("utf-8", "replace")[:500])
        os.replace(staged_cipher, final_cipher)
        age_sha, age_bytes = hash_path(final_cipher)
        pending = {"index": index, "path": path, "plain_bytes": len(raw),
                   "plain_sha256": plain_sha, "age_bytes": age_bytes,
                   "age_sha256": age_sha}
        progress = {"schema": SCHEMA, "status": "UPLOADING_NOT_READY",
                    "repo": REPO, "branch": branch,
                    "source_bytes": state["source_bytes"],
                    "source_sha256": state["source_sha256"],
                    "recipient": recipient, "completed_parts": state["parts"],
                    "pending_part": pending}
        atomic_json(repo_dir / "TRANSFER_PROGRESS.json", progress)
        run(["git", "add", "-A", "--", "payload", "TRANSFER_PROGRESS.json"], cwd=repo_dir)
        run(["git", "-c", "gc.auto=0", "commit", "-m",
             f"Upload encrypted TSC48 checkpoint part {index:06d}"], cwd=repo_dir)
        commit = run(["git", "rev-parse", "HEAD"], cwd=repo_dir)
        blob = run(["git", "rev-parse", f"HEAD:{path}"], cwd=repo_dir)
        run(["git", "push", "origin", f"HEAD:refs/heads/{branch}"], cwd=repo_dir,
            capture=False)
        evidence = remote_verify(commit, path, blob, age_sha, age_bytes)
        record = {**pending, **evidence}
        state["parts"].append(record)
        atomic_json(state_path, state)
        final_cipher.unlink()
        return {"status": "REMOTE_VERIFIED", **record,
                "local_ciphertext_removed": True,
                "local_git_object_bytes": git_object_bytes(repo_dir / ".git"),
                "free_bytes": shutil.disk_usage(repo_dir).free}
    except BaseException:
        # Preserve any staged ciphertext and local commit for inspection/retry.
        raise
    finally:
        staged_cipher.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    init = subs.add_parser("init")
    init.add_argument("--source", type=Path, required=True)
    init.add_argument("--state", type=Path, required=True)
    init.add_argument("--recipient", required=True)
    init.add_argument("--expected-sha256", required=True)
    init.add_argument("--chunk-mib", type=int, default=32)
    push = subs.add_parser("push-one")
    push.add_argument("--source", type=Path, required=True)
    push.add_argument("--repo-dir", type=Path, required=True)
    push.add_argument("--state", type=Path, required=True)
    push.add_argument("--branch", required=True)
    push.add_argument("--recipient", required=True)
    push.add_argument("--age", default="age")
    push.add_argument("--max-object-bytes", type=int, default=12_000_000_000)
    push.add_argument("--min-free-bytes", type=int, default=1_000_000_000)
    args = parser.parse_args()
    if args.command == "init":
        if args.chunk_mib != 32:
            raise SystemExit("CHUNK_SIZE_MUST_BE_32_MIB")
        if args.state.exists():
            raise SystemExit("STATE_ALREADY_EXISTS")
        result = init_state(args.source, args.state, args.recipient,
                            args.expected_sha256, args.chunk_mib * MIB)
    else:
        result = push_one(args.source, args.repo_dir, args.state, args.branch,
                          args.recipient, args.age, args.max_object_bytes,
                          args.min_free_bytes)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
