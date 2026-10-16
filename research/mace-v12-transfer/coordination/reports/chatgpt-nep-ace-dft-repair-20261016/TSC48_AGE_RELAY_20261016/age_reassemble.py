#!/usr/bin/env python3
"""Stream independently age-encrypted checkpoint chunks from immutable Git commits.

The manifest contains only ciphertext locations and per-chunk/whole-file hashes.
Each ciphertext is downloaded and verified, decrypted with the receiver's local
identity, checked against the plaintext chunk hash, then appended to a temporary
output. The identity is passed to the age process and is never read or copied here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from urllib.parse import quote
from urllib.request import Request, urlopen

MIB = 1024 * 1024
SCHEMA = "h2-age-gpw-relay-v1"
REPO = "linkengde/-"
HEX64 = re.compile(r"^[0-9a-f]{64}$")
HEX40 = re.compile(r"^[0-9a-f]{40}$")
SAFE_PATH = re.compile(r"^[A-Za-z0-9._/-]+$")


def sha256_file(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(MIB), b""):
            digest.update(block)
            size += len(block)
    return digest.hexdigest(), size


def validate_manifest(data: dict) -> list[dict]:
    if data.get("schema") != SCHEMA or data.get("repo") != REPO:
        raise ValueError("MANIFEST_SCHEMA_OR_REPO_MISMATCH")
    if data.get("encryption") != "age-v1-per-chunk":
        raise ValueError("ENCRYPTION_FORMAT_MISMATCH")
    total = data.get("source_bytes")
    source_sha = data.get("source_sha256")
    parts = data.get("parts")
    if not isinstance(total, int) or total <= 0 or not isinstance(parts, list) or not parts:
        raise ValueError("MANIFEST_TOTAL_OR_PARTS_INVALID")
    if not isinstance(source_sha, str) or not HEX64.fullmatch(source_sha):
        raise ValueError("MANIFEST_SOURCE_SHA_INVALID")
    plain_total = 0
    for index, part in enumerate(parts):
        if not isinstance(part, dict) or part.get("index") != index:
            raise ValueError("PART_ORDER_INVALID")
        commit = part.get("commit")
        path = part.get("path")
        if not isinstance(commit, str) or not HEX40.fullmatch(commit):
            raise ValueError("PART_COMMIT_INVALID")
        if (not isinstance(path, str) or not SAFE_PATH.fullmatch(path)
                or path.startswith("/") or ".." in Path(path).parts):
            raise ValueError("PART_PATH_INVALID")
        for name in ("plain_sha256", "age_sha256"):
            if not isinstance(part.get(name), str) or not HEX64.fullmatch(part[name]):
                raise ValueError("PART_SHA_INVALID: " + name)
        for name in ("plain_bytes", "age_bytes"):
            if not isinstance(part.get(name), int) or part[name] <= 0:
                raise ValueError("PART_SIZE_INVALID: " + name)
        if part["plain_bytes"] > 32 * MIB:
            raise ValueError("PLAINTEXT_PART_EXCEEDS_32_MIB")
        plain_total += part["plain_bytes"]
    if plain_total != total:
        raise ValueError("PART_TOTAL_BYTES_MISMATCH")
    return parts


def remote_url(commit: str, path: str) -> str:
    return "https://raw.githubusercontent.com/" + REPO + "/" + commit + "/" + quote(path, safe="/")


def download_verified(url: str, destination: Path, expected_sha: str,
                      expected_bytes: int, timeout: int = 120) -> None:
    digest = hashlib.sha256()
    size = 0
    request = Request(url, headers={"User-Agent": "TSC48-age-relay/1"})
    with urlopen(request, timeout=timeout) as response, destination.open("xb") as out:
        while True:
            block = response.read(MIB)
            if not block:
                break
            size += len(block)
            if size > expected_bytes:
                raise ValueError("REMOTE_PART_EXCEEDS_MANIFEST_SIZE")
            digest.update(block)
            out.write(block)
        out.flush()
        os.fsync(out.fileno())
    if size != expected_bytes or digest.hexdigest() != expected_sha:
        destination.unlink(missing_ok=True)
        raise ValueError("REMOTE_CIPHERTEXT_SIZE_OR_SHA_MISMATCH")


def decrypt_verified(age_binary: str, identity: Path, encrypted: Path,
                     plaintext: Path, expected_sha: str,
                     expected_bytes: int) -> None:
    digest = hashlib.sha256()
    size = 0
    with encrypted.open("rb") as src, plaintext.open("xb") as dst:
        process = subprocess.Popen([age_binary, "-d", "-i", str(identity)],
                                   stdin=src, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        assert process.stdout is not None
        try:
            for block in iter(lambda: process.stdout.read(MIB), b""):
                size += len(block)
                if size > expected_bytes:
                    process.kill()
                    raise ValueError("DECRYPTED_PART_EXCEEDS_MANIFEST_SIZE")
                digest.update(block)
                dst.write(block)
            stderr = process.communicate()[1]
            if process.returncode != 0:
                raise RuntimeError("AGE_DECRYPT_FAILED: " + stderr.decode("utf-8", "replace")[:500])
            dst.flush()
            os.fsync(dst.fileno())
        except BaseException:
            if process.poll() is None:
                process.kill()
            process.wait()
            plaintext.unlink(missing_ok=True)
            raise
    if size != expected_bytes or digest.hexdigest() != expected_sha:
        plaintext.unlink(missing_ok=True)
        raise ValueError("DECRYPTED_PART_SIZE_OR_SHA_MISMATCH")


def reassemble(manifest_path: Path, identity: Path, destination: Path,
               age_binary: str = "age") -> dict:
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    parts = validate_manifest(data)
    if not identity.is_file() or identity.is_symlink():
        raise ValueError("AGE_IDENTITY_NOT_A_REGULAR_FILE")
    if destination.exists() or destination.is_symlink():
        raise ValueError("DESTINATION_ALREADY_EXISTS")
    if not destination.parent.is_dir():
        raise ValueError("DESTINATION_PARENT_NOT_FOUND")
    writing = destination.with_name(destination.name + ".writing")
    if writing.exists() or writing.is_symlink():
        raise ValueError("STALE_WRITING_FILE_REQUIRES_MANUAL_REVIEW")
    whole = hashlib.sha256()
    total = 0
    try:
        with tempfile.TemporaryDirectory(prefix="tsc48-age-", dir=destination.parent) as tmp:
            temp = Path(tmp)
            with writing.open("xb") as output:
                for part in parts:
                    encrypted = temp / f"part_{part['index']:06d}.age"
                    plaintext = temp / f"part_{part['index']:06d}.plain"
                    download_verified(remote_url(part["commit"], part["path"]),
                                      encrypted, part["age_sha256"], part["age_bytes"])
                    decrypt_verified(age_binary, identity, encrypted, plaintext,
                                     part["plain_sha256"], part["plain_bytes"])
                    with plaintext.open("rb") as src:
                        while True:
                            block = src.read(MIB)
                            if not block:
                                break
                            output.write(block)
                            whole.update(block)
                            total += len(block)
                    encrypted.unlink()
                    plaintext.unlink()
                output.flush()
                os.fsync(output.fileno())
        if total != data["source_bytes"] or whole.hexdigest() != data["source_sha256"]:
            raise ValueError("REASSEMBLED_SOURCE_SIZE_OR_SHA_MISMATCH")
        os.replace(writing, destination)
        dir_fd = os.open(destination.parent, os.O_RDONLY)
        try:
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)
    except BaseException:
        # Keep a partial .writing file as evidence; never present it as a valid GPW.
        raise
    return {"status": "PASS", "path": str(destination), "bytes": total,
            "sha256": whole.hexdigest(), "parts": len(parts),
            "encryption": "age-v1-per-chunk"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--identity", type=Path, required=True,
                        help="H2-local age identity; never upload or transmit this file")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--age", default="age")
    args = parser.parse_args()
    result = reassemble(args.manifest, args.identity, args.output, args.age)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
