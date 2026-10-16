#!/usr/bin/env python3
"""Offline, reproducible GPW chunk split / hash verification / reconstruction.

This file does not run DFT, connect to GitHub, delete original GPW, or push.
A separate, authenticated transfer process is required. Designed for
independent H windows. Uses constant-size I/O, SHA-256 for every chunk and
for the fully reconstructed file, plus atomic output replacement.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re

SCHEMA = "tsc48-gpw-chunks-v1"
MIB = 1024 * 1024
CHUNK_RE = re.compile(r"^part_[0-9]{6}\.bin$")


def hash_file(path: Path, bufsize: int = MIB) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(bufsize), b""):
            h.update(block)
    return h.hexdigest()


def save_json_atomically(path: Path, data: dict) -> None:
    tmp = path.with_name(path.name + ".writing")
    with tmp.open("xb") as f:
        f.write((json.dumps(data, sort_keys=True, indent=2) + "\n").encode("utf-8"))
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def split_file(source: Path, outdir: Path, chunk_mib: int = 32) -> dict:
    if not source.is_file() or source.is_symlink():
        raise ValueError("SOURCE_NOT_REGULAR_FILE")
    if outdir.exists():
        raise ValueError("OUTPUT_DIRECTORY_ALREADY_EXISTS")
    if not (1 <= chunk_mib <= 40):
        raise ValueError("CHUNK_MIB_MUST_BE_1_TO_40")
    size = source.stat().st_size
    if size == 0:
        raise ValueError("SOURCE_EMPTY")
    outdir.mkdir(parents=True, exist_ok=False)
    main = hashlib.sha256()
    parts = []
    total = 0
    with source.open("rb") as incoming:
        while True:
            payload = incoming.read(chunk_mib * MIB)
            if not payload:
                break
            index = len(parts)
            if index >= 1000000:
                raise ValueError("TOO_MANY_CHUNKS")
            name = f"part_{index:06d}.bin"
            target = outdir / name
            with target.open("xb") as f:
                f.write(payload)
                f.flush()
                os.fsync(f.fileno())
            main.update(payload)
            parts.append({"file": name, "bytes": len(payload),
                          "sha256": hashlib.sha256(payload).hexdigest()})
            total += len(payload)
    if total != size:
        raise RuntimeError("SOURCE_SIZE_CHANGED_DURING_SPLIT")
    manifest = {"schema": SCHEMA, "source_name": source.name,
                "original_bytes": total, "original_sha256": main.hexdigest(),
                "chunk_mib": chunk_mib, "parts": parts}
    save_json_atomically(outdir / "manifest.json", manifest)
    return manifest


def load_manifest(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != SCHEMA:
        raise ValueError("MANIFEST_SCHEMA_MISMATCH")
    parts = data.get("parts")
    if not isinstance(parts, list) or not parts:
        raise ValueError("MANIFEST_EMPTY")
    if not isinstance(data.get("original_bytes"), int) or data["original_bytes"] <= 0:
        raise ValueError("MANIFEST_SIZE_BAD")
    if not re.fullmatch("[0-9a-f]{64}", str(data.get("original_sha256"))):
        raise ValueError("MANIFEST_SHA_BAD")
    names = []
    sizes = 0
    for idx, part in enumerate(parts):
        name = part.get("file") if isinstance(part, dict) else None
        if not isinstance(name, str) or not CHUNK_RE.fullmatch(name):
            raise ValueError("BAD_PART_NAME")
        if name != f"part_{idx:06d}.bin":
            raise ValueError("PART_ORDER_MISMATCH")
        if not isinstance(part.get("bytes"), int) or part["bytes"] <= 0:
            raise ValueError("BAD_PART_SIZE")
        if not re.fullmatch("[0-9a-f]{64}", str(part.get("sha256"))):
            raise ValueError("BAD_PART_SHA")
        names.append(name)
        sizes += part["bytes"]
    if len(names) != len(set(names)) or sizes != data["original_bytes"]:
        raise ValueError("MANIFEST_TOTAL_MISMATCH")
    return data


def verify_parts(manifest_path: Path) -> dict:
    manifest = load_manifest(manifest_path)
    directory = manifest_path.parent
    for part in manifest["parts"]:
        f = directory / part["file"]
        if f.is_symlink() or not f.is_file():
            raise ValueError("MISSING_OR_SYMLINK_PART: " + part["file"])
        if f.stat().st_size != part["bytes"]:
            raise ValueError("PART_SIZE_MISMATCH: " + part["file"])
        if hash_file(f) != part["sha256"]:
            raise ValueError("PART_SHA_MISMATCH: " + part["file"])
    return manifest


def join_file(manifest_path: Path, destination: Path) -> dict:
    manifest = load_manifest(manifest_path)
    if destination.exists() or destination.is_symlink():
        raise ValueError("DESTINATION_ALREADY_EXISTS")
    if not destination.parent.is_dir():
        raise ValueError("DESTINATION_PARENT_NOT_FOUND")
    temp = destination.with_name(destination.name + ".writing")
    if temp.exists() or temp.is_symlink():
        raise ValueError("STALE_OUTPUT_WRITING_FILE")
    directory = manifest_path.parent
    assembled = hashlib.sha256()
    count = 0
    try:
        with temp.open("xb") as target:
            for part in manifest["parts"]:
                file = directory / part["file"]
                if not file.is_file() or file.is_symlink():
                    raise ValueError("MISSING_OR_SYMLINK_PART: " + part["file"])
                if file.stat().st_size != part["bytes"]:
                    raise ValueError("PART_SIZE_MISMATCH: " + part["file"])
                local = hashlib.sha256()
                subtotal = 0
                with file.open("rb") as incoming:
                    for block in iter(lambda: incoming.read(MIB), b""):
                        local.update(block)
                        assembled.update(block)
                        target.write(block)
                        subtotal += len(block)
                if subtotal != part["bytes"] or local.hexdigest() != part["sha256"]:
                    raise ValueError("PART_CONTENT_MISMATCH: " + part["file"])
                count += subtotal
            target.flush()
            os.fsync(target.fileno())
        if count != manifest["original_bytes"] or assembled.hexdigest() != manifest["original_sha256"]:
            raise ValueError("REASSEMBLED_WHOLE_FILE_SHA_MISMATCH")
        os.replace(temp, destination)
    except BaseException:
        if temp.is_file():
            temp.unlink()
        raise
    return {"path": str(destination), "bytes": count,
            "sha256": assembled.hexdigest(), "parts": len(manifest["parts"])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    a = subs.add_parser("split")
    a.add_argument("source", type=Path)
    a.add_argument("output_directory", type=Path)
    a.add_argument("--chunk-mib", type=int, default=32)
    a = subs.add_parser("verify")
    a.add_argument("manifest", type=Path)
    a = subs.add_parser("join")
    a.add_argument("manifest", type=Path)
    a.add_argument("destination", type=Path)
    args = parser.parse_args()
    if args.command == "split":
        result = split_file(args.source, args.output_directory, args.chunk_mib)
        print(json.dumps({"parts": len(result["parts"]),
                          "original_bytes": result["original_bytes"],
                          "original_sha256": result["original_sha256"]}))
    elif args.command == "verify":
        result = verify_parts(args.manifest)
        print(json.dumps({"status": "PASS", "parts": len(result["parts"]),
                          "original_sha256": result["original_sha256"]}))
    else:
        print(json.dumps(join_file(args.manifest, args.destination)))


if __name__ == "__main__":
    main()
