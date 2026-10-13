#!/usr/bin/env python3
"""Reassemble and verify the Stage67 300x300 cold model package."""
from __future__ import annotations
import hashlib
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parent
PARTS = ROOT / "parts"
OUTPUT = ROOT / "output"
ARCHIVE_NAME = "stage67_300x300_model_upload_20261013.zip"
ARCHIVE_SHA256 = "2ad27a42fd84a6cb85114155213e7d09940fce4e01db02c7a843142644d9730c"
MODEL_NAME = "stage67_model_150x300x300_300K_1x3x3.data"
MODEL_SHA256 = "36cbc94ec2aa4f7b13677d260fed4b55be7236583638be43b9f5a61e89a22cc1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    part_files = sorted(PARTS.glob("model_archive.part[0-9][0-9][0-9]"))
    if not part_files:
        raise SystemExit("No model archive parts found under parts/.")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    archive = OUTPUT / ARCHIVE_NAME
    with archive.open("wb") as out:
        for part in part_files:
            with part.open("rb") as src:
                for chunk in iter(lambda: src.read(1024 * 1024), b""):
                    out.write(chunk)
    got_archive = sha256(archive)
    if got_archive != ARCHIVE_SHA256:
        raise SystemExit(f"Archive SHA256 mismatch: {got_archive}")
    with zipfile.ZipFile(archive) as zf:
        bad = zf.testzip()
        if bad:
            raise SystemExit(f"Corrupt ZIP member: {bad}")
        model_data = zf.read(MODEL_NAME)
        got_model = hashlib.sha256(model_data).hexdigest()
        if got_model != MODEL_SHA256:
            raise SystemExit(f"Model SHA256 mismatch: {got_model}")
        target = OUTPUT / MODEL_NAME
        target.write_bytes(model_data)
        for name in ("README.md", "GEOMETRY_VALIDATION.md", "SHA256SUMS.txt"):
            (OUTPUT / name).write_bytes(zf.read(name))
    print(f"MODEL_READY = YES")
    print(f"MODEL_PATH = {target}")
    print(f"MODEL_SHA256 = {got_model}")
    print(f"ARCHIVE_SHA256 = {got_archive}")


if __name__ == "__main__":
    main()