#!/usr/bin/env python3
"""Write SHA256SUMS.txt for this directory, excluding the manifest itself."""

from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "SHA256SUMS.txt"


def main():
    entries = []
    for path in sorted(ROOT.rglob("*")):
        if path == OUTPUT or "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        if path.is_file():
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            entries.append(f"{digest}  {path.relative_to(ROOT).as_posix()}")
    OUTPUT.write_text("\n".join(entries) + "\n", encoding="utf-8")
    print(f"wrote {len(entries)} hashes to {OUTPUT}")


if __name__ == "__main__":
    main()
