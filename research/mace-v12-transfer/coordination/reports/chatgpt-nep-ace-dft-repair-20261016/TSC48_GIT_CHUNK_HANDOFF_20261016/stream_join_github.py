#!/usr/bin/env python3
"""Stream one verified GitHub raw-file chunk at a time to a reconstructed GPW.

For PUBLIC GitHub repositories only. This is not a substitute for confidentiality,
restart validation or disk-capacity qualification. No checkout or local chunk copy.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
from urllib.parse import quote, urlparse
from urllib.request import urlopen
from chunk_transfer import load_manifest, MIB


def stream_join(manifest_path: Path, raw_base_url: str, destination: Path,
                opener=urlopen) -> dict:
    url = urlparse(raw_base_url)
    if url.scheme != "https" or url.hostname != "raw.githubusercontent.com" or url.query or url.fragment:
        raise ValueError("ONLY_HTTPS_RAW_GITHUB_SUPPORTED")
    if not raw_base_url.endswith("/"):
        raise ValueError("BASE_URL_MUST_END_IN_SLASH")
    manifest = load_manifest(manifest_path)
    if destination.exists() or destination.is_symlink():
        raise ValueError("DESTINATION_ALREADY_EXISTS")
    temp = destination.with_name(destination.name + ".writing")
    if temp.exists() or temp.is_symlink():
        raise ValueError("STALE_OUTPUT_WRITING_FILE")
    total = 0
    whole = hashlib.sha256()
    try:
        with temp.open("xb") as out:
            for part in manifest["parts"]:
                h = hashlib.sha256()
                size = 0
                remote_url = raw_base_url + quote(part["file"], safe="")
                with opener(remote_url) as incoming:
                    while True:
                        data = incoming.read(MIB)
                        if not data:
                            break
                        size += len(data)
                        if size > part["bytes"]:
                            raise ValueError("REMOTE_PART_TOO_LARGE: " + part["file"])
                        h.update(data)
                        whole.update(data)
                        out.write(data)
                if size != part["bytes"] or h.hexdigest() != part["sha256"]:
                    raise ValueError("REMOTE_PART_HASH_MISMATCH: " + part["file"])
                total += size
            out.flush()
            os.fsync(out.fileno())
        if total != manifest["original_bytes"] or whole.hexdigest() != manifest["original_sha256"]:
            raise ValueError("REMOTE_WHOLE_FILE_HASH_MISMATCH")
        os.replace(temp, destination)
    except BaseException:
        if temp.is_file():
            temp.unlink()
        raise
    return {"status": "SHA_MATCH", "bytes": total, "sha256": whole.hexdigest(),
            "parts": len(manifest["parts"]), "destination": str(destination)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("manifest", type=Path)
    p.add_argument("raw_base_url")
    p.add_argument("destination", type=Path)
    a = p.parse_args()
    print(json.dumps(stream_join(a.manifest, a.raw_base_url, a.destination)))


if __name__ == "__main__":
    main()
