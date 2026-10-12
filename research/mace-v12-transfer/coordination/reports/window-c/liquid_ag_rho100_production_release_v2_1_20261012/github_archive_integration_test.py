#!/usr/bin/env python3
"""Push and fresh-fetch a synthetic archive on the real D GitHub branch.

The payload is fabricated and is never an MD checkpoint or scientific result.
"""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import archive_backend_v2 as backend


def main() -> int:
    repo = Path(subprocess.check_output(["git", "-C", str(HERE), "rev-parse", "--show-toplevel"], text=True).strip())
    branch = "cloud-d-b1-dft"
    remote = "origin"
    head_check = backend.verify_branch_head(repo, remote, branch)
    if not head_check["pass"]:
        raise RuntimeError("D local and GitHub branch heads differ before the synthetic test")
    base = head_check["local_head"]
    archive_relpath = (
        HERE.relative_to(repo).as_posix()
        + f"/github_synthetic_archive_test/head_{base[:12]}"
    )
    with tempfile.TemporaryDirectory(prefix="ag64-real-github-synthetic-") as temp_name:
        temp = Path(temp_name)
        source = temp / "synthetic_bundle"
        restore = temp / "restored_bundle"
        source.mkdir()
        checkpoint = {
            "test_only": True,
            "global_step": 77,
            "positions_A": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            "momenta_amu_A_per_fs": [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],
            "rng_states": {"synthetic": {"algorithm": "fixture", "state": 987654321}},
        }
        (source / "checkpoint.json").write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (source / "trajectory.synthetic.bin").write_bytes(b"synthetic trajectory bytes; not an MD trajectory\n")
        (source / "manifest.json").write_text(json.dumps({
            "test_only": True,
            "scientific_data": False,
            "purpose": "GitHub upload, fresh fetch, restore and SHA validation",
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (source / "README.md").write_text(
            "Synthetic archive integration fixture only. Contains no real positions, trajectory, potential, or simulation output.\n",
            encoding="utf-8",
        )
        result = backend.archive_push_fresh_restore(
            repo_root=repo,
            source_dir=source,
            archive_relpath=archive_relpath,
            branch=branch,
            remote=remote,
            expected_remote_head=base,
            restore_dir=restore,
            commit_message="Test synthetic Ag64 segment archive fresh-fetch and restore",
        )
        if not backend.segment_archive_passed(result):
            raise RuntimeError("synthetic archive result did not pass every required archive gate")
        restored = json.loads((restore / "checkpoint.json").read_text(encoding="utf-8"))
        if restored != checkpoint:
            raise RuntimeError("fresh-fetch checkpoint fixture differs from the uploaded payload")
        result["source_head"] = base
        result["archive_path"] = archive_relpath
        result["synthetic_payload_only"] = True
        result["restored_checkpoint_payload_match"] = True
        print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
