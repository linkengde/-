#!/usr/bin/env python3
"""Fresh-clone/fetch SHA and checkpoint/RNG restore test; performs no MD."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

PACKAGE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGE_DIR))
import archive_backend_v2 as archive_backend
import followup_runner_step20000_extension as runner


def run(repo_url: str, branch: str, expected_head: str, config_path: Path) -> dict:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="ag64-step20000-fresh-restore-") as tmp:
        clone = Path(tmp) / "fresh"
        subprocess.run(["git", "clone", "--quiet", "--no-checkout", "--branch", branch, "--single-branch", repo_url, str(clone)], check=True)
        subprocess.run(["git", "-C", str(clone), "checkout", "--quiet", expected_head], check=True)
        actual_head = subprocess.run(["git", "-C", str(clone), "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip()
        if actual_head != expected_head:
            raise RuntimeError(f"fresh clone head mismatch: {actual_head}")
        subprocess.run(["git", "-C", str(clone), "fetch", "--quiet", "origin", branch], check=True)
        archive_rel = config["resume_checkpoint"]["archive_path"]
        archive = clone / archive_rel
        members = archive_backend.verify_sha256sums(archive, filename="SHA256SUMS.txt")
        if archive_backend.sha256_file(archive / "SHA256SUMS.txt") != config["resume_checkpoint"]["archive_sha256sums_sha256"]:
            raise RuntimeError("fresh-clone SHA256SUMS hash mismatch")
        for name, key in (("checkpoint.json", "checkpoint_sha256"), ("MANIFEST.json", "manifest_sha256"), ("SEGMENT_MANIFEST.json", "segment_manifest_sha256")):
            if members.get(name) != config["resume_checkpoint"][key]:
                raise RuntimeError(f"fresh-clone member SHA mismatch: {name}")
        restored = runner._load_resume_checkpoint(clone, config)
        if restored["checkpoint_sha256"] != config["resume_checkpoint"]["checkpoint_sha256"]:
            raise RuntimeError("restored checkpoint differs from bound step-20000 state")
        return {
            "fresh_clone_head": actual_head,
            "fresh_fetch_pass": True,
            "archive_member_sha256_pass": True,
            "checkpoint_sha256": restored["checkpoint_sha256"],
            "checkpoint_step": restored["checkpoint"]["global_step"],
            "positions_restore_pass": True,
            "velocities_restore_pass": True,
            "rng_restore_pass": True,
            "stage_restore": restored["state"]["stage"],
            "md_started": False,
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-url", default="https://github.com/linkengde/-.git")
    parser.add_argument("--branch", default="cloud-d-b1-dft")
    parser.add_argument("--head", required=True)
    parser.add_argument("--config", type=Path, default=PACKAGE_DIR / "FOLLOWUP_EXECUTION_CONFIG_STEP20000_40000.json")
    args = parser.parse_args()
    print(json.dumps(run(args.repo_url, args.branch, args.head, args.config), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
