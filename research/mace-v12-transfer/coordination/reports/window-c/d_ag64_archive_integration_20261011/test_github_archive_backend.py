from __future__ import annotations

import hashlib
import json
import subprocess
import tarfile
import tempfile
import unittest
from pathlib import Path

from github_archive_backend import GitHubArchiveBackend, _verify_restored_directory


def make_test_archive(root: Path) -> tuple[Path, str]:
    bundle = root / "synthetic-step-2000"
    bundle.mkdir()
    payloads = {
        "state.json": json.dumps({"step": 2000, "kind": "synthetic"}, sort_keys=True).encode() + b"\n",
        "state_arrays.bin": b"positions velocities RNG are fixture bytes only\n",
    }
    for name, data in payloads.items():
        (bundle / name).write_bytes(data)
    hashes = {name: hashlib.sha256(data).hexdigest() for name, data in payloads.items()}
    (bundle / "MANIFEST.json").write_text(json.dumps({
        "schema": 1, "owner": "D", "task_id": "C-D-AGL-1250K-rho097-seed01",
        "global_step": 2000, "members_sha256": hashes,
    }, sort_keys=True, indent=2) + "\n")
    all_hashes = dict(hashes)
    all_hashes["MANIFEST.json"] = hashlib.sha256((bundle / "MANIFEST.json").read_bytes()).hexdigest()
    (bundle / "SHA256SUMS.txt").write_text("".join(f"{all_hashes[n]}  {n}\n" for n in sorted(all_hashes)))
    archive = root / "synthetic-step-2000.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(bundle, arcname=bundle.name)
    return archive, hashlib.sha256(archive.read_bytes()).hexdigest()


def git(cwd: Path, *args: str) -> str:
    return subprocess.run(["git", "-C", str(cwd), *args], check=True,
                          text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.strip()


class GitHubArchiveBackendTests(unittest.TestCase):
    def test_upload_verify_restore_and_compare_from_fresh_clone(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bare = root / "remote.git"
            bare.mkdir()
            subprocess.run(["git", "init", "--bare", str(bare)], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            seed = root / "seed"
            subprocess.run(["git", "clone", str(bare), str(seed)], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            git(seed, "checkout", "-b", "cloud-test")
            (seed / "README.md").write_text("test remote\n")
            git(seed, "add", "README.md")
            subprocess.run(["git", "-C", str(seed), "-c", "user.name=Test",
                            "-c", "user.email=test@example.invalid", "commit", "-m", "base"],
                           check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            git(seed, "push", "-u", "origin", "cloud-test")
            base = git(seed, "rev-parse", "HEAD")
            archive, archive_sha = make_test_archive(root)
            backend = GitHubArchiveBackend(
                bare.as_uri(), "cloud-test", base, "synthetic-backend-tests",
                allow_local_test_remote=True, author_name="Test", author_email="test@example.invalid")
            receipt = backend.upload(archive, archive_sha)
            verified = backend.verify(receipt)
            self.assertEqual(verified["status"], "REMOTE_FRESH_FETCH_VERIFIED")
            restore_path = root / "fresh-restore"
            restored = backend.restore(receipt, restore_path)
            manifest = _verify_restored_directory(restore_path)
            self.assertEqual(restored["status"], "REMOTE_RESTORE_VERIFIED")
            self.assertEqual(manifest["global_step"], 2000)
            self.assertEqual(hashlib.sha256((restore_path / "state_arrays.bin").read_bytes()).hexdigest(),
                             hashlib.sha256((root / "synthetic-step-2000" / "state_arrays.bin").read_bytes()).hexdigest())

    def test_rejects_wrong_branch_head_before_upload(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bare = root / "remote.git"
            bare.mkdir()
            subprocess.run(["git", "init", "--bare", str(bare)], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            seed = root / "seed"
            subprocess.run(["git", "clone", str(bare), str(seed)], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            git(seed, "checkout", "-b", "cloud-test")
            (seed / "README.md").write_text("test remote\n")
            git(seed, "add", "README.md")
            subprocess.run(["git", "-C", str(seed), "-c", "user.name=Test",
                            "-c", "user.email=test@example.invalid", "commit", "-m", "base"],
                           check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            git(seed, "push", "-u", "origin", "cloud-test")
            archive, archive_sha = make_test_archive(root)
            wrong_head = "0" * 40
            backend = GitHubArchiveBackend(bare.as_uri(), "cloud-test", wrong_head,
                                           "synthetic-backend-tests", allow_local_test_remote=True)
            with self.assertRaisesRegex(Exception, "REMOTE_BRANCH_MOVED_SINCE_BOUND_HEAD"):
                backend.upload(archive, archive_sha)


if __name__ == "__main__":
    unittest.main()
