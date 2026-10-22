#!/usr/bin/env python3
"""Fail-closed owner, role, method, and hash checks for V19 dev labels."""
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
MANIFEST = HERE / "input_manifest.json"
TASKS = REPO / "research/mace-v12-transfer/coordination/tasks/window-a.json"
sys.path.insert(0, str(REPO / "research/mace-v12-transfer/coordination"))
import sync_tasks


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load(label, launch=False):
    manifest = json.loads(MANIFEST.read_text())
    if launch:
        assert manifest["launch_enabled"] is True, "Manifest launch disabled"
    record = next(row for row in manifest["records"] if row["label"] == label)
    if launch:
        assert record["launch_enabled"] is True, "Record launch disabled"
    assert record["role"] == "independent_development"
    assert record["A_lineage_review"] == "APPROVED_NOT_REUSED_V17_DEVELOPMENT_FAMILY"
    input_path = (REPO / record["input"]).resolve()
    assert input_path.is_relative_to(REPO) and input_path.is_file()
    assert sha(input_path) == record["input_sha256"], "Input SHA256 mismatch"
    assert record["source_CIF_sha256"] == manifest["source_CIF_sha256"]
    assert record["method"] == manifest["method"]
    assert record["pbc"] == manifest["pbc"]
    source_path = (REPO / manifest["source_CIF"]).resolve()
    audit_path = (REPO / manifest["geometry_audit"]).resolve()
    assert source_path.is_relative_to(REPO) and sha(source_path) == manifest["source_CIF_sha256"]
    assert audit_path.is_relative_to(REPO) and sha(audit_path) == manifest["geometry_audit_sha256"]
    assert json.loads(audit_path.read_text())["status"] == "PASS_GEOMETRY_AND_ROLE_SCREEN_PENDING_EXTERNAL_B_AUDIT"
    task = json.loads(TASKS.read_text())
    assert task.get("owner_instance") == sync_tasks.identity(), "A task belongs to another instance"
    assert label in task.get("labels", {}), "Label is outside the A assignment"
    assert label in task.get("label_result_directories", {})
    if launch:
        assert task["labels"][label] in ("queued", "running")
    return manifest, record


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: entry.py LABEL")
    print(json.dumps(load(sys.argv[1], launch=True)[1], indent=2))
