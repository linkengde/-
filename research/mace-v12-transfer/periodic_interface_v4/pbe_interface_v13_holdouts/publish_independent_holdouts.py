#!/usr/bin/env python3
"""Verify and publish A's completed v13 holdouts without force-pushing."""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
TRANSFER = ROOT.parents[1]
REPO = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True
).strip())
REL_TRANSFER = TRANSFER.relative_to(REPO)
TASK_REL = REL_TRANSFER / "coordination/tasks/window-a.json"
LOG_REL = REL_TRANSFER / "WORK_LOG_2026-10.md"
MANIFEST_REL = REL_TRANSFER / "file_manifest.json"
PROGRESS_REL = REL_TRANSFER / "coordination/progress/window-a.json"
RESULT_REL = ROOT.relative_to(REPO)
LABELS = (
    "AgC_lateral_registry_holdout_v13",
    "AgSi_lateral_registry_holdout_v13",
)
MAX_FILE_BYTES = 100_000_000


def git(*args, check=True):
    result = subprocess.run(
        ["git", *args], cwd=REPO, text=True, capture_output=True
    )
    if check and result.returncode:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip())
    return result


def read_json(path):
    return json.loads(path.read_text())


def hash_file(path):
    return sha256(path.read_bytes()).hexdigest()


def validate_archives():
    manifest_path = ROOT / "archive_manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError("archive_manifest.json is missing; run the archive verifier first.")
    archive = read_json(manifest_path)
    records = {row["label"]: row for row in archive.get("records", [])}
    if set(records) != set(LABELS):
        raise RuntimeError("Archive manifest must contain exactly the two assigned v13 labels.")

    energies = {}
    for label in LABELS:
        folder = ROOT / "calculations" / label
        summary_path = folder / "summary.json"
        progress_path = folder / "progress.json"
        verification_path = folder / "verification.json"
        output_path = folder / f"{label}_PW_PBE.extxyz"
        log_path = folder / "gpaw.log"
        for path in (summary_path, progress_path, verification_path, output_path, log_path):
            if not path.is_file():
                raise RuntimeError(f"Required archived file is missing: {path.relative_to(REPO)}")

        summary = read_json(summary_path)
        progress = read_json(progress_path)
        verification = read_json(verification_path)
        row = records[label]
        if summary.get("scf_converged") is not True or progress.get("status") != "complete":
            raise RuntimeError(f"{label}: SCF is not marked complete.")
        if verification.get("status") != "PASS" or not all(verification.get("checks", {}).values()):
            raise RuntimeError(f"{label}: archive verification is not PASS.")
        if row.get("input_sha256") != summary.get("source_sha256"):
            raise RuntimeError(f"{label}: input hash differs from the summary.")
        if row.get("extxyz_sha256") != hash_file(output_path):
            raise RuntimeError(f"{label}: extxyz hash differs from the archive manifest.")
        if row.get("summary_sha256") != hash_file(summary_path):
            raise RuntimeError(f"{label}: summary hash differs from the archive manifest.")
        if row.get("state_gpw_archived") is not False:
            raise RuntimeError(f"{label}: archive manifest claims a state.gpw was archived.")
        if f"Converged in {summary['scf_iterations']} steps" not in log_path.read_text(errors="replace"):
            raise RuntimeError(f"{label}: convergence confirmation is missing from gpaw.log.")
        energies[label] = (summary["scf_iterations"], summary["energy_eV_cell"])

    for path in (ROOT / "calculations").rglob("*"):
        if path.is_file() and (path.suffix == ".gpw" or path.stat().st_size >= MAX_FILE_BYTES):
            raise RuntimeError(f"Refusing to publish checkpoint or oversized file: {path.relative_to(REPO)}")
    return energies


def current_identity():
    return subprocess.check_output(
        [sys.executable, str(TRANSFER / "coordination/sync_tasks.py"), "identity"],
        cwd=REPO,
        text=True,
    ).strip()


def assert_owner(identity):
    task = read_json(REPO / TASK_REL)
    if task.get("owner_instance") != identity:
        raise RuntimeError("Local window-a task is not owned by this machine.")
    remote_text = git("show", f"origin/main:{TASK_REL.as_posix()}").stdout
    remote = json.loads(remote_text)
    if remote.get("owner_instance") != identity:
        raise RuntimeError("Remote window-a owner changed; refusing to publish or take over.")
    return task


def verify_changes_are_scoped():
    if git("diff", "--cached", "--name-only").stdout.strip():
        raise RuntimeError("The Git index already has staged changes; publish them separately first.")
    allowed = {
        TASK_REL.as_posix(),
        LOG_REL.as_posix(),
        MANIFEST_REL.as_posix(),
    }
    tracked_changes = set(git("diff", "--name-only").stdout.splitlines())
    unexpected = sorted(path for path in tracked_changes if path not in allowed)
    if unexpected:
        raise RuntimeError("Unrelated tracked changes are present: " + ", ".join(unexpected))
    untracked = git("ls-files", "--others", "--exclude-standard").stdout.splitlines()
    archive_manifest = (RESULT_REL / "archive_manifest.json").as_posix()
    progress_log = PROGRESS_REL.as_posix()
    outside = sorted(
        path for path in untracked
        if path not in (archive_manifest, progress_log)
        and not path.startswith(RESULT_REL.as_posix() + "/calculations/")
    )
    if outside:
        raise RuntimeError("Unrelated untracked files are present: " + ", ".join(outside))


def update_work_log(energies):
    path = REPO / LOG_REL
    text = path.read_text()
    marker = "### v13 independent holdout archive"
    if marker not in text:
        lines = [
            "",
            marker,
            "",
            "Both assigned lateral-registry v13 PW-PBE holdouts passed the archive checks and were published independently of window B.",
        ]
        for label in LABELS:
            iterations, energy = energies[label]
            lines.append(f"- `{label}`: SCF converged in {iterations} iterations; energy {energy:.10f} eV/cell.")
        lines.extend([
            "- Evidence: `periodic_interface_v4/pbe_interface_v13_holdouts/archive_manifest.json` and each label's `verification.json`.",
            "- These labels remain excluded from v13 training; score them before considering them for v14.",
            "",
        ])
        path.write_text(text.rstrip() + "\n" + "\n".join(lines))


def update_manifest():
    path = REPO / MANIFEST_REL
    manifest = read_json(path)
    candidates = [REPO / LOG_REL, ROOT / "archive_manifest.json"]
    candidates.extend(p for p in (ROOT / "calculations").rglob("*") if p.is_file())
    for candidate in candidates:
        if candidate.suffix == ".gpw" or candidate.stat().st_size >= MAX_FILE_BYTES:
            raise RuntimeError(f"Refusing to include checkpoint or oversized file: {candidate.relative_to(REPO)}")
        key = candidate.relative_to(TRANSFER).as_posix()
        manifest[key] = {"bytes": candidate.stat().st_size, "sha256": hash_file(candidate)}
    path.write_text(json.dumps(dict(sorted(manifest.items())), indent=2) + "\n")


def push_without_force(identity):
    for attempt in range(3):
        pushed = git("push", "origin", "HEAD:main", check=False)
        if pushed.returncode == 0:
            print(pushed.stdout + pushed.stderr)
            return
        git("fetch", "origin", "main")
        assert_owner(identity)
        rebased = git("rebase", "origin/main", check=False)
        if rebased.returncode:
            raise RuntimeError(
                "Push needs a rebase conflict resolution. Results and local commit are preserved; "
                "resolve the conflict on this A machine, then push HEAD:main.\n"
                + rebased.stderr.strip()
            )
    raise RuntimeError("Push was not confirmed after three ordinary attempts; keep this A machine and retry later.")


def main():
    energies = validate_archives()
    identity = current_identity()
    git("fetch", "origin", "main")
    task = assert_owner(identity)
    verify_changes_are_scoped()
    rebased = git("rebase", "origin/main", check=False)
    if rebased.returncode:
        raise RuntimeError("Cannot sync with main safely. Preserve outputs and resolve the rebase on this A machine.\n" + rebased.stderr.strip())

    task = read_json(REPO / TASK_REL)
    prepared = all(
        task["current_jobs"].get(label, "").startswith("complete; archive verification PASS")
        for label in LABELS
    )
    if prepared:
        push_without_force(identity)
        print("A v13 archive commit is already prepared; publication retry completed.")
        return

    task["current_jobs"][LABELS[0]] = "complete; archive verification PASS"
    task["current_jobs"][LABELS[1]] = "complete; archive verification PASS"
    task["next_action"] = (
        "Pull main and verify window-b's published acquisition labels; then score the v13 holdouts, "
        "run the assigned v14 holdouts, and continue the v14 data/train/evaluate cycle."
    )
    task["updated_utc"] = datetime.now(timezone.utc).isoformat()
    progress_path = REPO / PROGRESS_REL
    if progress_path.exists():
        progress_log = read_json(progress_path)
        if progress_log.get("owner_instance") != identity:
            raise RuntimeError("A progress log belongs to another machine; refusing to take it over.")
    else:
        progress_log = {"task_id": "window-a", "owner_instance": identity, "events": []}
    for label in LABELS:
        iterations, _energy = energies[label]
        summary = read_json(ROOT / "calculations" / label / "summary.json")
        event = {
            "job": label,
            "state": "completed",
            "iteration": iterations,
            "elapsed_s": summary.get("elapsed_s"),
            "note": "SCF and archive verification PASS; result included in this publication.",
            "updated_utc": task["updated_utc"],
        }
        progress_log["events"].append(event)
        task["live_progress"] = {
            "task_id": "window-a",
            "owner_instance": identity,
            **event,
        }
    progress_log["updated_utc"] = task["updated_utc"]
    progress_log["events"] = progress_log["events"][-100:]
    (REPO / TASK_REL).write_text(json.dumps(task, ensure_ascii=False, indent=2) + "\n")
    progress_path.parent.mkdir(parents=True, exist_ok=True)
    progress_path.write_text(json.dumps(progress_log, ensure_ascii=False, indent=2) + "\n")
    update_work_log(energies)
    update_manifest()

    git("add", "--", TASK_REL.as_posix(), LOG_REL.as_posix(), MANIFEST_REL.as_posix(),
        PROGRESS_REL.as_posix(),
        (RESULT_REL / "archive_manifest.json").as_posix(),
        (RESULT_REL / "calculations").as_posix())
    staged = git("diff", "--cached", "--name-only").stdout.splitlines()
    allowed_prefixes = (RESULT_REL.as_posix() + "/",)
    allowed_exact = {TASK_REL.as_posix(), LOG_REL.as_posix(), MANIFEST_REL.as_posix(), PROGRESS_REL.as_posix()}
    forbidden = [p for p in staged if p not in allowed_exact and not p.startswith(allowed_prefixes)]
    if forbidden:
        raise RuntimeError("Refusing to commit files outside the A v13 publication set: " + ", ".join(forbidden))
    if not staged:
        print("No new archive files to publish.")
        return
    git("commit", "-m", "Archive and publish A v13 holdout labels")
    push_without_force(identity)
    print("A v13 holdouts, work log, task status, and SHA256 manifest are published to main.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        raise SystemExit(str(exc))
