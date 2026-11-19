#!/usr/bin/env python3
"""Archive verified DFT labels and the exploratory ACE fit, then push the task branch."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "v24_c_response_run"
FIT_STATUS = ROOT / "fit_after_dft_status.json"
SEQ_STATUS = RUN / "sequence_status.json"
OUT_STATUS = ROOT / "archive_after_dft_status.json"
REPO = Path("/workspace/-")
BRANCH = "codex/mace-fit-ace-exploratory-20261119"
REL = Path("archive/mace-fit-ace-exploratory-20261119/followup_20261119/acquisition/target_acquisition_20261119")
ARCHIVE = REPO / REL
SELECTED = [
    "AgTi_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24",
    "AgC_JN1967_registryB_surfaceC_dz_m0p02_sigma0p10_k6x6x1_v24",
    "AgC_JN1967_registryB_surfaceC_dz_p0p02_sigma0p10_k6x6x1_v24",
]


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def write_status(value: dict) -> None:
    tmp = OUT_STATUS.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    tmp.replace(OUT_STATUS)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def git(*args: str, check: bool = True) -> str:
    result = subprocess.run(["git", "-C", str(REPO), *args], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            check=False)
    if check and result.returncode:
        raise RuntimeError(f"git {' '.join(args)} failed: {result.stdout.strip()}")
    return result.stdout.strip()


def copy_file(source: Path, destination: Path) -> None:
    if source.is_file():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def model_report(fit_summary: dict, fit_status: dict, seq_status: dict) -> str:
    lines = [
        "# DFT acquisition and ACE refit report",
        "",
        f"Updated: {utc()}",
        f"DFT sequence status: `{seq_status.get('status', 'unknown')}`",
        f"ACE pipeline status: `{fit_status.get('status', 'unknown')}`",
        "",
        "## Registered exploratory gates",
        "",
        "- Energy RMSE: ≤ 10 meV/atom on full-label train and periodic development.",
        "- Force-component RMSE: ≤ 0.1 eV/Å on full-label train and periodic development.",
        "- Registered Si hotspot center-force vector error: ≤ 0.5 eV/Å for training centers and the v11 development center.",
        "- The v11 development frame is not a blind test.",
        "",
        "## Verified DFT labels",
        "",
    ]
    for label in SELECTED:
        folder = RUN / "calculations" / label
        verification = read_json(folder / "verification.json")
        summary = read_json(folder / "summary.json")
        lines.append(f"- `{label}`: verification `{verification.get('status', 'missing')}`, SCF converged `{summary.get('scf_converged', False)}`.")
    if fit_summary:
        lines += ["", "## ACE candidates", ""]
        passing = fit_summary.get("passing_models", [])
        lines.append("Passing models: " + (", ".join(f"`{x}`" for x in passing) if passing else "none"))
        for model in fit_summary.get("models", []):
            lines.append(
                f"- `{model.get('model')}`: gate_pass={model.get('gate_pass')}; "
                f"periodic-dev energy={model.get('energy_rmse_meV_atom')} meV/atom; "
                f"force-component RMSE={model.get('force_component_rmse_eV_A')} eV/Å; "
                f"Si-center vector error={model.get('hotspot_center_vector_error_eV_A')} eV/Å."
            )
    else:
        lines += ["", "ACE fit summary not present; inspect the archived fit status and logs before using these labels."]
    lines += [
        "",
        "## Limits and next step",
        "",
        "This is an exploratory candidate. Passing these gates does not establish transfer to the 1M-atom MACE-derived QSGS trajectory. The four direct QSGS crops remain geometry-only until their own DFT center-force labels are calculated. If no candidate meets every registered gate, use the reported residuals to choose the next acquisition; do not treat a failed fit as a valid production potential.",
        "",
    ]
    return "\n".join(lines)


def archive_outputs() -> dict:
    seq = read_json(SEQ_STATUS)
    fit = read_json(FIT_STATUS)
    fit_summary_path = ROOT / "ace_acquisition_refits" / "fit_summary.json"
    fit_summary = read_json(fit_summary_path)

    verified_labels = []
    label_root = ARCHIVE / "v24_verified_labels"
    for label in SELECTED:
        source = RUN / "calculations" / label
        verification = read_json(source / "verification.json")
        summary = read_json(source / "summary.json")
        if verification.get("status") == "PASS" and summary.get("scf_converged") is True:
            dest = label_root / label
            for name in ["summary.json", "verification.json", f"{label}_PW_PBE.extxyz", "gpaw.log"]:
                copy_file(source / name, dest / name)
            verified_labels.append(label)

    (ARCHIVE / "v24_c_response_run" / "sequence_status.json").write_text(
        json.dumps(seq, indent=2) + "\n")
    copy_file(FIT_STATUS, ARCHIVE / "fit_pipeline_status.json")
    copy_file(Path(__file__), ARCHIVE / "archive_fit_after_dft.py")

    candidate_dir = ROOT / "fit_candidate_acquisition"
    if candidate_dir.is_dir():
        shutil.copytree(candidate_dir, ARCHIVE / "fit_candidate_acquisition", dirs_exist_ok=True)
    fit_dir = ROOT / "ace_acquisition_refits"
    if fit_dir.is_dir():
        shutil.copytree(fit_dir, ARCHIVE / "ace_acquisition_refits", dirs_exist_ok=True)
    logs_dir = ROOT / "fit_after_dft_logs"
    if logs_dir.is_dir():
        shutil.copytree(logs_dir, ARCHIVE / "fit_logs", dirs_exist_ok=True)

    report = model_report(fit_summary, fit, seq)
    (ARCHIVE / "FINAL_REPORT.md").write_text(report)
    status = {
        "status": "archived_and_pushed",
        "updated_utc": utc(),
        "dft_status": seq.get("status"),
        "fit_status": fit.get("status"),
        "verified_labels": verified_labels,
        "passing_models": fit_summary.get("passing_models", []),
        "archive_path": str(REL),
        "branch": BRANCH,
    }

    readme_path = ARCHIVE / "README.md"
    if readme_path.is_file():
        readme = readme_path.read_text()
        readme = readme.replace("No new fit result is claimed yet.", "The verified DFT labels and ACE refit are recorded in `v24_verified_labels/`, `fit_candidate_acquisition/`, `ace_acquisition_refits/`, and `FINAL_REPORT.md`.")
        readme = re.sub(r"As of [^,]+, sequence status is `[^`]+`\.", f"Latest sequence status: `{seq.get('status', 'unknown')}`; verified labels: {len(verified_labels)}/{len(SELECTED)}.", readme)
        readme_path.write_text(readme)

    manifest = ARCHIVE / "SHA256SUMS.txt"
    rows = []
    for path in sorted(ARCHIVE.rglob("*")):
        if path.is_file() and path != manifest:
            rows.append(f"{sha(path)}  {path.relative_to(ARCHIVE).as_posix()}")
    manifest.write_text("\n".join(rows) + "\n")
    status["archive_files"] = len(rows)
    write_status(status)
    return status


def push_archive() -> str:
    if git("branch", "--show-current") != BRANCH:
        raise RuntimeError(f"working tree is not on expected branch {BRANCH}; refusing to commit")
    staged_before = git("diff", "--cached", "--name-only")
    if staged_before:
        raise RuntimeError("staged user changes exist; refusing to include them in the archive commit")
    git("add", "--", str(REL))
    staged = git("diff", "--cached", "--name-only").splitlines()
    if any(not item.startswith(str(REL) + "/") for item in staged):
        raise RuntimeError("staging included a path outside the acquisition archive")
    if staged:
        git("commit", "-m", "Archive verified DFT labels and ACE acquisition refit")
    return git("push", "origin", BRANCH)


def main() -> None:
    state = {"status": "waiting_for_dft_and_fit", "created_utc": utc(), "automatic_retries": 0}
    write_status(state)
    terminal_fit = {"fit_complete_needs_review", "failed_stop_no_retry", "blocked_by_dft_stop"}
    terminal_dft = {"selected_three_archives_verified", "failed_stop_no_retry"}
    while True:
        seq = read_json(SEQ_STATUS)
        fit = read_json(FIT_STATUS)
        if seq.get("status") in terminal_dft and fit.get("status") in terminal_fit:
            break
        time.sleep(30)
    state.update({"status": "archiving", "started_utc": utc()})
    write_status(state)
    try:
        result = archive_outputs()
        output = push_archive()
        result.update({"status": "archived_and_pushed", "push_output": output, "finished_utc": utc()})
        write_status(result)
    except Exception as exc:
        state.update({"status": "failed_stop_no_retry", "error": str(exc), "finished_utc": utc()})
        write_status(state)
        raise


if __name__ == "__main__":
    main()
