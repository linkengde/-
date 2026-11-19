#!/usr/bin/env python3
"""After all selected v24 labels verify, prepare the candidate set and refit ACE once."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "v24_c_response_run"
STATUS = RUN / "sequence_status.json"
PIPELINE = ROOT / "fit_after_dft_status.json"
LOGDIR = ROOT / "fit_after_dft_logs"
PREPARE = ROOT / "prepare_acquisition_fit.py"
FIT = ROOT / "fit_ace_acquisition.py"
PY_ASE = "/workspace/.venvs/gpaw-diagnosis/bin/python"
PY_ACE = "/workspace/.venvs/ace/bin/python"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_status(value: dict) -> None:
    tmp = PIPELINE.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    tmp.replace(PIPELINE)


def run_stage(name: str, python: str, script: Path, status: dict) -> None:
    status.update({"stage": name, "stage_started_utc": utc()})
    write_status(status)
    with (LOGDIR / f"{name}.log").open("w") as log:
        proc = subprocess.run([python, str(script)], cwd=ROOT, stdout=log,
                              stderr=subprocess.STDOUT, check=False)
    status["stages"].append({"name": name, "return_code": proc.returncode,
                             "finished_utc": utc(), "log": str(LOGDIR / f"{name}.log")})
    if proc.returncode != 0:
        status.update({"status": "failed_stop_no_retry", "error": f"{name} returned {proc.returncode}",
                       "stage": name, "finished_utc": utc()})
        write_status(status)
        raise SystemExit(proc.returncode)


def main() -> None:
    LOGDIR.mkdir(parents=True, exist_ok=True)
    status = {"status": "waiting_for_verified_dft", "created_utc": utc(), "stage": "wait",
              "stages": [], "automatic_retries": 0,
              "note": "Only starts data preparation/fitting after the three registered v24 archives verify."}
    write_status(status)
    while True:
        try:
            dft_status = json.loads(STATUS.read_text())
        except (OSError, json.JSONDecodeError):
            time.sleep(20)
            continue
        if dft_status.get("status") == "selected_three_archives_verified":
            status["status"] = "running"
            status["dft_status_sha256"] = __import__("hashlib").sha256(STATUS.read_bytes()).hexdigest()
            write_status(status)
            run_stage("prepare_acquisition_fit", PY_ASE, PREPARE, status)
            run_stage("fit_ace_acquisition", PY_ACE, FIT, status)
            status.update({"status": "fit_complete_needs_review", "stage": None,
                           "finished_utc": utc(), "fit_summary": str(ROOT / "ace_acquisition_refits/fit_summary.json")})
            write_status(status)
            return
        if dft_status.get("status") == "failed_stop_no_retry":
            status.update({"status": "blocked_by_dft_stop", "stage": None,
                           "dft_status": str(STATUS), "finished_utc": utc()})
            write_status(status)
            return
        time.sleep(20)


if __name__ == "__main__":
    main()
