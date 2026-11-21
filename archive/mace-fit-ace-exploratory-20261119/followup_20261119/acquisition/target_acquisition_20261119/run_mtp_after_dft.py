#!/usr/bin/env python3
"""Wait for verified v24 archives, then prepare and fit independent MTP candidates."""
from __future__ import annotations

import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "v24_c_response_run"
STATUS_PATH = RUN / "sequence_status.json"
OUT_STATUS = ROOT / "mtp_after_dft_status.json"
LOGDIR = ROOT / "mtp_after_dft_logs"
PYTHON = "/workspace/.venvs/gpaw-diagnosis/bin/python"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def save(status: dict) -> None:
    temp = OUT_STATUS.with_suffix(".json.tmp")
    temp.write_text(json.dumps(status, indent=2) + "\n")
    temp.replace(OUT_STATUS)


def stage(name: str, script: Path, status: dict) -> None:
    status.update({"stage": name, "stage_started_utc": utc()})
    save(status)
    log_path = LOGDIR / f"{name}.log"
    with log_path.open("w") as log:
        proc = subprocess.run([PYTHON, str(script)], cwd=ROOT, stdout=log,
                              stderr=subprocess.STDOUT, check=False)
    status.setdefault("stages", []).append({"name": name, "return_code": proc.returncode,
                                             "finished_utc": utc(), "log": str(log_path)})
    if proc.returncode != 0:
        status.update({"status": "failed_stop_no_retry", "stage": name,
                       "error": f"{name} returned {proc.returncode}", "finished_utc": utc()})
        save(status)
        raise SystemExit(proc.returncode)
    save(status)


def main() -> None:
    LOGDIR.mkdir(parents=True, exist_ok=True)
    status = {"status": "waiting_for_verified_dft", "stage": "wait",
              "created_utc": utc(), "stages": [],
              "note": ("MTP uses verified periodic full labels plus v11 center-Si forces; "
                       "crop energies and noncenter forces carry zero weight.")}
    save(status)
    while True:
        try:
            dft = json.loads(STATUS_PATH.read_text())
        except (OSError, json.JSONDecodeError):
            time.sleep(30)
            continue
        if dft.get("status") == "selected_three_archives_verified":
            status["dft_sequence_status"] = dft
            status["status"] = "running"
            save(status)
            stage("prepare_mtp_acquisition", ROOT / "prepare_mtp_acquisition.py", status)
            stage("fit_mtp_acquisition", ROOT / "fit_mtp_acquisition.py", status)
            status.update({"status": "fit_complete_needs_review", "stage": None,
                           "finished_utc": utc(),
                           "fit_summary": str(ROOT / "mtp_acquisition_refits/fit_summary.json")})
            save(status)
            return
        if dft.get("status") == "failed_stop_no_retry":
            status.update({"status": "blocked_by_dft_stop", "stage": None,
                           "dft_sequence_status": dft, "finished_utc": utc()})
            save(status)
            return
        time.sleep(30)


if __name__ == "__main__":
    main()
