#!/usr/bin/env python3
"""Continue the authorized slab sequence after the active AgSi 6x6 run.

Sequence: verify AgSi 6x6, run the same-geometry 8x8 control, compare the two,
then acquire AgC and AgTi at 6x6 only if the mesh comparison passes. No fitting,
official dataset mutation, or automatic retries are performed.
"""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
PYTHON = "/workspace/.venvs/gpaw-mpi/bin/python"
CALC = ROOT / "calculations"
STATUS = CALC / "mesh_sequence_status.json"
LOG = CALC / "mesh_sequence.log"
FIRST = "AgSi_JN1967_registryB_k6x6x1_sigma0p10_preflight_v21"
SECOND = "AgSi_JN1967_registryB_k8x8x1_sigma0p10_control_v21"
CANDIDATES = [
    "AgC_JN1967_registryB_k6x6x1_sigma0p10_acquisition_v21",
    "AgTi_JN1967_registryB_k6x6x1_sigma0p10_acquisition_v21",
]

def update(stage: str, **fields) -> None:
    data = {"stage": stage, "updated_epoch_s": time.time(), **fields}
    STATUS.write_text(json.dumps(data, indent=2) + "\n")

def run(args: list[str], log) -> None:
    log.write("\n$ " + " ".join(args) + "\n")
    log.flush()
    subprocess.run(args, cwd=ROOT.parents[3], stdout=log, stderr=subprocess.STDOUT, check=True)

def main() -> None:
    CALC.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as log:
        started = time.time()
        update("waiting_for_AgSi_6x6", label=FIRST)
        while not (CALC / FIRST).is_dir():
            live = CALC / f".running_{FIRST}"
            if live.is_dir():
                time.sleep(30)
                continue
            attempts = list((CALC / "failed_attempts").glob(FIRST + "__attempt_*")) if (CALC / "failed_attempts").exists() else []
            failed = [path for path in attempts if path.stat().st_mtime >= started]
            if failed:
                update("halted_AgSi_6x6_failed_no_retry", failure_archive=str(failed[-1]))
                return
            time.sleep(30)
        update("verify_AgSi_6x6", label=FIRST)
        run([PYTHON, str(ROOT / "verify_preflight.py"), FIRST], log)
        update("run_AgSi_8x8_control", label=SECOND)
        run([PYTHON, str(ROOT / "launch_preflight.py"), SECOND], log)
        run([PYTHON, str(ROOT / "verify_preflight.py"), SECOND], log)
        update("compare_meshes", labels=[FIRST, SECOND])
        run([PYTHON, str(ROOT / "compare_kmesh.py")], log)
        comparison = json.loads((ROOT / "comparisons" / f"comparison_{FIRST}__{SECOND}.json").read_text())
        if comparison["status"] != "PASS_WITHIN_PROVISIONAL_BUDGET":
            update("halted_mesh_sensitivity_exceeds_budget", comparison_status=comparison["status"], comparison_file=str((ROOT / "comparisons" / f"comparison_{FIRST}__{SECOND}.json").relative_to(ROOT)))
            return
        for label in CANDIDATES:
            update("acquire_candidate", label=label, mesh_comparison="PASS_WITHIN_PROVISIONAL_BUDGET")
            run([PYTHON, str(ROOT / "launch_preflight.py"), label], log)
            run([PYTHON, str(ROOT / "verify_preflight.py"), label], log)
        update("acquisition_labels_complete_training_still_disabled", labels=[FIRST, *CANDIDATES], control=SECOND)

if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        update("halted_subprocess_failure_no_automatic_retry", return_code=exc.returncode, command=exc.cmd)
        raise
