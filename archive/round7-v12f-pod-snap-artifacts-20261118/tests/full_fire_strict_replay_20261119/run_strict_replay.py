#!/usr/bin/env python3
"""Repeat the 347,464-atom FIRE replay with strict lost-atom handling."""
from __future__ import annotations
import os
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
LMP = Path("/workspace/software/lammps-pod-snap-clean/build/lmp")
INPUT = HERE / "in.fire_repro_strict"
start = time.time()
env = os.environ.copy()
env["LD_LIBRARY_PATH"] = (
    "/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:"
    + str(LMP.parent) + ":" + env.get("LD_LIBRARY_PATH", "")
)
env["OMP_NUM_THREADS"] = "1"
cmd = ["/usr/bin/mpirun", "-np", "5", str(LMP), "-in", str(INPUT), "-log", "log.strict.lammps"]
result = subprocess.run(
    cmd, cwd=HERE, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    text=True, preexec_fn=lambda: os.nice(10),
)
(HERE / "stdout.strict.log").write_text(result.stdout)
(HERE / "strict_run_metadata.txt").write_text(
    f"Strict FIRE replay; np=5; input={INPUT.name}; dt not applicable; "
    f"elapsed_s={time.time()-start:.3f}; rc={result.returncode}\n"
)
print(result.stdout[-5000:])
raise SystemExit(result.returncode)
