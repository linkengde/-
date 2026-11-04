#!/usr/bin/env python3
"""Run automatic-curvature normal/tangential force-closure checks."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
CASES = HERE / "cases"
QCC = Path(os.environ.get("BASILISK_QCC", "/tmp/basilisk-build/src/qcc"))
BASILISK = Path(os.environ.get("BASILISK_SRC", str(QCC.parent)))
SOURCE = HERE / "p1_surface_force_closure.c"

env = os.environ.copy()
env["BASILISK"] = str(BASILISK)
env["PATH"] = f"{BASILISK}:{env.get('PATH', '')}"
build = subprocess.run([str(QCC), "-grid=quadtree", "-O2", SOURCE.name,
                        "-o", "p1_surface_force_closure", "-lm"],
                       cwd=HERE, env=env, text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, check=False)
(RESULTS / "compile_closure.log").write_text(build.stdout)
if build.returncode:
    raise SystemExit("closure test compile failed; see results/compile_closure.log")

rows = []
for level in (4, 5, 6):
    tag = f"circle_closure_Rdx{1 << level}"
    proc = subprocess.run([str(HERE / "p1_surface_force_closure"), str(level)],
                          cwd=HERE, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, check=False)
    (CASES / f"{tag}.log").write_text(proc.stdout)
    if proc.returncode:
        raise SystemExit(f"{tag} failed; see cases/{tag}.log")
    line = next((line for line in proc.stdout.splitlines() if line.startswith("RESULT ")), None)
    if not line:
        raise SystemExit(f"no RESULT row from {tag}")
    rows.append(dict(re.findall(r"([A-Za-z0-9_]+)=([^ ]+)", line[len("RESULT "):])) )

with (RESULTS / "p1_surface_force_closure.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
config_path = RESULTS / "p1_solver_configuration.json"
config = json.loads(config_path.read_text())
config["qcc_sha256"] = hashlib.sha256(QCC.read_bytes()).hexdigest()
config["basilisk_header_sha256"] = {
    name: hashlib.sha256((BASILISK / name).read_bytes()).hexdigest()
    for name in ("tension.h", "iforce.h", "curvature.h", "two-phase.h")
}
config["basilisk_release_or_darcs_revision"] = "not recorded in the installed build; executable and key-header hashes are recorded"
config["official_basilisk_marangoni_case"] = "https://basilisk.fr/src/test/marangoni.c; not present in the installed source mirror, so no code-to-code reproduction is claimed"
config["surface_force_closure"] = {
    "source_file": SOURCE.name,
    "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    "analytic_profile": "sigma(theta)=1+0.1*cos(theta), R=0.5",
    "normal_force": "sigma*kappa*n using Basilisk automatic height-function curvature",
    "tangential_force": "projected surface gradient grad_s(sigma); both forces integrated over PLIC interface segments",
    "expected_normal_resultant_x": "-pi*sigma1",
    "expected_tangential_resultant_x": "+pi*sigma1",
    "expected_total_resultant": "zero",
    "results_csv": "p1_surface_force_closure.csv",
}
config_path.write_text(json.dumps(config, indent=2) + "\n")
print(f"wrote {RESULTS / 'p1_surface_force_closure.csv'}")
