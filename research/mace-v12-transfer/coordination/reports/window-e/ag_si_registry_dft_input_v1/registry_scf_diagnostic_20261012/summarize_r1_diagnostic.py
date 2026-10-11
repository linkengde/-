#!/usr/bin/env python3
"""Extract the independent R1 diagnostic's SCF energy/residual history."""
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "R1_diagnostic_case/run_maxiter200_mixer_beta010"
LOG = RUN / "gpaw.log"
CSV = ROOT / "R1_DIAGNOSTIC_SCF_HISTORY.csv"
STATUS = ROOT / "R1_DIAGNOSTIC_CONVERGENCE_STATUS.json"

rows = []
for line in LOG.read_text(errors="replace").splitlines():
    if not line.startswith("|iter:"):
        continue
    fields = [field.strip() for field in line.strip().strip("|").split("|")]
    iter_match = re.fullmatch(r"iter:\s*(\d+)", fields[0])
    if len(fields) < 5 or not iter_match:
        raise SystemExit(f"Could not parse SCF row: {line}")
    def number(value):
        value = value.strip().removesuffix("c")
        return float(value) if re.fullmatch(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", value) else None
    rows.append({
        "iteration": int(iter_match.group(1)), "time": fields[1],
        "energy_eV": number(fields[2]), "eigst_log10_residual": number(fields[3]),
        "density_log10_residual": number(fields[4]),
        "energy_converged_marker": fields[2].endswith("c"),
        "eigenstate_converged_marker": fields[3].endswith("c"),
        "density_converged_marker": fields[4].endswith("c"),
    })

if not rows:
    raise SystemExit("No SCF rows found")
for i, row in enumerate(rows):
    row["delta_energy_eV"] = None if i == 0 else row["energy_eV"] - rows[i - 1]["energy_eV"]

with CSV.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

text = LOG.read_text(errors="replace")
converged = re.search(r"Converged in (\d+) steps", text)
maxiter_reached = len(rows) >= 200 and not converged
last = rows[-1]
status = {
    "case_id": "R1_DIAGNOSTIC_ONLY_MAXITER200_MIXER_BETA010",
    "run_log": str(LOG.relative_to(ROOT)),
    "history_csv": CSV.name,
    "iterations_recorded": len(rows),
    "convergence_status": "CONVERGED" if converged else ("MAXITER_REACHED_NOT_CONVERGED" if maxiter_reached else "IN_PROGRESS_OR_STATUS_UNRESOLVED"),
    "scf_steps_if_converged": int(converged.group(1)) if converged else None,
    "maxiter_reached": maxiter_reached,
    "last_iteration": last,
    "residual_columns": "GPAW eigst and dens columns are log10(error); c suffixes mark individual converged criteria.",
    "geometry_modified": False,
    "original_R1_R4_modified": False,
}
STATUS.write_text(json.dumps(status, indent=2) + "\n")
print(json.dumps(status, indent=2))
