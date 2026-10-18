#!/usr/bin/env python3
"""Write metrics and reports from the preserved partial heating trajectory."""
from pathlib import Path
import csv
import json
import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list

ROOT = Path(__file__).resolve().parent
frames = read(ROOT / "HEAT_COOL.extxyz", index=":", format="extxyz")
if not frames:
    raise RuntimeError("No heating frames found")
records = []
for atoms in frames:
    symbols = np.asarray(atoms.get_chemical_symbols())
    phase = np.asarray(atoms.arrays["phase"])
    norm = np.asarray(atoms.arrays["mace_force_norm_eV_A"], dtype=float)
    i, j, d = neighbor_list("ijd", atoms, 4.5)
    undirected = i < j
    tic = undirected & (((symbols[i] == "Ti") & (symbols[j] == "C")) |
                        ((symbols[i] == "C") & (symbols[j] == "Ti")))
    agag = undirected & (symbols[i] == "Ag") & (symbols[j] == "Ag")
    het = undirected & ((symbols[i] == "Ag") != (symbols[j] == "Ag"))
    info = atoms.info
    records.append({
        "step": int(info.get("step", -1)),
        "time_ps": float(info.get("time_ps", 0.0)),
        "T_system_K": float(info.get("T_system_K", float("nan"))),
        "T_hot_Ag_K": float(info.get("T_hot_Ag_K", float("nan"))),
        "T_Ag_K": float(info.get("T_Ag_K", float("nan"))),
        "T_TSC_K": float(info.get("T_TSC_K", float("nan"))),
        "energy_eV": float(info.get("energy_eV", float("nan"))),
        "Fmax_eV_A": float(np.max(norm)),
        "P95_eV_A": float(np.percentile(norm, 95)),
        "P99_eV_A": float(np.percentile(norm, 99)),
        "Fmax_Ag_eV_A": float(np.max(norm[symbols == "Ag"])),
        "Fmax_TSC_eV_A": float(np.max(norm[phase == "Ti3SiC2"])),
        "Fmax_atom_id": int(np.asarray(atoms.arrays["id"])[np.argmax(norm)]),
        "min_Ti_C_A": float(np.min(d[tic])) if np.any(tic) else float("nan"),
        "min_Ag_Ag_A": float(np.min(d[agag])) if np.any(agag) else float("nan"),
        "min_Ag_TSC_A": float(np.min(d[het])) if np.any(het) else float("nan"),
        "atoms_outside_z_margin": int(np.sum((atoms.positions[:, 2] < -3.0) |
                                               (atoms.positions[:, 2] > atoms.cell.lengths()[2] + 3.0))),
    })

with (ROOT / "HEATING_METRICS.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(records[0]))
    w.writeheader(); w.writerows(records)

result = {
    "stage": "heating",
    "requested_steps": 100,
    "completed_steps": records[-1]["step"],
    "dt_fs": 1.0,
    "simulated_time_ps": records[-1]["time_ps"],
    "wall_seconds_approx": 1710,
    "cpu_seconds_approx": 6194,
    "observed_process_RSS_peak_GiB_at_least": 16.2,
    "atom_count": len(frames[-1]),
    "final_frame": records[-1],
    "status": "stopped_due_to_short_Ag_Ag_contact_and_force_spike",
    "cooling_started": False,
    "trajectory": "HEAT_COOL.extxyz",
    "note": "The heating trajectory was stopped at saved step 90 (0.09 ps); no cooling segment was run. Timing is from process monitoring (approximately 28.5 min wall, 1 h 43 min aggregate CPU).",
}
(ROOT / "HEATING_REPORT.json").write_text(json.dumps(result, indent=2, allow_nan=True) + "\n")
combined = {
    "heating": result,
    "cooling": {"started": False, "reason": "heating was stopped after an anomalously short Ag-Ag contact correlated with a force spike"},
    "overall_status": "stopped_for_local_contact_anomaly_before_cooling",
}
(ROOT / "HEAT_COOL_REPORT.json").write_text(json.dumps(combined, indent=2, allow_nan=True) + "\n")
print(json.dumps(combined, allow_nan=True))
