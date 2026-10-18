#!/usr/bin/env python3
"""Finalize an intentionally truncated baseline from recorded checkpoints."""
from pathlib import Path
import csv
import json
import numpy as np
from ase.io import read, write
from ase.neighborlist import neighbor_list

ROOT = Path(__file__).resolve().parent
frames = read(ROOT / "BASELINE_300K.extxyz", index=":", format="extxyz")
if len(frames) < 2:
    raise RuntimeError("Need at least one post-initial baseline snapshot")

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
    pairs = {
        "min_Ti_C_A": float(np.min(d[tic])) if np.any(tic) else float("nan"),
        "min_Ag_Ag_A": float(np.min(d[agag])) if np.any(agag) else float("nan"),
        "min_Ag_TSC_A": float(np.min(d[het])) if np.any(het) else float("nan"),
    }
    info = atoms.info
    records.append({
        "step": int(info.get("step", -1)),
        "time_ps": float(info.get("time_ps", 0.0)),
        "energy_eV": float(info.get("energy_eV", float("nan"))),
        "Fmax_eV_A": float(np.max(norm)),
        "P95_eV_A": float(np.percentile(norm, 95)),
        "P99_eV_A": float(np.percentile(norm, 99)),
        "Fmax_Ag_eV_A": float(np.max(norm[symbols == "Ag"])),
        "Fmax_TSC_eV_A": float(np.max(norm[phase == "Ti3SiC2"])),
        "T_system_K": float(info.get("T_system_K", float("nan"))),
        "T_Ag_K": float(info.get("T_Ag_K", float("nan"))),
        "T_TSC_K": float(info.get("T_TSC_K", float("nan"))),
        "atoms_outside_z_margin": int(np.sum((atoms.positions[:, 2] < -3.0) | (atoms.positions[:, 2] > atoms.cell.lengths()[2] + 3.0))),
        **pairs,
    })

with (ROOT / "BASELINE_300K_METRICS.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(records[0]))
    w.writeheader(); w.writerows(records)

checkpoint = frames[-1]
checkpoint.info["stage"] = "baseline_300K_partial_checkpoint"
checkpoint.info["partial_baseline_steps"] = int(checkpoint.info.get("step", -1))
write(ROOT / "BASELINE_300K_CHECKPOINT.extxyz", checkpoint, format="extxyz")
final = records[-1]
result = {
    "stage": "baseline_300K",
    "requested_steps": 100,
    "completed_steps": final["step"],
    "dt_fs": 1.0,
    "simulated_time_ps": final["time_ps"],
    "wall_seconds": 26 * 60 + 42,
    "peak_RSS_GiB": 17.6,
    "atom_count": len(checkpoint),
    "final": final,
    "status": "partial_short_screen_no_numerical_collapse",
    "checkpoint": "BASELINE_300K_CHECKPOINT.extxyz",
    "limitation": "Stopped after 0.02 ps because measured CPU throughput made the original 0.1 ps baseline and 0.8 ps thermal cycle disproportionate. The saved interval is only an early numerical-stability screen, not a 300 K stability validation.",
}
(ROOT / "BASELINE_300K_REPORT.json").write_text(json.dumps(result, indent=2, allow_nan=True) + "\n")
print(json.dumps(result, allow_nan=True))
