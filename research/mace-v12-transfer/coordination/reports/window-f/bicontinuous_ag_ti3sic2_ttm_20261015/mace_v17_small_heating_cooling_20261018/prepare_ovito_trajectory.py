#!/usr/bin/env python3
"""Create a single OVITO-ready heat/cool timeline with per-frame scalar audits."""
from pathlib import Path
import csv
import numpy as np
from ase.io import read, write
from ase.neighborlist import neighbor_list
from analyze_ag_structure import hot_mask

ROOT = Path(__file__).resolve().parent
sources = [ROOT / "HEAT_COOL.extxyz", ROOT / "HEATING_DT0P2FS.extxyz", ROOT / "COOLING_DT0P2FS.extxyz"]
frames = []
for path in sources:
    if not path.exists():
        continue
    frames.extend(read(path, index=":", format="extxyz"))
frames.sort(key=lambda a: float(a.info.get("time_ps", 0.0)))
unique = []
for atoms in frames:
    t = float(atoms.info.get("time_ps", 0.0))
    if unique and t <= float(unique[-1].info.get("time_ps", 0.0)) + 1e-12:
        continue
    unique.append(atoms)
if not unique:
    raise RuntimeError("No usable trajectory frames were found")

hot_start = read(ROOT / "HEAT_COOL.extxyz", index=0, format="extxyz")
hot_ag_indices = np.flatnonzero(np.asarray(hot_start.get_chemical_symbols()) == "Ag")
fixed_hot_ids = np.asarray(hot_start.arrays["id"])[hot_ag_indices[hot_mask(hot_start, hot_ag_indices)]]

rows = []
for frame_id, atoms in enumerate(unique):
    symbols = np.asarray(atoms.get_chemical_symbols())
    phase = np.asarray(atoms.arrays["phase"])
    f = np.asarray(atoms.arrays["mace_force_norm_eV_A"], dtype=float)
    i, j, d = neighbor_list("ijd", atoms, 4.5)
    u = i < j
    tic = u & (((symbols[i] == "Ti") & (symbols[j] == "C")) |
               ((symbols[i] == "C") & (symbols[j] == "Ti")))
    agag = u & (symbols[i] == "Ag") & (symbols[j] == "Ag")
    cross = u & (phase[i] != phase[j])
    ag = phase == "Ag"
    tsc = phase == "Ti3SiC2"
    hot_indices = np.flatnonzero(np.isin(np.asarray(atoms.arrays["id"]), fixed_hot_ids))
    hot_v = atoms.get_velocities()[hot_indices]
    hot_m = atoms.get_masses()[hot_indices, None]
    hot_vcm = np.sum(hot_m * hot_v, axis=0) / np.sum(hot_m)
    hot_temp = float(np.sum(hot_m * (hot_v - hot_vcm) ** 2) /
                     (max(1, 3 * len(hot_indices) - 3) * 8.617333262145e-5))
    info = dict(atoms.info)
    info.update({
        "frame_id": int(frame_id),
        "Fmax_all_eV_A": float(np.max(f)),
        "P95_all_eV_A": float(np.percentile(f, 95)),
        "P99_all_eV_A": float(np.percentile(f, 99)),
        "Fmax_Ag_eV_A": float(np.max(f[ag])),
        "Fmax_TSC_eV_A": float(np.max(f[tsc])),
        "min_Ti_C_A": float(np.min(d[tic])) if np.any(tic) else float("nan"),
        "min_Ag_Ag_A": float(np.min(d[agag])) if np.any(agag) else float("nan"),
        "min_Ag_TSC_A": float(np.min(d[cross])) if np.any(cross) else float("nan"),
        "n_atoms": int(len(atoms)),
        "n_hot_Ag_fixed_id_set": int(len(hot_indices)),
        "T_hot_Ag_fixed_id_set_K": hot_temp,
    })
    atoms.info = info
    rows.append({k: info[k] for k in ["frame_id", "stage", "step", "time_ps", "T_system_K", "T_hot_Ag_fixed_id_set_K",
                                       "T_Ag_K", "T_TSC_K", "energy_eV", "Fmax_all_eV_A", "P95_all_eV_A",
                                       "P99_all_eV_A", "Fmax_Ag_eV_A", "Fmax_TSC_eV_A", "min_Ti_C_A",
                                       "min_Ag_Ag_A", "min_Ag_TSC_A", "n_atoms"] if k in info})

out = ROOT / "AG_TSC_MACE_V17_HEAT_COOL_OVITO.extxyz"
if out.exists(): out.unlink()
write(out, unique, format="extxyz")
with (ROOT / "AG_TSC_MACE_V17_HEAT_COOL_OVITO_METRICS.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print(f"wrote {out.name}: {len(unique)} frames, {len(unique[0])} atoms/frame")
