#!/usr/bin/env python3
"""Track closest contacts and the pair implicated in the heating force spike."""
from pathlib import Path
import csv
import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list

ROOT = Path(__file__).resolve().parent
rows = []
for atoms in read(ROOT / "HEAT_COOL.extxyz", index=":", format="extxyz"):
    ids = np.asarray(atoms.arrays["id"], dtype=int)
    symbols = np.asarray(atoms.get_chemical_symbols())
    phase = np.asarray(atoms.arrays["phase"])
    force = np.asarray(atoms.arrays["mace_force_norm_eV_A"], dtype=float)
    i, j, d = neighbor_list("ijd", atoms, 4.5)
    und = i < j
    masks = {
        "Ag_Ag": und & (symbols[i] == "Ag") & (symbols[j] == "Ag"),
        "Ag_TSC": und & (phase[i] != phase[j]),
        "Ti_C": und & (((symbols[i] == "Ti") & (symbols[j] == "C")) |
                       ((symbols[i] == "C") & (symbols[j] == "Ti"))),
    }
    row = {"step": int(atoms.info.get("step", -1)), "time_ps": float(atoms.info.get("time_ps", 0.0)),
           "T_hot_Ag_K": float(atoms.info.get("T_hot_Ag_K", float("nan"))),
           "Fmax_eV_A": float(force.max()), "Fmax_atom_id": int(ids[np.argmax(force)])}
    for name, mask in masks.items():
        if np.any(mask):
            k = np.flatnonzero(mask)[np.argmin(d[mask])]
            row[f"min_{name}_A"] = float(d[k])
            row[f"min_{name}_id1"] = int(ids[i[k]])
            row[f"min_{name}_element1"] = str(symbols[i[k]])
            row[f"min_{name}_id2"] = int(ids[j[k]])
            row[f"min_{name}_element2"] = str(symbols[j[k]])
        else:
            row[f"min_{name}_A"] = float("nan")
            row[f"min_{name}_id1"] = row[f"min_{name}_id2"] = -1
            row[f"min_{name}_element1"] = row[f"min_{name}_element2"] = ""
    pair = (((ids[i] == 1757) & (ids[j] == 1648)) |
            ((ids[i] == 1648) & (ids[j] == 1757)))
    row["Ag1757_Ag1648_A"] = float(d[pair][0]) if np.any(pair) else float("nan")
    rows.append(row)

with (ROOT / "SUSPECT_CONTACTS.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print(f"wrote {len(rows)} frames to SUSPECT_CONTACTS.csv")
