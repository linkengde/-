#!/usr/bin/env python3
"""Select spatially separated high-force centers for follow-up DFT design.

These are model-generated acquisition candidates only; no DFT labels are
created by this script. Type IDs follow the archived million-atom LAMMPS run:
1 Ti, 2 Si, 3 C, 4 Ag.
"""
import csv
import gzip
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DUMP = Path("/workspace/-/archive/mace-fit-ace-exploratory-20261119/million_atom_static_ace/force_hotspots.lammpstrj.gz")
OUT = ROOT / "label_design_center_shortlist.csv"
Lx, Ly, Lz = 299.635, 301.857484616085, 258.6115
ELEMENTS = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}
Z_BINS = 10
MAX_PER_ELEMENT_ZBIN = 3
MIN_SEPARATION_A = 8.0

rows = []
with gzip.open(DUMP, "rt") as f:
    for line in f:
        fields = line.split()
        if len(fields) != 9:
            continue
        try:
            atom_id, typ = int(fields[0]), int(fields[1])
            x, y, z, fmag = map(float, (fields[2], fields[3], fields[4], fields[5]))
            fx, fy, fz = map(float, fields[6:9])
        except ValueError:
            continue
        if typ not in ELEMENTS:
            continue
        zbin = min(Z_BINS - 1, max(0, int(z / Lz * Z_BINS)))
        rows.append((atom_id, typ, x, y, z, fmag, fx, fy, fz, zbin))

groups = defaultdict(list)
for row in rows:
    groups[(row[1], row[9])].append(row)

selected = []
for key in sorted(groups):
    candidates = sorted(groups[key], key=lambda r: (-r[5], r[0]))
    accepted = []
    for candidate in candidates:
        _, _, x, y, z, *_ = candidate
        separated = True
        for prior in accepted:
            dx = abs(x - prior[2]); dx = min(dx, Lx - dx)
            dy = abs(y - prior[3]); dy = min(dy, Ly - dy)
            dz = z - prior[4]
            if dx * dx + dy * dy + dz * dz < MIN_SEPARATION_A ** 2:
                separated = False
                break
        if separated:
            accepted.append(candidate)
        if len(accepted) >= MAX_PER_ELEMENT_ZBIN:
            break
    for rank, row in enumerate(accepted, 1):
        selected.append((key, rank, row))

fields = ["element", "z_bin", "within_bin_rank", "atom_id", "x_A", "y_A", "z_A",
          "model_force_eV_A", "fx_eV_A", "fy_eV_A", "fz_eV_A", "acquisition_status"]
with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for (typ, zbin), rank, row in sorted(selected, key=lambda x: (-x[2][5], x[2][0])):
        atom_id, _, x, y, z, fmag, fx, fy, fz, _ = row
        writer.writerow({
            "element": ELEMENTS[typ], "z_bin": zbin, "within_bin_rank": rank,
            "atom_id": atom_id, "x_A": f"{x:.8f}", "y_A": f"{y:.8f}", "z_A": f"{z:.8f}",
            "model_force_eV_A": f"{fmag:.8f}", "fx_eV_A": f"{fx:.8f}",
            "fy_eV_A": f"{fy:.8f}", "fz_eV_A": f"{fz:.8f}",
            "acquisition_status": "candidate center only; extract local environment and obtain DFT label",
        })
print(f"dumped_centers={len(rows)} selected_centers={len(selected)} output={OUT}")
