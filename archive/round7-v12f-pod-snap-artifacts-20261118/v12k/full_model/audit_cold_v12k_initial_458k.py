#!/usr/bin/env python3
"""Audit atom IDs, motion, forces, and periodic Ti-C contacts for v12k cold run."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
N = 458248
BOX_X = BOX_Y = 200.0
DATA = ROOT / "Ag_Ti3SiC2_bicontinuous_200x200x300.data"
DUMP = ROOT / "cold_v12k_initial_458k_1fs/end_0p01ps.dump"
LOG = ROOT / "cold_v12k_initial_458k_1fs/log.lammps"
OUT = ROOT / "cold_v12k_initial_458k_1fs/cold_screen_summary.json"

# Source data order: Ti, Si, C, Ag. The LAMMPS run remaps to Ag, C, Si, Ti.
lines = DATA.read_text().splitlines()
section = next(i for i, row in enumerate(lines) if row.strip().startswith("Atoms"))
source: dict[int, tuple[int, np.ndarray]] = {}
for row in lines[section + 1 :]:
    fields = row.split()
    if len(fields) >= 5 and fields[0].isdigit():
        source[int(fields[0])] = (int(fields[1]), np.asarray(fields[2:5], float))
    elif source:
        break

rows = DUMP.read_text().splitlines()
h = next(i for i, row in enumerate(rows) if row.startswith("ITEM: ATOMS "))
cols = rows[h].split()[2:]
index = {name: cols.index(name) for name in ("id", "type", "x", "y", "z", "fx", "fy", "fz")}
a = np.asarray([[float(value) for value in rows[j].split()] for j in range(h + 1, h + 1 + N)])
ids = a[:, index["id"]].astype(int)
assert len(source) == len(ids) == N and len(np.unique(ids)) == N
final_type = a[:, index["type"]].astype(int)
final_pos = a[:, [index["x"], index["y"], index["z"]]]
final_force = a[:, [index["fx"], index["fy"], index["fz"]]]
initial_type = np.asarray([source[i][0] for i in ids])
initial_pos = np.asarray([source[i][1] for i in ids])

# Neighbor search is periodic in x/y and nonperiodic in z.
def tic_neighbor_data(coords: np.ndarray, atom_types: np.ndarray, ti_t: int, c_t: int):
    ti_mask, c_mask = atom_types == ti_t, atom_types == c_t
    ti_pos, ti_ids = coords[ti_mask], ids[ti_mask]
    c_pos, c_ids = coords[c_mask], ids[c_mask]
    shifts = np.asarray([(x, y, 0.0) for x in (-BOX_X, 0.0, BOX_X) for y in (-BOX_Y, 0.0, BOX_Y)])
    image_pos = np.concatenate([ti_pos + shift for shift in shifts], axis=0)
    image_ids = np.tile(ti_ids, len(shifts))
    tree = cKDTree(image_pos)
    distances, nearest = tree.query(c_pos, k=1, workers=4)
    k = int(distances.argmin())
    all_neighbors = tree.query_ball_point(c_pos, 1.8, workers=4)
    pairs = set()
    for j, neighbors in enumerate(all_neighbors):
        for q in neighbors:
            pairs.add((int(image_ids[q]), int(c_ids[j])))
    return float(distances[k]), int(image_ids[nearest[k]]), int(c_ids[k]), len(pairs)

# LAMMPS type IDs: source Ti=1,C=3; remapped final Ti=4,C=2.
initial_min, initial_ti, initial_c, initial_short = tic_neighbor_data(initial_pos, initial_type, 1, 3)
final_min, final_ti, final_c, final_short = tic_neighbor_data(final_pos, final_type, 4, 2)
delta = final_pos - initial_pos
delta[:, 0] -= BOX_X * np.round(delta[:, 0] / BOX_X)
delta[:, 1] -= BOX_Y * np.round(delta[:, 1] / BOX_Y)
displacement = np.linalg.norm(delta, axis=1)
force_norm = np.linalg.norm(final_force, axis=1)
name = {1: "Ag", 2: "C", 3: "Si", 4: "Ti"}
thermo = []
for row in LOG.read_text().splitlines():
    fields = row.split()
    if len(fields) == 8 and fields[0].isdigit():
        try:
            thermo.append([float(value) for value in fields])
        except ValueError:
            pass
assert len(thermo) >= 2, "could not parse LAMMPS thermodynamic table"
start, end = thermo[0], thermo[-1]
summary = {
    "atoms_initial": len(source), "atoms_final": len(ids), "unique_ids_final": int(len(np.unique(ids))), "atom_loss": len(source) - len(ids),
    "dt_fs": 1.0, "steps": 10, "duration_ps": 0.01,
    "temperature_mobile_K_initial": start[3], "temperature_mobile_K_final": end[3],
    "potential_energy_eV_initial": start[4], "potential_energy_eV_final": end[4], "potential_energy_delta_eV": end[4] - start[4],
    "kinetic_energy_eV_initial": start[5], "kinetic_energy_eV_final": end[5], "kinetic_energy_delta_eV": end[5] - start[5],
    "total_energy_eV_initial": start[6], "total_energy_eV_final": end[6], "total_energy_delta_eV": end[6] - start[6],
    "Fmax_eV_A_initial": start[7], "Fmax_eV_A_final": end[7],
    "Fmax_by_element_eV_A_final": {name[t]: float(force_norm[final_type == t].max()) for t in np.unique(final_type)},
    "max_displacement_A": float(displacement.max()), "displacement_p99_A": float(np.quantile(displacement, 0.99)),
    "TiC_pairs_below_1p8A_initial": initial_short, "TiC_pairs_below_1p8A_final": final_short,
    "min_TiC_distance_A_initial": initial_min, "min_TiC_pair_ids_initial": [initial_ti, initial_c],
    "min_TiC_distance_A_final": final_min, "min_TiC_pair_ids_final": [final_ti, final_c],
    "role": "no-laser short cold stability screen; not heating approval",
}
OUT.write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
