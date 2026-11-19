#!/usr/bin/env python3
"""Audit full-model POD v12n cold trajectory from initial state through 50 fs."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
N = 458248
BOX_X = BOX_Y = 200.0
DATA = ROOT / "Ag_Ti3SiC2_bicontinuous_200x200x300.data"
DUMP = ROOT / "cold_v12n_continue_to_50fs/end_0p05ps.dump"
LOG0 = ROOT / "cold_v12n_initial_458k_1fs/log.lammps"
LOG50 = ROOT / "cold_v12n_continue_to_50fs/log.lammps"
OUT = ROOT / "cold_v12n_continue_to_50fs/cold_50fs_summary.json"

# Original source order: Ti, Si, C, Ag. LAMMPS potential order: Ag, C, Si, Ti.
source_lines = DATA.read_text().splitlines()
section = next(i for i, row in enumerate(source_lines) if row.strip().startswith("Atoms"))
source: dict[int, tuple[int, np.ndarray]] = {}
for row in source_lines[section + 1:]:
    fields = row.split()
    if len(fields) >= 5 and fields[0].isdigit():
        source[int(fields[0])] = (int(fields[1]), np.asarray(fields[2:5], float))
    elif source:
        break

rows = DUMP.read_text().splitlines()
h = next(i for i, row in enumerate(rows) if row.startswith("ITEM: ATOMS "))
cols = rows[h].split()[2:]
needed = ("id", "type", "x", "y", "z", "vx", "vy", "vz", "fx", "fy", "fz")
index = {name: cols.index(name) for name in needed}
a = np.asarray([[float(value) for value in rows[j].split()]
                for j in range(h + 1, h + 1 + N)])
ids = a[:, index["id"]].astype(int)
assert len(source) == len(ids) == N and len(np.unique(ids)) == N
final_type = a[:, index["type"]].astype(int)
final_pos = a[:, [index["x"], index["y"], index["z"]]]
final_force = a[:, [index["fx"], index["fy"], index["fz"]]]
initial_type = np.asarray([source[i][0] for i in ids])
initial_pos = np.asarray([source[i][1] for i in ids])

source_to_model_type = {1: 4, 2: 3, 3: 2, 4: 1}
assert np.array_equal(final_type, np.asarray([source_to_model_type[t] for t in initial_type]))

def min_cross_distance(coords, atom_types, id_a, id_b):
    pa, ia = coords[atom_types == id_a], ids[atom_types == id_a]
    pb, ib = coords[atom_types == id_b], ids[atom_types == id_b]
    shifts = np.asarray([(x, y, 0.0) for x in (-BOX_X, 0.0, BOX_X)
                         for y in (-BOX_Y, 0.0, BOX_Y)])
    image_pos = np.concatenate([pa + shift for shift in shifts], axis=0)
    image_ids = np.tile(ia, len(shifts))
    tree = cKDTree(image_pos)
    if id_a == id_b:
        dist, nearest = tree.query(pb, k=2, workers=4)
        dist, nearest = dist[:, 1], nearest[:, 1]
    else:
        dist, nearest = tree.query(pb, k=1, workers=4)
    k = int(dist.argmin())
    return float(dist[k]), int(image_ids[nearest[k]]), int(ib[k])

def tic_pairs(coords, atom_types):
    dmin, ti_id, c_id = min_cross_distance(coords, atom_types, 4, 2)
    ti_mask, c_mask = atom_types == 4, atom_types == 2
    ti_pos, ti_ids = coords[ti_mask], ids[ti_mask]
    c_pos, c_ids = coords[c_mask], ids[c_mask]
    shifts = np.asarray([(x, y, 0.0) for x in (-BOX_X, 0.0, BOX_X)
                         for y in (-BOX_Y, 0.0, BOX_Y)])
    image_pos = np.concatenate([ti_pos + shift for shift in shifts], axis=0)
    image_ids = np.tile(ti_ids, len(shifts))
    tree = cKDTree(image_pos)
    near = tree.query_ball_point(c_pos, 1.8, workers=4)
    pairs = {(int(image_ids[q]), int(c_ids[j])) for j, js in enumerate(near) for q in js}
    return dmin, ti_id, c_id, len(pairs)

def parse_thermo(path):
    out = []
    for row in path.read_text().splitlines():
        fields = row.split()
        if len(fields) == 8 and fields[0].isdigit():
            try:
                out.append([float(x) for x in fields])
            except ValueError:
                pass
    if not out:
        raise RuntimeError(f"No thermo rows parsed from {path}")
    return out

initial = parse_thermo(LOG0)[0]
final = parse_thermo(LOG50)[-1]
delta = final_pos - initial_pos
delta[:, 0] -= BOX_X * np.round(delta[:, 0] / BOX_X)
delta[:, 1] -= BOX_Y * np.round(delta[:, 1] / BOX_Y)
disp = np.linalg.norm(delta, axis=1)
fmag = np.linalg.norm(final_force, axis=1)
type_name = {1: "Ag", 2: "C", 3: "Si", 4: "Ti"}

tic_min, tic_ti, tic_c, tic_short = tic_pairs(final_pos, final_type)
cc_dist, cc_id_a, cc_id_b = min_cross_distance(final_pos, final_type, 2, 2)
# Track the specific pair that collapsed in v12m. Apply x/y periodic wrapping.
id_to_index = {int(atom_id): i for i, atom_id in enumerate(ids)}
i, j = id_to_index[128985], id_to_index[134363]
pair_delta = final_pos[i] - final_pos[j]
pair_delta[0] -= BOX_X * np.round(pair_delta[0] / BOX_X)
pair_delta[1] -= BOX_Y * np.round(pair_delta[1] / BOX_Y)
cc_target_dist = float(np.linalg.norm(pair_delta))

summary = {
    "model": "POD v12n",
    "atoms_initial": len(source), "atoms_final": len(ids),
    "unique_ids_final": int(len(np.unique(ids))), "atom_loss": len(source) - len(ids),
    "dt_fs": 1.0, "total_steps": 50, "duration_ps": 0.05,
    "temperature_mobile_K_initial": initial[3], "temperature_mobile_K_final": final[3],
    "potential_energy_eV_initial": initial[4], "potential_energy_eV_final": final[4],
    "potential_energy_delta_eV": final[4] - initial[4],
    "kinetic_energy_eV_initial": initial[5], "kinetic_energy_eV_final": final[5],
    "kinetic_energy_delta_eV": final[5] - initial[5],
    "total_energy_eV_initial": initial[6], "total_energy_eV_final": final[6],
    "total_energy_delta_eV": final[6] - initial[6],
    "Fmax_eV_A_final": float(fmag.max()), "F99_eV_A_final": float(np.quantile(fmag, .99)),
    "Fmax_by_element_eV_A_final": {type_name[int(t)]: float(fmag[final_type == t].max()) for t in np.unique(final_type)},
    "atoms_force_above_500_final": int((fmag > 500).sum()),
    "max_displacement_A_from_initial": float(disp.max()),
    "displacement_p99_A_from_initial": float(np.quantile(disp, .99)),
    "min_TiC_distance_A_final": tic_min, "min_TiC_pair_ids_final": [tic_ti, tic_c],
    "TiC_pairs_below_1p8A_final": tic_short,
    "min_CC_distance_A_final": cc_dist, "min_CC_pair_ids_final": [cc_id_a, cc_id_b],
    "tracked_v12m_collapse_pair_ids": [128985, 134363],
    "tracked_CC_pair_distance_A_final": cc_target_dist,
    "role": "no-laser 50 fs cold stability screen; not heating approval",
}
OUT.write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
