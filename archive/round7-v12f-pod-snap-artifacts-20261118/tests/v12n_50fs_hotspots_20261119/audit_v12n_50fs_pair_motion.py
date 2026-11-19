#!/workspace/venvs/dft-labels/bin/python
"""Audit radial motion and net-force tendency for selected 50 fs short contacts."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DUMP = HERE.parent / "cold_v12n_continue_to_50fs" / "end_0p05ps.dump"
PAIRS = {
    "Ti-C": [94075, 99595],
    "C-C": [38346, 33016],
}
MASS = {"Ag": 107.8682, "C": 12.011, "Si": 28.085, "Ti": 47.867}
SYMBOL = {1: "Ag", 2: "C", 3: "Si", 4: "Ti"}

lines = DUMP.read_text().splitlines()
h = next(i for i, row in enumerate(lines) if row.startswith("ITEM: ATOMS "))
n = int(lines[3])
columns = lines[h].split()[2:]
fields = ("id", "type", "x", "y", "z", "vx", "vy", "vz", "fx", "fy", "fz")
idx = {name: columns.index(name) for name in fields}
rows = np.loadtxt(lines[h + 1:h + 1 + n], dtype=float, ndmin=2)
ids = rows[:, idx["id"]].astype(np.int64)
types = rows[:, idx["type"]].astype(np.int32)
lookup = {int(atom_id): i for i, atom_id in enumerate(ids)}
xyz = rows[:, [idx[k] for k in ("x", "y", "z")]]
vel = rows[:, [idx[k] for k in ("vx", "vy", "vz")]]
force = rows[:, [idx[k] for k in ("fx", "fy", "fz")]]
results = []
for name, atom_ids in PAIRS.items():
    i, j = [lookup[a] for a in atom_ids]
    dr = xyz[j] - xyz[i]
    dr[:2] -= 200.0 * np.rint(dr[:2] / 200.0)
    distance = float(np.linalg.norm(dr))
    unit = dr / distance
    dv = vel[j] - vel[i]
    mi, mj = MASS[SYMBOL[int(types[i])]], MASS[SYMBOL[int(types[j])]]
    relative_net_force_per_mass = force[j] / mj - force[i] / mi
    results.append({
        "pair": name,
        "atom_ids": atom_ids,
        "distance_A": distance,
        "radial_relative_velocity_A_per_ps": float(np.dot(dv, unit)),
        "relative_net_force_projection_eV_A_per_amu": float(np.dot(relative_net_force_per_mass, unit)),
        "interpretation": "positive radial velocity and projection mean the pair is moving/accelerating apart",
        "note": "This uses full atomic net forces; it is not a pairwise bond-force decomposition.",
    })
out = {
    "model": "POD v12n",
    "source_dump_timestep": 40,
    "source_physical_time_fs": 50,
    "source_units": "LAMMPS metal; velocities in Angstrom/ps",
    "boundary": "x/y periodic, z nonperiodic",
    "pairs": results,
}
(HERE / "pair_motion_audit.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
