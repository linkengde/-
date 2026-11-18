"""Extract the residual Ti-C contact after the v12k SNAP 10 fs cold screen."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import write

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
PAIR = (276867, 276870)  # Ti-C, 1.38724 A after 10 fs with SNAP v12k
RADIUS, CELL = 5.5, 22.0
TYPE_TO_ELEMENT = {1: "Ag", 2: "C", 3: "Si", 4: "Ti"}
path = MODEL / "cold_snap_v12k_initial_458k_1fs/end_0p01ps.dump"
lines = path.read_text().splitlines()
h = next(i for i, row in enumerate(lines) if row.startswith("ITEM: ATOMS"))
cols = lines[h].split()[2:]
rows = np.loadtxt(lines[h + 1 :], dtype=float)
ix = {key: cols.index(key) for key in ("id", "type", "x", "y", "z", "fx", "fy", "fz")}
data = {int(row[ix["id"]]): row for row in rows}
if any(atom_id not in data for atom_id in PAIR):
    raise RuntimeError("Target pair is absent from v12k SNAP cold dump")
types = [int(data[atom_id][ix["type"]]) for atom_id in PAIR]
if types != [4, 2]:
    raise RuntimeError(f"Expected Ti-C type pair (4,2), got {types}")
pair_pos = np.asarray([[data[atom_id][ix["x"]], data[atom_id][ix["y"]], data[atom_id][ix["z"]]] for atom_id in PAIR])
delta = pair_pos[1] - pair_pos[0]
delta[:2] -= 200.0 * np.rint(delta[:2] / 200.0)
center = pair_pos[0] + delta / 2.0
ids = np.asarray(sorted(data), dtype=np.int64)
coords = np.asarray([[data[int(atom_id)][ix["x"]], data[int(atom_id)][ix["y"]], data[int(atom_id)][ix["z"]]] for atom_id in ids])
rel = coords - center
rel[:, :2] -= 200.0 * np.rint(rel[:, :2] / 200.0)
keep = np.linalg.norm(rel, axis=1) <= RADIUS
kept = ids[keep]
symbols = [TYPE_TO_ELEMENT[int(data[int(atom_id)][ix["type"]])] for atom_id in kept]
forces = np.asarray([[data[int(atom_id)][ix[key]] for key in ("fx", "fy", "fz")] for atom_id in kept])
target = np.isin(kept, PAIR).astype(np.int8)
atoms = Atoms(symbols=symbols, positions=rel[keep] + CELL / 2.0, cell=[CELL] * 3, pbc=True)
atoms.new_array("source_atom_id", kept.astype(np.int32))
atoms.new_array("SNAP_v12k_forces", forces)
atoms.new_array("is_DFT_force_target", target)
distance = float(np.linalg.norm(delta))
atoms.info.update({
    "config_type": "full200_v12k_tic_pair6_cold_10fs", "source_timestep": 10,
    "center_atom_id": PAIR[0], "center_species": "Ti", "partner_atom_id": PAIR[1], "partner_species": "C",
    "center_partner_distance_A": distance, "crop_radius_A": RADIUS, "crop_cell_A": CELL,
    "geometry_source": str(path), "source_pbc": "x/y periodic, z open; crop embedded in cubic vacuum box",
    "label_role": "v12l training candidate; Si01 center 104359 remains held out",
    "DFT_edge_force_note": "Only the center Ti-C pair forces will be used as DFT targets; crop edge rows excluded.",
    "SNAP_v12k_force_pair_eV_A": json.dumps(forces[target == 1].tolist()),
})
out = HERE / "candidates/tic_pair6_v12k_snap_cold_10fs.extxyz"
out.parent.mkdir(exist_ok=True)
write(out, atoms, format="extxyz")
print(json.dumps({"target_ids": PAIR, "elements": ["Ti", "C"], "pair_distance_A": distance, "crop_atoms": len(atoms), "DFT_target_rows": int(target.sum()), "output": str(out)}, indent=2))
