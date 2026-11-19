#!/workspace/venvs/dft-labels/bin/python
"""Extract the v12l POD Ag-C-C collapse precursor at relative step 18."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import write

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "cold_v12l_hotspot_trace/Ag_C_C_local_by_step.dump"
TARGET_IDS = [414900, 198354, 203868]  # Ag-C-C, in model type order
TYPE_TO_ELEMENT = {1: "Ag", 2: "C", 3: "Si", 4: "Ti"}
STEP = 18
RADIUS_A = 5.5
CELL_A = 22.0

lines = SOURCE.read_text().splitlines()
frame = None
i = 0
while i < len(lines):
    if not lines[i].startswith("ITEM: TIMESTEP"):
        i += 1
        continue
    step = int(lines[i + 1])
    n = int(lines[i + 3])
    h = i + 8
    cols = lines[h].split()[2:]
    rows = np.loadtxt(lines[h + 1 : h + 1 + n], dtype=float)
    if step == STEP:
        frame = (cols, rows)
        break
    i = h + 1 + n
if frame is None:
    raise RuntimeError(f"Could not find local trace step {STEP}")

cols, rows = frame
ix = {key: cols.index(key) for key in ("id", "type", "x", "y", "z", "fx", "fy", "fz")}
ids = rows[:, ix["id"]].astype(int)
types = rows[:, ix["type"]].astype(int)
xyz = rows[:, [ix[k] for k in ("x", "y", "z")]]
forces = rows[:, [ix[k] for k in ("fx", "fy", "fz")]]
lookup = {int(atom_id): j for j, atom_id in enumerate(ids)}
if any(atom_id not in lookup for atom_id in TARGET_IDS):
    raise RuntimeError("A target ID is absent from the saved local region")
expected = [1, 2, 2]
if [int(types[lookup[atom_id]]) for atom_id in TARGET_IDS] != expected:
    raise RuntimeError("Target element mapping is not Ag-C-C")

# Unwrap the three target atoms relative to the Ag center in periodic x/y.
anchor = xyz[lookup[TARGET_IDS[0]]]
unwrapped = []
for atom_id in TARGET_IDS:
    delta = xyz[lookup[atom_id]] - anchor
    delta[:2] -= 200.0 * np.rint(delta[:2] / 200.0)
    unwrapped.append(anchor + delta)
unwrapped = np.asarray(unwrapped)
center = unwrapped.mean(axis=0)
relative = xyz - center
relative[:, :2] -= 200.0 * np.rint(relative[:, :2] / 200.0)
keep = np.linalg.norm(relative, axis=1) <= RADIUS_A
kept_ids = ids[keep]
kept_types = types[keep]
symbols = [TYPE_TO_ELEMENT[int(t)] for t in kept_types]
target_mask = np.isin(kept_ids, TARGET_IDS).astype(np.int8)
if len(kept_ids) > 64 or not np.isin(TARGET_IDS, kept_ids).all() or int(target_mask.sum()) != 3:
    raise RuntimeError("Unexpected crop size or target-row count")

atoms = Atoms(
    symbols=symbols,
    positions=relative[keep] + CELL_A / 2.0,
    cell=[CELL_A] * 3,
    pbc=True,
)
atoms.new_array("source_atom_id", kept_ids.astype(np.int32))
atoms.new_array("POD_v12l_forces", forces[keep])
atoms.new_array("is_DFT_force_target", target_mask)
pairs = {
    f"{TARGET_IDS[a]}-{TARGET_IDS[b]}": float(np.linalg.norm(unwrapped[a] - unwrapped[b]))
    for a, b in ((0, 1), (0, 2), (1, 2))
}
atoms.info.update(
    {
        "config_type": "full200_v12l_POD_cold_Ag-C-C_collapse_precursor",
        "source_timestep": STEP,
        "cumulative_source_time_fs": 10 + STEP,
        "center_atom_ids": json.dumps(TARGET_IDS),
        "elements": json.dumps(["Ag", "C", "C"]),
        "center_pair_distances_A": json.dumps(pairs),
        "crop_radius_A": RADIUS_A,
        "crop_cell_A": CELL_A,
        "geometry_source": str(SOURCE),
        "source_pbc": "x/y periodic, z open; 5.5 A crop embedded in cubic vacuum box",
        "DFT_target_rows": 3,
        "label_role": "v12m training candidate; Si01 center 104359 remains held out",
    }
)
out = HERE / "candidates/agc_triad_v12l_t28fs_step18.extxyz"
out.parent.mkdir(exist_ok=True)
write(out, atoms, format="extxyz")
print(
    json.dumps(
        {
            "status": "EXTRACTED",
            "output": str(out),
            "natoms": len(atoms),
            "target_ids": TARGET_IDS,
            "target_rows": int(target_mask.sum()),
            "pair_distances_A": pairs,
        },
        indent=2,
    )
)
