#!/workspace/venvs/dft-labels/bin/python
"""Extract short-distance C-C and Ti-C environments from the POD v12n 50 fs state."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import write

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "cold_v12n_continue_to_50fs" / "end_0p05ps.dump"
BOX_X = BOX_Y = 200.0
RADIUS_A = 5.5
CELL_A = 22.0
PAIRS = {
    "tic_94075_99595": ([94075, 99595], [4, 2]),  # dump type: Ti=4, C=2
    "cc_38346_33016": ([38346, 33016], [2, 2]),
}
SYMBOL = {1: "Ag", 2: "C", 3: "Si", 4: "Ti"}


def read_dump(path: Path):
    lines = path.read_text().splitlines()
    h = next(i for i, row in enumerate(lines) if row.startswith("ITEM: ATOMS "))
    natoms = int(lines[3])
    columns = lines[h].split()[2:]
    names = ("id", "type", "x", "y", "z", "fx", "fy", "fz")
    index = {name: columns.index(name) for name in names}
    rows = np.loadtxt(lines[h + 1:h + 1 + natoms], dtype=float, ndmin=2)
    return {name: rows[:, index[name]] for name in names}


def main():
    data = read_dump(SOURCE)
    ids = data["id"].astype(np.int64)
    types = data["type"].astype(np.int32)
    xyz = np.column_stack((data["x"], data["y"], data["z"]))
    forces = np.column_stack((data["fx"], data["fy"], data["fz"]))
    by_id = {int(atom_id): i for i, atom_id in enumerate(ids)}
    if len(by_id) != len(ids):
        raise RuntimeError("Duplicate atom IDs in source dump")
    results = []
    for name, (target_ids, expected_types) in PAIRS.items():
        if any(atom_id not in by_id for atom_id in target_ids):
            raise RuntimeError(f"Missing target ID in source: {target_ids}")
        target_indices = [by_id[atom_id] for atom_id in target_ids]
        actual_types = [int(types[j]) for j in target_indices]
        if actual_types != expected_types:
            raise RuntimeError(f"Type mismatch for {target_ids}: {actual_types}")

        anchor = xyz[target_indices[0]].copy()
        target_xyz = []
        for j in target_indices:
            d = xyz[j] - anchor
            d[0] -= BOX_X * np.rint(d[0] / BOX_X)
            d[1] -= BOX_Y * np.rint(d[1] / BOX_Y)
            target_xyz.append(anchor + d)
        target_xyz = np.asarray(target_xyz)
        center = target_xyz.mean(axis=0)
        relative = xyz - center
        relative[:, 0] -= BOX_X * np.rint(relative[:, 0] / BOX_X)
        relative[:, 1] -= BOX_Y * np.rint(relative[:, 1] / BOX_Y)
        keep = np.linalg.norm(relative, axis=1) <= RADIUS_A
        kept_ids = ids[keep]
        if len(kept_ids) > 96 or not np.isin(target_ids, kept_ids).all():
            raise RuntimeError(f"Unexpected crop size/target omission for {name}: {keep.sum()}")

        target_mask = np.isin(kept_ids, target_ids).astype(np.int8)
        atoms = Atoms(
            symbols=[SYMBOL[int(t)] for t in types[keep]],
            positions=relative[keep] + CELL_A / 2,
            cell=[CELL_A] * 3,
            pbc=True,
        )
        atoms.new_array("source_atom_id", kept_ids.astype(np.int32))
        atoms.new_array("POD_v12n_forces", forces[keep])
        atoms.new_array("is_DFT_force_target", target_mask)
        pair_distance = float(np.linalg.norm(target_xyz[0] - target_xyz[1]))
        atoms.info.update({
            "config_type": f"full200_v12n_POD_cold_50fs_{name}",
            "source_timestep": 40,
            "source_time_fs": 50,
            "center_atom_ids": json.dumps(target_ids),
            "elements": json.dumps([SYMBOL[t] for t in expected_types]),
            "center_pair_distance_A": pair_distance,
            "crop_radius_A": RADIUS_A,
            "crop_cell_A": CELL_A,
            "geometry_source": str(SOURCE),
            "source_pbc": "x/y periodic, z open; local crop embedded in cubic vacuum box",
            "DFT_target_rows": len(target_ids),
            "label_role": "v12o training candidate; not a validation or holdout configuration",
        })
        out = HERE / "candidates" / f"{name}_t50fs.extxyz"
        out.parent.mkdir(parents=True, exist_ok=True)
        write(out, atoms, format="extxyz")
        record = {
            "candidate": name,
            "output": str(out),
            "natoms": len(atoms),
            "target_ids": target_ids,
            "target_types": [SYMBOL[t] for t in expected_types],
            "pair_distance_A": pair_distance,
            "radius_A": RADIUS_A,
            "source": str(SOURCE),
            "source_time_fs": 50,
            "max_model_force_eV_A": float(np.linalg.norm(forces[keep], axis=1).max()),
            "target_model_forces_eV_A": {
                str(atom_id): forces[by_id[atom_id]].tolist() for atom_id in target_ids
            },
        }
        results.append(record)
    manifest = {"status": "EXTRACTED", "model": "POD v12n", "pairs": results}
    (HERE / "extraction_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
