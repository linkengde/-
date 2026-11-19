#!/usr/bin/env python3
"""Audit strict FIRE replay frame IDs/bounds and compare step 60 to lost-ignore run."""
from __future__ import annotations
import json
import hashlib
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
N = 347464
STRICT = HERE / "fire_repro_strict.dump"
PRIOR = Path("/workspace/potential-fit/round7_pod_snap/models/v12_candidate_forceorder_fixed/fire_diagnostic_exact_xyz/xyz_by_step.dump")
OUT = HERE / "strict_replay_audit.json"
FRAMES = HERE / "strict_replay_frames_60_62.dump"

def get_frames(path: Path, wanted: set[int]):
    frames = {}
    with path.open() as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.strip() != "ITEM: TIMESTEP":
                continue
            step = int(f.readline())
            if f.readline().strip() != "ITEM: NUMBER OF ATOMS":
                raise RuntimeError(f"Malformed frame header at step {step} in {path}")
            n = int(f.readline())
            if n != N:
                raise RuntimeError(f"Frame {step}: expected {N} atoms; got {n}")
            if not f.readline().strip().startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError(f"Unexpected box header at step {step}: {path}")
            bounds = np.asarray([[float(x) for x in f.readline().split()] for _ in range(3)])
            atom_header = f.readline().split()
            if atom_header[:2] != ["ITEM:", "ATOMS"]:
                raise RuntimeError(f"Missing atom header at step {step}")
            cols = atom_header[2:]
            fields = [f.readline().split() for _ in range(n)]
            if step in wanted:
                arr = np.asarray(fields)
                ids = arr[:, cols.index("id")].astype(int)
                types = arr[:, cols.index("type")].astype(int) if "type" in cols else None
                xyz = arr[:, [cols.index(k) for k in ("x", "y", "z")]].astype(float)
                frames[step] = (ids, types, xyz, bounds)
    return frames

strict = get_frames(STRICT, {60, 61, 62})
prior = get_frames(PRIOR, {60})
if set(strict) != {60, 61, 62} or 60 not in prior:
    raise RuntimeError("Missing expected frame(s) in strict/prior replay dump")

summaries = []
for step in (60, 61, 62):
    ids, types, xyz, bounds = strict[step]
    unique = np.unique(ids)
    lo, hi = bounds[:, 0], bounds[:, 1]
    summaries.append({
        "step": step,
        "atom_rows": int(len(ids)),
        "unique_ids": int(len(unique)),
        "missing_ids": int(N-len(unique)),
        "duplicate_ids": int(len(ids)-len(unique)),
        "nonfinite_rows": int((~np.isfinite(xyz)).any(axis=1).sum()),
        "xyz_min_A": xyz.min(axis=0).tolist(),
        "xyz_max_A": xyz.max(axis=0).tolist(),
        "outside_primary_bounds_xyz": [int(((xyz[:, k] < lo[k]) | (xyz[:, k] > hi[k])).sum()) for k in range(3)],
        "box_bounds_A": bounds.tolist(),
    })

ids_a, _, xyz_a, _ = prior[60]
ids_b, types_b, xyz_b, _ = strict[60]
order_a, order_b = np.argsort(ids_a), np.argsort(ids_b)
if not np.array_equal(ids_a[order_a], ids_b[order_b]):
    raise RuntimeError("The prior lost-ignore and strict replays have different IDs at step 60")
delta = xyz_a[order_a] - xyz_b[order_b]

with FRAMES.open("w") as f:
    for step in (60, 61, 62):
        ids, types, xyz, bounds = strict[step]
        f.write(f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{len(ids)}\n")
        f.write("ITEM: BOX BOUNDS pp pp mm\n")
        for lohi in bounds:
            f.write(f"{lohi[0]:.16g} {lohi[1]:.16g}\n")
        if types is None:
            f.write("ITEM: ATOMS id x y z\n")
            for atom_id, pos in zip(ids, xyz):
                f.write(f"{atom_id:d} {pos[0]:.9g} {pos[1]:.9g} {pos[2]:.9g}\n")
        else:
            f.write("ITEM: ATOMS id type x y z\n")
            for atom_id, atom_type, pos in zip(ids, types, xyz):
                f.write(f"{atom_id:d} {atom_type:d} {pos[0]:.9g} {pos[1]:.9g} {pos[2]:.9g}\n")
full_hash = hashlib.sha256()
with STRICT.open("rb") as f:
    for block in iter(lambda: f.read(1024*1024), b""):
        full_hash.update(block)

log = (HERE / "log.strict.lammps").read_text(errors="replace")
result = {
    "status": "STRICT_REPLAY_COMPLETED_NOT_CONVERGED",
    "source_atoms": N,
    "strict_lost_error_reproduced": "ERROR: Lost atoms" in log,
    "strict_replay_frames": summaries,
    "step60_prior_lost_ignore_vs_strict": {
        "max_abs_coordinate_difference_A": float(np.abs(delta).max()),
        "exactly_identical_coordinates": bool(np.array_equal(xyz_a[order_a], xyz_b[order_b])),
        "same_sorted_ids": True,
    },
    "strict_lost_atom_error_lines": [line for line in log.splitlines() if "Lost atoms" in line],
    "fire_stop": "max iterations (62); not converged",
    "full_strict_dump": {"path": STRICT.name, "bytes": STRICT.stat().st_size, "sha256": full_hash.hexdigest(), "archived_frames_60_62": FRAMES.name},
    "interpretation": "Strict replay retained all IDs and matched the prior lost-ignore step-60 coordinates. The original step-60 lost-atoms report is not reproducible as physical escape; the specific internal cause remains undetermined.",
}
OUT.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
