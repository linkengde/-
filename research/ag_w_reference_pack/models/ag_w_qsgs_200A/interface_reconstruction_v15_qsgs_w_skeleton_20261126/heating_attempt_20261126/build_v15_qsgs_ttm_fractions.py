#!/usr/bin/env python3
"""Average V15's binary QSGS phase mask onto 5/10/20 A Newton TTM grids."""
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
MASK_FILE = HERE.parent / "atomic_candidate" / "qsgs_morphology.npz"
OUTPUT = HERE / "v15_qsgs_phase_fractions_5A_10A_20A.npz"

with np.load(MASK_FILE) as src:
    w = np.asarray(src["W_backbone"], dtype=bool)
    box = np.asarray(src["material_size_A"], dtype=float)
    grid = np.asarray(src["grid"], dtype=int)

voxel = box / grid
fields = {}
for g in (5, 10, 20):
    block = g / voxel
    if not np.allclose(block, np.rint(block), atol=1.0e-12):
        raise RuntimeError(f"{g} A grid does not align with QSGS voxels ({voxel.tolist()} A)")
    nblock = np.rint(block).astype(int)
    start_z = int(round(40.0 / voxel[2]))
    stop_z = int(round(200.0 / voxel[2]))
    if not np.isclose(start_z * voxel[2], 40.0) or not np.isclose(stop_z * voxel[2], 200.0):
        raise RuntimeError("Newton layer edges do not align with QSGS voxel boundaries")
    wz = w[:, :, start_z:stop_z]
    nx, ny, nz = grid[0] // nblock[0], grid[1] // nblock[1], (stop_z - start_z) // nblock[2]
    # Block average of the binary voxel mask gives geometric solid fractions.
    fw = wz.reshape(nx, nblock[0], ny, nblock[1], nz, nblock[2]).mean(axis=(1, 3, 5), dtype=np.float64)
    fa = 1.0 - fw
    fields[f"fraction_w_{g}A"] = fw
    fields[f"fraction_ag_{g}A"] = fa
    print(f"{g} A grid {fw.shape}, W fraction {fw.mean():.8f}, max combined fill {(fw+fa).max():.8f}")

np.savez_compressed(OUTPUT, **fields)
print(OUTPUT)
