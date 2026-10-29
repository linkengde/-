#!/usr/bin/env python3
"""Export one Ag70 continuum checkpoint as a LAMMPS-style voxel-particle dump.

Each record represents one native grid voxel. These are visualization glyphs,
not atomistic coordinates or a molecular-dynamics trajectory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


BASE = Path(__file__).resolve().parent
DEFAULT_GEOMETRY = BASE.parent / "qsgs_morphology_ag70vol.npz"
DEFAULT_CHECKPOINT = BASE / "checkpoints" / "caseA_ag70_1x_heat_18p636ps.npz"
DEFAULT_OUTPUT = BASE / "ag70_voxel_pseudoparticles_18p636ps.lammpstrj"
PHASE_NAMES = {
    1: "Ti3SiC2_scaffold",
    2: "Ag_solid",
    3: "Ag_partial_visual_bin_liquid_fraction_0_to_0p5",
    4: "Ag_liquid_visual_bin_liquid_fraction_ge_0p5",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    with np.load(args.geometry, allow_pickle=False) as data:
        scaffold = np.asarray(data["skeleton"], dtype=bool)
        size_A = np.asarray(data["size_A"], dtype=np.float64)
    with np.load(args.checkpoint, allow_pickle=False) as data:
        temperature = np.asarray(data["temperature_K"], dtype=np.float64)
        liquid_fraction = np.asarray(data["ag_liquid_fraction"], dtype=np.float64)
        metadata = json.loads(str(data["metadata_json"].item()))

    if scaffold.shape != temperature.shape or scaffold.shape != liquid_fraction.shape:
        raise ValueError(
            f"Geometry/checkpoint shape mismatch: {scaffold.shape}, "
            f"{temperature.shape}, {liquid_fraction.shape}"
        )
    if not np.all(np.isfinite(temperature)) or not np.all(np.isfinite(liquid_fraction)):
        raise ValueError("Checkpoint contains non-finite temperature or liquid fraction")
    if np.any((liquid_fraction < 0.0) | (liquid_fraction > 1.0)):
        raise ValueError("Ag liquid fraction falls outside [0, 1]")

    time_ps = float(metadata["time_ps"])
    spacing_A = size_A / np.asarray(scaffold.shape, dtype=np.float64)
    ag = ~scaffold
    # Match the categories used by render_atom_style_phase_3d.py.
    phase = np.ones(scaffold.shape, dtype=np.uint8)
    phase[ag & (liquid_fraction <= 1e-12)] = 2
    phase[ag & (liquid_fraction > 1e-12) & (liquid_fraction < 0.5)] = 3
    phase[ag & (liquid_fraction >= 0.5)] = 4

    flat_phase = phase.ravel(order="C")
    flat_temperature = temperature.ravel(order="C")
    flat_liquid = liquid_fraction.ravel(order="C")
    nx, ny, nz = scaffold.shape
    n_atoms = int(scaffold.size)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    step_index = int(round(time_ps * 1000.0))  # frame index expressed in fs

    with args.output.open("w", encoding="ascii", newline="\n") as out:
        out.write("ITEM: TIMESTEP\n")
        out.write(f"{step_index}\n")
        out.write("ITEM: NUMBER OF ATOMS\n")
        out.write(f"{n_atoms}\n")
        out.write("ITEM: BOX BOUNDS pp pp ff\n")
        for length in size_A:
            out.write(f"0.0000000000 {length:.10f}\n")
        out.write("ITEM: ATOMS id type x y z temperature liquid_fraction\n")

        # Stream chunks to keep memory low; coordinates are native voxel centres.
        chunk_size = 50_000
        for start in range(0, n_atoms, chunk_size):
            stop = min(start + chunk_size, n_atoms)
            ids = np.arange(start + 1, stop + 1, dtype=np.int64)
            linear = np.arange(start, stop, dtype=np.int64)
            iz = linear % nz
            iy = (linear // nz) % ny
            ix = linear // (ny * nz)
            x = (ix + 0.5) * spacing_A[0]
            y = (iy + 0.5) * spacing_A[1]
            z = (iz + 0.5) * spacing_A[2]
            block = np.column_stack((
                ids,
                flat_phase[start:stop],
                x,
                y,
                z,
                flat_temperature[start:stop],
                flat_liquid[start:stop],
            ))
            np.savetxt(
                out,
                block,
                fmt=("%d", "%d", "%.8f", "%.8f", "%.8f", "%.6f", "%.8f"),
                delimiter=" ",
            )

    counts = {str(code): int(np.count_nonzero(phase == code)) for code in PHASE_NAMES}
    manifest = {
        "format": "LAMMPS custom dump (one pseudo-particle per native continuum voxel)",
        "not_atomistic": True,
        "time_ps": time_ps,
        "dump_step_index_fs": step_index,
        "shape_xyz": [nx, ny, nz],
        "voxel_count": n_atoms,
        "box_size_A": size_A.tolist(),
        "voxel_spacing_A": spacing_A.tolist(),
        "boundary_labels": ["periodic", "periodic", "fixed_nonperiodic"],
        "phase_type_mapping": {str(k): v for k, v in PHASE_NAMES.items()},
        "phase_counts": counts,
        "liquid_visual_threshold": "Ag liquid_fraction >= 0.5; matches sphere render",
        "temperature_column_unit": "K",
        "coordinate_unit": "angstrom",
        "geometry_file": str(args.geometry.name),
        "geometry_sha256": sha256(args.geometry),
        "checkpoint_file": str(args.checkpoint.name),
        "checkpoint_sha256": sha256(args.checkpoint),
        "output_file": args.output.name,
        "output_sha256": sha256(args.output),
        "output_bytes": args.output.stat().st_size,
    }
    manifest_path = args.output.with_suffix(args.output.suffix + ".json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
