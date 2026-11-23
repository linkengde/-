#!/usr/bin/env python3
"""Regenerate the temporary LAMMPS table files used in overlap preconditioning."""
from pathlib import Path
import numpy as np


ROOT = Path(__file__).resolve().parent
N = 2401
R_MIN = 0.001
R_MAX = 2.4
STIFFNESS = 100.0


def write_table(path: Path, sections: list[tuple[str, float]]) -> None:
    radii = np.linspace(R_MIN, R_MAX, N)
    with path.open("w") as stream:
        stream.write("Temporary overlap wall functions; not physical pair potentials.\n")
        for name, cutoff in sections:
            stream.write(f"{name}\nN {N} R {R_MIN:.12g} {R_MAX:.12g}\n\n")
            for index, radius in enumerate(radii, 1):
                gap = max(cutoff - radius, 0.0)
                energy = 0.5 * STIFFNESS * gap**2
                force = STIFFNESS * gap if radius < cutoff else 0.0
                stream.write(f"{index} {radius:.12g} {energy:.12g} {force:.12g}\n")


if __name__ == "__main__":
    write_table(ROOT / "cross_wall.table", [("ZERO", 0.0), ("REPULSIVE", R_MAX)])
    write_table(ROOT / "cross_wall_v2.table", [("SELF", 2.2), ("REPULSIVE", R_MAX)])
