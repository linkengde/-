#!/usr/bin/env python3
"""Static Ag-mask/grid sensitivity check for the candidate TTM heat source.

This script does not generate a LAMMPS TTM field or run MD. It compares a
per-Ag-atom quadrature with finite-volume cell-centre occupancy estimates for
the current geometry. The bulk FCC Ag atomic volume is an explicit proxy.
"""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path

import numpy as np


LX = 299.635
LY = 301.857484616085
LZ = 158.6115
AG_LATTICE_A = 4.09
R_A = 30.0
LSKIN_A = 20.0
I0_EV_PS_A2 = 5.0
MOBILE_Z_MIN_A = 40.0


def read_newton_ag(path: Path, z_surface: float, z_min: float) -> np.ndarray:
    """Read type-4 atoms from the LAMMPS atomic-style data file."""
    opener = gzip.open if path.suffix == ".gz" else open
    xyz: list[tuple[float, float, float]] = []
    in_atoms = False
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped.startswith("Atoms"):
                in_atoms = True
                continue
            if not in_atoms:
                continue
            if not stripped:
                continue
            fields = stripped.split()
            if len(fields) < 5:
                continue
            try:
                atom_type = int(fields[1])
                x, y, z = map(float, fields[2:5])
            except ValueError:
                # Section headers after Atoms do not have numeric type fields.
                if fields[0] in {"Velocities", "Bonds", "Angles", "Masses"}:
                    break
                continue
            if atom_type == 4 and z_min <= z <= z_surface:
                xyz.append((x, y, z))
    if not xyz:
        raise ValueError(f"No Newton-layer Ag atoms found in {path}")
    return np.asarray(xyz, dtype=np.float64)


def gaussian_cell_average(lo: np.ndarray, hi: np.ndarray, center: float, width: float) -> np.ndarray:
    scale = width * math.sqrt(math.pi) / 2.0
    erf_hi = np.fromiter((math.erf((v - center) / width) for v in hi), dtype=np.float64)
    erf_lo = np.fromiter((math.erf((v - center) / width) for v in lo), dtype=np.float64)
    return scale * (erf_hi - erf_lo) / (hi - lo)


def power_density(xyz: np.ndarray, z_surface: float) -> np.ndarray:
    dx = (xyz[:, 0] - LX / 2.0 + LX / 2.0) % LX - LX / 2.0
    dy = (xyz[:, 1] - LY / 2.0 + LY / 2.0) % LY - LY / 2.0
    return (I0_EV_PS_A2 / LSKIN_A) * np.exp(
        -(dx * dx + dy * dy) / (R_A * R_A) - (z_surface - xyz[:, 2]) / LSKIN_A
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data",
        type=Path,
        default=Path(__file__).with_name("ag_ti3sic2_qsgs.data.gz"),
        help="LAMMPS atomic-style data file (plain text or .gz)",
    )
    parser.add_argument("--grid-candidates", type=float, nargs="+", default=[5, 8, 10, 12, 15, 20])
    parser.add_argument("--z-surface", type=float, default=LZ)
    parser.add_argument("--z-mobile-min", type=float, default=MOBILE_Z_MIN_A)
    args = parser.parse_args()

    xyz = read_newton_ag(args.data, args.z_surface, args.z_mobile_min)
    atom_volume = AG_LATTICE_A**3 / 4.0
    q_atom = power_density(xyz, args.z_surface)
    p_atom = atom_volume * float(q_atom.sum())
    p_continuum = (
        I0_EV_PS_A2
        * math.pi
        * R_A**2
        * (1.0 - math.exp(-args.z_surface / LSKIN_A))
    )

    print(f"input={args.data}")
    print(f"Newton-layer Ag atoms={len(xyz)}; Ag atomic-volume proxy={atom_volume:.6f} A^3")
    print(f"full-material continuum power={p_continuum:.6f} eV/ps")
    print(
        f"atom-quadrature Newton-Ag power={p_atom:.6f} eV/ps "
        f"({100.0 * p_atom / p_continuum:.4f}% of continuum)"
    )
    print("grid_A dimensions occupied_cells power_eV_per_ps error_vs_atom_pct fAg_gt1 max_fAg")

    for nominal_h in args.grid_candidates:
        lengths = np.asarray([LX, LY, args.z_surface], dtype=np.float64)
        n = np.ceil(lengths / nominal_h).astype(int)
        d = lengths / n
        ijk = np.minimum((xyz / d).astype(int), n - 1)
        key = (ijk[:, 0] * n[1] + ijk[:, 1]) * n[2] + ijk[:, 2]
        occupied, counts = np.unique(key, return_counts=True)
        ix = occupied // (n[1] * n[2])
        iy = (occupied // n[2]) % n[1]
        iz = occupied % n[2]

        xlo, xhi = ix * d[0], (ix + 1) * d[0]
        ylo, yhi = iy * d[1], (iy + 1) * d[1]
        zlo, zhi = iz * d[2], (iz + 1) * d[2]
        gx = gaussian_cell_average(xlo, xhi, LX / 2.0, R_A)
        gy = gaussian_cell_average(ylo, yhi, LY / 2.0, R_A)
        gz = LSKIN_A * (
            np.exp(-(args.z_surface - zhi) / LSKIN_A)
            - np.exp(-(args.z_surface - zlo) / LSKIN_A)
        ) / d[2]
        q_cell_average = (I0_EV_PS_A2 / LSKIN_A) * gx * gy * gz
        p_grid = atom_volume * float(np.dot(counts, q_cell_average))
        f_ag_proxy = counts * atom_volume / float(np.prod(d))
        max_f = float(f_ag_proxy.max())
        error_pct = 100.0 * (p_grid / p_atom - 1.0)
        dims = "x".join(str(v) for v in n)
        print(
            f"{nominal_h:g} {dims} {len(occupied)} {p_grid:.6f} "
            f"{error_pct:+.3f} {(f_ag_proxy > 1.0).sum()} {max_f:.3f}"
        )


if __name__ == "__main__":
    main()
