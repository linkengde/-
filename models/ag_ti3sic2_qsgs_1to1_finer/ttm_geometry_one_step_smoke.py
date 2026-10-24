#!/usr/bin/env python3
"""One reference-TTM step on the actual Ag geometry with synthetic 300 K v.

No LAMMPS, interatomic potential, or MD step is run. This checks full-size
atom-to-grid mapping, one 0.5 fs source/coupling update, and the equal/opposite
kinetic-energy handoff on reproducible Maxwell velocities.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from ag_ttm_reference import (
    AgParameters,
    AgTTMGrid,
    KB_EV_K,
    M_VV2E_METAL,
    build_ag_grid_fields,
    lattice_temperature_by_cell,
    rescale_lattice_kinetic,
)
from preflight_ttm_grid import (
    AG_LATTICE_A,
    LX,
    LY,
    LZ,
    power_density,
    read_newton_ag,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spacing", type=float, default=10.0, help="nominal grid spacing in Angstrom")
    parser.add_argument("--dt-fs", type=float, default=0.5, help="MD-equivalent operator step in femtoseconds")
    parser.add_argument("--seed", type=int, default=100)
    args = parser.parse_args()

    model_dir = Path(__file__).resolve().parent
    xyz = read_newton_ag(model_dir / "ag_ti3sic2_qsgs.data.gz", LZ, 40.0)
    q_atom = power_density(xyz, LZ)
    lengths = np.asarray((LX, LY, LZ), dtype=float)
    grid_shape = np.ceil(lengths / args.spacing).astype(int)
    spacing = lengths / grid_shape
    ag_volume = AG_LATTICE_A**3 / 4.0
    fields = build_ag_grid_fields(
        xyz, q_atom, tuple(grid_shape), tuple(spacing), ag_volume
    )
    parameters = AgParameters.from_csv(
        model_dir / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv",
        allow_low_temperature_extension=True,
    )

    rng = np.random.default_rng(args.seed)
    mass = np.full(len(xyz), 107.8682)
    sigma = math.sqrt(KB_EV_K * 300.0 / (mass[0] * M_VV2E_METAL))
    velocity = rng.normal(0.0, sigma, size=(len(xyz), 3))
    com = np.average(velocity, axis=0, weights=mass)
    velocity -= com
    kinetic = 0.5 * M_VV2E_METAL * np.sum(mass[:, None] * velocity**2)
    target_kinetic = 1.5 * (len(xyz) - 1) * KB_EV_K * 300.0
    velocity *= math.sqrt(target_kinetic / kinetic)

    lattice_temperature = lattice_temperature_by_cell(
        velocity, mass, fields.atom_map, smoothing_radius_cells=1
    )
    grid = AgTTMGrid(
        parameters,
        fields.fraction_ag,
        tuple(spacing),
        fields.source_power_eV_ps,
    )
    electron_before = float(grid.electron_energy_eV.sum())
    kinetic_before = 0.5 * M_VV2E_METAL * float(np.sum(mass[:, None] * velocity**2))
    accounting = grid.advance(args.dt_fs * 1.0e-3, lattice_temperature, heat_on=True)
    velocity_after, mapped_lattice_energy = rescale_lattice_kinetic(
        velocity,
        mass,
        fields.atom_map,
        accounting.lattice_energy_per_cell_eV,
    )
    electron_after = float(grid.electron_energy_eV.sum())
    kinetic_after = 0.5 * M_VV2E_METAL * float(np.sum(mass[:, None] * velocity_after**2))
    total_residual = (electron_after + kinetic_after) - (electron_before + kinetic_before)
    total_residual -= accounting.source_to_electron_eV

    print(f"MD steps: 0; LAMMPS calls: 0; synthetic velocity seed: {args.seed}")
    print(f"grid: {grid_shape[0]}x{grid_shape[1]}x{grid_shape[2]}, spacing_A={spacing}")
    print(f"active Ag cells: {np.count_nonzero(fields.fraction_ag)} / {fields.fraction_ag.size}")
    print(f"Newton-layer Ag atoms: {len(xyz)}; source power: {fields.source_power_eV_ps.sum():.9f} eV/ps")
    print(f"operator step: {accounting.elapsed_ps:.9g} ps; electron substeps: {accounting.substeps}")
    print(
        f"synthetic local lattice T range: "
        f"{np.nanmin(lattice_temperature[fields.fraction_ag > 0.0]):.3f}–"
        f"{np.nanmax(lattice_temperature[fields.fraction_ag > 0.0]):.3f} K"
    )
    print(f"electron T range: {grid.temperature_K[fields.fraction_ag > 0.0].min():.6f}–{grid.temperature_K.max():.6f} K")
    print(f"source energy: {accounting.source_to_electron_eV:.12g} eV")
    print(f"electron-to-lattice energy: {accounting.electron_to_lattice_eV:.12g} eV")
    print(f"mapped lattice kinetic-energy change: {mapped_lattice_energy.sum():.12g} eV")
    print(f"TTM ledger residual: {accounting.conservation_residual_eV:.3e} eV")
    print(f"electron+lattice residual after source subtraction: {total_residual:.3e} eV")
    if abs(total_residual) > max(abs(accounting.source_to_electron_eV) * 1.0e-10, 1.0e-12):
        raise SystemExit("Energy handoff did not pass the one-step smoke threshold")


if __name__ == "__main__":
    main()
