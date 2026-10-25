#!/usr/bin/env python3
"""Short real-MD energy-exchange test for the pure-Ag TTM reference operator.

This is a numerical bridge check, not a melt-pool simulation. A 4,000-atom
Ag EAM slab with p/p/open boundaries is advanced with LAMMPS NVE; the reference
electron grid receives the established lateral Gaussian/depth exponential
source and exchanges energy with the atoms. A matched no-TTM NVE run gives the
integration-drift baseline.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
sys.path.extend((
    "/workspace/.venvs/lammps/lib/python3.12/site-packages",
    str(MODEL),
))

from mpi4py import MPI
from lammps import lammps

from ag_ttm_reference import (
    AgParameters,
    AgTTMGrid,
    M_VV2E_METAL,
    cloud_in_cell_atom_map,
    lattice_temperature_by_cell,
    rescale_lattice_kinetic,
)


def make_lammps(log_path: Path, seed: int):
    lmp = lammps(cmdargs=["-screen", "none", "-log", str(log_path), "-sf", "omp"], comm=MPI.COMM_WORLD)
    commands = [
        "units metal",
        "atom_style atomic",
        "boundary p p f",
        "atom_modify sort 0 0",
        "lattice fcc 4.09",
        "region simbox block 0 10 0 10 -5 15 units lattice",
        "create_box 1 simbox",
        "region slab block 0 10 0 10 0 9.999 units lattice",
        "create_atoms 1 region slab",
        "mass 1 107.8682",
        "pair_style eam",
        f"pair_coeff * * {MODEL / 'potentials/Ag_u3.eam'}",
        "neighbor 2.0 bin",
        "neigh_modify delay 0 every 1 check yes",
        "timestep 0.0005",
        f"velocity all create 300.0 {seed} mom yes rot no dist gaussian",
        # The free surfaces need a short thermal relaxation before they can
        # serve as a clean energy-exchange reference.
        "fix equilibrate all nvt temp 300.0 300.0 0.1",
        "thermo 1000",
        "run 2000",
        "unfix equilibrate",
        "reset_timestep 0",
        "fix integrate all nve",
        "thermo 20",
        "thermo_style custom step time temp pe ke etotal",
        "thermo_modify flush yes",
        "run 0",
    ]
    for command in commands:
        lmp.command(command)
    return lmp


def atom_energy_and_temperature(lmp, velocity: np.ndarray, masses: np.ndarray):
    ke = 0.5 * M_VV2E_METAL * float(np.sum(masses[:, None] * velocity**2))
    pe = float(lmp.get_thermo("pe"))
    temperature = 2.0 * ke / (3.0 * len(masses) * 8.617333262145e-5)
    return pe + ke, pe, ke, temperature


def main() -> None:
    steps = int(os.environ.get("PURE_AG_TTM_STEPS", "100"))
    nominal_spacing_A = float(os.environ.get("PURE_AG_TTM_GRID_A", "10"))
    source_scale = float(os.environ.get("PURE_AG_TTM_SOURCE_SCALE", "1"))
    if steps <= 0 or nominal_spacing_A <= 0.0 or not np.isfinite(source_scale) or source_scale < 0.0:
        raise ValueError("Invalid pilot step count, source scale, or grid spacing")

    n_atoms = 4000
    seed = 190427
    dt_ps = 0.0005
    r_A = 30.0
    l_skin_A = 20.0
    i0_eV_ps_A2 = 5.0 * source_scale
    mass = np.full(n_atoms, 107.8682)
    grid_tag = f"{nominal_spacing_A:g}A"
    no_ttm_log = HERE / f"log.pure_ag_nve_reference_{grid_tag}"
    ttm_log = HERE / f"log.pure_ag_ttm_md_pilot_{grid_tag}"

    baseline = make_lammps(no_ttm_log, seed)
    baseline_start = float(baseline.get_thermo("etotal"))
    baseline.command(f"run {steps}")
    baseline_final = float(baseline.get_thermo("etotal"))
    nve_drift = baseline_final - baseline_start
    baseline.close()

    lmp = make_lammps(ttm_log, seed)
    nlocal = int(lmp.extract_setting("nlocal"))
    if nlocal != n_atoms:
        raise RuntimeError(f"Expected one MPI rank with {n_atoms} local atoms, got {nlocal}")
    x = lmp.numpy.extract_atom("x")
    v = lmp.numpy.extract_atom("v")
    xyz = np.array(x[:nlocal], copy=True)
    length = np.full(3, 10.0 * 4.09)
    grid_shape = np.ceil(length / nominal_spacing_A).astype(int)
    spacing = length / grid_shape
    dx = (xyz[:, 0] - length[0] / 2.0 + length[0] / 2.0) % length[0] - length[0] / 2.0
    dy = (xyz[:, 1] - length[1] / 2.0 + length[1] / 2.0) % length[1] - length[1] / 2.0
    q_atom = (i0_eV_ps_A2 / l_skin_A) * np.exp(
        -(dx * dx + dy * dy) / (r_A * r_A) - (length[2] - xyz[:, 2]) / l_skin_A
    )
    grid_shape_tuple = tuple(int(x) for x in grid_shape)
    spacing_tuple = tuple(float(x) for x in spacing)
    atom_map = cloud_in_cell_atom_map(xyz, grid_shape_tuple, spacing_tuple)
    # This is a dense, fully occupied pure-Ag reference slab. The QSGS helper
    # computes a bounded atom-volume fill field and is intended for porous
    # composite geometry; applying it to a nearly close-packed FCC solid
    # overfills individual cells. Here the material fraction is exactly one,
    # while the source is conservatively mapped from per-atom quadrature.
    fraction_ag = np.ones(grid_shape_tuple, dtype=float)
    atomic_volume_A3 = 4.09**3 / 4.0
    source_power = np.bincount(
        atom_map.cell_ids.ravel(),
        weights=np.repeat(atomic_volume_A3 * q_atom, atom_map.cell_ids.shape[1])
        * atom_map.weights.ravel(),
        minlength=int(np.prod(grid_shape_tuple)),
    ).reshape(grid_shape_tuple)
    expected_source_power = float(np.sum(atomic_volume_A3 * q_atom))
    if not np.isclose(source_power.sum(), expected_source_power, rtol=1.0e-12, atol=1.0e-12):
        raise RuntimeError("Pure-Ag CIC source mapping is not conservative")
    parameters = AgParameters.from_csv(
        MODEL / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv",
        allow_low_temperature_extension=True,
    )
    grid = AgTTMGrid(
        parameters,
        fraction_ag,
        tuple(float(x) for x in spacing),
        source_power,
        initial_temperature_K=300.0,
    )
    v_view = v[:nlocal]
    atom_energy_0, pe_0, ke_0, temp_0 = atom_energy_and_temperature(lmp, v_view, mass)
    electron_energy_0 = float(grid.electron_energy_eV.sum())
    history = []
    cumulative_transfer = 0.0
    cumulative_source = 0.0
    cumulative_ttm_residual = 0.0
    cumulative_electron_substeps = 0
    max_substeps = 0

    for step in range(1, steps + 1):
        lmp.command("run 1 pre no post no")
        v_md = np.array(v_view, copy=True)
        lattice_temperature = lattice_temperature_by_cell(
            v_md, mass, atom_map, smoothing_radius_cells=1
        )
        accounting = grid.advance(dt_ps, lattice_temperature, heat_on=True)
        cumulative_source += accounting.source_to_electron_eV
        cumulative_ttm_residual += accounting.conservation_residual_eV
        cumulative_electron_substeps += accounting.substeps
        v_new, applied = rescale_lattice_kinetic(
            v_md, mass, atom_map, accounting.lattice_energy_per_cell_eV
        )
        v_view[:] = v_new
        cumulative_transfer += float(applied.sum())
        max_substeps = max(max_substeps, accounting.substeps)

        if step == 1 or step % 20 == 0 or step == steps:
            atom_energy, pe, ke, temp = atom_energy_and_temperature(lmp, v_view, mass)
            combined_residual = (
                atom_energy + float(grid.electron_energy_eV.sum())
                - atom_energy_0 - electron_energy_0 - cumulative_source
            )
            history.append({
                "step": step,
                "time_ps": step * dt_ps,
                "atom_temperature_K": temp,
                "electron_temperature_mean_K": float(np.average(
                    grid.temperature_K[fraction_ag > 0.0],
                    weights=fraction_ag[fraction_ag > 0.0],
                )),
                "electron_temperature_max_K": float(grid.temperature_K.max()),
                "atom_energy_eV": atom_energy,
                "pe_eV": pe,
                "ke_eV": ke,
                "cumulative_electron_to_lattice_eV": cumulative_transfer,
                "cumulative_source_to_electron_eV": cumulative_source,
                "combined_energy_residual_eV": combined_residual,
                "ttm_grid_energy_residual_eV": cumulative_ttm_residual,
                "max_electron_substeps_per_md_step": max_substeps,
            })

    atom_energy_final, pe_final, ke_final, temp_final = atom_energy_and_temperature(lmp, v_view, mass)
    electron_energy_final = float(grid.electron_energy_eV.sum())
    lmp.command("run 0 pre no post no")
    combined_residual = (
        atom_energy_final + electron_energy_final - atom_energy_0
        - electron_energy_0 - cumulative_source
    )
    excess_over_nve_drift = combined_residual - nve_drift
    energy_check_passed = (
        abs(excess_over_nve_drift) <= 1.0e-3
        and abs(cumulative_ttm_residual) <= 1.0e-9
    )
    result = {
        "status": "PURE_AG_EAM_REAL_MD_TTM_REFERENCE_COUPLING_PILOT",
        "purpose": "Validate source-to-electron and electron-to-lattice energy exchange on actual LAMMPS pure-Ag MD steps; not a melt-pool run.",
        "atom_model": "4000-atom FCC Ag slab; p/p/open boundaries; Ag_u3.eam; 0.5 fs NVE",
        "electron_model": "Ag-only grid initialized at 300 K; candidate Ag Ce/G/ke table; lateral exp(-r^2/R^2), depth exp(-d/l_skin) source; opt-in low-temperature closure enabled",
        "steps": steps,
        "duration_ps": steps * dt_ps,
        "pre_equilibration_ps": 1.0,
        "source_scale": source_scale,
        "electron_grid_shape": list(grid_shape_tuple),
        "initial_atom_temperature_K": temp_0,
        "initial_electron_temperature_K": 300.0,
        "final_atom_temperature_K": temp_final,
        "final_electron_temperature_mean_K": float(np.average(
            grid.temperature_K[fraction_ag > 0.0],
            weights=fraction_ag[fraction_ag > 0.0],
        )),
        "final_electron_temperature_max_K": float(grid.temperature_K.max()),
        "initial_atom_energy_eV": atom_energy_0,
        "final_atom_energy_eV": atom_energy_final,
        "initial_pe_eV": pe_0,
        "initial_ke_eV": ke_0,
        "final_pe_eV": pe_final,
        "final_ke_eV": ke_final,
        "initial_electron_energy_eV": electron_energy_0,
        "final_electron_energy_eV": electron_energy_final,
        "electron_to_lattice_energy_eV": cumulative_transfer,
        "source_peak_flux_eV_ps_A2": i0_eV_ps_A2,
        "source_R_A": r_A,
        "source_l_skin_A": l_skin_A,
        "source_power_eV_ps": float(source_power.sum()),
        "gross_combined_energy_change_eV": atom_energy_final + electron_energy_final - atom_energy_0 - electron_energy_0,
        "source_energy_eV": cumulative_source,
        "combined_energy_residual_eV": combined_residual,
        "ttm_grid_energy_residual_eV": cumulative_ttm_residual,
        "matched_no_ttm_nve_drift_eV": nve_drift,
        "residual_minus_no_ttm_drift_eV": excess_over_nve_drift,
        "energy_exchange_check_passed": energy_check_passed,
        "energy_exchange_check_tolerance_eV": 1.0e-3,
        "max_electron_substeps_per_md_step": max_substeps,
        "total_electron_substeps": cumulative_electron_substeps,
        "smallest_observed_electron_substep_fs": dt_ps * 1000.0 / max_substeps,
        "energy_history": history,
        "no_external_source": False,
        "no_production_or_composite_claim": True,
        "logs": {"ttm": ttm_log.name, "matched_nve": no_ttm_log.name},
        "limitations": [
            "The 40.9-A slab has two clean Ag free surfaces and a small periodic cross-section; it does not test QSGS confinement.",
            "The coupled solver remains a Python reference operator, not a compiled LAMMPS TTM fix.",
            "The electron parameters are candidate cross-paper fits; this short sourced run only checks numerical integration and energy handoff.",
        ],
    }
    (HERE / f"pure_ag_ttm_md_pilot_metrics_{grid_tag}.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)
    lmp.close()
    MPI.Finalize()


if __name__ == "__main__":
    main()
