#!/usr/bin/env python3
"""One bounded 0.25 ps no-source structural diagnostic for Model B."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.extend([
    "/workspace/.venvs/lammps/lib/python3.12/site-packages",
    "/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer",
])
from mpi4py import MPI
from lammps import lammps
from ag_ttm_reference import KB_EV_K, M_VV2E_METAL


ROOT = Path("/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer")
BDIR = ROOT / "potential_validation/geometry_constraint_B"
CANDIDATE = ROOT / "potential_validation/relaxed_backbone_candidate"


def ag_kinetic_state(lmp: lammps) -> tuple[float, float]:
    nlocal = int(lmp.extract_setting("nlocal"))
    atom_type = np.asarray(lmp.numpy.extract_atom("type")[:nlocal])
    xyz = np.asarray(lmp.numpy.extract_atom("x")[:nlocal])
    velocity = np.asarray(lmp.numpy.extract_atom("v")[:nlocal])
    force = np.asarray(lmp.numpy.extract_atom("f")[:nlocal])
    mobile_ag = (atom_type == 4) & (xyz[:, 2] >= 20.0)
    kinetic_eV = 0.5 * M_VV2E_METAL * 107.8682 * np.einsum(
        "ij,ij->", velocity[mobile_ag], velocity[mobile_ag]
    )
    temperature_K = 2.0 * kinetic_eV / (3.0 * int(mobile_ag.sum()) * KB_EV_K)
    max_force_eV_A = float(np.linalg.norm(force[mobile_ag], axis=1).max())
    return float(temperature_K), max_force_eV_A


def main() -> None:
    os.chdir(ROOT)
    geometry = json.loads((CANDIDATE / "geometry_report.json").read_text())
    lx, ly, lz = geometry["material_dimensions_A"]

    lmp = lammps(cmdargs=["-screen", "none", "-log", "none", "-sf", "omp"], comm=MPI.COMM_WORLD)
    commands = [
        "units metal",
        "atom_style atomic",
        "boundary p p f",
        "atom_modify sort 0 0",
        f"read_data {CANDIDATE}/ag_ti3sic2_qsgs.data.gz",
        f"read_dump {BDIR}/ag_minimized_100.dump.gz 100 x y z box no replace yes",
        f"include {BDIR}/potential.geometry_wca_fixed_skeleton.inc",
        "neighbor 2.0 bin",
        "neigh_modify delay 0 every 1 check yes",
        "region r_bottom block INF INF INF INF 0.0 20.0 units box",
        "group bottom region r_bottom",
        "group Ag type 4",
        "group skeleton type 1 2 3",
        "group mobile_Ag subtract Ag bottom",
        "fix freeze_bottom bottom setforce 0.0 0.0 0.0",
        "fix freeze_skeleton skeleton setforce 0.0 0.0 0.0",
        "mass 4 107.8682",
        "velocity mobile_Ag create 300.0 190427 mom yes rot no dist gaussian",
        "timestep 0.0005",
        "fix integrate mobile_Ag nve",
        "compute cna all cna/atom 3.35",
        "compute peatom all pe/atom",
        "thermo 100",
        "thermo_style custom step pe ke etotal",
        "thermo_modify flush yes",
        "run 0",
    ]
    for command in commands:
        lmp.command(command)

    start_dump = BDIR / "cold_diagnostic_start.dump"
    final_dump = BDIR / "cold_diagnostic_final.dump"
    dump_fields = "id type x y z c_cna c_peatom"
    lmp.command(f"write_dump mobile_Ag custom {start_dump} {dump_fields} modify sort id")
    temperature_K, max_force = ag_kinetic_state(lmp)
    print(
        "start "
        f"T_Ag={temperature_K:.6f} K "
        f"pe={lmp.get_thermo('pe'):.9f} eV "
        f"ke={lmp.get_thermo('ke'):.9f} eV "
        f"etotal={lmp.get_thermo('etotal'):.9f} eV "
        f"max_force={max_force:.9f} eV/A",
        flush=True,
    )

    for checkpoint in range(1, 6):
        lmp.command("run 100 pre yes post no")
        temperature_K, max_force = ag_kinetic_state(lmp)
        print(
            f"step={checkpoint * 100} time_ps={checkpoint * 0.05:.5f} "
            f"T_Ag={temperature_K:.6f} K "
            f"pe={lmp.get_thermo('pe'):.9f} eV "
            f"ke={lmp.get_thermo('ke'):.9f} eV "
            f"etotal={lmp.get_thermo('etotal'):.9f} eV "
            f"max_force={max_force:.9f} eV/A",
            flush=True,
        )

    lmp.command("run 0 pre no post no")
    lmp.command(f"write_dump mobile_Ag custom {final_dump} {dump_fields} modify sort id")
    temperature_K, max_force = ag_kinetic_state(lmp)
    print(
        "final "
        f"T_Ag={temperature_K:.6f} K "
        f"pe={lmp.get_thermo('pe'):.9f} eV "
        f"ke={lmp.get_thermo('ke'):.9f} eV "
        f"etotal={lmp.get_thermo('etotal'):.9f} eV "
        f"max_force={max_force:.9f} eV/A",
        flush=True,
    )
    print(f"snapshots={start_dump.name},{final_dump.name}; material_box_A={lx:.6f},{ly:.6f},{lz:.6f}", flush=True)
    lmp.close()
    MPI.Finalize()


if __name__ == "__main__":
    main()
