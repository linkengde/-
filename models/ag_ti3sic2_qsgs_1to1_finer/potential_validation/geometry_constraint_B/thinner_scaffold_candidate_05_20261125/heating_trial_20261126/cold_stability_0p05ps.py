#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""候选 05 的 300 K、0.05 ps 无热源稳定性检查。"""

import csv
import json
import os
from pathlib import Path

import numpy as np
from mpi4py import MPI
from lammps import lammps


ROOT = Path("/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer")
B_DIR = ROOT / "potential_validation/geometry_constraint_B"
CANDIDATE = B_DIR / "thinner_scaffold_candidate_05_20261125"
OUT = CANDIDATE / "heating_trial_20261126/cold_stability_0p05ps"
DATA = CANDIDATE / "ag_ti3sic2_qsgs.data"


def main():
    """只对 Ag 施加速度并积分 100 步，骨架和底部 Ag 保持固定。"""
    OUT.mkdir(parents=True, exist_ok=True)
    os.chdir(ROOT)
    simulation = lammps(
        cmdargs=["-screen", "none", "-log", "none", "-sf", "omp"],
        comm=MPI.COMM_WORLD,
    )
    setup = [
        "units metal",
        "atom_style atomic",
        "boundary p p f",
        "atom_modify sort 0 0",
        f"read_data {DATA}",
        f"include {B_DIR}/potential.geometry_wca_fixed_skeleton.inc",
        "neighbor 2.0 bin",
        "neigh_modify delay 0 every 1 check yes",
        "region r_bottom block INF INF INF INF 0.0 20.0 units box",
        "group bottom region r_bottom",
        "group Ag type 4",
        "group skeleton type 1 2 3",
        "group bottom_Ag intersect Ag bottom",
        "group mobile_Ag subtract Ag bottom_Ag",
        "fix freeze_skeleton skeleton setforce 0.0 0.0 0.0",
        "fix freeze_bottom_Ag bottom_Ag setforce 0.0 0.0 0.0",
        "mass 4 107.8682",
        "velocity mobile_Ag create 300.0 190427 mom yes rot no dist gaussian",
        "timestep 0.0005",
        "fix integrate mobile_Ag nve",
        "thermo 100",
        "thermo_style custom step pe ke etotal temp",
        "thermo_modify lost error flush yes",
        "run 0",
    ]
    for command in setup:
        simulation.command(command)

    nlocal = int(simulation.extract_setting("nlocal"))
    atom_type = np.asarray(simulation.numpy.extract_atom("type")[:nlocal], dtype=np.int32)
    position = np.asarray(simulation.numpy.extract_atom("x")[:nlocal], dtype=float)
    atom_ids = np.asarray(simulation.numpy.extract_atom("id")[:nlocal], dtype=np.int64)
    mobile = (atom_type == 4) & (position[:, 2] >= 20.0)
    force = simulation.numpy.extract_atom("f")
    velocity = simulation.numpy.extract_atom("v")
    initial_energy = float(simulation.get_thermo("etotal"))
    maximum_force = float(np.linalg.norm(force[:nlocal][mobile], axis=1).max())
    maximum_force_id = int(atom_ids[np.flatnonzero(mobile)[np.argmax(
        np.linalg.norm(force[:nlocal][mobile], axis=1)
    )]])

    metrics_path = OUT / "metrics.csv"
    with metrics_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["step", "time_ps", "Etot_eV", "delta_E_eV", "Fmax_mobile_Ag_eV_A"])
        writer.writerow([0, 0.0, initial_energy, 0.0, maximum_force])
        maximum_energy_deviation = 0.0
        completed = 0
        for step in range(1, 101):
            simulation.command("run 1 pre no post no")
            completed = step
            energy = float(simulation.get_thermo("etotal"))
            current_force = np.linalg.norm(force[:nlocal][mobile], axis=1)
            fmax = float(current_force.max())
            if fmax > maximum_force:
                maximum_force = fmax
                maximum_force_id = int(atom_ids[np.flatnonzero(mobile)[np.argmax(current_force)]])
            deviation = energy - initial_energy
            maximum_energy_deviation = max(maximum_energy_deviation, abs(deviation))
            writer.writerow([step, f"{step * 0.0005:.6f}", f"{energy:.12f}",
                             f"{deviation:.12g}", f"{fmax:.8f}"])
            if step % 20 == 0:
                stream.flush()
            if not np.isfinite(fmax) or fmax > 50.0:
                break

    simulation.command(
        f"write_dump all custom/gz {OUT / 'cold_end.dump.gz'} "
        "id type x y z vx vy vz fx fy fz modify sort id"
    )
    report = {
        "status": "COMPLETED" if completed == 100 and maximum_force < 50.0 else "FORCE_SAFETY_STOP",
        "stage": "cold NVE only; no laser or TTM source",
        "atoms": int(simulation.get_natoms()),
        "Ag_mobile": int(mobile.sum()),
        "duration_ps": completed * 0.0005,
        "timestep_fs": 0.5,
        "initial_Etotal_eV": initial_energy,
        "final_Etotal_eV": float(simulation.get_thermo("etotal")),
        "final_energy_drift_eV": float(simulation.get_thermo("etotal")) - initial_energy,
        "max_abs_energy_drift_eV": maximum_energy_deviation,
        "Fmax_mobile_Ag_eV_A": maximum_force,
        "Fmax_atom_id": maximum_force_id,
        "rigid_skeleton": "all Ti/Si/C fixed; Ag-skeleton WCA repulsion",
        "input_data": str(DATA),
    }
    (OUT / "cold_stability_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    simulation.close()
    MPI.Finalize()


if __name__ == "__main__":
    main()
