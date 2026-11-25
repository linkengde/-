#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""对新细骨架候选执行 300 K run 0 受力预检，不输入热量。"""

import json
import os
from pathlib import Path

import numpy as np
from mpi4py import MPI
from lammps import lammps


ROOT = Path("/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer")
CANDIDATE = ROOT / "potential_validation/geometry_constraint_B/thinner_scaffold_candidate_05_20261125"
OUTPUT = CANDIDATE / "heating_trial_20261126"
B_DIR = ROOT / "potential_validation/geometry_constraint_B"


def main():
    """装载结构与 WCA 几何模型，生成速度并记录未积分初态受力。"""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    os.chdir(ROOT)
    lmp = lammps(cmdargs=["-screen", "none", "-log", "none", "-sf", "omp"],
                 comm=MPI.COMM_WORLD)
    commands = [
        "units metal", "atom_style atomic", "boundary p p f", "atom_modify sort 0 0",
        f"read_data {CANDIDATE}/ag_ti3sic2_qsgs.data",
        f"include {B_DIR}/potential.geometry_wca_fixed_skeleton.inc",
        "neighbor 2.0 bin", "neigh_modify delay 0 every 1 check yes",
        "region r_bottom block INF INF INF INF 0.0 20.0 units box",
        "group bottom region r_bottom", "group Ag type 4", "group skeleton type 1 2 3",
        "group bottom_Ag intersect Ag bottom", "group mobile_Ag subtract Ag bottom_Ag",
        "fix freeze_skeleton skeleton setforce 0.0 0.0 0.0",
        "fix freeze_bottom_Ag bottom_Ag setforce 0.0 0.0 0.0",
        "mass 4 107.8682", "velocity mobile_Ag create 300.0 190427 mom yes rot no dist gaussian",
        "timestep 0.0005", "fix integrate mobile_Ag nve", "thermo 1",
        "thermo_style custom step pe ke etotal temp", "thermo_modify lost error flush yes",
        "run 0",
    ]
    for command in commands:
        lmp.command(command)

    nlocal = int(lmp.extract_setting("nlocal"))
    types = np.asarray(lmp.numpy.extract_atom("type")[:nlocal], dtype=np.int32)
    positions = np.asarray(lmp.numpy.extract_atom("x")[:nlocal], dtype=float)
    forces = np.asarray(lmp.numpy.extract_atom("f")[:nlocal], dtype=float)
    ids = np.asarray(lmp.numpy.extract_atom("id")[:nlocal], dtype=np.int64)
    mobile = (types == 4) & (positions[:, 2] >= 20.0)
    force_norm = np.linalg.norm(forces, axis=1)
    mobile_ids = np.flatnonzero(mobile)
    worst = mobile_ids[np.argmax(force_norm[mobile])]
    report = {
        "status": "COLD_RUN0_ONLY_NO_HEATING",
        "atoms": int(lmp.get_natoms()),
        "mobile_Ag": int(mobile.sum()),
        "Ag_temperature_K": float(lmp.get_thermo("temp")),
        "pe_eV": float(lmp.get_thermo("pe")),
        "ke_eV": float(lmp.get_thermo("ke")),
        "etotal_eV": float(lmp.get_thermo("etotal")),
        "Fmax_mobile_Ag_eV_A": float(force_norm[mobile].max()),
        "F99_mobile_Ag_eV_A": float(np.percentile(force_norm[mobile], 99)),
        "Fmax_atom_id": int(ids[worst]),
    }
    lmp.command(
        f"write_dump all custom/gz {OUTPUT}/cold_start_300K_run0.dump.gz "
        "id type x y z vx vy vz fx fy fz modify sort id"
    )
    (OUTPUT / "cold_preflight.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False), flush=True)
    lmp.close()
    MPI.Finalize()


if __name__ == "__main__":
    main()
