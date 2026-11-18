#!/usr/bin/env python3
"""Compare zero-step LAMMPS PACE outputs with PyACE on the same smoke-test cell."""
import json
from pathlib import Path

import numpy as np
from ase.io import read
from pyace import ACEBBasisSet, PyACECalculator


HERE = Path(__file__).resolve().parent
CHECKS = HERE / "lammps_candidate_checks"
MODELS = {
    "ACE_weight_1": HERE / "ace_kconv_refits" / "AgCSiTi_B1616_E10000_F1_a0.001_Ti_gap_weight1.ace.yaml",
    "ACE_weight_005": HERE / "ace_kconv_refits" / "AgCSiTi_B1616_E10000_F1_a0.001_Ti_gap_weight0.05.ace.yaml",
}


def dump_forces(path, natoms):
    lines = Path(path).read_text().splitlines()
    at = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    rows = [line.split() for line in lines[at + 1:at + 1 + natoms]]
    forces = np.zeros((natoms, 3))
    for row in rows:
        atom_id = int(row[0]) - 1
        forces[atom_id] = np.asarray(row[5:8], dtype=float)
    return forces


def thermo_energy(path):
    lines = Path(path).read_text().splitlines()
    for i, line in enumerate(lines):
        if line.split() == ["Step", "PotEng"]:
            return float(lines[i + 1].split()[1])
    raise ValueError(f"No thermo energy in {path}")


def main():
    summary = {}
    for dirname, model_path in MODELS.items():
        directory = CHECKS / dirname
        structure = read(directory / "data.lammps", format="lammps-data", atom_style="atomic",
                         Z_of_type={1: 47, 2: 6, 3: 14, 4: 22})
        structure.calc = PyACECalculator(ACEBBasisSet(str(model_path)))
        energy_pyace = float(structure.get_potential_energy())
        forces_pyace = np.asarray(structure.get_forces(), dtype=float)
        energy_lammps = thermo_energy(directory / "lammps.log")
        forces_lammps = dump_forces(directory / "forces.dump", len(structure))
        delta = forces_lammps - forces_pyace
        summary[dirname] = {
            "natoms": len(structure),
            "energy_pyace_eV": energy_pyace,
            "energy_lammps_eV": energy_lammps,
            "energy_abs_error_eV": abs(energy_lammps - energy_pyace),
            "force_component_rmse_eV_A": float(np.sqrt(np.mean(delta ** 2))),
            "force_vector_max_error_eV_A": float(np.linalg.norm(delta, axis=1).max()),
        }
    out = {"status": "zero-step PACE/PyACE consistency check", "models": summary}
    (CHECKS / "ace_consistency.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
