#!/usr/bin/env python3
"""Compare two-rank LAMMPS run 0 with ASE using the identical v11 MACE model."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
MODEL = ROOT / "checkpoints/MACE_periodic_v11_energyfocus_run-44.model"
WRAPPER = OUT / "mace_v11_lammps_mliap_forcecpu.pt"
EXTXYZ = Path("/workspace/.local/lammps-ttm/run_1ps_diagnostic_20261001/potential_screening/periodic_interface_v4/pbe_external_v9_holdouts/AgTi_d2p55/lammps_v9_serial_smoke_retry2/cluster_periodic.extxyz")
DUMP = OUT / "forces.run0.mpi.retry4.dump"
LOG = OUT / "log.run0.mpi.retry4"

dump_lines = DUMP.read_text().splitlines()
header = next(i for i, line in enumerate(dump_lines) if line.startswith("ITEM: ATOMS"))
columns = dump_lines[header].split()[2:]
rows = [line.split() for line in dump_lines[header + 1:] if line.strip() and not line.startswith("ITEM:")]
dump = {int(row[columns.index("id")]): row for row in rows}
if len(dump) != 49:
    raise SystemExit(f"Expected 49 dumped atoms, found {len(dump)}")
log_text = LOG.read_text()
match = re.search(r"^\s*0\s+49\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", log_text, re.M)
if not match:
    raise SystemExit("Could not parse step-0 PE/Fmax from LAMMPS log")
lmp_energy, lmp_fmax = map(float, match.groups())

atoms = read(EXTXYZ)
# The LAMMPS data converter assigned sequential IDs in this extxyz row order.
ids = np.arange(1, len(atoms) + 1, dtype=int)
symbols = atoms.get_chemical_symbols()
z_to_type = {"C": 1, "Si": 2, "Ti": 3, "Ag": 4}
if len(atoms) != 49 or set(ids.tolist()) != set(dump):
    raise SystemExit("ASE and LAMMPS atom IDs do not match")
for i, atom_id in enumerate(ids):
    row = dump[int(atom_id)]
    if int(row[columns.index("type")]) != z_to_type[symbols[i]]:
        raise SystemExit(f"Type mapping mismatch for atom id {atom_id}")

torch.set_default_dtype(torch.float64)
model = torch.load(MODEL, map_location="cpu", weights_only=False).double().eval()
atoms.calc = MACECalculator(models=[model], device="cpu", default_dtype="float64")
ase_energy = float(atoms.get_potential_energy())
ase_forces = np.asarray(atoms.get_forces())
lmp_forces = np.zeros_like(ase_forces)
lmp_positions = np.zeros_like(atoms.positions)
for i, atom_id in enumerate(ids):
    row = dump[int(atom_id)]
    lmp_positions[i] = [float(row[columns.index(k)]) for k in ("x", "y", "z")]
    lmp_forces[i] = [float(row[columns.index(k)]) for k in ("fx", "fy", "fz")]
delta = lmp_forces - ase_forces
result = {
    "mpi_ranks": 2,
    "atoms": len(atoms),
    "step": 0,
    "lammps_energy_eV": lmp_energy,
    "ase_energy_eV": ase_energy,
    "abs_energy_error_eV": abs(lmp_energy - ase_energy),
    "lammps_Fmax_eV_A": lmp_fmax,
    "ase_Fmax_eV_A": float(np.linalg.norm(ase_forces, axis=1).max()),
    "position_max_abs_error_A": float(np.abs(lmp_positions - atoms.positions).max()),
    "force_component_max_abs_error_eV_A": float(np.abs(delta).max()),
    "force_component_RMSE_eV_A": float(np.sqrt(np.mean(delta**2))),
    "force_vector_max_abs_error_eV_A": float(np.linalg.norm(delta, axis=1).max()),
    "force_vector_RMSE_eV_A": float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
    "dump": str(DUMP),
    "log": str(LOG),
    "model_sha256": hashlib.sha256(MODEL.read_bytes()).hexdigest(),
    "mliap_wrapper_sha256": hashlib.sha256(WRAPPER.read_bytes()).hexdigest(),
    "warnings": re.findall(r"^WARNING:.*$", log_text, re.M),
}
(OUT / "run0_ase_parity.retry4.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
