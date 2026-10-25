#!/usr/bin/env python3
"""Extract the TSCBOP-relaxed orthorhombic Ti3SiC2 parent basis."""
from __future__ import annotations
import hashlib
import json
import re
from pathlib import Path
import numpy as np
from ase.build import make_supercell
from ase.io import read

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1]
REPO = HERE.parents[3]
CIF = REPO / "research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif"
BULK_DIR = MODEL / "potential_validation/ti3sic2_bulk_tscbop"
INITIAL_DATA = BULK_DIR / "ti3sic2_bulk_3x3x2.data"
RELAXED_DATA = BULK_DIR / "zero_pressure_minimized.data"
POTENTIAL = MODEL / "potentials/TSCBOP.tersoff"
REPEATS_PROBE = np.array([3, 3, 2], dtype=int)
REPEATS_MODEL = [98, 57, 9]
TYPE_TO_Z = {1: 22, 2: 14, 3: 6}


def box_lo(data_path: Path) -> np.ndarray:
    lines = data_path.read_text().splitlines()[:20]
    result = []
    for axis in "xyz":
        line = next(line for line in lines if f"{axis}lo {axis}hi" in line)
        result.append(float(line.split()[0]))
    return np.asarray(result, dtype=float)


bulk = read(CIF)
base = make_supercell(bulk, [[1, 0, 0], [1, 2, 0], [0, 0, 1]])
initial = read(INITIAL_DATA, format="lammps-data", atom_style="atomic", units="metal", Z_of_type=TYPE_TO_Z)
relaxed = read(RELAXED_DATA, format="lammps-data", atom_style="atomic", units="metal", Z_of_type=TYPE_TO_Z)
if len(base) != 24 or len(initial) != 432 or len(relaxed) != 432:
    raise RuntimeError("Unexpected Ti3SiC2 basis or probe supercell size")

# Match atoms by their LAMMPS IDs, not by file-row ordering.
order_i = np.argsort(initial.arrays["id"])
order_r = np.argsort(relaxed.arrays["id"])
if not np.array_equal(initial.arrays["id"][order_i], relaxed.arrays["id"][order_r]):
    raise RuntimeError("Probe atom IDs changed during box relaxation")
initial_pos = initial.positions[order_i]
relaxed_pos = relaxed.positions[order_r]
initial_z = initial.numbers[order_i]
relaxed_z = relaxed.numbers[order_r]
if not np.array_equal(initial_z, relaxed_z):
    raise RuntimeError("Probe elements changed during relaxation")

initial_lengths = initial.cell.lengths()
relaxed_lengths = relaxed.cell.lengths()
origin = box_lo(RELAXED_DATA)
initial_frac = np.mod(initial_pos / initial_lengths, 1.0)
relaxed_frac = np.mod((relaxed_pos - origin) / relaxed_lengths, 1.0)
source_basis_frac = base.get_scaled_positions(wrap=True)
source_basis_z = base.numbers
basis_initial = np.mod(initial_frac * REPEATS_PROBE, 1.0)
basis_relaxed = np.mod(relaxed_frac * REPEATS_PROBE, 1.0)

basis_ids = np.empty(len(initial), dtype=int)
for i, (frac, z) in enumerate(zip(basis_initial, initial_z)):
    delta = (source_basis_frac - frac + 0.5) % 1.0 - 0.5
    delta[source_basis_z != z] = 100.0
    basis_ids[i] = int(np.argmin(np.sum(delta * delta, axis=1)))
if not np.all(np.bincount(basis_ids, minlength=len(base)) == np.prod(REPEATS_PROBE)):
    raise RuntimeError("Could not map exactly 18 repetitions onto each parent basis site")

relaxed_basis = []
for basis_id in range(len(base)):
    samples = basis_relaxed[basis_ids == basis_id]
    # Circular mean preserves fractional coordinates across periodic seams.
    mean = np.angle(np.exp(2j * np.pi * samples).mean(axis=0)) / (2 * np.pi)
    relaxed_basis.append(np.mod(mean, 1.0))
relaxed_basis = np.asarray(relaxed_basis)
unit_cell_lengths = relaxed_lengths / REPEATS_PROBE
basis_delta = (relaxed_basis - source_basis_frac + 0.5) % 1.0 - 0.5
basis_delta_A = basis_delta * unit_cell_lengths

result = {
    "status": "TSCBOP_ZERO_STRESS_PARENT_BASIS_DERIVED_FOR_SEPARATE_QSGS_CANDIDATE",
    "cif_sha256": hashlib.sha256(CIF.read_bytes()).hexdigest(),
    "tscbop_sha256": hashlib.sha256(POTENTIAL.read_bytes()).hexdigest(),
    "probe_relaxed_data_sha256": hashlib.sha256(RELAXED_DATA.read_bytes()).hexdigest(),
    "probe_supercell_repeats": REPEATS_PROBE.tolist(),
    "model_repeats": REPEATS_MODEL,
    "source_unit_cell_lengths_A": [float(x) for x in base.cell.lengths()],
    "relaxed_unit_cell_lengths_A": [float(x) for x in unit_cell_lengths],
    "model_material_dimensions_A": [float(x) for x in unit_cell_lengths * np.asarray(REPEATS_MODEL)],
    "cell_scale_percent": [float(x) for x in (unit_cell_lengths / base.cell.lengths() - 1.0) * 100],
    "source_basis_symbols": base.get_chemical_symbols(),
    "source_basis_fractional": source_basis_frac.tolist(),
    "relaxed_basis_fractional": relaxed_basis.tolist(),
    "parent_basis_relaxation_displacement_A": basis_delta_A.tolist(),
    "maximum_parent_basis_displacement_A": float(np.max(np.linalg.norm(basis_delta_A, axis=1))),
    "mean_parent_basis_displacement_A": float(np.mean(np.linalg.norm(basis_delta_A, axis=1))),
    "mapping_count_per_basis_site": [int(x) for x in np.bincount(basis_ids, minlength=len(base))],
}
(HERE / "relaxed_parent_basis.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({
    "unit_cell_lengths_A": result["relaxed_unit_cell_lengths_A"],
    "model_dimensions_A": result["model_material_dimensions_A"],
    "cell_scale_percent": result["cell_scale_percent"],
    "max_basis_displacement_A": result["maximum_parent_basis_displacement_A"],
    "mean_basis_displacement_A": result["mean_parent_basis_displacement_A"],
    "output": str(HERE / "relaxed_parent_basis.json"),
}, indent=2))
