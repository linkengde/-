#!/workspace/.venvs/gpaw-mpi/bin/python
"""Build one paired tangential Ti-displacement diagnostic from archived REF."""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import ase
import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[5]
SOURCE_ROOT = REPO / "research/mace-v12-transfer/coordination/reports/window-e/tsc_surface_dev_prismatic_11m20_20261021"
SOURCE_TASK_ID = "TSC120_DEV_PRISM_11m20_COD9009647_10L_REF"
TASK_ID = "TSC120_DIAG_PRISM_11m20_COD9009647_10L_RESP_TIX020"
DISPLACEMENT_A = 0.02
BOTTOM_ID = 4
TOP_ID = 113


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometry_sha(atoms) -> str:
    h = hashlib.sha256()
    h.update(" ".join(atoms.get_chemical_symbols()).encode())
    h.update(np.asarray(atoms.cell, dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(np.round(atoms.positions, 8), dtype="<f8").tobytes(order="C"))
    h.update(np.asarray(atoms.pbc, dtype=np.uint8).tobytes())
    h.update(np.asarray(atoms.arrays["source_atom_id"], dtype="<i8").tobytes())
    return h.hexdigest()


def minimum_distance(atoms) -> tuple[float, int, int]:
    distances = atoms.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    i, j = np.unravel_index(np.argmin(distances), distances.shape)
    return float(distances[i, j]), int(i + 1), int(j + 1)


def inversion_symmetry_error(atoms) -> float:
    positions = np.asarray(atoms.positions)
    fractional = atoms.get_scaled_positions(wrap=True)
    symbols = atoms.get_chemical_symbols()
    zlo, zhi = float(positions[:, 2].min()), float(positions[:, 2].max())
    errors = []
    for i, symbol in enumerate(symbols):
        target_xy = (-fractional[i, :2]) % 1.0
        target_z = zlo + zhi - positions[i, 2]
        candidates = []
        for j, other in enumerate(symbols):
            if other != symbol:
                continue
            dxy = (fractional[j, :2] - target_xy + 0.5) % 1.0 - 0.5
            drxy = dxy @ atoms.cell[:2]
            dz = positions[j, 2] - target_z
            candidates.append(float(np.sqrt(np.dot(drxy, drxy) + dz * dz)))
        errors.append(min(candidates))
    return max(errors)


if ase.__version__ != "3.29.0":
    raise RuntimeError(f"Expected ASE 3.29.0, got {ase.__version__}")
source_path = SOURCE_ROOT / "inputs" / f"{SOURCE_TASK_ID}.extxyz"
expected_source_sha = "fd5ff5c1cae1490bcbde9ffd08f329c7ec8cac687555f29568d725f4230cc5c3"
if sha256(source_path) != expected_source_sha:
    raise RuntimeError("Archived reference input SHA mismatch")

ref = read(source_path, format="extxyz")
if len(ref) != 120 or Counter(ref.get_chemical_symbols()) != {"Ti": 60, "Si": 20, "C": 40}:
    raise RuntimeError("Unexpected source atom count or composition")
if not np.array_equal(ref.arrays.get("source_atom_id"), np.arange(1, 121)):
    raise RuntimeError("Source atom ID map is missing or noncontiguous")
if ref.pbc.tolist() != [True, True, True]:
    raise RuntimeError("Unexpected PBC")

response = ref.copy()
response.info = dict(ref.info)
response.info.update({
    "structure_id": TASK_ID,
    "source_family": "surface_prismatic_11m20_cod9009647_10layer",
    "candidate_role": "SAME_SOURCE_SURFACE_RESPONSE_DIAGNOSTIC_ONLY",
    "diagnostic_scope": "in-plane tangential response; not independent crystallographic-source validation",
    "displacement_definition": "inversion-paired terminal Ti atoms displaced oppositely along the first in-plane cell vector",
    "displacement_magnitude_angstrom": DISPLACEMENT_A,
    "displaced_source_atom_ids": [BOTTOM_ID, TOP_ID],
    "source_reference_input_sha256": expected_source_sha,
})
direction = np.array(ref.cell[0], dtype=float, copy=True)
direction /= np.linalg.norm(direction)
response.positions[BOTTOM_ID - 1] += DISPLACEMENT_A * direction
response.positions[TOP_ID - 1] -= DISPLACEMENT_A * direction

delta = response.positions - ref.positions
expected_delta = np.zeros_like(delta)
expected_delta[BOTTOM_ID - 1] = DISPLACEMENT_A * direction
expected_delta[TOP_ID - 1] = -DISPLACEMENT_A * direction
if not np.allclose(delta, expected_delta, atol=1e-12, rtol=0):
    raise RuntimeError("Unexpected coordinate changes outside the designated Ti pair")
if not np.array_equal(response.cell.array, ref.cell.array) or not np.array_equal(response.pbc, ref.pbc):
    raise RuntimeError("Cell or PBC changed")
minimum, min_i, min_j = minimum_distance(response)
if minimum < 1.5:
    raise RuntimeError(f"New geometry has a short contact: {minimum:.6f} A at {min_i}/{min_j}")
symmetry_error = inversion_symmetry_error(response)
if symmetry_error > 1e-5:
    raise RuntimeError(f"Inversion symmetry failed: {symmetry_error:.6g} A")

input_path = ROOT / "inputs" / f"{TASK_ID}.extxyz"
write(input_path, response, format="extxyz")
serialized = read(input_path, format="extxyz")
input_sha = sha256(input_path)
geom_sha = geometry_sha(serialized)
if len(serialized) != len(ref) or serialized.get_chemical_symbols() != ref.get_chemical_symbols():
    raise RuntimeError("Serialized input changed species ordering")

source_tasks = json.loads((SOURCE_ROOT / "TASKS.json").read_text())
source_commit = source_tasks["source_git_commit_at_preparation"]
task = {
    "structure_id": TASK_ID,
    "input_path": f"inputs/{TASK_ID}.extxyz",
    "input_sha256": input_sha,
    "geometry_sha256": geom_sha,
    "atom_count": len(serialized),
    "composition": dict(sorted(Counter(serialized.get_chemical_symbols()).items())),
    "source_family": "surface_prismatic_11m20_cod9009647_10layer",
    "role_candidate": "SAME_SOURCE_SURFACE_RESPONSE_DIAGNOSTIC_ONLY",
    "displaced_source_atom_ids": [BOTTOM_ID, TOP_ID],
    "displacement_cartesian_angstrom": {
        str(BOTTOM_ID): (DISPLACEMENT_A * direction).tolist(),
        str(TOP_ID): (-DISPLACEMENT_A * direction).tolist(),
    },
}
tasks = {
    "schema_version": "1.0",
    "source_cif": source_tasks["source_cif"],
    "source_cif_sha256": source_tasks["source_cif_sha256"],
    "source_git_commit_at_preparation": source_commit,
    "source_reference_input": str(source_path.relative_to(REPO)),
    "source_reference_input_sha256": expected_source_sha,
    "role_policy": "diagnostic candidate only; TRAIN/DEV/TEST are unchanged",
    "tasks": [task],
}
(ROOT / "TASKS.json").write_text(json.dumps(tasks, indent=2, sort_keys=True) + "\n")

config_path = ROOT / "SCIENTIFIC_CONFIG.json"
config = json.loads(config_path.read_text())
config["scope"] = "One fixed-geometry PBE/PAW static point for tangential Ti-displacement response on the COD 9009647 (11-20) slab; reuse archived same-geometry reference force label for paired response analysis."
config["role"] = "SAME_SOURCE_DIAGNOSTIC_ONLY; no TRAIN/DEV/TEST assignment"
config["authorization"] = {
    "authorized_by_current_user": True,
    "scope": "Build and calculate only this one deterministic fixed-geometry response point; preserve the source reference, do not change any data role, and do not run relaxation, MD, fitting, or other DFT points.",
    "c_approval_required": False,
}
(ROOT / "SCIENTIFIC_CONFIG.json").write_text(json.dumps(config, indent=2) + "\n")

audit = {
    "status": "GEOMETRY_AUDIT_PASS_SAME_SOURCE_DIAGNOSTIC_ONLY",
    "ASE_version": ase.__version__,
    "source_cif_sha256": tasks["source_cif_sha256"],
    "source_reference_input_sha256": expected_source_sha,
    "input_sha256": input_sha,
    "geometry_sha256": geom_sha,
    "atom_count": len(serialized),
    "composition": dict(sorted(Counter(serialized.get_chemical_symbols()).items())),
    "source_atom_ids_preserved": True,
    "cell_unchanged": True,
    "pbc": serialized.pbc.tolist(),
    "displacement_source_atom_ids": [BOTTOM_ID, TOP_ID],
    "displacement_magnitude_angstrom": DISPLACEMENT_A,
    "displacement_direction_unit_vector": direction.tolist(),
    "minimum_distance_angstrom": minimum,
    "minimum_distance_source_atom_ids": [min_i, min_j],
    "inversion_symmetry_max_error_angstrom": symmetry_error,
    "data_roles_changed": False,
    "independent_crystallographic_source": False,
}
(ROOT / "STRUCTURE_AUDIT.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
print(json.dumps(audit, indent=2))
