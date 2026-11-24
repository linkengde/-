#!/usr/bin/env python3
"""Compare a completed V14 Ag-liquid DFT candidate with its exact V13 frame.

Reports force differences; it does not impose an arbitrary acceptance threshold.
The candidate remains diagnostic until its DFT protocol and geometry are reviewed.
"""
import json
from pathlib import Path

import numpy as np
from ase.io import read
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "v14_stage2_p0_ag_dft_diagnostics" / "ag_liquid_precollapse"
manifest = json.loads((OUT / "candidate_manifest.json").read_text())
label_path = OUT / "label.extxyz"
label_meta_path = OUT / "label.json"
if not label_path.is_file() or not label_meta_path.is_file():
    raise SystemExit("DFT candidate label is not complete yet")
meta = json.loads(label_meta_path.read_text())
if meta.get("run_status") != "complete_candidate_only" or not meta.get("scf_converged"):
    raise SystemExit(f"Candidate DFT label is not converged: {meta.get('run_status')}")

atoms = read(label_path, format="extxyz")
f_dft = np.asarray(atoms.get_forces(), dtype=float)

dump_path = Path(manifest.get("source_trajectory_relative", manifest["source_trajectory"]))
if not dump_path.is_absolute():
    dump_path = (OUT / dump_path).resolve()
target_step = int(manifest["source_step"])
lines = dump_path.read_text().splitlines()
frame = None
i = 0
while i < len(lines):
    if lines[i].strip() != "ITEM: TIMESTEP":
        i += 1
        continue
    step = int(lines[i + 1])
    i += 2
    if lines[i].strip() != "ITEM: NUMBER OF ATOMS":
        raise RuntimeError("Unexpected LAMMPS dump structure")
    n = int(lines[i + 1]); i += 2
    box_header = lines[i].split(); i += 1
    bounds = np.array([[float(v) for v in lines[i + j].split()[:2]] for j in range(3)])
    i += 3
    atom_header = lines[i].split()[2:]; i += 1
    rows = np.array([[float(v) for v in lines[i + j].split()] for j in range(n)])
    i += n
    if step == target_step:
        col = {name: k for k, name in enumerate(atom_header)}
        needed = ["id", "type", "x", "y", "z", "fx", "fy", "fz"]
        missing = [name for name in needed if name not in col]
        if missing:
            raise RuntimeError(f"Missing source dump columns: {missing}")
        order = np.argsort(rows[:, col["id"]].astype(int))
        frame = {
            "step": step,
            "ids": rows[order, col["id"]].astype(int),
            "types": rows[order, col["type"]].astype(int),
            "positions": rows[order][:, [col["x"], col["y"], col["z"]]],
            "forces": rows[order][:, [col["fx"], col["fy"], col["fz"]]],
            "lengths": bounds[:, 1] - bounds[:, 0],
        }
        break
if frame is None:
    raise RuntimeError(f"Step {target_step} not found in {dump_path}")

if len(atoms) != len(frame["ids"]):
    raise RuntimeError(f"Atom count mismatch: DFT {len(atoms)} vs source {len(frame['ids'])}")

# Map by exact geometry under periodic boundaries, so atom ordering cannot silently
# corrupt the force comparison.
delta = atoms.positions[:, None, :] - frame["positions"][None, :, :]
lengths = frame["lengths"]
delta -= np.round(delta / lengths[None, None, :]) * lengths[None, None, :]
dist = np.linalg.norm(delta, axis=2)
rows_dft, cols_src = linear_sum_assignment(dist)
max_pos_error = float(dist[rows_dft, cols_src].max())
if max_pos_error > 1e-5:
    raise RuntimeError(f"Could not map DFT atoms to original frame; max MIC mismatch {max_pos_error:g} Å")

f_ref = frame["forces"][cols_src]
err = f_dft[rows_dft] - f_ref
norm = np.linalg.norm(err, axis=1)
summary = {
    "candidate_id": manifest["candidate_id"],
    "scf_converged": bool(meta["scf_converged"]),
    "scf_iterations": meta.get("scf_iterations"),
    "source_step": target_step,
    "source_time_ps": manifest["source_time_ps"],
    "natoms": len(atoms),
    "max_position_mapping_error_A": max_pos_error,
    "force_rmse_component_eV_A": float(np.sqrt(np.mean(err ** 2))),
    "force_rmse_vector_norm_eV_A": float(np.sqrt(np.mean(np.sum(err ** 2, axis=1)))),
    "force_error_norm_percentiles_eV_A": {
        str(p): float(np.percentile(norm, p)) for p in [50, 90, 95, 99, 100]
    },
    "interpretation": "Diagnostic comparison only; no universal force-error cutoff is applied.",
    "top_force_error_atoms": [],
}
id_to_dft = {int(frame["ids"][j]): int(rows_dft[k]) for k, j in enumerate(cols_src)}
pair_ids = [int(x) for x in manifest.get("minimum_pair_atom_indices", [])]
if len(pair_ids) == 2 and all(atom_id in id_to_dft for atom_id in pair_ids):
    atom_i, atom_j = [id_to_dft[x] for x in pair_ids]
    rij = atoms.positions[atom_j] - atoms.positions[atom_i]
    rij -= np.round(rij / lengths) * lengths
    distance = float(np.linalg.norm(rij))
    rhat = rij / distance
    fi_v13 = frame["forces"][np.where(frame["ids"] == pair_ids[0])[0][0]]
    fj_v13 = frame["forces"][np.where(frame["ids"] == pair_ids[1])[0][0]]
    summary["identified_short_pair_radial_force"] = {
        "lammps_ids": pair_ids,
        "distance_A": distance,
        "dft_repulsive_outward_force_on_j_eV_A": float(np.dot(f_dft[atom_j], rhat)),
        "dft_repulsive_outward_force_on_i_eV_A": float(np.dot(f_dft[atom_i], -rhat)),
        "v13_repulsive_outward_force_on_j_eV_A": float(np.dot(fj_v13, rhat)),
        "v13_repulsive_outward_force_on_i_eV_A": float(np.dot(fi_v13, -rhat)),
    }
for atom_idx in np.argsort(norm)[-5:][::-1]:
    j = cols_src[atom_idx]
    summary["top_force_error_atoms"].append({
        "dft_atom_index_1based": int(atom_idx + 1),
        "source_lammps_id": int(frame["ids"][j]),
        "element": atoms[atom_idx].symbol,
        "position_A": atoms.positions[atom_idx].tolist(),
        "force_dft_eV_A": f_dft[atom_idx].tolist(),
        "force_v13_eV_A": f_ref[atom_idx].tolist(),
        "force_error_norm_eV_A": float(norm[atom_idx]),
    })
out_path = OUT / "dft_v13_force_comparison.json"
out_path.write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
