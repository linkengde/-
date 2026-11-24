#!/usr/bin/env python3
"""Audit a full-size spring-free Ag-W EFS minimization for v5 gates."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
PACK = HERE.parents[2]
MODEL = PACK / "models/ag_w_qsgs_200A"
MASS_W = 183.84
MASS_AG = 107.8682
BOX_XY = 200.0
MATERIAL_Z = 200.0
FLOORS = {"Ag_W": 2.54, "W_W": 2.60, "Ag_Ag": 2.48}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_dump(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        lines = stream.read().splitlines()
    atom_header = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    columns = lines[atom_header].split()[2:]
    index = {name: columns.index(name) for name in ("id", "type", "x", "y", "z", "fx", "fy", "fz")}
    rows = np.loadtxt(lines[atom_header + 1 :], dtype=float, ndmin=2)
    ids = rows[:, index["id"]].astype(np.int64)
    order = np.argsort(ids)
    types = rows[order, index["type"]].astype(np.int8)
    positions = rows[order][:, [index["x"], index["y"], index["z"]]].copy()
    forces = rows[order][:, [index["fx"], index["fy"], index["fz"]]].copy()
    return ids[order], types, positions, forces


def minimizer_status(path: Path, max_iterations: int) -> dict[str, object]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as stream:
        text = stream.read()
    matches = re.findall(r"Stopping criterion = ([^\n]+)", text)
    iterations = re.findall(r"Iterations, force evaluations =\s*(\d+)", text)
    reason = matches[-1].strip() if matches else "MISSING"
    iters = int(iterations[-1]) if iterations else None
    return {
        "stopping_criterion": reason,
        "iterations": iters,
        "iteration_cap": max_iterations,
        "converged": bool("force tolerance" in reason.lower() or "energy tolerance" in reason.lower()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--iteration-cap", type=int, default=10000)
    parser.add_argument("--fixed-z", type=float, default=0.0,
                        help="Exclude atoms in the fixed bottom layer from force gates.")
    args = parser.parse_args()
    dump_path, data_path, log_path, output_path = (p.resolve() for p in (args.dump, args.data, args.log, args.output))
    ids, types, positions, forces = read_dump(dump_path)
    positions[:, :2] %= BOX_XY
    if len(ids) != len(np.unique(ids)):
        raise RuntimeError("Dump contains duplicate atom IDs")

    n_w = int(np.count_nonzero(types == 1))
    n_ag = int(np.count_nonzero(types == 2))
    mass_ratio = n_ag * MASS_AG / (n_w * MASS_W)
    morphology = np.load(MODEL / "qsgs_morphology.npz")
    mask_w = morphology["W_backbone"].astype(bool)
    volume = BOX_XY * BOX_XY * MATERIAL_Z
    volume_w = float(mask_w.mean() * volume)
    volume_ag = volume - volume_w
    sys.path.insert(0, str(PACK / "scripts"))
    import build_ag_w_qsgs as builder

    w_pool, _ = builder.make_lattice_sites("W")
    ag_pool, _ = builder.make_lattice_sites("Ag")
    density_w = n_w / volume_w
    density_ag = n_ag / volume_ag
    density_w_target = len(w_pool) / volume
    density_ag_target = len(ag_pool) / volume
    density_error = {
        "W": density_w / density_w_target - 1.0,
        "Ag": density_ag / density_ag_target - 1.0,
    }

    z_out = int(np.count_nonzero((positions[:, 2] < -1e-8) | (positions[:, 2] > MATERIAL_Z + 1e-8)))
    phase_positions = {"W": positions[types == 1], "Ag": positions[types == 2]}
    distances: dict[str, float] = {}
    below: dict[str, int] = {}
    for element, key in (("W", "W_W"), ("Ag", "Ag_Ag")):
        tree = cKDTree(phase_positions[element], boxsize=[BOX_XY, BOX_XY, 1000.0])
        nearest = tree.query(phase_positions[element], k=2, workers=-1)[0][:, 1]
        distances[key] = float(nearest.min())
        below[key] = int(np.count_nonzero(nearest < FLOORS[key]) // 2)
    cross = cKDTree(phase_positions["W"], boxsize=[BOX_XY, BOX_XY, 1000.0]).query(
        phase_positions["Ag"], workers=-1
    )[0]
    distances["Ag_W"] = float(cross.min())
    below["Ag_W"] = int(np.count_nonzero(cross < FLOORS["Ag_W"]))

    moving = positions[:, 2] > args.fixed_z + 1e-8
    force_magnitude = np.linalg.norm(forces[moving], axis=1)
    topology = {
        element: builder.topology(phase_positions[element], element)
        for element in ("W", "Ag")
    }
    minimizer = minimizer_status(log_path, args.iteration_cap)
    gates = {
        "mass_ratio_within_0.1_percent": abs(mass_ratio - 1.0) <= 1.0e-3,
        "phase_number_density_within_1_percent": all(abs(v) <= 0.01 for v in density_error.values()),
        "largest_component_share_at_least_98_percent": all(
            topology[e]["largest_component_share"] >= 0.98 for e in ("W", "Ag")
        ),
        "both_phases_xy_winding_rank_2": all(topology[e]["xy_winding_rank"] == 2 for e in ("W", "Ag")),
        "both_phases_span_z": all(topology[e]["largest_component_spans_z"] for e in ("W", "Ag")),
        "no_frozen_distance_screen_violations": all(value == 0 for value in below.values()),
        "no_atoms_outside_material_z_0_200_A": z_out == 0,
        "maximum_force_at_most_0.1_eV_A": float(force_magnitude.max()) <= 0.1,
            "rms_force_at_most_0.02_eV_A": float(np.sqrt(np.mean(force_magnitude**2))) <= 0.02,
        "spring_free_efs_minimizer_converged": bool(minimizer["converged"]),
        "atom_count_matches_input": len(ids) == 482601,
    }
    passed = all(gates.values())
    report = {
        "title": "Ag-W QSGS v5 spring-free Dai EFS static audit",
        "status": "STATIC_EFS_PASS" if passed else "STATIC_EFS_FAIL",
        "md_start_gate": "PASS" if passed else "FAIL",
        "thermal_dynamics_started": False,
        "potential": "Dai 2006 EFS converted to LAMMPS eam/fs",
        "boundary": "p p f",
        "fixed_atoms_during_static_minimization": 0,
        "input": {"dump": dump_path.name, "data": data_path.name, "data_sha256": sha256(data_path)},
        "counts": {"W": n_w, "Ag": n_ag, "total": int(len(ids))},
        "mass_ratio_Ag_over_W": mass_ratio,
        "phase_number_density_A^-3": {"W": density_w, "Ag": density_ag},
        "phase_number_density_relative_error": density_error,
        "minimum_distances_A": distances,
        "pairs_below_frozen_screen": below,
        "frozen_distance_screens_A": FLOORS,
        "position_bounds_A": {
            "x": [float(positions[:, 0].min()), float(positions[:, 0].max())],
            "y": [float(positions[:, 1].min()), float(positions[:, 1].max())],
            "z": [float(positions[:, 2].min()), float(positions[:, 2].max())],
        },
        "atoms_outside_material_z_0_200_A": z_out,
        "forces_eV_A": {
            "maximum_vector": float(force_magnitude.max()),
            "rms_vector": float(np.sqrt(np.mean(force_magnitude**2))),
            "p99_vector": float(np.quantile(force_magnitude, 0.99)),
            "moving_atom_count": int(np.count_nonzero(moving)),
            "fixed_bottom_z_A": args.fixed_z,
        },
        "minimizer": minimizer,
        "topology": topology,
        "gates": gates,
        "interpretation": (
            "A static EFS pass is only permission to begin a small no-source 300 K stability pilot."
            if passed else "Do not start 300 K MD or TTM from this structure. Preserve the failed audit and rebuild."
        ),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "md_start_gate": report["md_start_gate"],
                      "gates": gates, "minimum_distances_A": distances,
                      "forces_eV_A": report["forces_eV_A"], "minimizer": minimizer}, indent=2))


if __name__ == "__main__":
    main()
