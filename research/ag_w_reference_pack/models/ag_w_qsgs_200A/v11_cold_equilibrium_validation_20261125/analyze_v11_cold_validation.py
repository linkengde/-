#!/usr/bin/env python3
"""Summarize the V11 force-relaxation and staged cold-equilibrium diagnostics."""
from __future__ import annotations

import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


HERE = Path(__file__).resolve().parent
MASS_U = {1: 183.84, 2: 107.8682}  # LAMMPS type 1=W, 2=Ag
KB_EV_K = 8.617333262145e-5
M_VV2E_METAL = 1.036426965e-4
LXY_A = 200.0
BASELINE_COMMIT = "40d6b5676c447c204026346957fee05714561817"
REPO = HERE.parents[4]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_dump(path: Path) -> dict[str, np.ndarray]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        lines = stream.read().splitlines()
    atom_header = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS"))
    columns = lines[atom_header].split()[2:]
    values = np.loadtxt(lines[atom_header + 1 :], dtype=float)
    if values.ndim == 1:
        values = values[None, :]
    return {name: values[:, i] for i, name in enumerate(columns)}


def fixed_ids_from_include(path: Path) -> np.ndarray:
    ids: list[int] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[:3] == ["group", "fixed", "id"]:
            ids.extend(int(item) for item in fields[3:])
    return np.asarray(ids, dtype=np.int64)


def phase_temperatures(data: dict[str, np.ndarray], fixed_ids: np.ndarray) -> dict:
    atom_ids = data["id"].astype(np.int64)
    types = data["type"].astype(int)
    mobile = ~np.isin(atom_ids, fixed_ids)
    velocity = np.column_stack([data["vx"], data["vy"], data["vz"]])
    output = {}
    groups = {
        "mobile": mobile,
        "Ag": (types == 2) & mobile,
        "W": (types == 1) & mobile,
    }
    for label, mask in groups.items():
        masses = np.asarray([MASS_U[int(t)] for t in types[mask]])
        group_velocity = velocity[mask]
        com_velocity = np.average(group_velocity, axis=0, weights=masses)
        peculiar_velocity = group_velocity - com_velocity
        dof = 3 * len(masses) - 3
        temperature = M_VV2E_METAL * np.sum(masses[:, None] * peculiar_velocity**2) / (dof * KB_EV_K)
        output[label] = {
            "atoms": int(mask.sum()),
            "temperature_K": float(temperature),
            "com_velocity_A_per_ps": com_velocity.tolist(),
        }
    return output


def geometry_and_forces(data: dict[str, np.ndarray], fixed_ids: np.ndarray) -> dict:
    atom_ids = data["id"].astype(np.int64)
    types = data["type"].astype(int)
    mobile = ~np.isin(atom_ids, fixed_ids)
    xyz = np.column_stack([data["x"], data["y"], data["z"]])
    xyz[:, :2] %= LXY_A
    xyz[:, 2] += 1.0  # leave z nonperiodic using a large dummy period
    forces = np.column_stack([data["fx"], data["fy"], data["fz"]])
    force_norm = np.linalg.norm(forces, axis=1)
    trees = {t: cKDTree(xyz[types == t], boxsize=[LXY_A, LXY_A, 1.0e6]) for t in (1, 2)}
    output = {}
    for atom_type, element in ((2, "Ag"), (1, "W")):
        indices = np.where(types == atom_type)[0]
        same, _ = trees[atom_type].query(xyz[indices], k=2)
        cross, _ = trees[3 - atom_type].query(xyz[indices], k=1)
        phase_forces = force_norm[indices[mobile[indices]]]
        output[element] = {
            "atom_count": int(len(indices)),
            "same_phase_nearest_min_A": float(same[:, 1].min()),
            "same_phase_nearest_p01_A": float(np.quantile(same[:, 1], 0.01)),
            "cross_phase_nearest_min_A": float(cross.min()),
            "cross_phase_nearest_p01_A": float(np.quantile(cross, 0.01)),
            "mobile_force_max_eV_A": float(phase_forces.max()),
            "mobile_force_rms_eV_A": float(np.sqrt(np.mean(phase_forces**2))),
            "mobile_force_p99_eV_A": float(np.quantile(phase_forces, 0.99)),
            "mobile_atoms_force_gt_2_eV_A": int(np.count_nonzero(phase_forces > 2.0)),
            "z_min_A": float(data["z"][indices].min()),
            "z_max_A": float(data["z"][indices].max()),
        }
    return output


def read_series(path: Path) -> np.ndarray:
    return np.loadtxt(path, comments="#", dtype=float)


def nvt_summary(name: str, series: np.ndarray) -> dict:
    # Each archived run has 200 NVE steps after its NVT segment.
    nvt = series[series[:, 0] <= series[-1, 0] - 200]
    tail = nvt[nvt[:, 0] >= nvt[-1, 0] - 400]
    # With samples every 20 steps and dt=0.00025 ps, a 0.1 ps window
    # contains 20 intervals (21 samples), not merely five intervals.
    start_0p1ps = tail[0]
    tail_slope = np.polyfit(tail[:, 0] * 0.00025, tail[:, 8], 1)[0]
    return {
        "stage": name,
        "start_step": int(nvt[0, 0]),
        "end_step": int(nvt[-1, 0]),
        "start_T_mobile_Ag_W_K": nvt[0, 1:4].tolist(),
        "end_T_mobile_Ag_W_K": nvt[-1, 1:4].tolist(),
        "end_W_minus_Ag_K": float(nvt[-1, 3] - nvt[-1, 2]),
        "last_0p1ps_mobile_dT_K": float(nvt[-1, 1] - start_0p1ps[1]),
        "last_0p1ps_W_minus_Ag_change_K": float(
            (nvt[-1, 3] - nvt[-1, 2]) - (start_0p1ps[3] - start_0p1ps[2])
        ),
        "last_0p1ps_potential_energy_change_eV": float(nvt[-1, 8] - start_0p1ps[8]),
        "last_0p1ps_potential_energy_slope_eV_ps": float(tail_slope),
    }


def precise_nve_rows(path: Path) -> list[list[float]]:
    rows = []
    pattern = re.compile(r"^\s*(6600|6800)\s+")
    for line in path.read_text().splitlines():
        if pattern.match(line):
            rows.append([float(item) for item in line.split()])
    if len(rows) != 2:
        raise RuntimeError(f"Expected precise NVE start/end rows, found {len(rows)}")
    return rows


def main() -> None:
    fixed_ids = fixed_ids_from_include(HERE / "v11_cold_fixed_ids.include")
    final = read_dump(HERE / "v11_nve_precision_check.dump.gz")
    relaxed = read_dump(HERE / "v11_force_relaxed.dump.gz")
    nve_rows = precise_nve_rows(HERE / "v11_nve_precision_check.lammps.log")

    final_ids = final["id"].astype(np.int64)
    relaxed_ids = relaxed["id"].astype(np.int64)
    final_order = np.argsort(final_ids)
    relaxed_order = np.argsort(relaxed_ids)
    if not np.array_equal(final_ids[final_order], relaxed_ids[relaxed_order]):
        raise RuntimeError("Atom IDs changed between force-relaxed candidate and final NVE dump")
    displacement = np.column_stack(
        [final[k][final_order] - relaxed[k][relaxed_order] for k in ("x", "y", "z")]
    )
    displacement[:, :2] -= LXY_A * np.rint(displacement[:, :2] / LXY_A)
    displacement_norm = np.linalg.norm(displacement, axis=1)
    mobile = ~np.isin(final_ids[final_order], fixed_ids)
    final_types = final["type"][final_order].astype(int)

    nve_start, nve_end = nve_rows
    # Thermo columns: step, time, Tmobile, TAg, TW, PE, KE, total E, max force, RMS force.
    nve = {
        "start_step": int(nve_start[0]),
        "end_step": int(nve_end[0]),
        "duration_ps": float(nve_end[1] - nve_start[1]),
        "total_energy_start_eV": float(nve_start[7]),
        "total_energy_end_eV": float(nve_end[7]),
        "total_energy_drift_eV": float(nve_end[7] - nve_start[7]),
        "potential_energy_change_eV": float(nve_end[5] - nve_start[5]),
        "kinetic_energy_change_eV": float(nve_end[6] - nve_start[6]),
        "start_T_mobile_Ag_W_K": nve_start[2:5],
        "end_T_mobile_Ag_W_K": nve_end[2:5],
        "end_W_minus_Ag_K": float(nve_end[4] - nve_end[3]),
        "end_mobile_force_max_eV_A": float(nve_end[8]),
        "end_mobile_force_rms_eV_A": float(nve_end[9]),
    }
    report = {
        "title": "Ag-W V11 cold equilibrium validation diagnostic",
        "status": "COLD_BASELINE_FAIL",
        "scope": "Independent V11-only candidate; no laser, no TTM, no V5-V10 files modified.",
        "baseline_commit": BASELINE_COMMIT,
        "lammps_version": "22 Jul 2025 - Update 4",
        "input_sha256": {
            "v11_generated_geometry": sha256_file(
                REPO / "research/ag_w_reference_pack/models/ag_w_qsgs_200A/v11_generated_data/agw_v11_gyroid.data.gz"
            ),
            "Dai_2006_AgW_EFS": sha256_file(
                REPO / "research/ag_w_reference_pack/potential_files/derived/Dai_2006_EFS_AgW.eam.fs"
            ),
        },
        "composition": {
            "atoms": int(len(final_ids)),
            "fixed_ids": int(len(fixed_ids)),
            "mobile_atoms": int(len(final_ids) - len(fixed_ids)),
            "Ag_atoms": int(np.count_nonzero(final["type"] == 2)),
            "W_atoms": int(np.count_nonzero(final["type"] == 1)),
            "Ag_over_W_mass_ratio": float(291274 * 107.8682 / (170906 * 183.84)),
        },
        "force_relaxation": {
            "stopping_reason": "maximum iteration count (2000/2000), not force convergence",
            "potential_energy_start_eV": -2347019.88742023,
            "potential_energy_end_eV": -2357564.75665815,
            "energy_change_eV": -10544.86923792,
            "mobile_force_max_eV_A": 1.2388023,
            "mobile_force_rms_eV_A": 0.0099953913,
            "Ag_force_rms_eV_A": 0.0078860349,
            "W_force_rms_eV_A": 0.012826995,
        },
        "velocity_temperature_check": phase_temperatures(final, fixed_ids),
        "nvt_stages": [
            nvt_summary("initial_0p5ps", read_series(HERE / "v11_cold_staged_timeseries.dat")),
            nvt_summary("continuation_0p5ps", read_series(HERE / "v11_cold_continue_timeseries.dat")),
            nvt_summary("continuation2_0p5ps", read_series(HERE / "v11_cold_continue2_timeseries.dat")),
        ],
        "precise_nve_check": nve,
        "final_geometry_and_forces": geometry_and_forces(final, fixed_ids),
        "mobile_displacement_from_force_relaxed_candidate_A": {
            "mean": float(displacement_norm[mobile].mean()),
            "p99": float(np.quantile(displacement_norm[mobile], 0.99)),
            "maximum": float(displacement_norm[mobile].max()),
            "Ag_mean": float(displacement_norm[mobile & (final_types == 2)].mean()),
            "W_mean": float(displacement_norm[mobile & (final_types == 1)].mean()),
        },
        "interpretation": [
            "Initial Ag/W Maxwell temperatures were both approximately 300 K; type/mass mapping is consistent, and one shared NVT group/target was used.",
            "A 2000-step stricter static relaxation substantially reduced RMS forces, but stopped at the iteration cap and retained a 1.24 eV/A maximum force.",
            "The Ag/W split narrowed at low temperature, then increased to about 20.7 K by the high-precision NVE endpoint near 286 K.",
            "The 0.05 ps NVE total-energy drift is about 4.2e-4 eV, while potential energy fell by about 83.8 eV and kinetic energy rose by the same amount; energy conservation therefore does not imply a settled structure.",
            "Force hotspots and short same/cross-phase contacts remain, so this is not an accepted 300 K baseline.",
        ],
    }
    (HERE / "v11_cold_diagnostics_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
