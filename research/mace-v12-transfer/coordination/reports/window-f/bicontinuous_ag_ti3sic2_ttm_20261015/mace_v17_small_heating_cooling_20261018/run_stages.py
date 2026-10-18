#!/usr/bin/env python3
"""Run V17-only relaxation, 300 K baseline, and local heat/cool stages."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import resource
import time
import numpy as np
import torch
from ase import units
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import FIRE
from ase.neighborlist import neighbor_list
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parent
REPO = Path.cwd()
MODEL = REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_03/models/MACE_periodic_v17_provisional_screening.model"
EXPECTED_MODEL_SHA = "559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3"
ORIGINAL = REPO / "research/mace-v12-transfer/coordination/reports/window-f/bicontinuous_ag_ti3sic2_ttm_20261015/exploratory_melt_cool_20261018/input/NEW_MODEL_ATOMS.data"
EXPECTED_ORIGINAL_SHA = "e2b62a68c46e23037370a6ca20951aa0774edea2d728d3c3aadbcba4b3755674"
DT_FS = float(os.environ.get("MACE_SMALL_DT_FS", "1.0"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def attach_model(atoms):
    if sha(MODEL) != EXPECTED_MODEL_SHA:
        raise RuntimeError("V17 weight SHA mismatch; refusing to load another model")
    if sha(ORIGINAL) != EXPECTED_ORIGINAL_SHA:
        raise RuntimeError("Original 39,486-atom source SHA mismatch")
    support = np.asarray(atoms.arrays["support_fixed"], dtype=bool)
    atoms.set_constraint(FixAtoms(indices=np.flatnonzero(support)))
    torch.set_num_threads(4)
    atoms.calc = MACECalculator(model_paths=str(MODEL), device="cpu", default_dtype="float64", enable_cueq=False)


def phase_mask(atoms, phase):
    return np.asarray(atoms.arrays["phase"]) == phase


def min_pair_metrics(atoms):
    phase = np.asarray(atoms.arrays["phase"])
    symbols = np.asarray(atoms.get_chemical_symbols())
    i, j, d = neighbor_list("ijd", atoms, 4.5)
    m = i < j
    out = {}
    specs = {"Ti_C": ("Ti", "C", None), "Ag_Ag": ("Ag", "Ag", "Ag"),
             "Ag_TSC": (None, None, "cross")}
    for name, (s1, s2, mode) in specs.items():
        if mode == "cross":
            keep = phase[i] != phase[j]
        elif s1 == s2:
            keep = (symbols[i] == s1) & (symbols[j] == s2)
        else:
            keep = ((symbols[i] == s1) & (symbols[j] == s2)) | ((symbols[i] == s2) & (symbols[j] == s1))
        keep &= m
        out["min_" + name + "_A"] = float(d[keep].min()) if np.any(keep) else float("nan")
    return out


def forces_and_metrics(atoms):
    forces = np.asarray(atoms.get_forces(), dtype=float)
    norms = np.linalg.norm(forces, axis=1)
    if not np.all(np.isfinite(forces)) or not np.all(np.isfinite(atoms.positions)):
        raise FloatingPointError("NaN/Inf in coordinates or MACE forces")
    row = {
        "energy_eV": float(atoms.get_potential_energy()),
        "Fmax_eV_A": float(norms.max()),
        "P95_eV_A": float(np.percentile(norms, 95)),
        "P99_eV_A": float(np.percentile(norms, 99)),
        "Fmax_Ag_eV_A": float(norms[phase_mask(atoms, "Ag")].max()),
        "Fmax_TSC_eV_A": float(norms[phase_mask(atoms, "Ti3SiC2")].max()),
        **min_pair_metrics(atoms),
    }
    return forces, norms, row


def set_observation(atoms, forces, norms, stage, step, time_ps, extra=None):
    atoms.arrays["mace_force_eV_A"] = np.asarray(forces, dtype=float).copy()
    atoms.arrays["mace_force_norm_eV_A"] = np.asarray(norms, dtype=float).copy()
    atoms.info["stage"] = stage
    atoms.info["step"] = int(step)
    atoms.info["time_ps"] = float(time_ps)
    atoms.info["potential"] = "MACE_V17_provisional"
    if extra:
        for k, v in extra.items():
            atoms.info[k] = v


def write_snapshot(path, atoms, forces, norms, stage, step, time_ps, extra=None):
    set_observation(atoms, forces, norms, stage, step, time_ps, extra)
    write(path, atoms, format="extxyz", append=Path(path).exists())


def rss_gib():
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0 / 1024.0)


def static_stage():
    atoms = read(ROOT / "SMALL_MODEL_INITIAL.extxyz", format="extxyz")
    attach_model(atoms)
    t0 = time.perf_counter()
    forces, norms, row = forces_and_metrics(atoms)
    elapsed = time.perf_counter() - t0
    symbols = np.asarray(atoms.get_chemical_symbols())
    phase = np.asarray(atoms.arrays["phase"])
    rows = []
    for i in range(len(atoms)):
        rows.append({"id": int(atoms.arrays["id"][i]), "element": symbols[i], "phase": phase[i],
                     "z_A": float(atoms.positions[i, 2]), "force_norm_eV_A": float(norms[i]),
                     "fx_eV_A": float(forces[i, 0]), "fy_eV_A": float(forces[i, 1]), "fz_eV_A": float(forces[i, 2]),
                     "support_fixed": int(atoms.arrays["support_fixed"][i])})
    with (ROOT / "INITIAL_STATIC_FORCES.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    set_observation(atoms, forces, norms, "initial_static", 0, 0.0, row)
    write(ROOT / "SMALL_MODEL_STATIC.extxyz", atoms, format="extxyz")
    result = {"model": str(MODEL.relative_to(REPO)), "model_sha256": sha(MODEL), "device": "CPU", "dtype": "float64",
              "MACE_version": "0.3.16", "torch_version": torch.__version__, "torch_threads": torch.get_num_threads(),
              "atom_count": len(atoms), "original_structure_sha256": sha(ORIGINAL), **row,
              "force_energy_wall_seconds": elapsed, "peak_RSS_GiB": rss_gib(), "status": "finite_static_force_computed"}
    (ROOT / "INITIAL_STATIC_REPORT.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


def relax_stage(steps=60):
    atoms = read(ROOT / "SMALL_MODEL_INITIAL.extxyz", format="extxyz")
    attach_model(atoms)
    forces, norms, initial = forces_and_metrics(atoms)
    start = time.perf_counter()
    traj = ROOT / "RELAXATION.extxyz"
    if traj.exists(): traj.unlink()
    csv_path = ROOT / "RELAXATION_METRICS.csv"
    records = []
    dyn = FIRE(atoms, dt=0.10, maxstep=0.06, downhill_check=False, logfile=None)
    for step in range(steps + 1):
        forces, norms, row = forces_and_metrics(atoms)
        row.update({"step": step, "time_wall_s": time.perf_counter() - start})
        records.append(row)
        if step % 5 == 0 or step == steps:
            write_snapshot(traj, atoms, forces, norms, "relaxation", step, 0.0)
            print("relax", step, "Fmax", row["Fmax_eV_A"], "PE", row["energy_eV"], flush=True)
        if row["Fmax_eV_A"] > 100.0:
            raise RuntimeError(f"Relaxation force runaway at step {step}: {row['Fmax_eV_A']}")
        if row["min_Ti_C_A"] < 1.45 or row["min_Ag_TSC_A"] < 1.70 or row["min_Ag_Ag_A"] < 1.70:
            raise RuntimeError(f"Severe short contact during relaxation at step {step}: {row}")
        if step == steps or row["Fmax_eV_A"] < 0.35:
            break
        # Reuse the force vector just evaluated above. This avoids a second
        # calculator request per minimizer step on the CPU-only environment.
        dyn.step(forces)
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
    forces, norms, final = forces_and_metrics(atoms)
    set_observation(atoms, forces, norms, "relaxed_checkpoint", records[-1]["step"], 0.0, final)
    write(ROOT / "RELAXED_STRUCTURE.extxyz", atoms, format="extxyz")
    result = {"steps": int(records[-1]["step"]), "initial": initial, "final": final,
              "wall_seconds": time.perf_counter() - start, "peak_RSS_GiB": rss_gib(),
              "atom_count_before": 1926, "atom_count_after": len(atoms), "fixed_support_atoms": int(np.sum(atoms.arrays["support_fixed"]))}
    result["status"] = "relaxation_short_pass" if final["Fmax_eV_A"] < 10.0 and final["min_Ti_C_A"] > 1.45 and final["min_Ag_TSC_A"] > 1.70 else "review_needed"
    (ROOT / "RELAXATION_REPORT.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


def regional_temperature(atoms, mask):
    ids = np.flatnonzero(mask)
    if len(ids) < 2: return float("nan")
    v = atoms.get_velocities()[ids]
    m = atoms.get_masses()[ids, None]
    vcm = np.sum(m * v, axis=0) / np.sum(m)
    dof = max(1, 3 * len(ids) - 3)
    return float(np.sum(m * (v - vcm) ** 2) / (dof * units.kB))


def base_run(atoms, nsteps, dt_fs, stage, trajpath, metrics_path, temperature_K=None,
             friction_per_fs=0.01, hot_mask=None, sink_mask=None, hot_temperature_K=1800.0,
             sample_every=10, seed=1919, append_existing=False, time_offset_ps=0.0):
    rng = np.random.default_rng(seed)
    support = np.asarray(atoms.arrays["support_fixed"], dtype=bool)
    if temperature_K is not None:
        MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_K, rng=rng)
        v = atoms.get_velocities(); v[support] = 0.0; atoms.set_velocities(v)
        t0 = time.perf_counter()
        dyn = Langevin(atoms, timestep=dt_fs * units.fs, temperature_K=temperature_K,
                       friction=friction_per_fs / units.fs, fixcm=False, rng=rng)
        # A global baseline thermostat is allowed only in the baseline stage.
        temp_array = None; fric_array = None
    else:
        t0 = time.perf_counter()
        dyn = None
        # Region-dependent Langevin bath: NVE elsewhere; hot source and cold
        # substrate sink are the only stochastic reservoirs during heat/cool.
        is_cooling = stage == "cooling"
        temp_array = np.full((len(atoms), 1), 300.0, dtype=float)
        fric_array = np.zeros((len(atoms), 1), dtype=float)
        fric_array[sink_mask] = friction_per_fs / units.fs
        if not is_cooling:
            temp_array[hot_mask] = hot_temperature_K
            fric_array[hot_mask] = friction_per_fs / units.fs
        dyn = Langevin(atoms, timestep=dt_fs * units.fs, temperature_K=temp_array,
                       friction=fric_array, fixcm=False, rng=rng)
    if Path(trajpath).exists() and not append_existing: Path(trajpath).unlink()
    records = []
    start = time.perf_counter()
    for step in range(nsteps + 1):
        forces, norms, row = forces_and_metrics(atoms)
        row.update({"stage": stage, "step": step, "stage_time_ps": step * dt_fs / 1000.0,
                    "time_ps": time_offset_ps + step * dt_fs / 1000.0,
                    "T_system_K": float(atoms.get_temperature()), "time_wall_s": time.perf_counter() - start,
                    "n_hot_Ag": int(np.sum(hot_mask)) if hot_mask is not None else 0,
                    "T_hot_Ag_K": regional_temperature(atoms, hot_mask) if hot_mask is not None else float("nan"),
                    "T_Ag_K": regional_temperature(atoms, phase_mask(atoms, "Ag")),
                    "T_TSC_K": regional_temperature(atoms, phase_mask(atoms, "Ti3SiC2")),
                    "T_sink_K": regional_temperature(atoms, sink_mask) if sink_mask is not None else float("nan"),
                    "atoms_outside_z_margin": int(np.sum((atoms.positions[:, 2] < -3.0) | (atoms.positions[:, 2] > atoms.cell.lengths()[2] + 3.0)))} )
        records.append(row)
        if step % sample_every == 0 or step == nsteps:
            # Store the complete scalar snapshot audit in extxyz. In particular,
            # per-frame phase forces and minimum contacts must not inherit stale
            # values from the input structure in OVITO.
            info = {k: v for k, v in row.items() if isinstance(v, (str, int, float, bool))}
            write_snapshot(trajpath, atoms, forces, norms, stage, step, row["time_ps"], info)
            print(stage, step, "t_ps", row["time_ps"], "T", row["T_system_K"], "Thot", row["T_hot_Ag_K"],
                  "FAg", row["Fmax_Ag_eV_A"], "FTSC", row["Fmax_TSC_eV_A"], "dmin", row["min_Ti_C_A"], flush=True)
        # Stop only for clear numerical collapse, not for a merely nonzero force.
        if row["Fmax_eV_A"] > 100.0 or row["min_Ti_C_A"] < 1.45 or row["min_Ag_TSC_A"] < 1.70 or row["min_Ag_Ag_A"] < 1.70:
            row["early_stop_reason"] = "severe_force_or_short_contact"
            break
        if row["atoms_outside_z_margin"] > 0:
            row["early_stop_reason"] = "atom_escaped_open_z_boundary"
            break
        if step == nsteps: break
        dyn.run(1)
    with Path(metrics_path).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
    forces, norms, final = forces_and_metrics(atoms)
    set_observation(atoms, forces, norms, stage + "_checkpoint", records[-1]["step"], records[-1]["time_ps"], final)
    checkpoint = ROOT / (stage.upper() + "_CHECKPOINT.extxyz")
    write(checkpoint, atoms, format="extxyz")
    result = {"stage": stage, "requested_steps": int(nsteps), "completed_steps": int(records[-1]["step"]),
              "dt_fs": float(dt_fs), "simulated_time_ps": float(records[-1]["time_ps"]),
              "wall_seconds": float(time.perf_counter() - start), "peak_RSS_GiB": rss_gib(),
              "atom_count": len(atoms), "final": final, "first": records[0], "last": records[-1],
              "status": records[-1].get("early_stop_reason", "completed_requested_steps"), "checkpoint": checkpoint.name}
    (ROOT / (stage.upper() + "_REPORT.json")).write_text(json.dumps(result, indent=2, allow_nan=True) + "\n")
    print(json.dumps(result, allow_nan=True), flush=True)
    return result, atoms


def baseline_stage(steps=100):
    atoms = read(ROOT / "RELAXED_STRUCTURE.extxyz", format="extxyz")
    attach_model(atoms)
    result, _ = base_run(atoms, steps, 1.0, "baseline_300K", ROOT / "BASELINE_300K.extxyz",
                         ROOT / "BASELINE_300K_METRICS.csv", temperature_K=300.0,
                         friction_per_fs=0.01, sample_every=10, seed=20261018)


def lateral_distance_to_center(atoms, center_frac=(0.5, 0.5)):
    cell = np.asarray(atoms.cell.array)
    frac = atoms.positions @ np.linalg.inv(cell)
    duv = frac[:, :2] - np.asarray(center_frac)[None, :]
    # Search nearby image shifts; this gives the true shortest distance for the
    # oblique hexagonal two-dimensional cell.
    best = np.full(len(atoms), np.inf)
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            dv = (duv + np.array([i, j]))[:, 0, None] * cell[0, :2] + (duv + np.array([i, j]))[:, 1, None] * cell[1, :2]
            best = np.minimum(best, np.linalg.norm(dv, axis=1))
    return best


def heatcool_stage(heat_steps=300, cool_steps=500):
    atoms = read(ROOT / "BASELINE_300K_CHECKPOINT.extxyz", format="extxyz")
    attach_model(atoms)
    n_tsc = int(np.sum(phase_mask(atoms, "Ti3SiC2")))
    n_ag = int(np.sum(phase_mask(atoms, "Ag")))
    audit = json.loads((ROOT / "MODEL_GEOMETRY_AUDIT.json").read_text())
    z_ag0 = audit["material_z_max_A"] - (6 - 1) * audit["Ag111_interplane_spacing_A"]
    layer = np.asarray(atoms.arrays["Ag111_layer"], dtype=int)
    radial = lateral_distance_to_center(atoms)
    hot = phase_mask(atoms, "Ag") & (layer >= 3) & (radial <= 8.0)
    sink = phase_mask(atoms, "Ti3SiC2") & (atoms.positions[:, 2] > 0.05) & (atoms.positions[:, 2] <= 4.50)
    support = np.asarray(atoms.arrays["support_fixed"], dtype=bool)
    sink &= ~support
    if hot.sum() < 20 or sink.sum() < 50:
        raise RuntimeError(f"Unexpected region sizes: hotAg={hot.sum()}, sink={sink.sum()}")
    run_dir = ROOT / "HEAT_COOL.extxyz"
    if run_dir.exists(): run_dir.unlink()
    heat, atoms = base_run(atoms, heat_steps, 1.0, "heating", run_dir,
                           ROOT / "HEATING_METRICS.csv", temperature_K=None,
                           friction_per_fs=0.02, hot_mask=hot, sink_mask=sink,
                           hot_temperature_K=1900.0, sample_every=10, seed=20261019)
    # Same trajectory file: append a phase separator frame, then turn the local
    # Ag thermostat off. The 300 K substrate sink remains active as the cooling path.
    cool, atoms = base_run(atoms, cool_steps, 1.0, "cooling", run_dir,
                           ROOT / "COOLING_METRICS.csv", temperature_K=None,
                           friction_per_fs=0.02, hot_mask=hot, sink_mask=sink,
                           hot_temperature_K=1900.0, sample_every=10, seed=20261020,
                           append_existing=True, time_offset_ps=heat_steps / 1000.0)
    combined = {"heating": heat, "cooling": cool, "hot_source_atom_count": int(hot.sum()),
                "hot_source_definition": "Ag within 8 A lateral radius of cell center and in top 3 Ag(111) planes; Langevin target 1900 K, friction 0.02 fs^-1",
                "sink_atom_count": int(sink.sum()),
                "sink_definition": "mobile Ti3SiC2 atoms at 0 < z <= 4.50 A, Langevin target 300 K; bottom z=0 Ti plane fixed",
                "other_atoms": "unthermostatted NVE; no global thermostat in heating/cooling",
                "note": "Cooling is heat removal through a local fixed-temperature substrate sink, not an isolated-system quench.",
                "pbc": [True, True, False], "dt_fs": 1.0, "model_sha256": sha(MODEL)}
    (ROOT / "HEAT_COOL_REPORT.json").write_text(json.dumps(combined, indent=2, allow_nan=True) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["static", "relax", "baseline", "heatcool"])
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--heat-steps", type=int, default=300)
    ap.add_argument("--cool-steps", type=int, default=500)
    args = ap.parse_args()
    if args.stage == "static": static_stage()
    elif args.stage == "relax": relax_stage(args.steps or 60)
    elif args.stage == "baseline": baseline_stage(args.steps or 100)
    elif args.stage == "heatcool": heatcool_stage(args.heat_steps, args.cool_steps)


if __name__ == "__main__":
    main()
