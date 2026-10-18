#!/usr/bin/env python3
"""MACE V17-only reduced-slab baseline and local heat/cool dynamics."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import resource
import time
import numpy as np
import torch
from ase import units
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.neighborlist import neighbor_list
from mace.calculators import MACECalculator
from scipy.special import sph_harm_y
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parent
REPO = Path.cwd()
MODEL = REPO / "research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_03/models/MACE_periodic_v17_provisional_screening.model"
ORIGINAL = REPO / "research/mace-v12-transfer/coordination/reports/window-f/bicontinuous_ag_ti3sic2_ttm_20261015/exploratory_melt_cool_20261018/input/NEW_MODEL_ATOMS.data"
EXPECTED_MODEL = "559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3"
EXPECTED_ORIGINAL = "e2b62a68c46e23037370a6ca20951aa0774edea2d728d3c3aadbcba4b3755674"
DT_FS = 0.5
SAMPLE_STEPS = 10  # 0.005 ps per record
Q6_CUTOFF = 3.45


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def attach(atoms):
    if sha(MODEL) != EXPECTED_MODEL or sha(ORIGINAL) != EXPECTED_ORIGINAL:
        raise RuntimeError("V17 weight or original 39,486-atom source SHA mismatch")
    support = np.asarray(atoms.arrays["support_fixed"], bool)
    atoms.set_constraint(FixAtoms(indices=np.flatnonzero(support)))
    torch.set_num_threads(4)
    atoms.calc = MACECalculator(model_paths=str(MODEL), device="cpu", default_dtype="float64", enable_cueq=False)
    return atoms.calc


def p95(x):
    return float(np.percentile(x, 95)) if len(x) else float("nan")


def distances(atoms):
    sy = np.asarray(atoms.get_chemical_symbols())
    ph = atoms.arrays["phase"].astype(str)
    i, j, d = neighbor_list("ijd", atoms, 4.5)
    pairs = i < j
    tic = pairs & (((sy[i] == "Ti") & (sy[j] == "C")) | ((sy[i] == "C") & (sy[j] == "Ti")))
    agag = pairs & (sy[i] == "Ag") & (sy[j] == "Ag")
    cross = pairs & (ph[i] != ph[j])
    def minimum(mask):
        return float(d[mask].min()) if np.any(mask) else float("nan")
    return {"min_Ti_C_A": minimum(tic), "min_Ag_Ag_A": minimum(agag), "min_Ag_TSC_A": minimum(cross)}


def graph_fraction(atoms, mask, cutoff):
    node_ids = np.flatnonzero(mask)
    i, j = neighbor_list("ij", atoms, cutoff)
    take = mask[i] & mask[j] & (i < j)
    lookup = np.full(len(atoms), -1, int)
    lookup[node_ids] = np.arange(len(node_ids))
    u, v = lookup[i[take]], lookup[j[take]]
    adj = coo_matrix((np.ones(2 * len(u)), (np.r_[u, v], np.r_[v, u])), shape=(len(node_ids), len(node_ids))).tocsr()
    _, labels = connected_components(adj, directed=False)
    largest = int(np.bincount(labels).max()) if len(labels) else 0
    return largest, float(largest / len(node_ids)) if len(node_ids) else 0.0


def ag_structure(atoms):
    sy = np.asarray(atoms.get_chemical_symbols())
    ids = np.flatnonzero(sy == "Ag")
    i, j, vec = neighbor_list("ijD", atoms, Q6_CUTOFF)
    pair = (sy[i] == "Ag") & (sy[j] == "Ag")
    i, j, vec = i[pair], j[pair], vec[pair]
    deg = np.bincount(i, minlength=len(atoms))
    qlm = np.zeros((len(atoms), 13), complex)
    if len(i):
        r = np.linalg.norm(vec, axis=1)
        theta = np.arccos(np.clip(vec[:, 2] / r, -1.0, 1.0))
        phi = np.arctan2(vec[:, 1], vec[:, 0])
        for k, m in enumerate(range(-6, 7)):
            y = sph_harm_y(6, m, theta, phi)
            qlm[:, k] = np.bincount(i, weights=y.real, minlength=len(atoms)) + 1j * np.bincount(i, weights=y.imag, minlength=len(atoms))
        valid = deg > 0
        qlm[valid] /= deg[valid, None]
    q6 = np.zeros(len(atoms))
    q6[ids] = np.sqrt(4 * np.pi / 13 * np.sum(np.abs(qlm[ids]) ** 2, axis=1))
    return ids, deg[ids], q6[ids], i, j


def largest_lowq_component(atoms, selected_ids, i, j, threshold=0.35):
    if not len(selected_ids):
        return 0
    allowed = np.zeros(len(atoms), bool)
    allowed[selected_ids] = True
    keep = allowed[i] & allowed[j]
    index = np.full(len(atoms), -1, int)
    index[selected_ids] = np.arange(len(selected_ids))
    u, v = index[i[keep]], index[j[keep]]
    adj = coo_matrix((np.ones(2 * len(u)), (np.r_[u, v], np.r_[v, u])), shape=(len(selected_ids), len(selected_ids))).tocsr()
    _, labels = connected_components(adj, directed=False)
    return int(np.bincount(labels).max()) if len(labels) else 0


def region_temp(atoms, mask):
    ids = np.flatnonzero(mask)
    if len(ids) < 2:
        return float("nan")
    vel = atoms.get_velocities()[ids]
    mass = atoms.get_masses()[ids, None]
    vcm = np.sum(mass * vel, axis=0) / np.sum(mass)
    dof = max(1, 3 * len(ids) - 3)
    return float(np.sum(mass * (vel - vcm) ** 2) / (dof * units.kB))


def hot_rmsd(atoms, reference, mask):
    cell = np.asarray(atoms.cell.array)
    frac = atoms.positions @ np.linalg.inv(cell)
    frac0 = reference.positions @ np.linalg.inv(np.asarray(reference.cell.array))
    delta = frac - frac0
    delta[:, :2] -= np.rint(delta[:, :2])
    xyz = delta[:, 0, None] * cell[0] + delta[:, 1, None] * cell[1]
    xyz[:, 2] = atoms.positions[:, 2] - reference.positions[:, 2]
    disp = xyz[mask]
    if len(disp) == 0:
        return float("nan")
    disp -= disp.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum(disp * disp, axis=1))))


def metrics(atoms, reference, stage, step, stage_time_ps, global_offset_ps, start_wall):
    force = np.asarray(atoms.get_forces(), float)
    mobile = ~np.asarray(atoms.arrays["support_fixed"], bool)
    fmag = np.linalg.norm(force, axis=1)
    if not np.all(np.isfinite(force)) or not np.all(np.isfinite(atoms.positions)):
        raise FloatingPointError("NaN/Inf coordinates or MACE forces")
    agids, cn, q6, ii, jj = ag_structure(atoms)
    hotmask = np.asarray(atoms.arrays["heat_region_initial"], bool)
    hotids = agids[hotmask[agids]]
    hotmap = hotmask[agids]
    hotlow = agids[hotmap & (q6 < 0.35)]
    largest = largest_lowq_component(atoms, hotlow, ii, jj)
    tsc = atoms.arrays["phase"].astype(str) == "Ti3SiC2"
    tsc_largest, tsc_fraction = graph_fraction(atoms, tsc, 3.25)
    posz = atoms.positions[:, 2]
    outside = int(np.sum((posz < -3.0) | (posz > atoms.cell.lengths()[2] + 3.0)))
    row = {
        "stage": stage, "step": int(step), "stage_time_ps": float(stage_time_ps),
        "time_ps": float(global_offset_ps + stage_time_ps), "wall_seconds": float(time.perf_counter() - start_wall),
        "energy_eV": float(atoms.get_potential_energy()), "T_all_K": float(atoms.get_temperature()),
        "T_Ag_K": region_temp(atoms, atoms.arrays["phase"].astype(str) == "Ag"),
        "T_TSC_K": region_temp(atoms, tsc), "T_hot_Ag_K": region_temp(atoms, hotmask),
        "T_cold_sink_K": region_temp(atoms, np.asarray(atoms.arrays["cold_sink"], bool)),
        "Fmax_mobile_eV_A": float(np.max(fmag[mobile])), "F95_mobile_eV_A": p95(fmag[mobile]),
        "Fmax_Ag_eV_A": float(np.max(fmag[atoms.arrays["phase"].astype(str) == "Ag"])),
        "n_hot_Ag": int(len(hotids)), "hot_Ag_median_CN_3p45": float(np.median(cn[hotmap])) if np.any(hotmap) else float("nan"),
        "hot_Ag_median_q6": float(np.median(q6[hotmap])) if np.any(hotmap) else float("nan"),
        "hot_Ag_low_q6_count": int(len(hotlow)), "hot_Ag_low_q6_fraction": float(len(hotlow) / len(hotids)) if len(hotids) else float("nan"),
        "hot_low_q6_largest_cluster_atoms": int(largest), "hot_Ag_RMSD_from_baseline_A": hot_rmsd(atoms, reference, hotmask),
        "TSC_largest_component_fraction_cut3p25": tsc_fraction, "TSC_largest_component_atoms": tsc_largest,
        "atoms_outside_open_z_margin": outside, **distances(atoms),
    }
    atoms.arrays["mace_force_eV_A"] = force.copy()
    atoms.arrays["mace_force_norm_eV_A"] = fmag.copy()
    atoms.info.update({k: v for k, v in row.items() if isinstance(v, (str, int, float, bool))})
    atoms.info["potential"] = "MACE_V17_provisional_single_potential"
    return row


def append_row(path, row):
    exists = Path(path).exists()
    with Path(path).open("a" if exists else "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row))
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def snapshot(path, atoms):
    write(path, atoms, format="extxyz", append=Path(path).exists())


def checkpoint(path, atoms):
    write(path, atoms, format="extxyz")


def crash_guard(row):
    if row["Fmax_mobile_eV_A"] > 100 or row["min_Ti_C_A"] < 1.45 or row["min_Ag_TSC_A"] < 1.70 or row["min_Ag_Ag_A"] < 1.70:
        return "severe_force_or_short_contact"
    if row["atoms_outside_open_z_margin"]:
        return "atom_escaped_open_z_boundary"
    return None


def initialize(steps):
    path = ROOT / "REDUCED_MODEL.extxyz"
    if not path.exists():
        raise FileNotFoundError("Run prepare_reduced_model.py first")
    if (ROOT / "BASELINE_300K.extxyz").exists():
        raise FileExistsError("Baseline exists; refusing to overwrite prior reduced-run output")
    atoms = read(path)
    calc = attach(atoms)
    static_start = time.perf_counter()
    row0 = metrics(atoms, atoms.copy(), "initial_static", 0, 0.0, 0.0, static_start)
    static_wall = time.perf_counter() - static_start
    rmax = None
    for model in getattr(calc, "models", []):
        if hasattr(model, "r_max"):
            rmax = float(model.r_max)
            break
    write(ROOT / "INITIAL_STATIC.extxyz", atoms, format="extxyz")
    static_report = {**row0, "atom_count": len(atoms), "static_wall_seconds": static_wall,
                     "peak_RSS_GiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024,
                     "MACE_version": "0.3.16", "torch_version": torch.__version__, "dtype": "float64", "device": "CPU", "threads": torch.get_num_threads(),
                     "MACE_r_max_A_if_exposed": rmax, "MACE_model_sha256": sha(MODEL), "status": "finite_static_force_computed"}
    (ROOT / "INITIAL_STATIC_REPORT.json").write_text(json.dumps(static_report, indent=2, allow_nan=True) + "\n")
    print("static", json.dumps(static_report, allow_nan=True), flush=True)
    if row0["Fmax_mobile_eV_A"] > 100 or row0["min_Ti_C_A"] < 1.45 or row0["min_Ag_TSC_A"] < 1.70 or row0["min_Ag_Ag_A"] < 1.70:
        raise RuntimeError("Initial reduced model failed severe-force/overlap safety screen")
    rng = np.random.default_rng(20261020)
    MaxwellBoltzmannDistribution(atoms, temperature_K=300.0, rng=rng)
    v = atoms.get_velocities()
    v[np.asarray(atoms.arrays["support_fixed"], bool)] = 0.0
    atoms.set_velocities(v)
    reference = atoms.copy()
    t0 = time.perf_counter()
    dyn = Langevin(atoms, timestep=DT_FS * units.fs, temperature_K=300.0, friction=0.01 / units.fs,
                   fixcm=False, rng=np.random.default_rng(20261021))
    metrics_path = ROOT / "BASELINE_300K_METRICS.csv"
    force_row = metrics(atoms, reference, "baseline_300K", 0, 0.0, 0.0, t0)
    atoms.info["baseline_hot_q6_median"] = force_row["hot_Ag_median_q6"]
    atoms.info["baseline_reference_time_ps"] = 0.0
    write(ROOT / "BASELINE_300K.extxyz", atoms, format="extxyz")
    append_row(metrics_path, force_row)
    for step in range(min(SAMPLE_STEPS, steps), steps + 1, SAMPLE_STEPS):
        dyn.run(SAMPLE_STEPS)
        row = metrics(atoms, reference, "baseline_300K", step, step * DT_FS / 1000.0, 0.0, t0)
        append_row(metrics_path, row)
        snapshot(ROOT / "BASELINE_300K.extxyz", atoms)
        checkpoint(ROOT / "BASELINE_CHECKPOINT.extxyz", atoms)
        stop = crash_guard(row)
        print("baseline", step, json.dumps(row, allow_nan=True), flush=True)
        if stop:
            (ROOT / "BASELINE_STATUS.txt").write_text(stop + "\n")
            raise RuntimeError("300 K baseline stopped: " + stop)
    report = {"atom_count": len(atoms), "dt_fs": DT_FS, "completed_steps": steps // SAMPLE_STEPS * SAMPLE_STEPS,
              "physical_time_ps": (steps // SAMPLE_STEPS * SAMPLE_STEPS) * DT_FS / 1000.0,
              "wall_seconds": time.perf_counter() - t0, "peak_RSS_GiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024,
              "status": "short_baseline_completed_not_long_term_validation"}
    (ROOT / "BASELINE_REPORT.json").write_text(json.dumps(report, indent=2) + "\n")
    print("baseline_report", json.dumps(report), flush=True)


def run_thermal(stage, steps, seed, sink_z_max=None):
    if stage not in ("heating", "cooling"):
        raise ValueError(stage)
    cp = ROOT / ("BASELINE_CHECKPOINT.extxyz" if stage == "heating" and not (ROOT / "HEATING_CHECKPOINT.extxyz").exists() else f"{stage.upper()}_CHECKPOINT.extxyz")
    if not cp.exists():
        raise FileNotFoundError(f"Missing start checkpoint {cp.name}; initialize 300 K first")
    atoms = read(cp)
    attach(atoms)
    if sink_z_max is not None:
        phase_now = atoms.arrays["phase"].astype(str)
        support_now = np.asarray(atoms.arrays["support_fixed"], bool)
        new_sink = (phase_now == "Ti3SiC2") & ~support_now & (atoms.positions[:, 2] > 0.05) & (atoms.positions[:, 2] <= sink_z_max)
        if new_sink.sum() < 100:
            raise RuntimeError(f"Requested sink z <= {sink_z_max} A has only {new_sink.sum()} mobile atoms")
        atoms.arrays["cold_sink"] = new_sink.astype(np.int32)
        atoms.info["cold_sink_definition_run"] = f"mobile TSC atoms 0 < z <= {sink_z_max} A; bottom Ti support plane fixed"
    if stage == "heating":
        elapsed = float(atoms.info.get("heating_elapsed_ps", 0.0))
        global0 = float(atoms.info.get("time_ps", 0.0))
        heat = np.asarray(atoms.arrays["heat_region_initial"], bool)
        sink = np.asarray(atoms.arrays["cold_sink"], bool)
        temp = np.full((len(atoms), 1), 300.0)
        fric = np.zeros((len(atoms), 1))
        temp[heat] = 1900.0
        fric[heat | sink] = 0.02 / units.fs
        heat_temp = 1900.0
    else:
        elapsed = float(atoms.info.get("cooling_elapsed_ps", 0.0))
        global0 = float(atoms.info.get("time_ps", 0.0))
        heat = np.asarray(atoms.arrays["heat_region_initial"], bool)
        sink = np.asarray(atoms.arrays["cold_sink"], bool)
        temp = np.full((len(atoms), 1), 300.0)
        fric = np.zeros((len(atoms), 1))
        fric[sink] = 0.02 / units.fs
        heat_temp = 300.0
    if heat.sum() < 50 or sink.sum() < 100:
        raise RuntimeError(f"invalid reservoirs: heat={heat.sum()} sink={sink.sum()}")
    dyn = Langevin(atoms, timestep=DT_FS * units.fs, temperature_K=temp, friction=fric,
                   fixcm=False, rng=np.random.default_rng(seed))
    reference = read(ROOT / "BASELINE_300K.extxyz", index=0)
    path = ROOT / f"{stage.upper()}.extxyz"
    metricpath = ROOT / f"{stage.upper()}_METRICS.csv"
    t0 = time.perf_counter()
    ncompleted = 0
    melt_run = int(atoms.info.get("heating_melt_consecutive", 0))
    for local_step in range(SAMPLE_STEPS, steps + 1, SAMPLE_STEPS):
        dyn.run(SAMPLE_STEPS)
        elapsed += SAMPLE_STEPS * DT_FS / 1000.0
        row = metrics(atoms, reference, stage, int(round(elapsed * 1000.0 / DT_FS)), elapsed, global0, t0)
        row["thermostat_target_hot_K"] = heat_temp
        row["hot_source_atom_count"] = int(heat.sum())
        row["cold_sink_atom_count"] = int(sink.sum())
        append_row(metricpath, row)
        atoms.info["time_ps"] = global0 + elapsed
        if stage == "heating":
            atoms.info["heating_elapsed_ps"] = elapsed
            hot_disorder = (row["T_hot_Ag_K"] >= 1235.0 and row["hot_Ag_RMSD_from_baseline_A"] >= 1.5 and
                            row["hot_Ag_median_q6"] < 0.70 * float(atoms.info.get("baseline_hot_q6_median", 1.0)) and
                            row["hot_Ag_low_q6_fraction"] >= 0.5 and row["hot_low_q6_largest_cluster_atoms"] >= 20)
            if hot_disorder:
                melt_run += 1
            else:
                melt_run = 0
            atoms.info["heating_melt_consecutive"] = melt_run
            if melt_run >= 3:
                atoms.info["melting_established"] = True
        else:
            atoms.info["cooling_elapsed_ps"] = elapsed
        snapshot(path, atoms)
        checkpoint(ROOT / f"{stage.upper()}_CHECKPOINT.extxyz", atoms)
        ncompleted += SAMPLE_STEPS
        stop = crash_guard(row)
        print(stage, row["step"], "elapsed_ps", elapsed, "Thot", row["T_hot_Ag_K"], "q6", row["hot_Ag_median_q6"],
              "lowq6", row["hot_Ag_low_q6_fraction"], "cluster", row["hot_low_q6_largest_cluster_atoms"],
              "RMSD", row["hot_Ag_RMSD_from_baseline_A"], "Fmax", row["Fmax_mobile_eV_A"],
              "TSCconn", row["TSC_largest_component_fraction_cut3p25"], flush=True)
        if stop:
            atoms.info["early_stop_reason"] = stop
            checkpoint(ROOT / f"{stage.upper()}_CHECKPOINT.extxyz", atoms)
            print("EARLY_STOP", stage, stop, flush=True)
            break
        if stage == "heating" and melt_run >= 3:
            print("MELTING_ESTABLISHED candidate criteria met for 3 consecutive samples; heating segment stopped", flush=True)
            break
    report_path = ROOT / f"{stage.upper()}_REPORT.json"
    summary = {"stage": stage, "model_atoms": len(atoms), "dt_fs": DT_FS, "completed_steps_this_call": ncompleted,
               "stage_elapsed_ps_total": elapsed, "wall_seconds_this_call": time.perf_counter() - t0,
               "peak_RSS_GiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024,
               "hot_source_atom_count": int(heat.sum()), "cold_sink_atom_count": int(sink.sum()),
               "thermostat": "Langevin 1900 K on fixed-label top-two Ag planes locally; 300 K TSC cold sink; remaining atoms NVE" if stage == "heating" else "Ag heat source off; only mobile lower TSC sink Langevin 300 K; remaining atoms NVE",
               "status": str(atoms.info.get("early_stop_reason", "completed_requested_steps")),
               "melting_established_from_three_consecutive_multimetric_samples": bool(atoms.info.get("melting_established", False))}
    report_path.write_text(json.dumps(summary, indent=2) + "\n")
    print("thermal_report", json.dumps(summary), flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("stage", choices=["initialize", "heating", "cooling"])
    p.add_argument("--steps", type=int, default=200)
    p.add_argument("--seed", type=int, default=20261022)
    p.add_argument("--sink-z-max", type=float, default=None, help="Override the TSC Langevin sink depth, in A, for this continued thermal segment")
    a = p.parse_args()
    if a.stage == "initialize":
        initialize(a.steps)
    else:
        run_thermal(a.stage, a.steps, a.seed, a.sink_z_max)


if __name__ == "__main__":
    main()
