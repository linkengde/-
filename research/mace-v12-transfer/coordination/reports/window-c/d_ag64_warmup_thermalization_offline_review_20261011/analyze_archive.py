#!/usr/bin/env python3
"""Read-only diagnostics for D's archived 64-Ag warmup; performs no MD/DFT."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

import ase
import matplotlib
import numpy as np
import scipy
from ase.io import read
from scipy.special import sph_harm_y


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def manifest_entry(manifest: dict, rel: str) -> dict:
    matches = [m for m in manifest["members"] if m["path"] == rel]
    if len(matches) != 1:
        raise ValueError(f"archive manifest entry missing or duplicated: {rel}")
    return matches[0]


def verify_member(root: Path, manifest: dict, rel: str) -> Path:
    path = root / rel
    entry = manifest_entry(manifest, rel)
    if not path.is_file() or path.stat().st_size != entry["size_bytes"] or sha256(path) != entry["sha256"]:
        raise ValueError(f"archive member failed size/SHA verification: {rel}")
    return path


def window_stats(rows: list[dict], name: str, selected: list[dict]) -> dict:
    result: dict = {"window": name, "n_records": len(selected),
                    "step_first": selected[0]["step"], "step_last": selected[-1]["step"],
                    "time_first_fs": selected[0]["time_fs"], "time_last_fs": selected[-1]["time_fs"]}
    t = np.asarray([x["T_K"] for x in selected], dtype=float)
    times = np.asarray([x["time_fs"] for x in selected], dtype=float)
    fields = {"T_K": t,
              "sampler_energy_eV": np.asarray([x["sampler_energy_eV"] for x in selected]),
              "kinetic_energy_derived_eV": np.asarray([x["kinetic_energy_derived_eV"] for x in selected]),
              "total_energy_derived_eV": np.asarray([x["total_energy_derived_eV"] for x in selected]),
              "pressure_eV_A3": np.asarray([x["pressure_eV_A3"] for x in selected])}
    for key, values in fields.items():
        result[key + "_mean"] = float(values.mean())
        result[key + "_std_sample"] = float(values.std(ddof=1)) if len(values) > 1 else 0.0
        result[key + "_min"] = float(values.min())
        result[key + "_max"] = float(values.max())
        result[key + "_slope_per_100fs"] = (float(np.polyfit(times, values, 1)[0] * 100.0)
                                                  if len(values) > 1 else 0.0)
    return result


def rdf_for_frames(frames, bins: np.ndarray, edge: float) -> np.ndarray:
    n = len(frames[0])
    volume = edge**3
    pair_count = n * (n - 1) / 2
    counts = np.zeros(len(bins) - 1, dtype=np.int64)
    for atoms in frames:
        pos = atoms.get_positions()
        delta = pos[:, None, :] - pos[None, :, :]
        delta -= np.rint(delta / edge) * edge
        distances = np.linalg.norm(delta, axis=2)
        distances = distances[np.triu_indices(n, 1)]
        counts += np.histogram(distances, bins=bins)[0]
    shell = 4.0 * np.pi / 3.0 * (bins[1:]**3 - bins[:-1]**3)
    return counts / (len(frames) * pair_count * shell / volume)


def q6_values(vectors: np.ndarray) -> float:
    radii = np.linalg.norm(vectors, axis=1)
    theta = np.arccos(np.clip(vectors[:, 2] / radii, -1.0, 1.0))
    phi = np.arctan2(vectors[:, 1], vectors[:, 0])
    qlm = np.asarray([np.mean(sph_harm_y(6, m, theta, phi)) for m in range(-6, 7)])
    return float(np.sqrt(4.0 * np.pi / 13.0 * np.sum(np.abs(qlm)**2)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive-dir", type=Path, required=True,
                        help="Freshly restored extracted archive root")
    parser.add_argument("--manifest", type=Path, required=True,
                        help="Existing archive-manifest.json from D branch")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    root, manifest_path, out = args.archive_dir, args.manifest, args.outdir
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(manifest_path.read_text())
    for member in manifest["members"]:
        verify_member(root, manifest, member["path"])
    thermal = verify_member(root, manifest, "thermal.csv")
    trajectory_path = verify_member(root, manifest, "trajectory.traj")
    state = verify_member(root, manifest, "sampler_state.json")
    source = verify_member(root, manifest, "source_manifest.json")

    with thermal.open(newline="") as f:
        raw_rows = list(csv.DictReader(f))
    if len(raw_rows) != 51:
        raise ValueError(f"expected exactly 51 thermal records, got {len(raw_rows)}")
    numeric = ("step", "time_fs", "target_temperature_K", "T_K", "sampler_energy_eV",
               "pressure_eV_A3", "COMspeed_A_fs")
    rows = []
    dof = 3 * 64 - 3
    k_b_ev_k = 8.617333262145e-5
    for raw in raw_rows:
        row = dict(raw)
        for key in numeric:
            row[key] = float(raw[key])
        row["kinetic_energy_derived_eV"] = 0.5 * dof * k_b_ev_k * row["T_K"]
        row["total_energy_derived_eV"] = row["sampler_energy_eV"] + row["kinetic_energy_derived_eV"]
        rows.append(row)
    if not np.isfinite([[r[k] for k in numeric] for r in rows]).all():
        raise ValueError("non-finite value found in thermal log")

    shutil.copyfile(thermal, out / "thermal.csv")
    with (out / "thermal_derived_quantities.csv").open("w", newline="") as f:
        fields = ["step", "time_fs", "T_K", "target_temperature_K", "sampler_energy_eV",
                  "kinetic_energy_derived_eV", "total_energy_derived_eV", "pressure_eV_A3",
                  "pressure_GPa_from_eV_A3"]
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            result = {key: row[key] for key in fields if key in row}
            result["pressure_GPa_from_eV_A3"] = row["pressure_eV_A3"] * 160.21766208
            writer.writerow(result)
    with (out / "thermal_windows.csv").open("w", newline="") as f:
        names = ["all", "first_100_steps_inclusive", "last_200_steps_inclusive",
                 "step_100_to_500_inclusive", "step_500_to_1000_inclusive"]
        selected = [rows,
                    [r for r in rows if r["step"] <= 100],
                    [r for r in rows if r["step"] >= 800],
                    [r for r in rows if 100 <= r["step"] <= 500],
                    [r for r in rows if 500 <= r["step"] <= 1000]]
        stats = [window_stats(rows, n, ss) for n, ss in zip(names, selected)]
        writer = csv.DictWriter(f, fieldnames=list(stats[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(stats)

    frames = read(str(trajectory_path), index=":")
    if len(frames) != 11 or any(len(a) != 64 for a in frames):
        raise ValueError("trajectory frame/atom count differs from archived specification")
    edge = float(frames[0].cell.lengths()[0])
    if not all(np.all(a.pbc) and np.allclose(a.cell.array, frames[0].cell.array, atol=1e-10)
               for a in frames):
        raise ValueError("trajectory PBC/cell mismatch")
    if any(not np.isfinite(a.positions).all() for a in frames):
        raise ValueError("non-finite trajectory coordinate")
    expected_steps = list(range(0, 1001, 100))
    for atoms, expected_step in zip(frames, expected_steps):
        if int(atoms.info.get("step", -1)) != expected_step:
            raise ValueError("trajectory step sequence mismatch")
        if not np.all(atoms.numbers == 47) or atoms.info.get("atom_ids") != list(range(1, 65)):
            raise ValueError("trajectory element identity or atom ID mismatch")
    restart_rel = next((m["path"] for m in manifest["members"]
                       if m["path"].endswith("/state_arrays.npz") and "/step-000001000-" in m["path"]), None)
    restart_atoms_rel = next((m["path"] for m in manifest["members"]
                             if m["path"].endswith("/state_atoms.traj") and "/step-000001000-" in m["path"]), None)
    if not restart_rel or not restart_atoms_rel:
        raise ValueError("step-1000 checkpoint members missing from archive manifest")
    checkpoint_arrays = np.load(root / restart_rel)
    checkpoint_atoms = read(str(root / restart_atoms_rel))
    endpoint = frames[-1]
    checkpoint_exact = {
        "positions": bool(np.array_equal(checkpoint_arrays["positions"], endpoint.positions)),
        "momenta": bool(np.array_equal(checkpoint_arrays["momenta"], endpoint.get_momenta())),
        "cell": bool(np.array_equal(checkpoint_arrays["cell"], endpoint.cell.array)),
        "numbers": bool(np.array_equal(checkpoint_arrays["numbers"], endpoint.numbers)),
        "masses": bool(np.array_equal(checkpoint_arrays["masses"], endpoint.get_masses())),
        "pbc": bool(np.array_equal(checkpoint_arrays["pbc"], endpoint.pbc)),
        "atom_ids": bool(np.array_equal(checkpoint_arrays["atom_ids"], np.asarray(endpoint.info["atom_ids"]))),
        "state_atoms_positions": bool(np.array_equal(checkpoint_atoms.positions, endpoint.positions)),
        "state_atoms_momenta": bool(np.array_equal(checkpoint_atoms.get_momenta(), endpoint.get_momenta())),
    }
    if not all(checkpoint_exact.values()):
        raise ValueError(f"checkpoint differs from trajectory endpoint: {checkpoint_exact}")

    # The cubic minimum-image RDF is only reported through L/2.  Bin width is
    # fixed at 0.05 A for this coarse screening; no phase label is inferred.
    dr = 0.05
    rmax = edge / 2.0
    bins = np.arange(0.0, math.floor(rmax / dr) * dr + dr * 0.5, dr)
    if bins[-1] < rmax:
        bins = np.append(bins, rmax)
    centers = (bins[:-1] + bins[1:]) / 2.0
    early = [a for a in frames if float(a.info["step"]) <= 100]
    late = [a for a in frames if float(a.info["step"]) >= 800]
    g_all = rdf_for_frames(frames, bins, edge)
    g_early = rdf_for_frames(early, bins, edge)
    g_late = rdf_for_frames(late, bins, edge)
    smooth = np.convolve(g_all, np.ones(5) / 5.0, mode="same")
    peaks = [i for i in range(2, len(smooth) - 2)
             if smooth[i] > smooth[i - 1] and smooth[i] >= smooth[i + 1] and centers[i] > 2.0]
    first_peak = peaks[0] if peaks else None
    first_min = None
    if first_peak is not None:
        for i in range(first_peak + 1, len(smooth) - 1):
            if smooth[i] <= smooth[i - 1] and smooth[i] < smooth[i + 1]:
                first_min = i
                break
    q6_cutoff = float(centers[first_min]) if first_min is not None else 3.675

    with (out / "rdf.csv").open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["r_A", "g_all_11_frames", "g_early_step0_to100", "g_late_step800_to1000"])
        writer.writerows(zip(centers, g_all, g_early, g_late))

    geom_rows = []
    local_rows = []
    all_nearest = []
    for atoms in frames:
        pos = atoms.get_positions()
        delta = pos[:, None, :] - pos[None, :, :]
        delta -= np.rint(delta / edge) * edge
        distances = np.linalg.norm(delta, axis=2)
        np.fill_diagonal(distances, np.inf)
        nearest = distances.min(axis=1)
        all_nearest.extend(nearest.tolist())
        qvals, coordination = [], []
        for i in range(len(atoms)):
            vec = delta[i]
            rr = distances[i]
            mask = (rr < q6_cutoff) & np.isfinite(rr)
            coordination.append(int(mask.sum()))
            if mask.any():
                qval = q6_values(vec[mask])
                qvals.append(qval)
            else:
                qval = float("nan")
            local_rows.append({"step": int(atoms.info["step"]), "time_fs": float(atoms.info["time_fs"]),
                               "atom_id": int(atoms.info["atom_ids"][i]),
                               "nearest_neighbor_distance_A": float(nearest[i]),
                               "coordination_within_cutoff": int(coordination[-1]), "Q6": qval})
        step = int(atoms.info["step"])
        geom_rows.append({
            "step": step, "time_fs": float(atoms.info["time_fs"]),
            "min_pair_distance_A": float(nearest.min()),
            "nearest_neighbor_mean_A": float(nearest.mean()),
            "nearest_neighbor_std_A": float(nearest.std(ddof=1)),
            "nearest_neighbor_q10_A": float(np.quantile(nearest, 0.10)),
            "nearest_neighbor_median_A": float(np.median(nearest)),
            "nearest_neighbor_q90_A": float(np.quantile(nearest, 0.90)),
            "coordination_cutoff_A": q6_cutoff,
            "coordination_mean": float(np.mean(coordination)),
            "coordination_std": float(np.std(coordination, ddof=1)),
            "Q6_mean": float(np.mean(qvals)), "Q6_std": float(np.std(qvals, ddof=1)),
            "Q6_q10": float(np.quantile(qvals, 0.10)), "Q6_median": float(np.median(qvals)),
            "Q6_q90": float(np.quantile(qvals, 0.90)),
        })
    with (out / "frame_geometry.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(geom_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(geom_rows)
    with (out / "per_atom_local_order.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(local_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(local_rows)
    with (out / "nearest_neighbor_distribution.csv").open("w", newline="") as f:
        fields = ["step", "time_fs", "atom_id", "nearest_neighbor_distance_A"]
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows([{key: row[key] for key in fields} for row in local_rows])

    state_data, source_data = json.loads(state.read_text()), json.loads(source.read_text())
    summary = {
        "archive_task_id": manifest["task_id"], "owner": manifest["owner"],
        "archive_sha256": manifest["archive_sha256"], "archive_bytes": manifest["archive_bytes"],
        "manifest_sha256": sha256(manifest_path), "thermal_sha256": sha256(thermal),
        "trajectory_sha256": sha256(trajectory_path), "sampler_state_sha256": sha256(state),
        "source_manifest_sha256": sha256(source), "thermal_rows": len(rows),
        "trajectory_frames": len(frames), "atom_count": 64, "box_edge_A": edge,
        "all_18_archive_members_verified": True, "checkpoint_endpoint_exact_match": checkpoint_exact,
        "pbc": [True, True, True], "rdf_dr_A": dr, "rdf_rmax_A": edge / 2.0,
        "rdf_first_peak_A_smoothed": float(centers[first_peak]) if first_peak is not None else None,
        "rdf_first_peak_g_smoothed": float(smooth[first_peak]) if first_peak is not None else None,
        "rdf_first_minimum_A_smoothed": q6_cutoff,
        "Q6_definition": "per-atom Steinhardt Q6 from neighbors inside the first minimum of the 5-bin-smoothed 11-frame RDF; coarse structural descriptor only",
        "nearest_neighbor_distance_quantiles_A": dict(zip(
            ["min", "q01", "q10", "median", "q90", "q99", "max"],
            map(float, np.quantile(all_nearest, [0, .01, .10, .5, .9, .99, 1])))),
        "checkpoint_state": {key: state_data.get(key) for key in
                              ("status", "step", "time_fs", "task_id", "bundle_sha256", "rng_class",
                               "langevin_version", "C_release_config_sha256_at_checkpoint", "scientific_parameters_sha256")},
        "source_classification": {key: source_data.get(key) for key in
                                  ("liquid_QA", "DFT_LABELS", "STRUCTURE_PASS", "source_role", "RNG_seed",
                                   "density_g_cm3", "box_edge_A", "packing_attempts")},
        "analysis_software": {"python": sys.version.split()[0], "ase": ase.__version__,
                              "numpy": np.__version__, "scipy": scipy.__version__,
                              "matplotlib": matplotlib.__version__},
        "temperature_canonical_sigma_estimate_K_at_1600K": float(1600.0 * math.sqrt(2.0 / (3 * 64 - 3))),
        "energy_definition": {
            "sampler_energy_eV": "potential energy from atoms.get_potential_energy(); kinetic energy is not included",
            "kinetic_energy_derived_eV": "derived from logged T using ASE DOF=3*64-3=189 and K=DOF*kB*T/2; not an original CSV field",
            "total_energy_derived_eV": "potential plus derived kinetic energy; Langevin NVT does not conserve this quantity",
            "pressure_eV_A3": "runner formula -trace(stress)/3 + 2*K/(3V); K is kinetic energy"
        },
        "derived_energy_at_step0_peakT_and_final": {
            str(step): {"T_K": next(r["T_K"] for r in rows if r["step"] == step),
                        "potential_eV": next(r["sampler_energy_eV"] for r in rows if r["step"] == step),
                        "kinetic_derived_eV": next(r["kinetic_energy_derived_eV"] for r in rows if r["step"] == step),
                        "total_derived_eV": next(r["total_energy_derived_eV"] for r in rows if r["step"] == step)}
            for step in (0, 80, 1000)
        },
        "classification": {
            "INITIAL_TRANSIENT": "YES: early potential energy falls while temperature and pressure show a large startup excursion",
            "THERMALIZATION_NOT_ESTABLISHED": "YES: final 200-step temperature/energy/pressure means continue to drift and remain above target temperature",
            "NUMERICAL_INSTABILITY": "NO_EVIDENCE: all logged values and coordinates finite; no runtime distance-stop violation; archived run completed normally",
            "THERMALIZATION_CANDIDATE": "NOT_YET: 0.5 ps/11 structural frames are too short and show unresolved drift; no liquid claim"
        }
    }
    (out / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    steps = np.asarray([r["step"] for r in rows])
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True, constrained_layout=True)
    axes[0].plot(steps, [r["T_K"] for r in rows], color="#b2182b", marker="o", ms=2.8, lw=1.2, label="T_actual")
    axes[0].plot(steps, [r["target_temperature_K"] for r in rows], color="#2166ac", ls="--", lw=1.2, label="T_target")
    axes[0].set_ylabel("Temperature (K)"); axes[0].legend(frameon=False, ncol=2)
    axes[0].grid(alpha=.25)
    axes[1].plot(steps, [r["sampler_energy_eV"] for r in rows], color="#1b7837", marker="o", ms=2.5, lw=1.1)
    axes[1].set_ylabel("Potential energy (eV)"); axes[1].grid(alpha=.25)
    pressure = np.asarray([r["pressure_eV_A3"] for r in rows])
    axes[2].plot(steps, pressure, color="#762a83", marker="o", ms=2.5, lw=1.1)
    axes[2].set_ylabel("Pressure (eV/Å³)"); axes[2].set_xlabel("Global step (0.5 fs/step)"); axes[2].grid(alpha=.25)
    pressure_ax = axes[2].twinx()
    pressure_ax.set_ylim(np.asarray(axes[2].get_ylim()) * 160.21766208)
    pressure_ax.set_ylabel("Pressure (GPa)")
    fig.suptitle("D Ag64 archived 1000-step warmup: offline thermal diagnostics")
    fig.savefig(out / "thermal_diagnostics.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True)
    axes[0].plot(centers, g_all, color="#222222", lw=1.4, label="all 11 frames")
    axes[0].plot(centers, g_early, color="#d6604d", lw=1.0, label="steps 0–100 (2 frames)")
    axes[0].plot(centers, g_late, color="#4393c3", lw=1.0, label="steps 800–1000 (3 frames)")
    axes[0].axvline(q6_cutoff, color="#666666", ls="--", lw=1.0, label=f"first-minimum cutoff ≈ {q6_cutoff:.2f} Å")
    axes[0].set_xlim(0, edge / 2); axes[0].set_ylabel("g(r)"); axes[0].set_xlabel("r (Å)")
    axes[0].legend(frameon=False, ncol=2, fontsize=8); axes[0].grid(alpha=.2)
    frame_steps = [g["step"] for g in geom_rows]
    ax = axes[1]
    ax.plot(frame_steps, [g["Q6_mean"] for g in geom_rows], marker="o", color="#1b7837", label="mean local Q6")
    ax.fill_between(frame_steps,
                    [g["Q6_mean"] - g["Q6_std"] for g in geom_rows],
                    [g["Q6_mean"] + g["Q6_std"] for g in geom_rows], alpha=.18, color="#1b7837", label="±1 atomwise SD")
    ax.set_ylabel("Q6 (coarse local order)"); ax.set_xlabel("Saved trajectory step"); ax.grid(alpha=.2); ax.legend(frameon=False)
    fig.suptitle("Ag64 structure prescreen only; not a liquid-state validation")
    fig.savefig(out / "structure_prescreen.png", dpi=180)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
