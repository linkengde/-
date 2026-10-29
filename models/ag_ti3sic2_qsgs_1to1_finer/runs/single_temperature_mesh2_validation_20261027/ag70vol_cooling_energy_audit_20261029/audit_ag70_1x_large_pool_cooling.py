#!/usr/bin/env python3
"""Audit energy, phase fraction, and cluster evolution in the Ag70 1x cooling run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np


AUDIT_DIR = Path(__file__).resolve().parent
RUN_DIR = AUDIT_DIR.parent
PROJECT_DIR = RUN_DIR.parents[1]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
HEAT_OFF_CHECKPOINT = RUN_DIR / "ag70vol_1x_continue_to_2300K_20261028/checkpoints/caseA_ag70_1x_final_18p75ps.npz"
COOLING_DIR = RUN_DIR / "ag70vol_cooling_from18p75_20261029"
EXTENSION_DIR = RUN_DIR / "ag70vol_cooling_extend_from19p25_20261029"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid  # noqa: E402
from qsgs_liquid_tracking import LiquidClusterTracker  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_metadata(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as archive:
        return json.loads(str(archive["metadata_json"].item()))


def latest_checkpoint() -> Path:
    candidates = list((EXTENSION_DIR / "checkpoints").glob("caseA_ag70_1x_cool_*.npz"))
    if not candidates:
        raise FileNotFoundError("No completed cooling checkpoints found in extension directory")
    valid: list[tuple[float, Path]] = []
    for candidate in candidates:
        try:
            valid.append((float(load_metadata(candidate)["time_ps"]), candidate))
        except Exception:
            continue
    if not valid:
        raise ValueError("No readable complete cooling checkpoint found")
    return max(valid, key=lambda item: item[0])[1]


def checkpoint_metrics(model: QSGSEnthalpyGrid, path: Path, label: str) -> dict[str, float | str | int]:
    tracker = LiquidClusterTracker(model.spacing_A, model.size_A, liquid_threshold=0.5)
    restored = model.load_checkpoint(path, tracker)
    with np.load(path, allow_pickle=False) as archive:
        enthalpy = np.asarray(archive["enthalpy_density_J_m3"], dtype=np.float64)
        temperature = np.asarray(archive["temperature_K"], dtype=np.float64)
        liquid = np.asarray(archive["ag_liquid_fraction"], dtype=np.float64)
        track_labels = np.asarray(archive["tracker_track_labels"], dtype=np.int32)
        metadata = json.loads(str(archive["metadata_json"].item()))

    if not np.all(np.isfinite(enthalpy)) or not np.all(np.isfinite(temperature)) or not np.all(np.isfinite(liquid)):
        raise FloatingPointError(f"Non-finite data in {path}")
    if np.any((liquid < 0.0) | (liquid > 1.0)):
        raise ValueError(f"Liquid fraction outside [0,1] in {path}")
    ag = model.ag_fraction
    ag_volume_weight = float(ag.sum())
    total_energy = float(enthalpy.sum() * model.cell_volume_m3)
    injected_energy = float(metadata["injected_energy_J"])
    energy_error = total_energy - injected_energy
    relative_error = energy_error / injected_energy if injected_energy else 0.0
    labels = track_labels.ravel()
    positive = labels[labels > 0]
    largest_voxels = int(np.bincount(positive).max()) if positive.size else 0
    largest_cluster_nm3 = largest_voxels * model.cell_volume_m3 * 1.0e27
    return {
        "checkpoint_label": label,
        "time_ps": float(metadata["time_ps"]),
        "checkpoint": str(path.relative_to(RUN_DIR)),
        "checkpoint_sha256": sha256(path),
        "total_enthalpy_J": total_energy,
        "cumulative_source_energy_J": injected_energy,
        "energy_error_J": energy_error,
        "relative_energy_error": relative_error,
        "peak_temperature_K": float(temperature.max()),
        "mean_ag_temperature_K": float(np.sum(temperature * ag) / ag_volume_weight),
        "ag_liquid_volume_fraction": float(np.sum(liquid * ag) / ag_volume_weight),
        "largest_cluster_voxels": largest_voxels,
        "largest_cluster_volume_nm3": largest_cluster_nm3,
        "restored_energy_error_J": float(restored["energy_balance_error_J"]),
    }


def collect_trace() -> list[dict[str, float]]:
    files = [COOLING_DIR / "cooling.csv"]
    files.extend(sorted((EXTENSION_DIR / "segments").glob("cooling_*.csv")))
    by_time: dict[float, dict[str, float]] = {}
    for path in files:
        if not path.is_file():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                t = round(float(row["time_ps"]), 9)
                by_time[t] = {
                    "time_ps": t,
                    "peak_temperature_K": float(row["max_temperature_K"]),
                    "mean_ag_temperature_K": float(row["mean_ag_temperature_K"]),
                    "ag_liquid_volume_fraction": float(row["ag_liquid_volume_fraction"]),
                    "largest_cluster_volume_nm3": float(row["largest_liquid_cluster_volume_A3"]) / 1000.0,
                    "total_enthalpy_J": float(row["total_energy_J"]),
                    "cumulative_source_energy_J": float(row["injected_energy_J"]),
                    "energy_error_J": float(row["energy_balance_error_J"]),
                }
    trace = [by_time[t] for t in sorted(by_time)]
    if not trace:
        raise FileNotFoundError("No cooling time-series samples found")
    return trace


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_trace(path: Path, trace: list[dict[str, float]], tm_K: float, f_eq: float | None) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-ag70-cooling-audit")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp/cache-ag70-cooling-audit")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    times = np.asarray([row["time_ps"] for row in trace])
    cooling_elapsed = times - 18.75
    peak = np.asarray([row["peak_temperature_K"] for row in trace])
    mean_ag = np.asarray([row["mean_ag_temperature_K"] for row in trace])
    liquid = np.asarray([row["ag_liquid_volume_fraction"] for row in trace]) * 100.0
    cluster = np.asarray([row["largest_cluster_volume_nm3"] for row in trace])

    fig, axes = plt.subplots(3, 1, figsize=(10.5, 9.0), sharex=True)
    axes[0].plot(cooling_elapsed, peak, label="Peak temperature", color="#c44536")
    axes[0].plot(cooling_elapsed, mean_ag, label="Ag-weighted mean temperature", color="#286f9b")
    axes[0].axhline(tm_K, color="#444444", linestyle="--", linewidth=1.0, label=f"Ag Tm = {tm_K:.2f} K")
    axes[0].set_ylabel("Temperature (K)")
    axes[0].legend(loc="best", frameon=True)
    axes[0].grid(alpha=0.25)

    axes[1].plot(cooling_elapsed, liquid, color="#e28e00")
    if f_eq is not None:
        axes[1].axhline(f_eq * 100.0, color="#555555", linestyle="--", linewidth=1.0,
                        label=f"Energy-based uniform-Tm estimate: {f_eq * 100.0:.2f}%")
    axes[1].set_ylabel("Ag liquid fraction (%)")
    axes[1].legend(loc="best", frameon=True)
    axes[1].grid(alpha=0.25)

    axes[2].plot(cooling_elapsed, cluster, color="#6b4c9a")
    axes[2].set_ylabel("Largest tracked cluster (nm³)")
    axes[2].set_xlabel("Time after heat-off (ps)")
    axes[2].grid(alpha=0.25)

    fig.suptitle("Ag70 1× Case A source-off evolution (insulated, zero-storage scaffold)")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--latest-checkpoint", type=Path, default=None,
                        help="Latest complete checkpoint; defaults to the newest readable extension checkpoint")
    args = parser.parse_args()
    latest = (args.latest_checkpoint or latest_checkpoint()).resolve()
    required = [
        ("heat_off_18p75ps", HEAT_OFF_CHECKPOINT),
        ("cool_19p25ps", COOLING_DIR / "checkpoints/caseA_ag70_1x_cool_final_19p25ps.npz"),
        ("cool_19p27ps", EXTENSION_DIR / "checkpoints/caseA_ag70_1x_cool_19p27ps.npz"),
        ("latest_complete", latest),
    ]
    for _, path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    model = QSGSEnthalpyGrid(None, "A", native_geometry_npz=GEOMETRY,
                             source_geometry_npz=GEOMETRY, mesh_refinement=1)
    geometry_ag_fraction = float(model.ag_fraction.sum() / model.ag_fraction.size)
    ag_volume_m3 = float(model.ag_fraction.sum() * model.cell_volume_m3)
    ag_mass_kg = ag_volume_m3 * model.ag.reference_density_kg_m3
    e_s_J = float(model._solid_at_tm_J_m3.sum() * model.cell_volume_m3)
    latent_capacity_J = float(model._latent_volumetric_J_m3.sum() * model.cell_volume_m3)
    e_l_J = e_s_J + latent_capacity_J
    sample_metrics = [checkpoint_metrics(model, path, label) for label, path in required]
    e_off_J = float(sample_metrics[0]["total_enthalpy_J"])
    if e_off_J < e_s_J:
        eq_case = "E_off < E_s: all-solid equilibrium is energetically possible"
        f_eq = None
    elif e_off_J <= e_l_J:
        f_eq = (e_off_J - e_s_J) / (e_l_J - e_s_J)
        eq_case = "E_s <= E_off <= E_l: residual liquid remains at uniform-Tm equilibrium"
    else:
        f_eq = None
        eq_case = "E_off > E_l: total energy exceeds the all-liquid-at-Tm enthalpy"

    trace = collect_trace()
    injected_values = np.asarray([row["cumulative_source_energy_J"] for row in trace])
    source_energy_constant = bool(np.allclose(injected_values, injected_values[0], rtol=0.0, atol=1e-28))
    out_metrics = AUDIT_DIR / "cooling_energy_checkpoint_metrics.csv"
    out_trace = AUDIT_DIR / "cooling_trace_18p75ps_to_latest.csv"
    out_plot = AUDIT_DIR / "cooling_trace_18p75ps_to_latest.png"
    out_json = AUDIT_DIR / "cooling_energy_audit_results.json"
    write_csv(out_metrics, sample_metrics)
    write_csv(out_trace, trace)
    plot_trace(out_plot, trace, model.melting_temperature_K, f_eq)

    result = {
        "model": "Ag70vol Case A 1x, single-temperature enthalpy, zero scaffold storage/conductivity, insulated external boundaries",
        "geometry_file": str(GEOMETRY.relative_to(RUN_DIR)),
        "geometry_sha256": sha256(GEOMETRY),
        "thermophysical_table_sha256": sha256(PHASE_DIR / "thermophysical_candidates.csv"),
        "thermophysical_config_sha256": sha256(PHASE_DIR / "thermophysical_candidates.json"),
        "latest_checkpoint": str(latest.relative_to(RUN_DIR)),
        "latest_checkpoint_sha256": sha256(latest),
        "latest_checkpoint_restore_validated": True,
        "time_heat_off_ps": 18.75,
        "latest_cooling_time_ps": float(sample_metrics[-1]["time_ps"]),
        "source_off_code_path": "cool_ag70_from_checkpoint.py invokes model.run(source_on=False)",
        "cumulative_source_energy_constant_during_cooling": source_energy_constant,
        "ag_volume_fraction_from_actual_grid": geometry_ag_fraction,
        "ag_volume_m3": ag_volume_m3,
        "ag_mass_kg_using_candidate_density": ag_mass_kg,
        "candidate_ag_density_kg_m3": model.ag.reference_density_kg_m3,
        "candidate_ag_melting_temperature_K": model.melting_temperature_K,
        "candidate_ag_latent_heat_J_kg": model.ag.latent_heat_J_kg,
        "candidate_ag_solid_cp_anchor_temperature_K": model.ag.solid_temperature_K.tolist(),
        "candidate_ag_solid_cp_anchor_J_kg_K": model.ag.solid_cp_J_kg_K.tolist(),
        "candidate_ag_liquid_cp_anchor_J_kg_K": model.ag.liquid_cp_J_kg_K.tolist(),
        "solid_specific_enthalpy_at_Tm_J_kg": model.ag.solid_h_at_melt_J_kg,
        "E_off_J": e_off_J,
        "E_s_all_ag_solid_at_Tm_J": e_s_J,
        "E_l_all_ag_liquid_at_Tm_J": e_l_J,
        "latent_capacity_J": latent_capacity_J,
        "E_off_minus_E_s_J": e_off_J - e_s_J,
        "E_l_minus_E_off_J": e_l_J - e_off_J,
        "energy_case": eq_case,
        "f_eq_if_uniform_at_Tm": f_eq,
        "source_off_energy_residual_latest_J": sample_metrics[-1]["energy_error_J"],
        "source_off_relative_energy_error_latest": sample_metrics[-1]["relative_energy_error"],
        "latest_sample_peak_temperature_K": sample_metrics[-1]["peak_temperature_K"],
        "latest_sample_mean_ag_temperature_K": sample_metrics[-1]["mean_ag_temperature_K"],
        "latest_sample_ag_liquid_volume_fraction": sample_metrics[-1]["ag_liquid_volume_fraction"],
        "latest_sample_largest_cluster_volume_nm3": sample_metrics[-1]["largest_cluster_volume_nm3"],
        "metric_checkpoints": sample_metrics,
        "trace_rows": len(trace),
        "trace_source_files": [str(COOLING_DIR / "cooling.csv"), str(EXTENSION_DIR / "segments")],
        "limitations": [
            "Thermophysical inputs remain candidate values and are not fully source-verified.",
            "Uniform-Tm f_eq is an energy-balance estimate, not proof that the finite-time transient reached equilibrium.",
            "No Ti3SiC2 heat storage/conduction or external heat rejection is present in Case A.",
            "Vaporization and hydrodynamic melt flow are not modeled.",
        ],
    }
    out_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
