#!/usr/bin/env python3
"""Continue the native-voxel QSGS Ag TTM heating run from a saved 50 fs state."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

MODEL_DIR = Path(__file__).resolve().parents[1]
PHASE_DIR = Path(__file__).resolve().parent
DIAG_DIR = PHASE_DIR / "nested_grid_diagnostics"
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))
if str(PHASE_DIR) not in sys.path:
    sys.path.insert(0, str(PHASE_DIR))

from qsgs_component_ttm_backend import (  # noqa: E402
    ANGSTROM3_TO_M3,
    EV_J,
    QSGSComponentTTMBackend,
)

DEFAULT_GEOMETRY = MODEL_DIR / "qsgs_morphology.npz"
DEFAULT_ELECTRON_TABLE = MODEL_DIR / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv"
DEFAULT_CHECKPOINT = DIAG_DIR / "qsgs_component_native_50fs_restart.npz"
DEFAULT_SCREEN_JSON = DIAG_DIR / "qsgs_component_full_domain_50fs_match.json"
DEFAULT_OUTPUT = DIAG_DIR / "qsgs_component_native_heating_to_0p1ps.json"
DEFAULT_LOG = DIAG_DIR / "qsgs_component_native_heating_to_0p1ps.log"


def _largest_liquid_cluster(backend: QSGSComponentTTMBackend, threshold: float = 0.5) -> dict:
    active = backend.lattice_liquid_fraction >= threshold
    count = int(np.count_nonzero(active))
    if count == 0:
        return {
            "threshold_liquid_fraction": threshold,
            "active_node_count": 0,
            "largest_cluster_node_count": 0,
            "largest_cluster_volume_A3": 0.0,
        }
    active_ids = np.flatnonzero(active)
    local_id = np.full(active.size, -1, dtype=np.int32)
    local_id[active_ids] = np.arange(count, dtype=np.int32)
    edge_i = backend.graph.edge_i
    edge_j = backend.graph.edge_j
    keep = active[edge_i] & active[edge_j]
    row = local_id[edge_i[keep]]
    col = local_id[edge_j[keep]]
    adjacency = coo_matrix(
        (np.ones(row.size * 2, dtype=np.uint8),
         (np.concatenate((row, col)), np.concatenate((col, row)))),
        shape=(count, count),
    ).tocsr()
    n_components, labels = connected_components(adjacency, directed=False, return_labels=True)
    sizes = np.bincount(labels, minlength=n_components)
    largest = int(sizes.max(initial=0))
    # Native graph nodes represent one source voxel each, but use node volumes
    # so this remains correct if the backend's native grouping changes.
    largest_labels = np.flatnonzero(sizes == largest)
    largest_mask = np.isin(labels, largest_labels)
    volume = float(backend.node_volume_A3[active_ids[largest_mask]].sum())
    return {
        "threshold_liquid_fraction": threshold,
        "active_node_count": count,
        "largest_cluster_node_count": largest,
        "largest_cluster_volume_A3": volume,
    }


def _restore_native_checkpoint(
    backend: QSGSComponentTTMBackend,
    checkpoint: Path,
    screen_json: Path,
    checkpoint_time_ps: float,
) -> dict:
    summary = json.loads(screen_json.read_text(encoding="utf-8"))
    run = summary["results_by_time_step"]["0.3125fs"]["native_voxel"]
    with np.load(checkpoint, allow_pickle=False) as state:
        ag_node_ids = backend.graph.component_id[backend.source_ag]
        n_nodes = backend.graph.voxel_count.size
        if ag_node_ids.size != n_nodes or np.unique(ag_node_ids).size != n_nodes:
            raise ValueError("Checkpoint restore requires one native graph node per Ag voxel")
        for name, destination in (
            ("electron_temperature_K", backend.electron_temperature_K),
            ("lattice_temperature_K", backend.lattice_temperature_K),
            ("ag_liquid_fraction", backend.lattice_liquid_fraction),
        ):
            values = np.asarray(state[name], dtype=float)
            if values.shape != ag_node_ids.shape:
                raise ValueError(f"Checkpoint field {name} has unexpected shape {values.shape}")
            destination[ag_node_ids] = values
        coupling_density = np.asarray(state["electron_to_lattice_energy_eV_A3"], dtype=float)
        backend.cumulative_coupling_eV[ag_node_ids] = (
            coupling_density * backend.node_volume_A3[ag_node_ids]
        )

    backend.time_ps = float(checkpoint_time_ps)
    backend.injected_energy_J = float(run["source_energy_J"])
    backend.coupling_energy_J = float(run["electron_to_lattice_cumulative_J"])
    backend.electron_diffusion_energy_J = float(run.get("electron_diffusion_energy_J", 0.0))
    backend.lattice_diffusion_energy_J = float(run.get("lattice_diffusion_energy_J", 0.0))
    checkpoint_total = backend.total_energy_J
    recorded_total = float(run["electron_energy_J"]) + float(run["lattice_enthalpy_J"])
    relative_restore_error = abs(checkpoint_total - recorded_total) / max(abs(recorded_total), 1.0e-30)
    if relative_restore_error > 1.0e-10:
        raise ValueError(f"Restored checkpoint energy differs from recorded state: {relative_restore_error:g}")
    residual = (
        checkpoint_total
        - backend.initial_total_energy_J
        - backend.injected_energy_J
        - backend.electron_diffusion_energy_J
        - backend.lattice_diffusion_energy_J
    )
    return {
        "checkpoint_time_ps": backend.time_ps,
        "restored_energy_relative_difference": relative_restore_error,
        "checkpoint_energy_residual_J": residual,
        "restored_source_energy_J": backend.injected_energy_J,
        "restored_cumulative_coupling_J": backend.coupling_energy_J,
        "max_electron_temperature_K": float(backend.electron_temperature_K.max()),
        "max_lattice_temperature_K": float(backend.lattice_temperature_K.max()),
        "ag_liquid_volume_fraction": float(
            np.dot(backend.node_volume_A3, backend.lattice_liquid_fraction)
            / backend.node_volume_A3.sum()
        ),
    }


def _sample(backend: QSGSComponentTTMBackend, label: str, dt_fs: float) -> dict:
    total = backend.total_energy_J
    residual = (
        total
        - backend.initial_total_energy_J
        - backend.injected_energy_J
        - backend.electron_diffusion_energy_J
        - backend.lattice_diffusion_energy_J
    )
    return {
        "label": label,
        "time_ps": backend.time_ps,
        "dt_fs": dt_fs,
        "electron_temperature_max_K": float(backend.electron_temperature_K.max()),
        "electron_temperature_mean_K": float(
            np.average(backend.electron_temperature_K, weights=backend.node_volume_A3)
        ),
        "lattice_temperature_max_K": float(backend.lattice_temperature_K.max()),
        "lattice_temperature_mean_K": float(
            np.average(backend.lattice_temperature_K, weights=backend.node_volume_A3)
        ),
        "ag_liquid_volume_fraction": float(
            np.dot(backend.node_volume_A3, backend.lattice_liquid_fraction)
            / backend.node_volume_A3.sum()
        ),
        "max_ag_liquid_fraction": float(backend.lattice_liquid_fraction.max()),
        "electron_energy_J": backend.electron_energy_J,
        "lattice_enthalpy_J": backend.lattice_energy_J,
        "source_energy_J": backend.injected_energy_J,
        "electron_to_lattice_cumulative_J": backend.coupling_energy_J,
        "electron_diffusion_net_J": backend.electron_diffusion_energy_J,
        "energy_balance_residual_J": residual,
        "energy_residual_relative_to_source": abs(residual)
        / max(abs(backend.injected_energy_J), 1.0e-30),
        "liquid_connectivity_at_f_l_ge_0p5": _largest_liquid_cluster(backend),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY)
    parser.add_argument("--electron-table", type=Path, default=DEFAULT_ELECTRON_TABLE)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--screen-json", type=Path, default=DEFAULT_SCREEN_JSON)
    parser.add_argument("--checkpoint-time-ps", type=float, default=0.05)
    parser.add_argument("--end-time-ps", type=float, default=0.1)
    parser.add_argument("--dt-fs", type=float, default=0.3125)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    args = parser.parse_args()
    if args.end_time_ps <= args.checkpoint_time_ps or args.dt_fs <= 0:
        parser.error("end time must exceed checkpoint time and time step must be positive")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    backend = QSGSComponentTTMBackend(
        args.geometry,
        args.electron_table,
        boundary_case="A",
        peak_flux_eV_ps_A2=5.0,
        gaussian_radius_A=30.0,
        absorption_depth_A=20.0,
        tile_xy=3,
        coarse_voxels=(1, 1, 1),
        linear_preconditioner="amg",
    )
    checkpoint = _restore_native_checkpoint(
        backend, args.checkpoint, args.screen_json, args.checkpoint_time_ps
    )
    dt_ps = args.dt_fs * 1.0e-3
    steps_float = (args.end_time_ps - args.checkpoint_time_ps) / dt_ps
    steps = int(round(steps_float))
    if not np.isclose(steps_float, steps, rtol=0.0, atol=1.0e-10):
        raise ValueError("extension duration must be an integer multiple of the time step")
    sample_stride = max(1, int(round(0.025 / dt_ps)))
    samples = [_sample(backend, "restored_50fs_checkpoint", args.dt_fs)]

    with args.log.open("w", encoding="utf-8", buffering=1) as log:
        log.write(json.dumps({"checkpoint": checkpoint}, sort_keys=True) + "\n")
        log.write(
            f"Starting native Ag voxel graph at t={backend.time_ps:.6g} ps; "
            f"continue heating to {args.end_time_ps:.6g} ps, dt={args.dt_fs:g} fs, "
            f"nodes={backend.graph.voxel_count.size:,}, P={backend.source_power_W:.12g} W\n"
        )
        for step in range(1, steps + 1):
            row = backend.advance(dt_ps, heat_on=True, coupling_on=True, electron_diffusion_on=True)
            if step % sample_stride == 0 or step == steps:
                sample = _sample(backend, "heated", args.dt_fs)
                samples.append(sample)
                log.write(json.dumps(sample, sort_keys=True) + "\n")
                log.write(
                    f"step={step}/{steps} t={backend.time_ps:.6g} ps "
                    f"Te,max={row['max_electron_temperature_K']:.3f} K "
                    f"Tl,max={row['max_lattice_temperature_K']:.6f} K "
                    f"liquid={row['ag_liquid_volume_fraction']:.8g} "
                    f"energy_residual={row['energy_balance_residual_J']:.4e} J\n"
                )

    backend.capture_snapshot(args.end_time_ps)
    snapshot_path = args.output.with_name("qsgs_component_native_heating_to_0p1ps_fields.npz")
    np.savez_compressed(snapshot_path, **backend.field_snapshots[args.end_time_ps])
    result = {
        "status": "NATIVE_QSGS_AG_HEATING_PILOT_COMPLETED",
        "baseline_commit": "3feb5bc",
        "started_from_50fs_native_voxel_checkpoint": True,
        "checkpoint": str(args.checkpoint),
        "checkpoint_validation": checkpoint,
        "geometry": args.geometry.name,
        "boundary_conditions": "periodic x/y; insulated z; boundary A; Ti3SiC2 is a geometric mask only",
        "graph": {
            "coarse_voxels_xyz": [1, 1, 1],
            "nodes": int(backend.graph.voxel_count.size),
            "edges": int(backend.graph.edge_i.size),
            "source_power_W": backend.source_power_W,
        },
        "heating": {
            "peak_flux_eV_ps_A2": 5.0,
            "gaussian_radius_A": 30.0,
            "gaussian_definition": "exp(-r^2/R^2)",
            "absorption_depth_A": 20.0,
            "dt_fs": args.dt_fs,
            "checkpoint_time_ps": args.checkpoint_time_ps,
            "end_time_ps": args.end_time_ps,
            "heat_source_remained_on": True,
            "source_energy_added_during_extension_J": backend.injected_energy_J
            - checkpoint["restored_source_energy_J"],
        },
        "samples": samples,
        "final_energy_balance": {
            "residual_J": samples[-1]["energy_balance_residual_J"],
            "residual_relative_to_source": samples[-1]["energy_residual_relative_to_source"],
            "electron_diffusion_net_J": backend.electron_diffusion_energy_J,
            "lattice_diffusion_net_J": backend.lattice_diffusion_energy_J,
        },
        "final_field_snapshot": snapshot_path.name,
        "runtime_seconds": time.perf_counter() - started,
        "interpretation": "Heating pilot only; determine melt-pool formation from Ag liquid fraction/connectivity, not electron temperature alone.",
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["samples"][-1], indent=2))
    print(f"Wrote {args.output}")
    print(f"Wrote {snapshot_path}")


if __name__ == "__main__":
    main()
