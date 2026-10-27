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


def _periodic_interval(indices: np.ndarray, n_cells: int, spacing_A: float) -> dict:
    occupied = np.unique(np.asarray(indices, dtype=np.int64))
    if occupied.size == 0:
        return {"start_index": None, "end_index": None, "extent_A": 0.0, "wraps_boundary": False}
    if occupied.size == n_cells:
        start, span = 0, n_cells
    else:
        gaps = (np.roll(occupied, -1) - occupied) % n_cells
        empty = gaps - 1
        gap_index = int(np.argmax(empty))
        start = int(occupied[(gap_index + 1) % occupied.size])
        span = int(n_cells - empty[gap_index])
    end = int((start + span - 1) % n_cells)
    return {
        "start_index": start,
        "end_index": end,
        "extent_A": float(span * spacing_A),
        "wraps_boundary": bool(start + span > n_cells),
    }


def _cluster_spatial_metrics(backend: QSGSComponentTTMBackend, node_ids: np.ndarray) -> dict:
    if node_ids.size == 0:
        return {
            "bbox_voxel_xyz": None,
            "bbox_A": None,
            "extent_fraction_xyz": None,
            "center_A": None,
            "source_center_A": None,
            "distance_to_source_center_A": None,
        }
    if not np.all(backend.graph.voxel_count == 1):
        raise ValueError("Cluster coordinates require the native one-node-per-Ag-voxel graph")
    coords = backend.graph.block_index[node_ids]
    shape = np.asarray(backend.source_ag.shape, dtype=int)
    spacing = np.asarray(backend.spacing_A, dtype=float)
    intervals = [_periodic_interval(coords[:, a], int(shape[a]), float(spacing[a])) for a in range(2)]
    intervals.append({
        "start_index": int(coords[:, 2].min()),
        "end_index": int(coords[:, 2].max()),
        "extent_A": float((coords[:, 2].max() - coords[:, 2].min() + 1) * spacing[2]),
        "wraps_boundary": False,
    })
    bbox_voxel = {}
    bbox_A = {}
    extents = []
    centers = []
    for axis, name in enumerate(("x", "y", "z")):
        row = intervals[axis]
        start = row["start_index"]
        end = row["end_index"]
        extent = row["extent_A"]
        bbox_voxel[name] = {"start": start, "end": end, "wraps_periodic_boundary": row["wraps_boundary"]}
        start_A = float(start * spacing[axis])
        bbox_A[name] = {
            "start": start_A,
            "end_unwrapped": float(start_A + extent),
            "extent": extent,
            "wraps_periodic_boundary": row["wraps_boundary"],
        }
        extents.append(float(extent))
        if axis < 2:
            angles = 2.0 * np.pi * (coords[:, axis] + 0.5) / shape[axis]
            vector = np.mean(np.exp(1j * angles))
            center_index = float((np.angle(vector) % (2.0 * np.pi)) * shape[axis] / (2.0 * np.pi))
            centers.append(float(center_index * spacing[axis]))
        else:
            centers.append(float((coords[:, axis].mean() + 0.5) * spacing[axis]))
    source_coords = backend.graph.block_index
    source_weights = np.asarray(backend.source_node_W, dtype=float)
    source_center = []
    for axis in range(3):
        if axis < 2:
            angles = 2.0 * np.pi * (source_coords[:, axis] + 0.5) / shape[axis]
            vector = np.average(np.exp(1j * angles), weights=source_weights)
            center_index = float((np.angle(vector) % (2.0 * np.pi)) * shape[axis] / (2.0 * np.pi))
            source_center.append(float(center_index * spacing[axis]))
        else:
            source_center.append(float(np.average(source_coords[:, axis] + 0.5, weights=source_weights) * spacing[axis]))
    displacement = np.asarray(centers) - np.asarray(source_center)
    for axis in range(2):
        displacement[axis] = (displacement[axis] + 0.5 * backend.size_A[axis]) % backend.size_A[axis] - 0.5 * backend.size_A[axis]
    return {
        "bbox_voxel_xyz": bbox_voxel,
        "bbox_A": bbox_A,
        "extent_fraction_xyz": [float(extents[a] / backend.size_A[a]) for a in range(3)],
        "center_A": centers,
        "source_center_A": source_center,
        "distance_to_source_center_A": float(np.linalg.norm(displacement)),
    }


def _largest_liquid_cluster(
    backend: QSGSComponentTTMBackend, threshold: float = 0.5
) -> tuple[dict, np.ndarray]:
    active = backend.lattice_liquid_fraction >= threshold
    count = int(np.count_nonzero(active))
    if count == 0:
        return ({
            "threshold_liquid_fraction": threshold,
            "active_node_count": 0,
            "active_threshold_liquid_volume_A3": 0.0,
            "largest_cluster_node_count": 0,
            "largest_cluster_volume_A3": 0.0,
            "largest_cluster_liquid_fraction_weighted_volume_A3": 0.0,
            "equal_size_largest_cluster_count": 0,
            "largest_cluster_is_tied": False,
            "largest_cluster_source_power_fraction": 0.0,
            **_cluster_spatial_metrics(backend, np.empty(0, dtype=np.int64)),
        }, np.empty(0, dtype=np.int64))
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
    # Select exactly one largest component. np.argmax deterministically chooses
    # the first label in a tie; never sum disconnected equal-size components.
    largest_label = int(np.argmax(sizes))
    cluster_ids = np.sort(active_ids[labels == largest_label])
    volume = float(backend.node_volume_A3[cluster_ids].sum())
    liquid_volume = float(np.dot(
        backend.node_volume_A3[cluster_ids], backend.lattice_liquid_fraction[cluster_ids]
    ))
    metadata = {
        "threshold_liquid_fraction": threshold,
        "active_node_count": count,
        "active_threshold_liquid_volume_A3": float(backend.node_volume_A3[active].sum()),
        "largest_cluster_node_count": largest,
        "largest_cluster_volume_A3": volume,
        "largest_cluster_liquid_fraction_weighted_volume_A3": liquid_volume,
        "equal_size_largest_cluster_count": int(np.count_nonzero(sizes == largest)),
        "largest_cluster_is_tied": bool(np.count_nonzero(sizes == largest) > 1),
        "largest_cluster_source_power_fraction": float(
            backend.source_node_W[cluster_ids].sum() / max(backend.source_power_W, 1.0e-300)
        ),
        **_cluster_spatial_metrics(backend, cluster_ids),
    }
    return metadata, cluster_ids


def _restore_native_checkpoint(
    backend: QSGSComponentTTMBackend,
    checkpoint: Path,
    checkpoint_json: Path,
    checkpoint_time_ps: float,
) -> dict:
    summary = json.loads(checkpoint_json.read_text(encoding="utf-8"))
    if "samples" in summary:
        matches = [
            row for row in summary["samples"]
            if abs(float(row["time_ps"]) - checkpoint_time_ps) <= 1.0e-9
        ]
        if not matches:
            raise ValueError(
                f"No saved sample at {checkpoint_time_ps:g} ps in {checkpoint_json}"
            )
        run = matches[-1]
    else:
        # Original 50 fs matched-screen format.
        if not np.isclose(checkpoint_time_ps, 0.05, rtol=0.0, atol=1.0e-9):
            raise ValueError("The matched-screen checkpoint is valid only at 0.05 ps")
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
            if values.shape == backend.source_ag.shape:
                values = values[backend.source_ag]
            if values.shape != ag_node_ids.shape:
                raise ValueError(f"Checkpoint field {name} has unexpected shape {values.shape}")
            destination[ag_node_ids] = values
        if "electron_to_lattice_energy_eV_A3" in state:
            coupling_density = np.asarray(state["electron_to_lattice_energy_eV_A3"], dtype=float)
            coupling_units = "eV_A3"
        elif "electron_to_lattice_cumulative_energy_density_J_m3" in state:
            coupling_density = np.asarray(
                state["electron_to_lattice_cumulative_energy_density_J_m3"], dtype=float
            )
            coupling_units = "J_m3"
        else:
            raise ValueError("Checkpoint is missing cumulative electron-lattice coupling field")
        if coupling_density.shape == backend.source_ag.shape:
            coupling_density = coupling_density[backend.source_ag]
        if coupling_density.shape != ag_node_ids.shape:
            raise ValueError("Checkpoint cumulative coupling field has unexpected shape")
        if coupling_units == "J_m3":
            coupling_density = coupling_density * ANGSTROM3_TO_M3 / EV_J
        backend.cumulative_coupling_eV[ag_node_ids] = (
            coupling_density * backend.node_volume_A3[ag_node_ids]
        )

    backend.time_ps = float(checkpoint_time_ps)
    backend.injected_energy_J = float(run["source_energy_J"])
    backend.coupling_energy_J = float(run["electron_to_lattice_cumulative_J"])
    backend.electron_diffusion_energy_J = float(
        run.get("electron_diffusion_energy_J", run.get("electron_diffusion_net_J", 0.0))
    )
    if "lattice_diffusion_net_J" in run:
        backend.lattice_diffusion_energy_J = float(run["lattice_diffusion_net_J"])
    else:
        backend.lattice_diffusion_energy_J = float(
            summary.get("final_energy_balance", {}).get("lattice_diffusion_net_J", 0.0)
        )
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
    liquid_volume = float(np.dot(backend.node_volume_A3, backend.lattice_liquid_fraction))
    liquid_cluster, cluster_ids = _largest_liquid_cluster(backend)
    sample = {
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
            liquid_volume / backend.node_volume_A3.sum()
        ),
        "ag_liquid_volume_A3": liquid_volume,
        "initial_melting_detected": bool(np.any(backend.lattice_liquid_fraction > 0.0)),
        "max_ag_liquid_fraction": float(backend.lattice_liquid_fraction.max()),
        "minimum_pool_nodes": 50,
        "minimum_pool_volume_A3": float(50 * backend.voxel_volume_A3),
        "electron_energy_J": backend.electron_energy_J,
        "lattice_enthalpy_J": backend.lattice_energy_J,
        "source_energy_J": backend.injected_energy_J,
        "electron_to_lattice_cumulative_J": backend.coupling_energy_J,
        "electron_diffusion_net_J": backend.electron_diffusion_energy_J,
        "energy_balance_residual_J": residual,
        "energy_residual_relative_to_source": abs(residual)
        / max(abs(backend.injected_energy_J), 1.0e-30),
        "liquid_connectivity_at_f_l_ge_0p5": liquid_cluster,
        "_largest_cluster_node_ids": cluster_ids,
    }
    sample["largest_cluster_volume_fraction_of_total_ag"] = float(
        liquid_cluster["largest_cluster_volume_A3"] / backend.node_volume_A3.sum()
    )
    return sample


def _add_cluster_overlap(sample: dict, previous_ids: np.ndarray | None) -> np.ndarray:
    current_ids = np.asarray(sample.pop("_largest_cluster_node_ids"), dtype=np.int64)
    if previous_ids is None or previous_ids.size == 0 or current_ids.size == 0:
        common = 0
        jaccard = None
        previous_fraction = None
        current_fraction = None
    else:
        common = int(np.intersect1d(previous_ids, current_ids, assume_unique=True).size)
        jaccard = float(common / (previous_ids.size + current_ids.size - common))
        previous_fraction = float(common / previous_ids.size)
        current_fraction = float(common / current_ids.size)
    sample["largest_cluster_overlap_with_previous_sample"] = {
        "shared_node_count": common,
        "jaccard_fraction": jaccard,
        "fraction_of_previous_cluster": previous_fraction,
        "fraction_of_current_cluster": current_fraction,
    }
    return current_ids


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY)
    parser.add_argument("--electron-table", type=Path, default=DEFAULT_ELECTRON_TABLE)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--screen-json", type=Path, default=DEFAULT_SCREEN_JSON)
    parser.add_argument(
        "--checkpoint-json",
        type=Path,
        help="Summary JSON matching --checkpoint; defaults to the original 50 fs screen JSON",
    )
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
        backend,
        args.checkpoint,
        args.checkpoint_json or args.screen_json,
        args.checkpoint_time_ps,
    )
    dt_ps = args.dt_fs * 1.0e-3
    steps_float = (args.end_time_ps - args.checkpoint_time_ps) / dt_ps
    steps = int(round(steps_float))
    if not np.isclose(steps_float, steps, rtol=0.0, atol=1.0e-10):
        raise ValueError("extension duration must be an integer multiple of the time step")
    sample_stride = max(1, int(round(0.025 / dt_ps)))
    first_sample = _sample(backend, f"restored_{args.checkpoint_time_ps:g}ps_checkpoint", args.dt_fs)
    previous_cluster_ids = _add_cluster_overlap(first_sample, None)
    samples = [first_sample]

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
                previous_cluster_ids = _add_cluster_overlap(sample, previous_cluster_ids)
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
    snapshot_path = args.output.with_name(f"{args.output.stem}_fields.npz")
    np.savez_compressed(snapshot_path, **backend.field_snapshots[args.end_time_ps])
    result = {
        "status": "NATIVE_QSGS_AG_HEATING_PILOT_COMPLETED",
        "baseline_commit": "3feb5bc",
        "started_from_50fs_native_voxel_checkpoint": bool(
            np.isclose(args.checkpoint_time_ps, 0.05, rtol=0.0, atol=1.0e-9)
        ),
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
