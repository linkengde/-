#!/usr/bin/env python3
"""Transient energy/convergence check for the component-lumped QSGS graph."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu

from qsgs_component_coarsegrain import build_component_graph
from qsgs_transport_benchmark import ROI_CASES, SOURCE


HERE = Path(__file__).resolve().parent
INTRINSIC_K = 1.0  # Manufactured benchmark value.
CE = 1.0  # Manufactured heat capacity per unit volume.
COMPONENT_GROUPS = {
    "1x1x2": (1, 1, 2),
    "1x2x1": (1, 2, 1),
    "2x1x1": (2, 1, 1),
}
FINAL_TIME_PS = 10.0
TIME_STEPS_PS = (0.5, 0.25)
ROIS = ("regression_1", "low_fill_through", "isolated_islands_1")


def native_voxel_laplacian(source_ag: np.ndarray, size_A: np.ndarray):
    ag = np.asarray(source_ag, dtype=bool)
    size = np.asarray(size_A, dtype=float)
    pitch = size / np.asarray(ag.shape)
    cell_volume = float(np.prod(pitch))
    ids = np.full(ag.shape, -1, dtype=np.int64)
    ids[ag] = np.arange(int(ag.sum()))
    edge_i, edge_j, edge_g = [], [], []
    for axis in range(3):
        low = [slice(None)] * 3
        high = [slice(None)] * 3
        low[axis], high[axis] = slice(0, -1), slice(1, None)
        i = ids[tuple(low)].ravel()
        j = ids[tuple(high)].ravel()
        keep = (i >= 0) & (j >= 0)
        edge_i.append(i[keep])
        edge_j.append(j[keep])
        edge_g.append(
            np.full(int(keep.sum()), INTRINSIC_K * cell_volume / pitch[axis] ** 2)
        )
    return _assemble_laplacian(int(ag.sum()), np.concatenate(edge_i), np.concatenate(edge_j), np.concatenate(edge_g)), cell_volume


def _assemble_laplacian(n: int, i: np.ndarray, j: np.ndarray, g: np.ndarray):
    diagonal = np.bincount(i, weights=g, minlength=n).astype(float)
    diagonal += np.bincount(j, weights=g, minlength=n)
    rows = np.concatenate((np.arange(n), i, j))
    cols = np.concatenate((np.arange(n), j, i))
    values = np.concatenate((diagonal, -g, -g))
    return coo_matrix((values, (rows, cols)), shape=(n, n)).tocsr()


def aggregate_component_temperature(
    component_id: np.ndarray,
    ag_mask: np.ndarray,
    temperature_K: np.ndarray,
    component_count: int,
) -> np.ndarray:
    """Volume-average native Ag voxel temperatures onto component nodes."""
    ids = np.asarray(component_id, dtype=np.int64)
    mask = np.asarray(ag_mask, dtype=bool)
    temperature = np.asarray(temperature_K, dtype=float)
    if ids.shape != mask.shape or temperature.shape != mask.shape:
        raise ValueError("component_id, ag_mask and temperature_K must have matching shapes")
    active_ids = ids[mask]
    if np.any(active_ids < 0) or np.any(active_ids >= component_count):
        raise ValueError("Every Ag voxel must map to a valid component node")
    return np.bincount(
        active_ids,
        weights=temperature[mask],
        minlength=component_count,
    ) / np.bincount(active_ids, minlength=component_count)


def implicit_history(laplacian, capacity: np.ndarray, initial: np.ndarray, dt_ps: float, steps: int):
    mass = diags(capacity)
    matrix = (mass + dt_ps * laplacian).tocsc()
    factor = splu(matrix)
    rhs_scale = capacity
    temperature = np.asarray(initial, dtype=float).copy()
    initial_energy = float(np.dot(capacity, temperature))
    for _ in range(steps):
        temperature = factor.solve(rhs_scale * temperature)
    final_energy = float(np.dot(capacity, temperature))
    return temperature, abs(final_energy - initial_energy) / max(abs(initial_energy), 1.0e-30)


def run_roi(case: dict, source_ag: np.ndarray, source_size_A: np.ndarray) -> dict:
    start, shape = case["start"], case["shape"]
    roi = source_ag[tuple(slice(s, s + n) for s, n in zip(start, shape))]
    size = source_size_A * np.asarray(shape) / np.asarray(source_ag.shape)
    pitch = size / np.asarray(shape)
    voxel_laplacian, voxel_volume = native_voxel_laplacian(roi, size)

    x = (np.arange(roi.shape[0]) + 0.5) * pitch[0]
    y = (np.arange(roi.shape[1]) + 0.5) * pitch[1]
    z = (np.arange(roi.shape[2]) + 0.5) * pitch[2]
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    hot = 500.0 * np.exp(
        -((xx - size[0] / 2.0) ** 2 + (yy - size[1] / 2.0) ** 2 + (zz - size[2] / 2.0) ** 2)
        / (2.0 * 12.0**2)
    )
    voxel_initial_field = 300.0 + hot
    voxel_initial_active = voxel_initial_field[roi]
    voxel_capacity = np.full(int(roi.sum()), CE * voxel_volume)
    component_models = {}
    for name, groups in COMPONENT_GROUPS.items():
        graph = build_component_graph(roi, size, groups, intrinsic_k=INTRINSIC_K)
        node_capacity = graph.voxel_count.astype(float) * CE * voxel_volume
        node_initial = aggregate_component_temperature(
            graph.component_id,
            roi,
            voxel_initial_field,
            graph.voxel_count.size,
        )
        laplacian = _assemble_laplacian(
            int(graph.voxel_count.size), graph.edge_i, graph.edge_j, graph.edge_conductance
        )
        component_models[name] = {
            "graph": graph,
            "capacity": node_capacity,
            "initial": node_initial,
            "laplacian": laplacian,
        }
    voxel_results = {}
    component_results = {name: {} for name in COMPONENT_GROUPS}
    for dt_ps in TIME_STEPS_PS:
        steps = round(FINAL_TIME_PS / dt_ps)
        voxel_final, voxel_energy_error = implicit_history(
            voxel_laplacian, voxel_capacity, voxel_initial_active, dt_ps, steps
        )
        voxel_field = np.full(roi.shape, np.nan)
        voxel_field[roi] = voxel_final
        voxel_results[str(dt_ps)] = {
            "energy_relative_drift": voxel_energy_error,
            "mean_temperature_K": float(np.average(voxel_final, weights=voxel_capacity)),
            "peak_temperature_K": float(voxel_final.max()),
        }
        for name, model in component_models.items():
            graph = model["graph"]
            capacity = model["capacity"]
            node_id = graph.component_id[roi]
            component_final, component_energy_error = implicit_history(
                model["laplacian"], capacity, model["initial"], dt_ps, steps
            )
            component_reference = np.bincount(
                node_id[node_id >= 0],
                weights=voxel_field[roi][node_id >= 0],
                minlength=graph.voxel_count.size,
            ) / graph.voxel_count
            difference = component_final - component_reference
            rms = float(np.sqrt(np.average(difference**2, weights=capacity)))
            component_results[name][str(dt_ps)] = {
                "energy_relative_drift": component_energy_error,
                "mean_temperature_K": float(np.average(component_final, weights=capacity)),
                "peak_temperature_K": float(component_final.max()),
                "component_vs_native_mean_max_abs_K": float(np.max(np.abs(difference))),
                "component_vs_native_mean_weighted_rms_K": rms,
            }
    component_graph_metrics = {
        name: {
            "coarse_voxels": list(COMPONENT_GROUPS[name]),
            "component_nodes": int(model["graph"].voxel_count.size),
            "component_capacity_voxel_count": int(model["graph"].voxel_count.sum()),
        }
        for name, model in component_models.items()
    }
    time_refinement = {
        "native_peak_delta_0p5_vs_0p25_K": abs(
            voxel_results["0.5"]["peak_temperature_K"]
            - voxel_results["0.25"]["peak_temperature_K"]
        ),
        "component_peak_delta_0p5_vs_0p25_K": {
            name: abs(
                values["0.5"]["peak_temperature_K"]
                - values["0.25"]["peak_temperature_K"]
            )
            for name, values in component_results.items()
        },
    }
    return {
        "motif": case["motif"],
        "roi_start_voxels_xyz": list(start),
        "roi_shape_voxels_xyz": list(shape),
        "ag_volume_fraction": float(roi.mean()),
        "native_ag_voxels": int(roi.sum()),
        "component_graphs": component_graph_metrics,
        "native_voxel_reference": voxel_results,
        "component_lumped": component_results,
        "time_refinement": time_refinement,
    }


def main() -> None:
    with np.load(SOURCE, allow_pickle=False) as data:
        source_ag = ~np.asarray(data["skeleton"], dtype=bool)
        source_size_A = np.asarray(data["size_A"], dtype=float)
    cases = {case["name"]: case for case in ROI_CASES}
    report = {
        "status": "QSGS_COMPONENT_GRAPH_TRANSIENT_AUDIT",
        "source_geometry": SOURCE.name,
        "intrinsic_k": INTRINSIC_K,
        "Ce_per_volume": CE,
        "coarse_voxels_per_component_block": {
            name: list(groups) for name, groups in COMPONENT_GROUPS.items()
        },
        "final_time_ps": FINAL_TIME_PS,
        "time_steps_ps": list(TIME_STEPS_PS),
        "initial_condition": "300 K plus a 500 K Gaussian excess centred in each window, sigma=12 Å; insulating ROI exterior; no source or interface physics",
        "interpretation": "manufactured constant-property transient; verifies capacity aggregation, energy conservation, and dynamic coarse-model error only",
        "regions": {
            name: run_roi(cases[name], source_ag, source_size_A) for name in ROIS
        },
    }
    output = HERE / "nested_grid_diagnostics" / "qsgs_component_transient_benchmark.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
