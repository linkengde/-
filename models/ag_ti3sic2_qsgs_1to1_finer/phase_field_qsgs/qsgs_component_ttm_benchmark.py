#!/usr/bin/env python3
"""Local 20 fs QSGS two-temperature comparison on native and component graphs.

The graph retains actual native-voxel Ag connectivity. Electron properties and
Ag enthalpy use the repository's candidate tables. Each ROI is an insulated
local diagnostic with a Gaussian/exponential source integrated over the binary
Ag mask; it is not a production-domain melt calculation.
"""

from __future__ import annotations

import json
import argparse
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import cg

MODEL_DIR = Path(__file__).resolve().parents[1]
PHASE_DIR = Path(__file__).resolve().parent
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))
if str(PHASE_DIR) not in sys.path:
    sys.path.insert(0, str(PHASE_DIR))

from ag_ttm_reference import AgParameters  # noqa: E402
from qsgs_component_coarsegrain import ComponentGraph, build_component_graph  # noqa: E402
from qsgs_source_integration import integrate_qsgs_source_power_W  # noqa: E402
from qsgs_transport_benchmark import ROI_CASES, SOURCE  # noqa: E402
from thermophysical_model import AgEnthalpy  # noqa: E402


EV_J = 1.602176634e-19
ANGSTROM3_TO_M3 = 1.0e-30
ELECTRON_TABLE = MODEL_DIR / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv"
GROUPS = {"native_voxel": (1, 1, 1), "component_2x1x1": (2, 1, 1)}
ROIS = ("regression_1", "low_fill_through", "isolated_islands_1")
TIME_STEPS_PS = (0.0003125, 0.00015625)  # 0.3125 and 0.15625 fs.
PEAK_FLUX_EV_PS_A2 = 5.0
RADIUS_A = 30.0
SKIN_DEPTH_A = 20.0


def _laplacian(graph: ComponentGraph, edge_conductance: np.ndarray):
    n = graph.voxel_count.size
    i, j = graph.edge_i, graph.edge_j
    g = np.asarray(edge_conductance, dtype=float)
    diagonal = np.bincount(i, weights=g, minlength=n)
    diagonal += np.bincount(j, weights=g, minlength=n)
    rows = np.concatenate((np.arange(n), i, j))
    cols = np.concatenate((np.arange(n), j, i))
    values = np.concatenate((diagonal, -g, -g))
    return coo_matrix((values, (rows, cols)), shape=(n, n)).tocsr()


def _ue_slope(parameters: AgParameters, temperature_K: np.ndarray) -> np.ndarray:
    table_t = parameters.temperature_K
    slopes = np.diff(parameters.ue_eV_A3) / np.diff(table_t)
    idx = np.clip(np.searchsorted(table_t, temperature_K, side="right") - 1, 0, len(slopes) - 1)
    return slopes[idx]


def _implicit_electron_diffusion(
    graph: ComponentGraph,
    node_volume_A3: np.ndarray,
    temperature_K: np.ndarray,
    dt_ps: float,
    parameters: AgParameters,
) -> np.ndarray:
    """Backward-Euler diffusion; conductivity is frozen at the step start."""
    if graph.edge_i.size == 0:
        return temperature_K.copy()
    t0 = temperature_K
    edge_ke = 2.0 * parameters.ke(t0[graph.edge_i]) * parameters.ke(t0[graph.edge_j]) / (
        parameters.ke(t0[graph.edge_i]) + parameters.ke(t0[graph.edge_j])
    )
    operator = _laplacian(graph, graph.edge_conductance * edge_ke)
    e0 = node_volume_A3 * parameters.ue(t0)
    temperature = t0.copy()
    for _ in range(12):
        residual = node_volume_A3 * (parameters.ue(temperature) - parameters.ue(t0)) + dt_ps * (
            operator @ temperature
        )
        if float(np.max(np.abs(residual))) < max(1.0e-13, float(np.max(np.abs(e0))) * 2.0e-10):
            return temperature
        jacobian = diags(node_volume_A3 * _ue_slope(parameters, temperature)) + dt_ps * operator
        diagonal = np.maximum(jacobian.diagonal(), 1.0e-30)
        delta, info = cg(
            jacobian,
            -residual,
            M=diags(1.0 / diagonal),
            rtol=1.0e-10,
            atol=1.0e-14,
            maxiter=min(3000, max(200, 5 * len(temperature))),
        )
        if info != 0:
            raise RuntimeError(f"Electron diffusion CG did not converge (info={info})")
        old_norm = float(np.linalg.norm(residual))
        alpha = 1.0
        accepted = False
        while alpha >= 1.0 / 128.0:
            trial = temperature + alpha * delta
            if np.any(trial < parameters.Tmin) or np.any(trial > parameters.Tmax):
                alpha *= 0.5
                continue
            trial_residual = node_volume_A3 * (parameters.ue(trial) - parameters.ue(t0)) + dt_ps * (
                operator @ trial
            )
            if float(np.linalg.norm(trial_residual)) < old_norm:
                temperature = trial
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            raise RuntimeError("Electron diffusion Newton line search failed")
    raise RuntimeError("Electron diffusion Newton solve did not converge")


def _lattice_h_vol_eV_A3(
    temperature_K: np.ndarray,
    liquid_fraction: np.ndarray,
    lattice: AgEnthalpy,
) -> np.ndarray:
    h_mass = lattice.specific_enthalpy_J_kg(temperature_K, liquid_fraction)
    return lattice.reference_density_kg_m3 * np.asarray(h_mass) * ANGSTROM3_TO_M3 / EV_J


def _state_from_lattice_h_vol(
    h_vol_eV_A3: np.ndarray,
    lattice: AgEnthalpy,
) -> tuple[np.ndarray, np.ndarray]:
    h_mass = np.asarray(h_vol_eV_A3) * EV_J / (
        lattice.reference_density_kg_m3 * ANGSTROM3_TO_M3
    )
    state = lattice.state_from_specific_enthalpy(h_mass)
    return np.asarray(state.temperature_K), np.asarray(state.liquid_fraction)


def _implicit_local_coupling(
    node_volume_A3: np.ndarray,
    electron_temperature_K: np.ndarray,
    lattice_temperature_K: np.ndarray,
    lattice_fraction: np.ndarray,
    dt_ps: float,
    parameters: AgParameters,
    lattice: AgEnthalpy,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Implicit, cellwise energy-conserving electron/lattice exchange."""
    te0 = electron_temperature_K
    tl0 = lattice_temperature_K
    h0 = _lattice_h_vol_eV_A3(tl0, lattice_fraction, lattice)
    total_density = parameters.ue(te0) + h0
    low, high = np.minimum(te0, tl0).copy(), np.maximum(te0, tl0).copy()

    def residual(te: np.ndarray):
        he = total_density - parameters.ue(te)
        tl, fl = _state_from_lattice_h_vol(he, lattice)
        return parameters.ue(te) - parameters.ue(te0) + dt_ps * parameters.g(te) * (te - tl), tl, fl

    f_low = residual(low)[0]
    f_high = residual(high)[0]
    if np.any(f_low > 2.0e-12) or np.any(f_high < -2.0e-12):
        raise RuntimeError("Electron-lattice coupling root is not bracketed")
    for _ in range(48):
        mid = 0.5 * (low + high)
        f_mid = residual(mid)[0]
        positive = f_mid > 0.0
        high = np.where(positive, mid, high)
        low = np.where(positive, low, mid)
    te1 = 0.5 * (low + high)
    h1 = total_density - parameters.ue(te1)
    tl1, fl1 = _state_from_lattice_h_vol(h1, lattice)
    transfer_eV = node_volume_A3 * (parameters.ue(te0) - parameters.ue(te1))
    return te1, tl1, fl1, float(transfer_eV.sum())


def run_graph(
    graph: ComponentGraph,
    roi_ag: np.ndarray,
    source_power_W: np.ndarray,
    voxel_volume_A3: float,
    dt_ps: float,
    final_time_ps: float,
    parameters: AgParameters,
    lattice: AgEnthalpy,
) -> dict:
    node_volume = graph.voxel_count.astype(float) * voxel_volume_A3
    node_source_W = np.bincount(
        graph.component_id[roi_ag],
        weights=source_power_W[roi_ag],
        minlength=graph.voxel_count.size,
    )
    nsteps = round(final_time_ps / dt_ps)
    if not np.isclose(nsteps * dt_ps, final_time_ps, rtol=0.0, atol=1.0e-14):
        raise ValueError("final_time_ps must be divisible by dt_ps")
    te = np.full(graph.voxel_count.size, parameters.temperature_K[0], dtype=float)
    tl = np.full_like(te, 300.0)
    fl = np.zeros_like(te)
    injected_eV = 0.0
    transferred_eV = 0.0
    initial_energy_eV = float(np.dot(node_volume, parameters.ue(te)) + np.dot(node_volume, _lattice_h_vol_eV_A3(tl, fl, lattice)))
    for _ in range(nsteps):
        half_source_eV = node_source_W * (0.5 * dt_ps * 1.0e-12 / EV_J)
        te = parameters.temperature_from_ue(parameters.ue(te) + half_source_eV / node_volume)
        injected_eV += float(half_source_eV.sum())
        te = _implicit_electron_diffusion(graph, node_volume, te, 0.5 * dt_ps, parameters)
        te, tl, fl, transfer = _implicit_local_coupling(
            node_volume, te, tl, fl, dt_ps, parameters, lattice
        )
        transferred_eV += transfer
        te = _implicit_electron_diffusion(graph, node_volume, te, 0.5 * dt_ps, parameters)
        te = parameters.temperature_from_ue(parameters.ue(te) + half_source_eV / node_volume)
        injected_eV += float(half_source_eV.sum())
    electron_eV = node_volume * parameters.ue(te)
    lattice_eV = node_volume * _lattice_h_vol_eV_A3(tl, fl, lattice)
    final_energy_eV = float(electron_eV.sum() + lattice_eV.sum())
    residual_eV = final_energy_eV - initial_energy_eV - injected_eV
    return {
        "dt_ps": dt_ps,
        "steps": nsteps,
        "electron_temperature_K": te,
        "lattice_temperature_K": tl,
        "lattice_liquid_fraction": fl,
        "electron_energy_eV": electron_eV,
        "lattice_enthalpy_eV": lattice_eV,
        "source_power_W": float(source_power_W.sum()),
        "injected_energy_eV": injected_eV,
        "electron_to_lattice_energy_eV": transferred_eV,
        "electron_energy_final_eV": float(electron_eV.sum()),
        "lattice_enthalpy_final_eV": float(lattice_eV.sum()),
        "energy_residual_eV": residual_eV,
        "energy_residual_relative": abs(residual_eV) / max(abs(injected_eV), 1.0e-30),
        "max_electron_temperature_K": float(te.max()),
        "mean_electron_temperature_K": float(np.average(te, weights=node_volume)),
        "max_lattice_temperature_K": float(tl.max()),
        "mean_lattice_temperature_K": float(np.average(tl, weights=node_volume)),
        "max_liquid_fraction": float(fl.max()),
        "liquid_ag_volume_fraction": float(np.dot(node_volume, fl) / node_volume.sum()),
    }


def _compare_component_to_native(
    native_graph: ComponentGraph,
    native_result: dict,
    coarse_graph: ComponentGraph,
    coarse_result: dict,
    roi_ag: np.ndarray,
    voxel_volume_A3: float,
) -> dict[str, float]:
    native_ids = native_graph.component_id[roi_ag]
    coarse_ids = coarse_graph.component_id[roi_ag]
    volume_weights = np.full(int(roi_ag.sum()), voxel_volume_A3)
    comparisons: dict[str, float] = {}
    comparisons: dict[str, float] = {}
    for field in ("electron_temperature_K", "lattice_temperature_K"):
        native = native_result[field][native_ids]
        coarse = coarse_result[field][coarse_ids]
        difference = coarse - native
        comparisons[f"{field}_weighted_rms"] = float(
            np.sqrt(np.average(difference**2, weights=volume_weights))
        )
        comparisons[f"{field}_absolute_max"] = float(np.max(np.abs(difference)))
    native_node_volume = native_graph.voxel_count.astype(float) * voxel_volume_A3
    coarse_node_volume = coarse_graph.voxel_count.astype(float) * voxel_volume_A3
    for field in ("electron_energy_eV", "lattice_enthalpy_eV"):
        native_density = native_result[field] / native_node_volume
        coarse_density = coarse_result[field] / coarse_node_volume
        native = native_density[native_ids]
        coarse = coarse_density[coarse_ids]
        difference = coarse - native
        comparisons[f"{field}_density_weighted_rms_eV_A3"] = float(
            np.sqrt(np.average(difference**2, weights=volume_weights))
        )
        scale = float(np.average(np.abs(native), weights=volume_weights))
        comparisons[f"{field}_density_relative_L1"] = float(
            np.average(np.abs(difference), weights=volume_weights) / max(scale, 1.0e-30)
        )
    return comparisons


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration-ps", type=float, default=0.02)
    args = parser.parse_args()
    final_time_ps = float(args.duration_ps)
    if not np.isfinite(final_time_ps) or final_time_ps <= 0.0:
        raise ValueError("duration-ps must be finite and positive")
    with np.load(SOURCE, allow_pickle=False) as data:
        full_ag = ~np.asarray(data["skeleton"], dtype=bool)
        full_size_A = np.asarray(data["size_A"], dtype=float)
    parameters = AgParameters.from_csv(ELECTRON_TABLE, allow_low_temperature_extension=False)
    lattice = AgEnthalpy()
    cases = {case["name"]: case for case in ROI_CASES}
    result = {
        "status": "LOCAL_QSGS_COMPONENT_TTM_20FS_DIAGNOSTIC",
        "source_geometry": SOURCE.name,
        "candidate_electron_table": ELECTRON_TABLE.name,
        "ag_enthalpy_table": lattice.table_path.name,
        "rois": list(ROIS),
        "component_group": [2, 1, 1],
        "final_time_ps": final_time_ps,
        "time_steps_ps": list(TIME_STEPS_PS),
        "heating": {
            "peak_flux_eV_ps_A2": PEAK_FLUX_EV_PS_A2,
            "gaussian_radius_A": RADIUS_A,
            "gaussian_definition": "exp(-r^2/R^2)",
            "absorption_depth_A": SKIN_DEPTH_A,
            "source_mapping": "exact separable integration over binary Ag voxels in each local ROI; Gaussian centred on ROI xy and exponential decay measured from ROI top",
        },
        "boundary_conditions": "local ROI exterior insulated; actual voxel-face Ag connectivity; no scaffold heat capacity or interface transport; Ag lattice conduction omitted (Case A)",
        "interpretation": "candidate-property continuum dual-temperature numerical comparison only; not full-domain melt prediction or validation of QSGS heat-sink physics",
        "regions": {},
    }
    for name in ROIS:
        case = cases[name]
        start, shape = case["start"], case["shape"]
        roi = full_ag[tuple(slice(s, s + n) for s, n in zip(start, shape))]
        roi_size = full_size_A * np.asarray(shape) / np.asarray(full_ag.shape)
        spacing = roi_size / np.asarray(shape)
        voxel_volume = float(np.prod(spacing))
        peak_flux_W_m2 = PEAK_FLUX_EV_PS_A2 * EV_J / 1.0e-12 / 1.0e-20
        source_power = integrate_qsgs_source_power_W(
            roi, roi.shape, roi_size, peak_flux_W_m2, RADIUS_A, SKIN_DEPTH_A
        )
        native_graph = build_component_graph(roi, roi_size, (1, 1, 1), intrinsic_k=1.0)
        coarse_graph = build_component_graph(roi, roi_size, (2, 1, 1), intrinsic_k=1.0)
        region = {
            "motif": case["motif"],
            "roi_start_voxels_xyz": list(start),
            "roi_shape_voxels_xyz": list(shape),
            "roi_size_A": roi_size.tolist(),
            "source_voxel_spacing_A": spacing.tolist(),
            "ag_volume_fraction": float(roi.mean()),
            "native_ag_voxels": int(roi.sum()),
            "component_nodes": int(coarse_graph.voxel_count.size),
            "component_node_fraction_of_native_voxels": float(coarse_graph.voxel_count.size / roi.sum()),
            "preserved_ag_voxel_count": int(coarse_graph.voxel_count.sum()),
            "source_power_W": float(source_power.sum()),
            "time_integrations": {},
        }
        for dt_ps in TIME_STEPS_PS:
            native_fields = run_graph(
                native_graph, roi, source_power, voxel_volume, dt_ps,
                final_time_ps, parameters, lattice,
            )
            coarse_fields = run_graph(
                coarse_graph, roi, source_power, voxel_volume, dt_ps,
                final_time_ps, parameters, lattice,
            )
            native_summary = {k: v for k, v in native_fields.items() if not isinstance(v, np.ndarray)}
            coarse_summary = {k: v for k, v in coarse_fields.items() if not isinstance(v, np.ndarray)}
            region["time_integrations"][str(dt_ps)] = {
                "native_voxel": native_summary,
                "component_2x1x1": coarse_summary,
                "component_vs_native_field_error": _compare_component_to_native(
                    native_graph, native_fields, coarse_graph, coarse_fields, roi, voxel_volume
                ),
            }
        result["regions"][name] = region
    duration_fs = final_time_ps * 1000.0
    duration_label = f"{duration_fs:g}fs".replace(".", "p")
    output = PHASE_DIR / "nested_grid_diagnostics" / f"qsgs_component_ttm_{duration_label}_benchmark.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
