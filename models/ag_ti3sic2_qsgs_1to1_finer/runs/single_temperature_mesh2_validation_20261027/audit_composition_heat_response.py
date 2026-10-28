#!/usr/bin/env python3
"""Read-only, matched-0.2-ps source/topology audit for Ag 41.7 vs 70 vol%.

This script initializes the unchanged Case A operators and inspects saved
geometries/checkpoints. It never advances time or changes checkpoint files.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT / "phase_field_qsgs"))

from qsgs_enthalpy_solver import QSGSEnthalpyGrid, ThermalState  # noqa: E402


CASES = {
    "ag_41p6654_volpct": {
        "geometry": PROJECT / "qsgs_morphology.npz",
        "checkpoint": PROJECT / "runs/single_temperature_restart_track_20261027/checkpoints/caseA_heat_0p200ps.npz",
    },
    "ag_69p9996_volpct": {
        "geometry": HERE / "qsgs_morphology_ag70vol.npz",
        "checkpoint": HERE / "ag70vol_melt_window/1x_0p2ps/checkpoints/caseA_ag70_0p2ps_0p2ps.npz",
    },
}


def analyze(label: str, geometry: Path, checkpoint: Path) -> dict[str, object]:
    model = QSGSEnthalpyGrid(
        None, "A", native_geometry_npz=geometry, source_geometry_npz=geometry
    )
    with np.load(geometry, allow_pickle=False) as archive:
        ag = ~np.asarray(archive["skeleton"], dtype=bool)
        size = np.asarray(archive["size_A"], dtype=float)
    with np.load(checkpoint, allow_pickle=False) as archive:
        temperature = np.asarray(archive["temperature_K"], dtype=float)
        liquid = np.asarray(archive["ag_liquid_fraction"], dtype=float)
        metadata = json.loads(str(archive["metadata_json"].item()))

    if temperature.shape != ag.shape or liquid.shape != ag.shape:
        raise ValueError(f"{label}: checkpoint fields do not match geometry")
    if not np.all(np.isfinite(temperature)) or not np.all(np.isfinite(liquid)):
        raise ValueError(f"{label}: non-finite saved field")

    shape = ag.shape
    x = (np.arange(shape[0]) + 0.5) * model.spacing_A[0]
    y = (np.arange(shape[1]) + 0.5) * model.spacing_A[1]
    z = (np.arange(shape[2]) + 0.5) * model.spacing_A[2]
    cx, cy = model.source_center_xy_A
    radius = np.broadcast_to(
        np.sqrt((x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2),
        shape,
    )
    depth = np.broadcast_to(size[2] - z[None, None, :], shape)
    source_power = model.source_cell_power_W
    total_source_power = float(source_power.sum())
    if total_source_power <= 0.0:
        raise ValueError(f"{label}: integrated source power is not positive")

    face_pairs: list[int] = []
    degree = np.zeros(shape, dtype=np.int8)
    for axis in (0, 1):
        face_pairs.append(int(np.count_nonzero(ag & np.roll(ag, -1, axis=axis))))
        degree += np.roll(ag, 1, axis=axis)
        degree += np.roll(ag, -1, axis=axis)
    face_pairs.append(int(np.count_nonzero(ag[:, :, :-1] & ag[:, :, 1:])))
    degree[:, :, 1:] += ag[:, :, :-1]
    degree[:, :, :-1] += ag[:, :, 1:]

    total_ag_capacity = float(
        model._solid_capacity_at_300_J_m3_K.sum() * model.cell_volume_m3
    )
    region = (radius < 30.0) & (depth < 20.0)
    source_region_power = float(source_power[region].sum())
    source_region_capacity = float(
        model._solid_capacity_at_300_J_m3_K[region].sum() * model.cell_volume_m3
    )
    initial_dtemp_dt = np.divide(
        model.source_W_m3,
        model._solid_capacity_at_300_J_m3_K,
        out=np.zeros(shape, dtype=float),
        where=model._solid_capacity_at_300_J_m3_K > 0.0,
    )
    hottest = np.unravel_index(int(np.argmax(temperature)), shape)
    state = ThermalState(temperature, liquid)
    conductivity = model.conductivity_W_m_K(state)
    diffusion = model.diffusion_power_density_W_m3(temperature, conductivity)
    local_source = float(source_power[hottest])
    local_diffusion = float(diffusion[hottest] * model.cell_volume_m3)

    neighbors: list[dict[str, object]] = []
    i, j, k = (int(value) for value in hottest)
    for axis, step in ((0, -1), (0, 1), (1, -1), (1, 1), (2, -1), (2, 1)):
        index = [i, j, k]
        index[axis] += step
        if axis < 2:
            index[axis] %= shape[axis]
        elif index[axis] < 0 or index[axis] >= shape[axis]:
            continue
        neighbor = tuple(index)
        neighbors.append(
            {
                "index_xyz": [int(value) for value in neighbor],
                "is_Ag": bool(ag[neighbor]),
                "temperature_K": float(temperature[neighbor]),
            }
        )

    ag_model = model.ag
    alpha_300 = ag_model.conductivity_W_m_K(300.0) / (
        ag_model.reference_density_kg_m3 * ag_model.cp_solid_J_kg_K(300.0)
    )
    return {
        "label": label,
        "geometry_path": str(geometry.relative_to(PROJECT)),
        "geometry_sha256": model._sha256_file(geometry),
        "checkpoint_path": str(checkpoint.relative_to(PROJECT)),
        "checkpoint_time_ps": float(metadata["time_ps"]),
        "checkpoint_identity_sha256": metadata["identity_sha256"],
        "grid_shape_xyz": list(shape),
        "domain_size_A": size.tolist(),
        "ag_volume_fraction": float(ag.mean()),
        "case": "A; zero scaffold conductivity and heat capacity; periodic x/y; insulated z",
        "source": {
            "peak_flux_W_m2": model.peak_flux_W_m2,
            "gaussian_radius_A": model.gaussian_radius_A,
            "absorption_depth_A": model.absorption_depth_A,
            "center_xy_A": list(model.source_center_xy_A),
            "integration_method": model.source_integration_method,
            "actual_power_W": total_source_power,
            "injected_energy_at_checkpoint_J": float(metadata["injected_energy_J"]),
            "power_weighted_mean_radius_A": float((source_power * radius).sum() / total_source_power),
            "power_weighted_mean_depth_A": float((source_power * depth).sum() / total_source_power),
            "power_fraction_radius_lt_30_A": float(source_power[radius < 30.0].sum() / total_source_power),
            "power_fraction_depth_lt_20_A": float(source_power[depth < 20.0].sum() / total_source_power),
            "power_fraction_radius_lt_30_and_depth_lt_20_A": float(source_power[region].sum() / total_source_power),
            "maximum_cell_power_W": float(source_power.max()),
            "maximum_source_density_W_m3": float(model.source_W_m3.max()),
            "source_weighted_initial_dTdt_K_per_ps": float(
                (source_power * initial_dtemp_dt).sum() / total_source_power * 1.0e-12
            ),
            "maximum_initial_dTdt_K_per_ps": float(initial_dtemp_dt.max() * 1.0e-12),
        },
        "case_A_sensible_heat_capacity_at_300_J_K": total_ag_capacity,
        "source_power_over_total_Ag_capacity_K_per_ps": total_source_power / total_ag_capacity * 1.0e-12,
        "saved_state_at_0p2ps": {
            "maximum_temperature_K": float(temperature.max()),
            "mean_Ag_temperature_K": float(temperature[ag].mean()),
            "Ag_liquid_volume_fraction": float((liquid * ag).sum() / ag.sum()),
            "maximum_liquid_fraction": float(liquid.max()),
            "hottest_voxel_index_xyz": [int(value) for value in hottest],
            "hottest_voxel_position_A": ((np.asarray(hottest) + 0.5) * model.spacing_A).tolist(),
            "hottest_voxel_radius_from_source_A": float(radius[hottest]),
            "hottest_voxel_depth_from_top_A": float(depth[hottest]),
            "hottest_voxel_Ag_face_neighbors": int(degree[hottest]),
            "hottest_voxel_source_power_W": local_source,
            "hottest_voxel_instantaneous_diffusion_power_W": local_diffusion,
            "hottest_voxel_diffusion_to_source_ratio": local_diffusion / local_source if local_source else None,
            "hottest_voxel_temperature_K": float(temperature[hottest]),
            "hottest_voxel_liquid_fraction": float(liquid[hottest]),
            "hottest_voxel_neighbors": neighbors,
        },
        "Ag_face_network_proxy": {
            "mean_face_neighbors_per_Ag_voxel": float(degree[ag].mean()),
            "source_power_weighted_face_neighbors": float((source_power * degree).sum() / total_source_power),
            "Ag_Ag_face_pairs_xyz_periodic_xy_open_z": face_pairs,
            "Ag_fraction_in_source_region_radius_lt_30_A_depth_lt_20_A": float(ag[region].mean()),
            "interpretation": "Voxel adjacency proxy only; not a homogenized/effective conductivity calculation.",
        },
        "source_region_radius_lt_30_A_depth_lt_20_A": {
            "source_power_W": source_region_power,
            "candidate_Ag_capacity_at_300_J_K": source_region_capacity,
            "source_power_over_region_capacity_K_per_ps": (
                source_region_power / source_region_capacity * 1.0e-12
                if source_region_capacity > 0.0 else None
            ),
        },
        "candidate_Ag_thermal_scale": {
            "density_kg_m3": ag_model.reference_density_kg_m3,
            "cp_at_300_J_kg_K": ag_model.cp_solid_J_kg_K(300.0),
            "k_at_300_W_m_K": ag_model.conductivity_W_m_K(300.0),
            "alpha_at_300_m2_s": alpha_300,
            "sqrt_2_alpha_t_at_0p2ps_A": float(np.sqrt(2.0 * alpha_300 * 0.2e-12) * 1.0e10),
            "Tm_K": ag_model.melting_temperature_K,
            "status": "candidate properties; not source-verified",
        },
        "transient_steps_performed": 0,
    }


def main() -> None:
    result = {
        "status": "READ_ONLY_MATCHED_TIME_COMPOSITION_RESPONSE_AUDIT",
        "question": "Why the 41.6654 vol% and 69.9996 vol% geometries differ at the same 0.2 ps input",
        "interpretation_limit": "Different QSGS realizations/compositions are not a controlled micro-versus-nano or composition-only experiment.",
        "cases": {
            name: analyze(name, paths["geometry"], paths["checkpoint"])
            for name, paths in CASES.items()
        },
    }
    output = HERE / "ag41_vs_ag70_static_diagnostics.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
