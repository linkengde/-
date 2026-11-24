#!/usr/bin/env python3
"""Diagnose v6 restart boundary leakage and conservative Ag/W CIC remapping."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

try:
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import maximum_flow
except ImportError:  # the production TTM mapper itself only needs NumPy
    coo_matrix = maximum_flow = None

from agw_ttm_reference import (
    AG_ATOMIC_VOLUME_A3,
    W_ATOMIC_VOLUME_A3,
    AgWMappingError,
    build_ag_grid_fields,
    build_agw_grid_fields_joint,
    cloud_in_cell_atom_map,
    gaussian_ag_source_power,
    map_agw_newton_layer,
    _raw_cic_volume_power,
)

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
V6 = MODEL / "interface_reconstruction_v6_exploratory"
OUTPUT = HERE / "agw_ttm_mapping_v7_report.json"
TRANSFER_LIMIT_A = 20.0


def read_dump(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("ITEM: ATOMS"):
                columns = line.split()[2:]
                break
        else:
            raise ValueError(f"No ATOMS section in {path}")
        table = np.loadtxt(stream)
    index = {name: columns.index(name) for name in ("id", "type", "x", "y", "z")}
    order = np.argsort(table[:, index["id"]].astype(np.int64))
    atom_ids = table[order, index["id"]].astype(np.int64)
    types = table[order, index["type"]].astype(np.int32)
    xyz = table[order][:, [index["x"], index["y"], index["z"]]].astype(float)
    return atom_ids, types, xyz


def read_lammps_data_atoms(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read sorted ID/type/xyz from an atomic-style LAMMPS data file."""
    rows: list[tuple[int, int, float, float, float]] = []
    in_atoms = False
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not in_atoms:
                if stripped.startswith("Atoms"):
                    in_atoms = True
                continue
            if not stripped:
                if rows:
                    break
                continue
            fields = stripped.split()
            try:
                rows.append((int(fields[0]), int(fields[1]), float(fields[2]), float(fields[3]), float(fields[4])))
            except (ValueError, IndexError):
                if rows:
                    break
                raise ValueError(f"Could not parse an Atoms record in {path}: {stripped}")
    if not rows:
        raise ValueError(f"No Atoms records found in {path}")
    table = np.asarray(rows, dtype=float)
    order = np.argsort(table[:, 0].astype(np.int64))
    return (
        table[order, 0].astype(np.int64),
        table[order, 1].astype(np.int32),
        table[order, 2:5].astype(float),
    )


def count_by_phase(mask: np.ndarray, types: np.ndarray) -> dict:
    return {
        "total": int(mask.sum()),
        "Ag": int(np.count_nonzero(mask & (types == 2))),
        "W": int(np.count_nonzero(mask & (types == 1))),
    }


def clipped_edge_weight(xyz: np.ndarray, mask: np.ndarray, shape: tuple[int, int, int], spacing: tuple[float, float, float]) -> dict:
    selected = xyz[mask] - np.array([0.0, 0.0, 40.0])
    mapping = cloud_in_cell_atom_map(selected, shape, spacing)
    iz = mapping.cell_ids % shape[2]
    lower = np.sum(mapping.weights * (iz == 0), axis=1)
    upper = np.sum(mapping.weights * (iz == shape[2] - 1), axis=1)
    return {
        "atoms": int(mask.sum()),
        "atoms_with_weight_in_lower_z_cell": int(np.count_nonzero(lower > 1e-12)),
        "atoms_with_weight_in_upper_z_cell": int(np.count_nonzero(upper > 1e-12)),
        "sum_weight_lower_z_cell": float(lower.sum()),
        "sum_weight_upper_z_cell": float(upper.sum()),
        "max_single_atom_weight_lower": float(lower.max(initial=0.0)),
        "max_single_atom_weight_upper": float(upper.max(initial=0.0)),
    }


def run_synthetic_joint_check() -> dict:
    ag = np.repeat([[5.0, 5.0, 5.0]], 40, axis=0)
    w = np.repeat([[5.0, 5.0, 5.0]], 30, axis=0)
    q = np.full(len(ag), 0.1)
    mapped = build_agw_grid_fields_joint(ag, w, q, (4, 4, 4), (10.0, 10.0, 10.0), 20.0)
    d = mapped["diagnostics"]
    expected_ag = len(ag) * AG_ATOMIC_VOLUME_A3
    expected_w = len(w) * W_ATOMIC_VOLUME_A3
    expected_source = expected_ag * 0.1
    checks = {
        "combined_fill_at_most_one": float((mapped["fraction_ag"] + mapped["fraction_w"]).max()) <= 1.0 + 1e-10,
        "Ag_volume_conserved": abs(float(mapped["fraction_ag"].sum() * 1000.0) - expected_ag) <= 1e-8,
        "W_volume_conserved": abs(float(mapped["fraction_w"].sum() * 1000.0) - expected_w) <= 1e-8,
        "Ag_source_conserved": abs(float(mapped["source_power_ag_eV_ps"].sum()) - expected_source) <= 1e-10,
        "Ag_atom_map_rows_normalized": d["Ag_atom_map_max_row_sum_error"] <= 1e-12,
        "W_atom_map_rows_normalized": d["W_atom_map_max_row_sum_error"] <= 1e-12,
        "transfer_within_limit": d["transfer_max_observed_distance_A"] <= 20.0 + 1e-12,
    }
    if not all(checks.values()):
        raise RuntimeError(f"Synthetic joint-mapping check failed: {checks}")
    return {"status": "PASS", "checks": checks, "diagnostics": d}


def run_outside_filter_check() -> dict:
    xyz = np.array([
        [100.0, 100.0, 150.0],  # Ag inside
        [100.0, 100.0, 201.0],  # Ag outside the TTM material plane
        [100.0, 100.0, 160.0],  # W inside
    ])
    types = np.array([2, 2, 1], dtype=int)
    mapped = map_agw_newton_layer(
        xyz, types, (20, 20, 16), (10.0, 10.0, 10.0), newton_mask=np.ones(3, dtype=bool)
    )
    expected_source = float(gaussian_ag_source_power(xyz[[0]])[0] * AG_ATOMIC_VOLUME_A3)
    checks = {
        "outside_atom_excluded": mapped["excluded_outside_ttm_region"]["above_200A"] == 1,
        "only_one_Ag_mapped": len(mapped["atom_indices_ag"]) == 1,
        "Ag_source_has_no_above_surface_contribution": abs(mapped["source_power_on_state_eV_ps"] - expected_source) < max(1e-10, expected_source * 1e-12),
        "mapped_fill_at_most_one": mapped["max_total_fill_fraction"] <= 1.0 + 1e-10,
    }
    if not all(checks.values()):
        raise RuntimeError(f"Synthetic outside-z filter check failed: {checks}")
    return {"status": "PASS", "checks": checks, "excluded": mapped["excluded_outside_ttm_region"]}


def local_capacity_flow(raw_total: np.ndarray, spacing: tuple[float, float, float], max_distance_A: float) -> dict:
    """Check whether all overflow can reach local free volume, independent of donor order."""
    if coo_matrix is None or maximum_flow is None:
        return {"status": "SCIPY_UNAVAILABLE", "diagnostic_only": True}
    shape = tuple(int(v) for v in raw_total.shape)
    cell_volume = float(np.prod(spacing))
    excess = np.maximum(raw_total.ravel() - cell_volume, 0.0)
    free = np.maximum(cell_volume - raw_total.ravel(), 0.0)
    donors = np.flatnonzero(excess > cell_volume * 1e-12)
    receivers = np.flatnonzero(free > cell_volume * 1e-12)
    ncell = int(np.prod(shape))
    donor_node = np.full(ncell, -1, dtype=int)
    receiver_node = np.full(ncell, -1, dtype=int)
    source = 0
    donor_node[donors] = np.arange(1, len(donors) + 1)
    receiver_offset = 1 + len(donors)
    receiver_node[receivers] = receiver_offset + np.arange(len(receivers))
    sink = receiver_offset + len(receivers)
    scale = 1000  # integer max-flow units correspond to 0.001 A^3
    demand_units = np.rint(excess[donors] * scale).astype(np.int64)
    free_units = np.floor(free[receivers] * scale).astype(np.int64)
    total_demand_units = int(demand_units.sum())
    rows: list[int] = []
    cols: list[int] = []
    caps: list[int] = []
    supply_limit = total_demand_units
    nx, ny, nz = shape
    sx, sy, sz = spacing
    rx = min(nx // 2, int(np.ceil(max_distance_A / sx)))
    ry = min(ny // 2, int(np.ceil(max_distance_A / sy)))
    rz = min(nz - 1, int(np.ceil(max_distance_A / sz)))
    for donor, units in zip(donors, demand_units):
        dn = int(donor_node[donor])
        rows.append(source); cols.append(dn); caps.append(int(units))
        ix = donor // (ny * nz)
        iy = (donor // nz) % ny
        iz = donor % nz
        seen: set[int] = set()
        for dx in range(-rx, rx + 1):
            tx = (ix + dx) % nx
            ddx = min(abs(dx), nx - abs(dx)) * sx
            for dy in range(-ry, ry + 1):
                ty = (iy + dy) % ny
                ddy = min(abs(dy), ny - abs(dy)) * sy
                for dz in range(-rz, rz + 1):
                    tz = iz + dz
                    if tz < 0 or tz >= nz or (dx == 0 and dy == 0 and dz == 0):
                        continue
                    if float(np.sqrt(ddx * ddx + ddy * ddy + (dz * sz) ** 2)) > max_distance_A + 1e-12:
                        continue
                    target = int(np.ravel_multi_index((tx, ty, tz), shape))
                    rn = int(receiver_node[target])
                    if rn < 0 or target in seen:
                        continue
                    seen.add(target)
                    rows.append(dn); cols.append(rn); caps.append(supply_limit)
    for receiver, units in zip(receivers, free_units):
        rows.append(int(receiver_node[receiver])); cols.append(sink); caps.append(int(units))
    graph = coo_matrix(
        (np.asarray(caps, dtype=np.int64), (rows, cols)), shape=(sink + 1, sink + 1)
    ).tocsr()
    solved = maximum_flow(graph, source, sink, method="dinic")
    flow_units = int(solved.flow_value)
    return {
        "status": "PASS" if flow_units >= total_demand_units else "NO_LOCAL_CAPACITY_ASSIGNMENT",
        "diagnostic_only": True,
        "algorithm": "SciPy maximum flow on donor-to-recipient edges within the fixed transfer bound",
        "capacity_resolution_A3": 1.0 / scale,
        "max_transfer_distance_A": float(max_distance_A),
        "donor_cells": int(len(donors)),
        "recipient_cells": int(len(receivers)),
        "candidate_edges": int(graph.nnz),
        "overflow_demand_A3": float(excess.sum()),
        "global_free_capacity_A3": float(free.sum()),
        "locally_assignable_volume_A3": flow_units / scale,
        "locally_unassignable_overflow_A3": max(0.0, float(excess.sum()) - flow_units / scale),
    }


def main() -> None:
    manifest_path = V6 / "exploratory_ttm_state_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    data_path = HERE.parents[4] / manifest["data_file"]
    template_ids, template_types, template_xyz = read_lammps_data_atoms(data_path)
    start_ids, types, xyz_start = read_dump(V6 / "nvt_fixedids_start.dump.gz")
    end_ids, end_types, xyz_end = read_dump(V6 / "nvt_fixedids_end.dump.gz")
    final_ids, final_types, xyz_final = read_dump(V6 / "nve_fixedids_after_nvt_final.dump.gz")
    if not (np.array_equal(start_ids, end_ids) and np.array_equal(start_ids, final_ids)
            and np.array_equal(types, end_types) and np.array_equal(types, final_types)
            and np.array_equal(start_ids, template_ids) and np.array_equal(types, template_types)):
        raise RuntimeError("The fixed-ID continuation snapshots do not contain the same atom IDs/types")
    # The runner creates its LAMMPS groups from the original data file before
    # read_dump replaces coordinates. Reconstruct that static Newton group
    # from the same template, rather than inferring membership from a later z.
    original_newton = (template_xyz[:, 2] >= 40.0 - 1e-9) & (template_xyz[:, 2] <= 300.0 + 1e-9)
    crossings: dict[str, dict] = {}
    for label, coords in (("NVT_end", xyz_end), ("NVE_end", xyz_final)):
        below = original_newton & (coords[:, 2] < 40.0 - 1e-9)
        above = original_newton & (coords[:, 2] > 200.0 + 1e-9)
        outside = below | above
        template_inside = original_newton & (template_xyz[:, 2] >= 40.0 - 1e-9) & (template_xyz[:, 2] <= 200.0 + 1e-9)
        template_outside = original_newton & ~template_inside
        newly_outside = template_inside & outside
        remained_outside = template_outside & outside
        returned_inside = template_outside & ~outside
        crossings[label] = {
            "original_newton_candidates": count_by_phase(original_newton, types),
            "template_data_already_outside_40_to_200A": count_by_phase(template_outside, types),
            "below_40A": count_by_phase(below, types),
            "above_200A": count_by_phase(above, types),
            "outside_40_to_200A": count_by_phase(outside, types),
            "newly_outside_relative_to_template_data": count_by_phase(newly_outside, types),
            "template_outside_atoms_still_outside": count_by_phase(remained_outside, types),
            "template_outside_atoms_returned_inside": count_by_phase(returned_inside, types),
            "legacy_CIC_clipped_edge_contribution": {
                "below_40A": clipped_edge_weight(coords, below, (20, 20, 16), (10.0, 10.0, 10.0)),
                "above_200A": clipped_edge_weight(coords, above, (20, 20, 16), (10.0, 10.0, 10.0)),
            },
            "Ag_source_power_from_above_surface_if_mask_not_filtered_eV_ps": float(
                np.sum(gaussian_ag_source_power(coords[above & (types == 2)], allow_above_surface=True)
                       * AG_ATOMIC_VOLUME_A3)
            ),
        }

    newton_region = original_newton
    final_inside = original_newton & (xyz_final[:, 2] >= 40.0 - 1e-9) & (xyz_final[:, 2] <= 200.0 + 1e-9)
    active_ids = np.flatnonzero(final_inside)
    xyz_active = xyz_final[active_ids]
    type_active = types[active_ids]
    xyz_ag_abs = xyz_active[type_active == 2]
    xyz_w_abs = xyz_active[type_active == 1]
    q_ag = gaussian_ag_source_power(xyz_ag_abs)
    mask_diagnostics = {
        "frozen_NVE_restart_atoms": int(len(xyz_final)),
        "static_Newton_group_from_template_data_z40_to_300": count_by_phase(newton_region, types),
        "mapped_in_range_atoms": count_by_phase(final_inside, types),
        "static_Newton_group_atoms_outside_TTM_mapping_range_at_restart": count_by_phase(newton_region & ~final_inside, types),
        "static_Newton_group_atoms_above_200A_at_restart": count_by_phase(newton_region & (xyz_final[:, 2] > 200.0 + 1e-9), types),
        "template_Newton_group_atoms_already_outside_TTM_mapping_range": count_by_phase(newton_region & ~((template_xyz[:, 2] >= 40.0 - 1e-9) & (template_xyz[:, 2] <= 200.0 + 1e-9)), types),
        "upper_surface_z_A": 200.0,
        "material_mapping_z_range_A": [40.0, 200.0],
        "outside_atoms_are_filtered_before_CIC": True,
    }

    mesh_reports = []
    for h in (5.0, 10.0, 20.0):
        spacing = (h, h, h)
        shape = (int(round(200.0 / h)), int(round(200.0 / h)), int(round(160.0 / h)))
        origin = np.array([0.0, 0.0, 40.0])
        rel_ag, rel_w = xyz_ag_abs - origin, xyz_w_abs - origin
        _, raw_volume_ag, raw_source_power = _raw_cic_volume_power(
            rel_ag, q_ag, shape, spacing, AG_ATOMIC_VOLUME_A3
        )
        _, raw_volume_w, _ = _raw_cic_volume_power(
            rel_w, np.zeros(len(rel_w)), shape, spacing, W_ATOMIC_VOLUME_A3
        )
        cell_volume = float(np.prod(spacing))
        raw_combined = (raw_volume_ag + raw_volume_w) / cell_volume
        flow_diagnostic = local_capacity_flow(raw_volume_ag + raw_volume_w, spacing, TRANSFER_LIMIT_A)
        ag_single = build_ag_grid_fields(rel_ag, q_ag, shape, spacing, AG_ATOMIC_VOLUME_A3)
        w_single = build_ag_grid_fields(rel_w, np.zeros(len(rel_w)), shape, spacing, W_ATOMIC_VOLUME_A3)
        raw_ag = ag_single.raw_max_fill_fraction
        raw_w = w_single.raw_max_fill_fraction
        separate_combined = ag_single.fraction_ag + w_single.fraction_ag
        expected_source_power = float(np.dot(np.full(len(q_ag), AG_ATOMIC_VOLUME_A3), q_ag))
        try:
            joint = build_agw_grid_fields_joint(rel_ag, rel_w, q_ag, shape, spacing, TRANSFER_LIMIT_A)
            mesh_reports.append({
                "spacing_A": list(spacing),
                "shape": list(shape),
                "raw_CIC": {
                    "Ag_max_fill": raw_ag,
                    "W_max_fill": raw_w,
                    "combined_max_fill": float(raw_combined.max()),
                    "AgW_raw_overfilled_cells": int(np.count_nonzero(raw_combined > 1.0 + 1e-10)),
                    "Ag_raw_overfilled_cells": ag_single.raw_overfilled_cells,
                    "W_raw_overfilled_cells": w_single.raw_overfilled_cells,
                    "Ag_source_power_eV_ps": float(raw_source_power.sum()),
                    "expected_Ag_source_power_eV_ps": expected_source_power,
                },
                "local_capacity_flow": flow_diagnostic,
                "independent_phase_redistribution": {
                    "Ag_volume_moved_A3": ag_single.redistributed_ag_volume_A3,
                    "W_volume_moved_A3": w_single.redistributed_ag_volume_A3,
                    "combined_max_fill": float(separate_combined.max()),
                    "combined_cells_over_capacity": int(np.count_nonzero(separate_combined > 1.0 + 1e-10)),
                    "Ag_source_power_eV_ps": float(ag_single.source_power_eV_ps.sum()),
                },
                "joint_mapping": {"status": "PASS", "diagnostics": joint["diagnostics"]},
            })
        except AgWMappingError as exc:
            mesh_reports.append({
                "spacing_A": list(spacing),
                "shape": list(shape),
                "raw_CIC": {
                    "Ag_max_fill": raw_ag,
                    "W_max_fill": raw_w,
                    "combined_max_fill": float(raw_combined.max()),
                    "AgW_raw_overfilled_cells": int(np.count_nonzero(raw_combined > 1.0 + 1e-10)),
                    "Ag_raw_overfilled_cells": ag_single.raw_overfilled_cells,
                    "W_raw_overfilled_cells": w_single.raw_overfilled_cells,
                    "Ag_source_power_eV_ps": float(raw_source_power.sum()),
                    "expected_Ag_source_power_eV_ps": expected_source_power,
                },
                "local_capacity_flow": flow_diagnostic,
                "independent_phase_redistribution": {
                    "Ag_volume_moved_A3": ag_single.redistributed_ag_volume_A3,
                    "W_volume_moved_A3": w_single.redistributed_ag_volume_A3,
                    "combined_max_fill": float(separate_combined.max()),
                    "combined_cells_over_capacity": int(np.count_nonzero(separate_combined > 1.0 + 1e-10)),
                    "Ag_source_power_eV_ps": float(ag_single.source_power_eV_ps.sum()),
                },
                "joint_mapping": {
                    "status": "REJECTED_WITHIN_TRANSFER_LIMIT",
                    "max_transfer_distance_A": TRANSFER_LIMIT_A,
                    "reason": str(exc),
                },
            })

    result = {
        "title": "Ag-W v7 frozen-restart TTM mapping diagnosis",
        "status": "DIAGNOSIS_COMPLETE_JOINT_MAPPING_REJECTED_AT_20A_LIMIT",
        "source_restart": "interface_reconstruction_v6_exploratory/nve_fixedids_after_nvt_final.dump.gz",
        "template_data_file": manifest["data_file"],
        "template_data_sha256": manifest["data_sha256"],
        "source_restart_sha256": hashlib.sha256(
            (V6 / "nve_fixedids_after_nvt_final.dump.gz").read_bytes()
        ).hexdigest(),
        "grid_origin_A": [0.0, 0.0, 40.0],
        "transfer_limit_A": TRANSFER_LIMIT_A,
        "transfer_limit_basis": "Exploratory spatial-locality screen selected for this audit; not a measured physical length or a calibrated mapping parameter.",
        "out_of_range_crossing_diagnostics": crossings,
        "frozen_restart_mask_diagnostics": mask_diagnostics,
        "meshes": mesh_reports,
        "synthetic_joint_transfer_check": run_synthetic_joint_check(),
        "synthetic_outside_z_filter_check": run_outside_filter_check(),
        "decision": "Do not run source/no-source LAMMPS/TTM pair until a 5/10/20 A joint map passes the documented transfer-distance bound and conservation checks.",
        "limitations": [
            "The actual v6 geometry still fails the 20 A local-transfer limit; this does not relax the TTM phase-capacity guard.",
            "The v5 EFS static audit remains STATIC_EFS_FAIL.",
            "The joint redistribution is a conservative CIC representation, not a physical movement of atoms or an atomistic interface repair.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "meshes": result["meshes"],
                      "crossings": result["out_of_range_crossing_diagnostics"],
                      "synthetic_joint_transfer_check": result["synthetic_joint_transfer_check"]["status"],
                      "synthetic_outside_z_filter_check": result["synthetic_outside_z_filter_check"]["status"]}, indent=2))


if __name__ == "__main__":
    main()
