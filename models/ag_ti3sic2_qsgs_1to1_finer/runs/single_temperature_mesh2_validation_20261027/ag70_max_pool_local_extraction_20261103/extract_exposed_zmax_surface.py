#!/usr/bin/env python3
"""Extract the material and melt state exposed at the physical z_max face.

This is a surface-geometry audit only. It does not assign gas pressure,
temperature flux, contact angle, evaporation, or fluid motion.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np


ROOT = Path(__file__).resolve().parent
STUDY = ROOT.parent
SELECTION_PATH = ROOT / "TARGET_GEOMETRY_SELECTION.json"
MASKS_PATH = ROOT / "maximum_pool_full_domain_masks.npz"
RUN_NAME = "ag70vol_1x_regrown_64x64x50_caseB_4x_composite_sc2300K_20261102"
TOL = 1e-12


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def transition_edges(a: np.ndarray, b: np.ndarray, axis: int) -> int:
    """Count periodic in-plane cell edges shared by disjoint classes a and b."""
    b_next = np.roll(b, -1, axis=axis)
    a_next = np.roll(a, -1, axis=axis)
    return int(np.count_nonzero((a & b_next) | (b & a_next)))


def main() -> None:
    selection = json.loads(SELECTION_PATH.read_text(encoding="utf-8"))
    geometry_path = STUDY / RUN_NAME / "geometry/qsgs_morphology_ag70vol_regrown_64x64x50.npz"
    checkpoint_path = STUDY / RUN_NAME / "checkpoints/heating_complete_2300K.npz"
    geometry_hash = sha256(geometry_path)
    checkpoint_hash = sha256(checkpoint_path)
    if geometry_hash != selection["geometry"]["sha256"]:
        raise ValueError("QSGS geometry hash does not match the frozen target selection")
    if checkpoint_hash != selection["maximum_pool_checkpoint"]["sha256"]:
        raise ValueError("Checkpoint hash does not match the frozen target selection")

    with np.load(geometry_path) as geometry, np.load(checkpoint_path) as checkpoint, np.load(MASKS_PATH) as masks:
        skeleton = geometry["skeleton"].astype(bool)
        liquid_fraction = checkpoint["ag_liquid_fraction"]
        temperature = checkpoint["temperature_K"]
        pool = masks["maximum_pool_track_mask"].astype(bool)
        if skeleton.shape != (64, 64, 50) or liquid_fraction.shape != skeleton.shape or pool.shape != skeleton.shape:
            raise ValueError("Unexpected source array dimensions")
        size_A = np.asarray(geometry["size_A"], dtype=float)
        expected_size_nm = np.asarray(selection["geometry"]["domain_size_nm_xyz"], dtype=float)
        if not np.allclose(size_A / 10.0, expected_size_nm, rtol=0.0, atol=1e-10):
            raise ValueError("Geometry dimensions do not match the frozen target selection")
        dx_nm, dy_nm, dz_nm = size_A / np.asarray(skeleton.shape, dtype=float) / 10.0
        exposed_skeleton = skeleton[:, :, -1]
        exposed_ag = ~exposed_skeleton
        exposed_fraction = liquid_fraction[:, :, -1]
        exposed_pool = pool[:, :, -1]
        exposed_temperature = temperature[:, :, -1]

        solid_ag = exposed_ag & (exposed_fraction <= TOL)
        partial_ag = exposed_ag & (exposed_fraction > TOL) & (exposed_fraction < 1.0 - TOL)
        partial_below_half = partial_ag & (exposed_fraction < 0.5)
        partial_at_or_above_half = partial_ag & (exposed_fraction >= 0.5)
        fully_liquid_ag = exposed_ag & (exposed_fraction >= 1.0 - TOL)
        liquid_dominant_ag = exposed_ag & (exposed_fraction >= 0.5)

        phase = np.zeros(exposed_skeleton.shape, dtype=np.uint8)
        phase[solid_ag] = 1
        phase[partial_below_half] = 2
        phase[partial_at_or_above_half] = 3
        phase[fully_liquid_ag] = 4

        face_area_nm2 = dx_nm * dy_nm
        voxel_face_area_nm = {"x_neighbor_edge": dy_nm, "y_neighbor_edge": dx_nm}
        cell_count = int(phase.size)

        def count_area(mask: np.ndarray) -> dict:
            n = int(mask.sum())
            return {"cells": n, "area_nm2": n * face_area_nm2, "area_fraction_of_zmax": n / cell_count}

        interface_lengths = {}
        contact_edge_rows = []
        for label, liquid_class in (
            ("any_melt_fraction_gt_0", exposed_ag & (exposed_fraction > TOL)),
            ("liquid_dominant_f_l_ge_0p5", liquid_dominant_ag),
            ("fully_liquid_f_l_ge_1_minus_tol", fully_liquid_ag),
        ):
            nx = transition_edges(liquid_class, exposed_skeleton, axis=0)
            ny = transition_edges(liquid_class, exposed_skeleton, axis=1)
            interface_lengths[label] = {
                "liquid_skeleton_shared_edges_x_direction": nx,
                "liquid_skeleton_shared_edges_y_direction": ny,
                "grid_edge_contact_length_nm": nx * voxel_face_area_nm["x_neighbor_edge"] + ny * voxel_face_area_nm["y_neighbor_edge"],
            }
            for axis, edge_count in ((0, nx), (1, ny)):
                n_axis = phase.shape[axis]
                n_other = phase.shape[1 - axis]
                emitted = 0
                for i in range(n_axis):
                    i_next = (i + 1) % n_axis
                    for j in range(n_other):
                        a_idx = (i, j) if axis == 0 else (j, i)
                        b_idx = (i_next, j) if axis == 0 else (j, i_next)
                        if liquid_class[a_idx] and exposed_skeleton[b_idx]:
                            liquid_idx = a_idx
                        elif exposed_skeleton[a_idx] and liquid_class[b_idx]:
                            liquid_idx = b_idx
                        else:
                            continue
                        emitted += 1
                        if axis == 0:
                            x_edge = 0.0 if i == n_axis - 1 else (i + 1) * dx_nm
                            x0 = x1 = x_edge
                            y0, y1 = j * dy_nm, (j + 1) * dy_nm
                            seam = i == n_axis - 1
                            normal_axis = "x"
                        else:
                            y_edge = 0.0 if i == n_axis - 1 else (i + 1) * dy_nm
                            x0, x1 = j * dx_nm, (j + 1) * dx_nm
                            y0 = y1 = y_edge
                            seam = i == n_axis - 1
                            normal_axis = "y"
                        contact_edge_rows.append([
                            label, normal_axis, x0, y0, x1, y1,
                            bool(seam), float(exposed_fraction[liquid_idx]),
                        ])
                if emitted != edge_count:
                    raise AssertionError(f"contact-edge export mismatch for {label} along {normal_axis}")

        metrics = {
            "status": "EXPOSED_ZMAX_SURFACE_GEOMETRY_ONLY",
            "geometry_sha256": geometry_hash,
            "checkpoint_sha256": checkpoint_hash,
            "extraction_masks_sha256": sha256(MASKS_PATH),
            "checkpoint_time_ps": selection["maximum_pool_checkpoint"]["time_ps"],
            "grid_xyz": [64, 64, 50],
            "zmax_index": 49,
            "grid_spacing_nm_xyz": [dx_nm, dy_nm, dz_nm],
            "x_y_boundary_for_edge_counts": "periodic, matching the source QSGS geometry",
            "surface_area_per_cell_nm2": face_area_nm2,
            "zmax_cell_count": cell_count,
            "zmax_total_projected_area_nm2": cell_count * face_area_nm2,
            "exposed_classes": {
                "Ti3SiC2_skeleton": count_area(exposed_skeleton),
                "solid_Ag": count_area(solid_ag),
                "partially_melted_Ag_f_l_lt_0p5": count_area(partial_below_half),
                "partially_melted_Ag_f_l_ge_0p5": count_area(partial_at_or_above_half),
                "fully_liquid_Ag": count_area(fully_liquid_ag),
                "all_Ag_with_any_liquid_fraction": count_area(exposed_ag & (exposed_fraction > TOL)),
                "liquid_dominant_Ag_f_l_ge_0p5": count_area(liquid_dominant_ag),
                "maximum_tracked_pool_intersecting_zmax": count_area(exposed_pool),
            },
            "maximum_pool_zmax_temperature_K": {
                "min": float(exposed_temperature[exposed_pool].min()) if exposed_pool.any() else None,
                "max": float(exposed_temperature[exposed_pool].max()) if exposed_pool.any() else None,
                "mean": float(exposed_temperature[exposed_pool].mean()) if exposed_pool.any() else None,
            },
            "estimated_top_surface_liquid_skeleton_contact_lines": interface_lengths,
            "contact_edge_coordinates_csv": "exposed_zmax_liquid_skeleton_contact_edges.csv",
            "limitations": [
                "Top-face cells are voxelized exposed facets; contact lengths are grid-edge estimates, not atomistic contact-line lengths.",
                "Cells with 0 < liquid fraction < 1 contain unresolved solid/liquid coexistence; thresholds bracket the apparent liquid footprint.",
                "The exterior above zmax is treated only as gas for geometric classification; no gas state or physical boundary condition is solved.",
                "No pressure, velocity, flow rate, surface tension, contact angle, evaporation, or crater prediction is produced.",
            ],
        }
        (ROOT / "exposed_zmax_surface_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

        with (ROOT / "exposed_zmax_surface_cells.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(["ix", "iy", "x_center_nm", "y_center_nm", "surface_class", "liquid_fraction", "temperature_K", "is_maximum_tracked_pool"])
            names = ["Ti3SiC2_skeleton", "solid_Ag", "partial_Ag_f_l_lt_0p5", "partial_Ag_f_l_ge_0p5", "fully_liquid_Ag"]
            for ix in range(phase.shape[0]):
                for iy in range(phase.shape[1]):
                    writer.writerow([
                        ix, iy, (ix + 0.5) * dx_nm, (iy + 0.5) * dy_nm,
                        names[int(phase[ix, iy])], float(exposed_fraction[ix, iy]),
                        float(exposed_temperature[ix, iy]), bool(exposed_pool[ix, iy]),
                    ])

        with (ROOT / "exposed_zmax_liquid_skeleton_contact_edges.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow([
                "liquid_threshold", "boundary_normal_axis", "x0_nm", "y0_nm", "x1_nm", "y1_nm",
                "is_xy_periodic_seam", "liquid_fraction_on_Ag_side",
            ])
            writer.writerows(contact_edge_rows)

        colors = ["#18242c", "#8bb9cc", "#77ef35", "#fff200", "#f38a26"]
        cmap = ListedColormap(colors)
        fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
        ax.imshow(phase.T, origin="lower", interpolation="nearest", cmap=cmap, vmin=0, vmax=4,
                  extent=(0, phase.shape[0] * dx_nm, 0, phase.shape[1] * dy_nm), aspect="equal")
        if exposed_pool.any():
            ax.contour(np.arange(phase.shape[0]) * dx_nm + dx_nm / 2,
                       np.arange(phase.shape[1]) * dy_nm + dy_nm / 2,
                       exposed_pool.T.astype(float), levels=[0.5], colors=["white"], linewidths=0.65)
        labels = [
            ("Ti₃SiC₂ skeleton", colors[0]),
            ("Solid Ag", colors[1]),
            ("Partial Ag, fₗ < 0.5", colors[2]),
            ("Partial Ag, fₗ ≥ 0.5", colors[3]),
            ("Fully liquid Ag", colors[4]),
        ]
        handles = [plt.Line2D([0], [0], marker="s", linestyle="", markersize=9, markerfacecolor=c, markeredgecolor="none", label=l) for l, c in labels]
        ax.legend(handles=handles, loc="upper right", framealpha=0.96, fontsize=8)
        ax.set_xlabel("x (nm)")
        ax.set_ylabel("y (nm)")
        ax.set_title(f"Actual exposed z_max face | t = {metrics['checkpoint_time_ps']:.6f} ps\nWhite outline: tracked maximum liquid cluster footprint")
        fig.savefig(ROOT / "exposed_zmax_surface_map.png", dpi=220)
        plt.close(fig)

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
