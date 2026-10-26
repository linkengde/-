#!/usr/bin/env python3
"""Conservatively map 5 Å snapshot energies/enthalpies onto the 10 Å grid."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

MODEL_DIR = Path(__file__).resolve().parents[1]
PHASE_DIR = Path(__file__).resolve().parent
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))
if str(PHASE_DIR) not in sys.path:
    sys.path.insert(0, str(PHASE_DIR))

from ag_ttm_reference import AgParameters
from thermophysical_model import AgEnthalpy
from qsgs_enthalpy_solver import QSGSEnthalpyGrid
from qsgs_projection import _overlap_weights
from qsgs_ttm_enthalpy import field_snapshot_metrics


def conservative_average(field: np.ndarray, target_shape: tuple[int, int, int], size_A: np.ndarray) -> np.ndarray:
    source = np.asarray(field, dtype=float)
    weights = tuple(
        _overlap_weights(source.shape[axis], target_shape[axis], float(size_A[axis]))
        for axis in range(3)
    )
    return np.einsum("ai,bj,ck,ijk->abc", *weights, source, optimize=True)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    valid = (weights > 0.0) & np.isfinite(values)
    if not np.any(valid):
        return 0.0
    v, w = values[valid].ravel(), weights[valid].ravel()
    order = np.argsort(v)
    cumulative = np.cumsum(w[order])
    index = min(int(np.searchsorted(cumulative, q * cumulative[-1])), v.size - 1)
    return float(v[order[index]])


def stats_on_ag(values: np.ndarray, ag: np.ndarray, roi: tuple[slice, slice, slice]) -> dict[str, float]:
    v, w = np.asarray(values)[roi], np.asarray(ag)[roi]
    active = w > 1.0e-12
    if not np.any(active):
        return {"weighted_mean": 0.0, "weighted_q95": 0.0, "weighted_q99": 0.0, "maximum": 0.0}
    return {
        "weighted_mean": float(np.sum(v[active] * w[active]) / np.sum(w[active])),
        "weighted_q95": weighted_quantile(v, w, 0.95),
        "weighted_q99": weighted_quantile(v, w, 0.99),
        "weighted_q99.9": weighted_quantile(v, w, 0.999),
        "maximum": float(np.max(v[active])),
    }


def relative_l1(delta: np.ndarray, reference: np.ndarray, volume_m3: float) -> float:
    denominator = float(np.sum(np.abs(reference)) * volume_m3)
    numerator = float(np.sum(np.abs(delta)) * volume_m3)
    return numerator / denominator if denominator > 0.0 else 0.0


def weighted_relative_l1(
    delta: np.ndarray, reference: np.ndarray, weights: np.ndarray, volume_m3: float
) -> float:
    """Compare an intensive field only over Ag volume common to both projections."""
    denominator = float(np.sum(np.abs(reference) * weights) * volume_m3)
    numerator = float(np.sum(np.abs(delta) * weights) * volume_m3)
    return numerator / denominator if denominator > 0.0 else 0.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse-fields", type=Path, default=PHASE_DIR / "qsgs_localization_matched_3x3_10A_0p02ps_dt0p3125fs_fields.npz")
    parser.add_argument("--fine-fields", type=Path, default=PHASE_DIR / "qsgs_localization_matched_3x3_5A_0p02ps_dt0p3125fs_fields.npz")
    parser.add_argument("--output", type=Path, default=PHASE_DIR / "qsgs_localization_matched_0p02ps_grid_analysis.json")
    args = parser.parse_args()

    with np.load(args.coarse_fields, allow_pickle=False) as csrc, np.load(args.fine_fields, allow_pickle=False) as fsrc:
        coarse = {key: csrc[key] for key in csrc.files}
        fine = {key: fsrc[key] for key in fsrc.files}

    size_c, size_f = coarse["size_A"].astype(float), fine["size_A"].astype(float)
    shape_c, shape_f = coarse["ag_fraction"].shape, fine["ag_fraction"].shape
    if not np.allclose(size_c, size_f, rtol=0.0, atol=1.0e-8):
        raise ValueError("The compared grids must cover the same physical domain")
    coarse_dv = float(np.prod(coarse["spacing_A"]) * 1.0e-30)
    fine_dv = float(np.prod(fine["spacing_A"]) * 1.0e-30)
    ag_c = coarse["ag_fraction"].astype(float)
    ag_f = fine["ag_fraction"].astype(float)
    ag_f_mapped = conservative_average(ag_f, shape_c, size_c)
    ag_geometry_difference = ag_f_mapped - ag_c

    base_shape = tuple(int(v) for v in np.load(PHASE_DIR / "qsgs_phase_fractions_10A.npz", allow_pickle=False)["ag_fraction"].shape)
    tile = int(round(shape_c[0] / base_shape[0]))
    if tile != 3 or shape_c[1] != base_shape[1] * tile:
        raise ValueError("This analysis expects the centered 3×3 lateral repeat")
    start_x, start_y = (tile // 2) * base_shape[0], (tile // 2) * base_shape[1]
    central_roi = (slice(start_x, start_x + base_shape[0]), slice(start_y, start_y + base_shape[1]), slice(None))
    full_roi = (slice(None), slice(None), slice(None))

    table = MODEL_DIR / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv"
    parameters = AgParameters.from_csv(table, allow_low_temperature_extension=False)
    ag_enthalpy = AgEnthalpy()
    cell_volume_A3_c = float(np.prod(coarse["spacing_A"]))
    cell_volume_A3_f = float(np.prod(fine["spacing_A"]))
    times = sorted(
        float(key[2:].split("ps_", 1)[0])
        for key in coarse
        if key.startswith("t_") and key.endswith("ps_electron_temperature_K")
    )
    fine_times = sorted(
        float(key[2:].split("ps_", 1)[0])
        for key in fine
        if key.startswith("t_") and key.endswith("ps_electron_temperature_K")
    )
    if times != fine_times:
        raise ValueError("The matched grids must have identical snapshot times")

    results = []
    for time_ps in times:
        label = f"t_{time_ps:.9g}ps_"
        te_c = coarse[label + "electron_temperature_K"].astype(float)
        tl_c = coarse[label + "lattice_temperature_K"].astype(float)
        liquid_c = coarse[label + "ag_liquid_fraction"].astype(float)
        e_c = coarse[label + "electron_energy_density_eV_A3"].astype(float)
        h_c = coarse[label + "lattice_enthalpy_density_J_m3"].astype(float)
        g_c = coarse[label + "electron_to_lattice_cumulative_energy_density_J_m3"].astype(float)

        te_f = fine[label + "electron_temperature_K"].astype(float)
        tl_f = fine[label + "lattice_temperature_K"].astype(float)
        liquid_f = fine[label + "ag_liquid_fraction"].astype(float)
        e_f = fine[label + "electron_energy_density_eV_A3"].astype(float)
        h_f = fine[label + "lattice_enthalpy_density_J_m3"].astype(float)
        g_f = fine[label + "electron_to_lattice_cumulative_energy_density_J_m3"].astype(float)

        e_f_mapped = conservative_average(e_f, shape_c, size_c)
        h_f_mapped = conservative_average(h_f, shape_c, size_c)
        g_f_mapped = conservative_average(g_f, shape_c, size_c)
        common_ag_fraction = np.minimum(ag_c, ag_f_mapped)
        ue_c_specific = np.divide(
            e_c, ag_c, out=np.zeros_like(e_c), where=ag_c > 1.0e-12
        )
        ue_f_specific = np.divide(
            e_f_mapped, ag_f_mapped, out=np.zeros_like(e_f_mapped),
            where=ag_f_mapped > 1.0e-12,
        )
        h_c_specific = np.divide(
            h_c, ag_c, out=np.zeros_like(h_c), where=ag_c > 1.0e-12
        )
        h_f_specific = np.divide(
            h_f_mapped, ag_f_mapped, out=np.zeros_like(h_f_mapped),
            where=ag_f_mapped > 1.0e-12,
        )
        g_c_specific = np.divide(
            g_c, ag_c, out=np.zeros_like(g_c), where=ag_c > 1.0e-12
        )
        g_f_specific = np.divide(
            g_f_mapped, ag_f_mapped, out=np.zeros_like(g_f_mapped),
            where=ag_f_mapped > 1.0e-12,
        )
        ue_300 = float(parameters.ue(300.0))
        liquid_volume_f_mapped = conservative_average(ag_f * liquid_f, shape_c, size_c)
        liquid_fraction_f_mapped = np.divide(
            liquid_volume_f_mapped,
            ag_f_mapped,
            out=np.zeros_like(ag_f_mapped),
            where=ag_f_mapped > 1.0e-12,
        )
        ue_f_per_ag = np.divide(
            e_f_mapped, ag_f_mapped, out=np.zeros_like(ag_f_mapped), where=ag_f_mapped > 1.0e-12
        )
        te_f_mapped = parameters.temperature_from_ue(ue_f_per_ag)
        h_f_mapped_specific = np.divide(
            h_f_mapped,
            ag_f_mapped * ag_enthalpy.reference_density_kg_m3,
            out=np.zeros_like(ag_f_mapped),
            where=ag_f_mapped > 1.0e-12,
        )
        state_f_mapped = ag_enthalpy.state_from_specific_enthalpy(h_f_mapped_specific)

        fine_e_total_J = float(e_f.sum() * fine_dv * 1.602176634e-19 * 1.0e30)
        mapped_e_total_J = float(e_f_mapped.sum() * coarse_dv * 1.602176634e-19 * 1.0e30)
        fine_h_total_J = float(h_f.sum() * fine_dv)
        mapped_h_total_J = float(h_f_mapped.sum() * coarse_dv)
        fine_g_total_J = float(g_f.sum() * fine_dv)
        mapped_g_total_J = float(g_f_mapped.sum() * coarse_dv)

        entry: dict[str, object] = {"time_ps": time_ps}
        entry["raw_electron_temperature_10A_K"] = stats_on_ag(te_c, ag_c, full_roi)
        entry["raw_electron_temperature_5A_K"] = stats_on_ag(te_f, ag_f, (full_roi[0], full_roi[1], full_roi[2]))
        entry["mapped_electron_temperature_5A_to_10A_K"] = stats_on_ag(te_f_mapped, ag_f_mapped, full_roi)
        entry["raw_lattice_temperature_10A_K"] = stats_on_ag(tl_c, ag_c, full_roi)
        entry["raw_lattice_temperature_5A_K"] = stats_on_ag(tl_f, ag_f, full_roi)
        entry["ag_fraction_remap_max_abs_difference"] = float(np.max(np.abs(ag_geometry_difference)))
        entry["ag_fraction_remap_relative_L1_full"] = relative_l1(
            ag_geometry_difference, ag_c, coarse_dv
        )
        entry["ag_fraction_remap_relative_L1_central_tile"] = relative_l1(
            ag_geometry_difference[central_roi], ag_c[central_roi], coarse_dv
        )
        metrics_coarse = field_snapshot_metrics(
            time_ps,
            {"electron_temperature_K": te_c, "lattice_temperature_K": tl_c,
             "ag_liquid_fraction": liquid_c},
            ag_c,
            np.asarray(coarse["spacing_A"], dtype=float),
            size_c,
            roi_xyz=central_roi,
        )
        metrics_fine_mapped = field_snapshot_metrics(
            time_ps,
            {"electron_temperature_K": te_f_mapped,
             "lattice_temperature_K": state_f_mapped.temperature_K,
             "ag_liquid_fraction": liquid_fraction_f_mapped},
            ag_f_mapped,
            np.asarray(coarse["spacing_A"], dtype=float),
            size_c,
            roi_xyz=central_roi,
        )
        for field_name in (
            "electron_half_rise_hotspot_diameter_2rms_A",
            "electron_ag_volume_fraction_above_300_plus_5K",
            "electron_ag_weighted_q0.99_K",
            "electron_ag_weighted_q0.999_K",
            "lattice_peak_K",
        ):
            entry[f"central_tile_{field_name}_10A"] = metrics_coarse[field_name]
            entry[f"central_tile_{field_name}_5A_mapped"] = metrics_fine_mapped[field_name]
        for scope, roi in (("full", full_roi), ("central_tile", central_roi)):
            weights = np.minimum(ag_c[roi], ag_f_mapped[roi])
            te_difference = np.abs(te_f_mapped[roi] - te_c[roi])
            tl_difference = np.abs(state_f_mapped.temperature_K[roi] - tl_c[roi])
            entry[f"{scope}_electron_temperature_abs_difference_K"] = stats_on_ag(te_difference, weights, (slice(None),) * 3)
            entry[f"{scope}_lattice_temperature_abs_difference_K"] = stats_on_ag(tl_difference, weights, (slice(None),) * 3)
            entry[f"{scope}_electron_energy_density_relative_L1"] = relative_l1(
                (e_f_mapped - e_c)[roi], e_c[roi], coarse_dv
            )
            entry[f"{scope}_lattice_enthalpy_density_relative_L1"] = relative_l1(
                (h_f_mapped - h_c)[roi], h_c[roi], coarse_dv
            )
            entry[f"{scope}_coupling_map_relative_L1"] = relative_l1(
                (g_f_mapped - g_c)[roi], g_c[roi], coarse_dv
            )
            entry[f"{scope}_common_ag_excess_electron_energy_relative_L1"] = weighted_relative_l1(
                (ue_f_specific - ue_c_specific)[roi],
                (ue_c_specific - ue_300)[roi], weights, coarse_dv,
            )
            entry[f"{scope}_common_ag_lattice_enthalpy_relative_L1"] = weighted_relative_l1(
                (h_f_specific - h_c_specific)[roi], h_c_specific[roi], weights, coarse_dv
            )
            entry[f"{scope}_common_ag_coupling_map_relative_L1"] = weighted_relative_l1(
                (g_f_specific - g_c_specific)[roi], g_c_specific[roi], weights, coarse_dv
            )
            entry[f"{scope}_liquid_Ag_volume_fraction_10A"] = float(
                np.sum(ag_c[roi] * liquid_c[roi]) / np.sum(ag_c[roi])
            )
            entry[f"{scope}_liquid_Ag_volume_fraction_mapped_5A"] = float(
                np.sum(liquid_volume_f_mapped[roi]) / np.sum(ag_f_mapped[roi])
            )

        entry["electron_energy_5A_total_J"] = fine_e_total_J
        entry["electron_energy_5A_mapped_total_J"] = mapped_e_total_J
        entry["electron_energy_10A_total_J"] = float(e_c.sum() * coarse_dv * 1.602176634e-19 * 1.0e30)
        entry["lattice_enthalpy_5A_total_J"] = fine_h_total_J
        entry["lattice_enthalpy_5A_mapped_total_J"] = mapped_h_total_J
        entry["lattice_enthalpy_10A_total_J"] = float(h_c.sum() * coarse_dv)
        entry["coupling_transfer_5A_total_J"] = fine_g_total_J
        entry["coupling_transfer_5A_mapped_total_J"] = mapped_g_total_J
        entry["coupling_transfer_10A_total_J"] = float(g_c.sum() * coarse_dv)
        entry["local_lattice_enthalpy_vs_coupling_relative_L1_10A"] = relative_l1(
            h_c - g_c, h_c, coarse_dv
        )
        entry["local_lattice_enthalpy_vs_coupling_relative_L1_5A"] = relative_l1(
            h_f - g_f, h_f, fine_dv
        )
        entry["liquid_Ag_volume_fraction_10A"] = float(np.sum(ag_c * liquid_c) / ag_c.sum())
        entry["liquid_Ag_volume_fraction_mapped_5A"] = float(
            np.sum(liquid_volume_f_mapped) / ag_f_mapped.sum()
        )
        results.append(entry)

    summary = {
        "status": "MATCHED_3X3_GRID_COMPARISON_0P02PS",
        "mapping": "conservative overlap-volume averaging of energy density, enthalpy density, cumulative local coupling density, Ag fraction, and liquid Ag volume density from 5 Å onto 10 Å cells",
        "boundary_case": "A",
        "domain_size_A": size_c.tolist(),
        "grid_10A_shape_xyz": list(shape_c),
        "grid_5A_shape_xyz": list(shape_f),
        "heat_duration_ps": 0.02,
        "step_fs": 0.3125,
        "source_peak_flux_eV_ps_A2": 5.0,
        "gaussian_radius_A": 30.0,
        "absorption_depth_A": 20.0,
        "times": results,
        "limitations": [
            "This is a short continuum two-temperature/enthalpy localization screen, not a melting or solidification result.",
            "The 3×3 domain repeats the same QSGS cell and is not an independent morphology sample.",
            "The Ag and Ti3SiC2 boundary assumptions remain idealized; this screen uses adiabatic scaffold case A.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    csv_path = args.output.with_suffix(".csv")
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["time_ps", "electron_energy_10A_total_J", "electron_energy_5A_total_J", "lattice_enthalpy_10A_total_J", "lattice_enthalpy_5A_total_J", "coupling_transfer_10A_total_J", "coupling_transfer_5A_total_J", "local_lattice_enthalpy_vs_coupling_relative_L1_10A", "local_lattice_enthalpy_vs_coupling_relative_L1_5A", "raw_electron_peak_10A_K", "raw_electron_peak_5A_K", "mapped_fine_electron_peak_K", "full_electron_temperature_abs_difference_K", "central_tile_electron_temperature_abs_difference_K", "full_lattice_temperature_abs_difference_K", "central_tile_lattice_temperature_abs_difference_K", "full_electron_energy_density_relative_L1", "central_tile_electron_energy_density_relative_L1", "full_lattice_enthalpy_density_relative_L1", "central_tile_lattice_enthalpy_density_relative_L1", "full_coupling_map_relative_L1", "central_tile_coupling_map_relative_L1", "full_common_ag_excess_electron_energy_relative_L1", "central_tile_common_ag_excess_electron_energy_relative_L1", "full_common_ag_lattice_enthalpy_relative_L1", "central_tile_common_ag_lattice_enthalpy_relative_L1", "full_common_ag_coupling_map_relative_L1", "central_tile_common_ag_coupling_map_relative_L1", "liquid_Ag_volume_fraction_10A", "liquid_Ag_volume_fraction_mapped_5A"], lineterminator="\n")
        writer.writeheader()
        for entry in results:
            flat = {key: entry[key] for key in writer.fieldnames if key in entry}
            flat["raw_electron_peak_10A_K"] = entry["raw_electron_temperature_10A_K"]["maximum"]
            flat["raw_electron_peak_5A_K"] = entry["raw_electron_temperature_5A_K"]["maximum"]
            flat["mapped_fine_electron_peak_K"] = entry["mapped_electron_temperature_5A_to_10A_K"]["maximum"]
            for scope in ("full", "central_tile"):
                flat[f"{scope}_electron_temperature_abs_difference_K"] = entry[f"{scope}_electron_temperature_abs_difference_K"]["weighted_mean"]
                flat[f"{scope}_lattice_temperature_abs_difference_K"] = entry[f"{scope}_lattice_temperature_abs_difference_K"]["weighted_mean"]
            writer.writerow(flat)
    print(args.output)
    print(csv_path)


if __name__ == "__main__":
    main()
