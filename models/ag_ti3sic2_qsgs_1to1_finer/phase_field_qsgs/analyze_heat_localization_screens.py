"""Recompute localization metrics while excluding round-off-only Ag slivers."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
AG_FRACTION_TOL = 1.0e-12
CASES = (
    "base_volume_fraction",
    "base_exact_contacts",
    "3x3_volume_fraction",
    "3x3_exact_contacts",
    "5A_exact_contacts",
)


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    use = (weights >= AG_FRACTION_TOL) & np.isfinite(values)
    values, weights = values[use], weights[use]
    order = np.argsort(values)
    sorted_values, sorted_weights = values[order], weights[order]
    cumulative = np.cumsum(sorted_weights)
    return float(sorted_values[min(np.searchsorted(cumulative, quantile * cumulative[-1]), len(values) - 1)])


def analyze_case(stem: str) -> dict[str, object]:
    summary_path = HERE / f"qsgs_localization_{stem}.json"
    fields_path = HERE / f"qsgs_localization_{stem}_fields.npz"
    summary = json.loads(summary_path.read_text())
    with np.load(fields_path, allow_pickle=False) as data:
        ag = np.asarray(data["ag_fraction"], dtype=float)
        spacing = np.asarray(data["spacing_A"], dtype=float)
        size = np.asarray(data["size_A"], dtype=float)
        te = np.asarray(data["t_0.1ps_electron_temperature_K"], dtype=float)
        tl = np.asarray(data["t_0.1ps_lattice_temperature_K"], dtype=float)
        liquid = np.asarray(data["t_0.1ps_ag_liquid_fraction"], dtype=float)

    active = ag >= AG_FRACTION_TOL
    weight = ag[active]
    total = float(weight.sum())
    raw_index = np.unravel_index(np.argmax(te), te.shape)
    meaningful_peak = float(te[active].max())
    half_rise_threshold = 300.0 + 0.5 * max(meaningful_peak - 300.0, 0.0)
    hot = active & (te >= half_rise_threshold)
    hot_weight = float(ag[hot].sum())
    x = (np.arange(ag.shape[0]) + 0.5) * spacing[0]
    y = (np.arange(ag.shape[1]) + 0.5) * spacing[1]
    cx, cy = size[0] / 2.0, size[1] / 2.0
    xrel = (x - cx + size[0] / 2.0) % size[0] - size[0] / 2.0
    yrel = (y - cy + size[1] / 2.0) % size[1] - size[1] / 2.0
    r2 = np.broadcast_to(
        xrel[:, None, None] ** 2 + yrel[None, :, None] ** 2, ag.shape
    )
    radius_rms = float(np.sqrt(np.sum(ag[hot] * r2[hot]) / hot_weight)) if hot_weight > 0 else 0.0
    depth = size[2] - (np.arange(ag.shape[2]) + 0.5) * spacing[2]
    hot_depth = np.broadcast_to(depth[None, None, :], ag.shape)[hot]

    result: dict[str, object] = {
        "case": stem,
        "grid_shape_xyz": summary["grid_shape_xyz"],
        "domain_size_A": summary["domain_size_A"],
        "xy_tile_count": summary["xy_tile_count"],
        "source_power_grid_integral_W": summary["source_power_grid_integral_W"],
        "source_energy_total_J": summary["source_energy_total_J"],
        "energy_balance_residual_J": summary["energy_balance_residual_J"],
        "electron_to_lattice_cumulative_J": summary["electron_to_lattice_cumulative_J"],
        "raw_electron_temperature_max_K": float(te.max()),
        "raw_electron_max_ag_fraction": float(ag[raw_index]),
        "electron_temperature_max_K_ag_fraction_ge_1e-12": meaningful_peak,
        "electron_temperature_ag_weighted_mean_K": float(np.sum(ag[active] * te[active]) / total),
        "electron_temperature_ag_weighted_quantiles_K": {
            str(q): _weighted_quantile(te, ag, q) for q in (0.5, 0.9, 0.99, 0.999)
        },
        "electron_half_rise_threshold_K_filtered_peak": half_rise_threshold,
        "electron_half_rise_hotspot_ag_volume_fraction": hot_weight / total if total > 0 else 0.0,
        "electron_half_rise_hotspot_diameter_2rms_A": 2.0 * radius_rms,
        "electron_half_rise_hotspot_depth_extent_A": (
            float(hot_depth.max() - hot_depth.min()) if hot_depth.size else 0.0
        ),
        "electron_ag_volume_fraction_above_500K": float(ag[active & (te >= 500.0)].sum() / total),
        "electron_ag_volume_fraction_above_1000K": float(ag[active & (te >= 1000.0)].sum() / total),
        "lattice_temperature_max_K_ag_fraction_ge_1e-12": float(tl[active].max()),
        "lattice_temperature_ag_weighted_mean_K": float(np.sum(ag[active] * tl[active]) / total),
        "max_liquid_fraction_at_0p1ps": float(liquid[active].max()),
        "ag_liquid_volume_fraction_at_0p1ps": float(np.sum(ag[active] * liquid[active]) / total),
        "phase_fraction_filter": "Ag fractions below 1e-12 are round-off-only cells and excluded from reported extrema and weighted metrics",
        "excluded_roundoff_ag_fraction_of_total": float(ag[~active].sum() / ag.sum()),
    }
    if summary["xy_tile_count"] > 1:
        nx, ny = ag.shape[0] // summary["xy_tile_count"], ag.shape[1] // summary["xy_tile_count"]
        center = summary["xy_tile_count"] // 2
        roi = (slice(center * nx, (center + 1) * nx), slice(center * ny, (center + 1) * ny), slice(None))
        roi_ag, roi_te = ag[roi], te[roi]
        roi_active = roi_ag >= AG_FRACTION_TOL
        result["central_tile_electron_half_rise_ag_volume_fraction"] = (
            float(roi_ag[roi_active & (roi_te >= half_rise_threshold)].sum() / roi_ag[roi_active].sum())
        )
        result["central_tile_electron_ag_weighted_mean_K"] = float(
            np.sum(roi_ag[roi_active] * roi_te[roi_active]) / roi_ag[roi_active].sum()
        )
    return result


def main() -> None:
    report = {
        "status": "0.1_PS_HEAT_LOCALIZATION_SCREEN_REANALYZED",
        "method": "10 Å source-QSGS face-contact correction and 3x3 periodic repeat-domain comparison; 5 Å corrected-grid check",
        "excluded_geometrical_roundoff": "Ag fraction < 1e-12",
        "cases": [analyze_case(stem) for stem in CASES],
    }
    output = HERE / "qsgs_localization_screen_metrics.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
