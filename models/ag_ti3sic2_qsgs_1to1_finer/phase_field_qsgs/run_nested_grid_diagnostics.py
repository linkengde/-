#!/usr/bin/env python3
"""Run matched nested-grid source, conduction, and full-coupling screens."""

from __future__ import annotations

import csv
import gc
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/nested-grid-mplconfig")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/nested-grid-cache")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.environ["XDG_CACHE_HOME"], exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
MODEL_DIR = HERE.parent
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ag_ttm_reference import AgParameters
from thermophysical_model import AgEnthalpy
from analyze_matched_grid_fields import conservative_average, stats_on_ag, weighted_relative_l1
from qsgs_ttm_enthalpy import QSGSDoubleTemperatureGrid, field_snapshot_metrics, write_csv


OUT = HERE / "nested_grid_diagnostics"
SOURCE_GEOMETRY = MODEL_DIR / "qsgs_morphology.npz"
SCHEDULE_PS = (0.001, 0.005, 0.02)
DT_PS = 0.0003125
DURATION_PS = 0.02
TILE_XY = 3
MODE_FLAGS = {
    "source_only": {"electron_diffusion_on": False, "coupling_on": False},
    "electron_diffusion": {"electron_diffusion_on": True, "coupling_on": False},
    "full_coupling": {"electron_diffusion_on": True, "coupling_on": True},
}
GRID_FILES = {
    "10A_nested": ("qsgs_phase_fractions_10A_nested.npz", "qsgs_ag_face_contacts_10A_nested.npz"),
    "5A_nested": ("qsgs_phase_fractions_5A_nested.npz", "qsgs_ag_face_contacts_5A_nested.npz"),
}


def integrated_energy_l1(actual: np.ndarray, expected: np.ndarray, volume_m3: float) -> float:
    denominator = float(np.sum(np.abs(expected)) * volume_m3)
    numerator = float(np.sum(np.abs(actual - expected)) * volume_m3)
    return numerator / denominator if denominator > 0.0 else 0.0


def run_case(mode: str, grid: str) -> tuple[dict, QSGSDoubleTemperatureGrid]:
    projected_name, contact_name = GRID_FILES[grid]
    projected = HERE / projected_name
    contacts = HERE / contact_name
    model = QSGSDoubleTemperatureGrid(
        projected,
        "A",
        face_contacts_npz=contacts,
        tile_xy=TILE_XY,
        source_geometry_npz=SOURCE_GEOMETRY,
    )
    rows = model.run(
        DURATION_PS,
        DT_PS,
        heat_on=True,
        sample_interval_ps=0.001,
        snapshot_times_ps=SCHEDULE_PS,
        **MODE_FLAGS[mode],
    )
    stem = f"{mode}_{grid}_3x3_0p02ps_dt0p3125fs"
    write_csv(OUT / f"{stem}.csv", rows)
    base_shape = np.load(HERE / "qsgs_phase_fractions_10A.npz", allow_pickle=False)["ag_fraction"].shape
    center = TILE_XY // 2
    refinement = 2 if grid == "5A_nested" else 1
    roi_nx, roi_ny = base_shape[0] * refinement, base_shape[1] * refinement
    roi = (
        slice(center * roi_nx, (center + 1) * roi_nx),
        slice(center * roi_ny, (center + 1) * roi_ny),
        slice(None),
    )
    snapshot_metrics = []
    for time_ps in SCHEDULE_PS:
        fields = model.field_snapshots[time_ps]
        snapshot_metrics.append(
            field_snapshot_metrics(
                time_ps,
                fields,
                model.lattice.ag_fraction,
                model.lattice.spacing_A,
                model.lattice.size_A,
                roi_xyz=roi,
            )
        )
    result = {
        "mode": mode,
        "grid": grid,
        "projected_phase_grid_npz": projected.name,
        "face_contacts_npz": contacts.name,
        "shape_xyz_3x3": list(model.lattice.temperature_shape),
        "spacing_A_xyz": model.lattice.spacing_A.tolist(),
        "source_power_W": model.actual_source_power_W,
        "source_energy_J": model.injected_energy_J,
        "electron_energy_final_J": model.electron_energy_J,
        "lattice_enthalpy_final_J": model.lattice_energy_J,
        "electron_to_lattice_energy_J": model.coupling_energy_J,
        "electron_diffusion_net_J": model.electron_diffusion_energy_J,
        "energy_balance_residual_J": (
            model.total_energy_J - model.initial_total_energy_J - model.injected_energy_J
            - model.electron_diffusion_energy_J - model.lattice_diffusion_energy_J
        ),
        "electron_diffusion_enabled": MODE_FLAGS[mode]["electron_diffusion_on"],
        "electron_lattice_coupling_enabled": MODE_FLAGS[mode]["coupling_on"],
        "snapshot_metrics_central_unit_cell": snapshot_metrics,
        "time_history_csv": f"{stem}.csv",
    }
    return result, model


def compare_pair(coarse: QSGSDoubleTemperatureGrid, fine: QSGSDoubleTemperatureGrid,
                 parameters: AgParameters, time_ps: float) -> dict:
    shape = coarse.lattice.ag_fraction.shape
    size_A = coarse.lattice.size_A
    coarse_ag = coarse.lattice.ag_fraction
    fine_ag_mapped = conservative_average(fine.lattice.ag_fraction, shape, size_A)
    common_ag = np.minimum(coarse_ag, fine_ag_mapped)
    volume_m3 = coarse.lattice.cell_volume_m3
    cell_volume_A3 = float(np.prod(coarse.lattice.spacing_A))
    label = f"t_{time_ps:.9g}ps"
    coarse_fields = coarse.field_snapshots[time_ps]
    fine_fields = fine.field_snapshots[time_ps]

    e_f = conservative_average(
        fine_fields["electron_energy_density_eV_A3"], shape, size_A
    )
    h_f = conservative_average(
        fine_fields["lattice_enthalpy_density_J_m3"], shape, size_A
    )
    g_f = conservative_average(
        fine_fields["electron_to_lattice_cumulative_energy_density_J_m3"], shape, size_A
    )
    ue_coarse = np.divide(
        coarse_fields["electron_energy_density_eV_A3"], coarse_ag,
        out=np.zeros_like(coarse_ag), where=coarse_ag > 1.0e-12,
    )
    ue_fine = np.divide(
        e_f, fine_ag_mapped, out=np.zeros_like(fine_ag_mapped),
        where=fine_ag_mapped > 1.0e-12,
    )
    h_coarse = np.divide(
        coarse_fields["lattice_enthalpy_density_J_m3"], coarse_ag,
        out=np.zeros_like(coarse_ag), where=coarse_ag > 1.0e-12,
    )
    h_fine = np.divide(
        h_f, fine_ag_mapped, out=np.zeros_like(fine_ag_mapped),
        where=fine_ag_mapped > 1.0e-12,
    )
    g_coarse = np.divide(
        coarse_fields["electron_to_lattice_cumulative_energy_density_J_m3"], coarse_ag,
        out=np.zeros_like(coarse_ag), where=coarse_ag > 1.0e-12,
    )
    g_fine = np.divide(
        g_f, fine_ag_mapped, out=np.zeros_like(fine_ag_mapped),
        where=fine_ag_mapped > 1.0e-12,
    )
    te_f = parameters.temperature_from_ue(ue_fine)
    rho_ag = coarse.lattice.ag.reference_density_kg_m3
    state_f = coarse.lattice.ag.state_from_specific_enthalpy(
        np.divide(h_f, fine_ag_mapped * rho_ag,
                  out=np.zeros_like(fine_ag_mapped), where=fine_ag_mapped > 1.0e-12)
    )
    te_c = coarse_fields["electron_temperature_K"]
    tl_c = coarse_fields["lattice_temperature_K"]
    roi_shape = np.load(HERE / "qsgs_phase_fractions_10A.npz", allow_pickle=False)["ag_fraction"].shape
    center = TILE_XY // 2
    roi = (
        slice(center * roi_shape[0], (center + 1) * roi_shape[0]),
        slice(center * roi_shape[1], (center + 1) * roi_shape[1]),
        slice(None),
    )
    te_difference = stats_on_ag(np.abs(te_f - te_c), common_ag, roi)
    tl_difference = stats_on_ag(np.abs(state_f.temperature_K - tl_c), common_ag, roi)
    coarse_peak = float(np.max(te_c[coarse_ag > 0.0]))
    fine_peak = float(np.max(fine_fields["electron_temperature_K"][fine.lattice.ag_fraction > 0.0]))
    return {
        "time_ps": time_ps,
        "raw_electron_peak_10A_K": coarse_peak,
        "raw_electron_peak_5A_nested_K": fine_peak,
        "raw_electron_peak_abs_difference_K": abs(fine_peak - coarse_peak),
        "raw_electron_peak_relative_difference": abs(fine_peak - coarse_peak) / max(abs(coarse_peak), 1.0e-30),
        "raw_electron_peak_rise_relative_difference": abs(fine_peak - coarse_peak) / max(abs(coarse_peak - 300.0), 1.0e-30),
        "central_hotspot_diameter_10A_A": field_snapshot_metrics(
            time_ps, coarse_fields, coarse_ag, coarse.lattice.spacing_A, size_A, roi_xyz=roi
        )["electron_half_rise_hotspot_diameter_2rms_A"],
        "central_hotspot_diameter_5A_nested_mapped_A": field_snapshot_metrics(
            time_ps,
            {"electron_temperature_K": te_f,
             "lattice_temperature_K": state_f.temperature_K,
             "ag_liquid_fraction": state_f.liquid_fraction},
            fine_ag_mapped, coarse.lattice.spacing_A, size_A, roi_xyz=roi,
        )["electron_half_rise_hotspot_diameter_2rms_A"],
        "common_ag_weighted_electron_temperature_abs_difference_K": te_difference,
        "common_ag_weighted_lattice_temperature_abs_difference_K": tl_difference,
        "common_ag_excess_electron_energy_relative_L1": weighted_relative_l1(
            ue_fine - ue_coarse,
            ue_coarse - float(parameters.ue(300.0)),
            common_ag, volume_m3,
        ),
        "common_ag_lattice_enthalpy_relative_L1": weighted_relative_l1(
            h_fine - h_coarse, h_coarse, common_ag, volume_m3,
        ),
        "common_ag_coupling_energy_relative_L1": weighted_relative_l1(
            g_fine - g_coarse, g_coarse, common_ag, volume_m3,
        ),
        "mapped_electron_energy_conservation_relative_difference": (
            abs(float(e_f.sum() * cell_volume_A3) - float(
                fine_fields["electron_energy_density_eV_A3"].sum()
                * np.prod(fine.lattice.spacing_A)
            ))
            / max(abs(float(
                fine_fields["electron_energy_density_eV_A3"].sum()
                * np.prod(fine.lattice.spacing_A)
            )), 1.0e-30)
        ),
        "ag_fraction_mapping_relative_L1": integrated_energy_l1(
            fine_ag_mapped, coarse_ag, volume_m3
        ),
        "liquid_ag_fraction_10A": float(coarse.lattice._liquid_ag_volume_fraction()),
        "liquid_ag_fraction_5A_nested": float(fine.lattice._liquid_ag_volume_fraction()),
    }


def save_mode_plot(plot_data: dict[str, dict[str, dict]]) -> None:
    modes = list(MODE_FLAGS)
    fig, axes = plt.subplots(len(modes), 2, figsize=(10, 11), constrained_layout=True)
    vmax = max(
        float(grid_data["temperature_rise_K"].max())
        for pair in plot_data.values() for grid_data in pair.values()
    )
    image = None
    for row, mode in enumerate(modes):
        for col, grid in enumerate(("10A_nested", "5A_nested")):
            data = plot_data[mode][grid]
            ag = data["ag_surface_fraction"]
            temperature_rise = data["temperature_rise_K"]
            spacing = np.asarray(data["spacing_A_xyz"])
            size = np.asarray(data["unit_size_A"])
            x_centers = (np.arange(ag.shape[0]) + 0.5) * spacing[0] - size[0] / 2.0
            y_centers = (np.arange(ag.shape[1]) + 0.5) * spacing[1] - size[1] / 2.0
            extent = (x_centers[0] - spacing[0] / 2.0, x_centers[-1] + spacing[0] / 2.0,
                      y_centers[0] - spacing[1] / 2.0, y_centers[-1] + spacing[1] / 2.0)
            ax = axes[row, col]
            image = ax.imshow(
                temperature_rise.T,
                origin="lower", extent=extent, interpolation="nearest",
                cmap="inferno", vmin=0.0, vmax=vmax,
                alpha=np.clip(ag.T, 0.0, 1.0),
            )
            ax.set_title(f"{mode.replace('_', ' ')} — {'10 Å' if col == 0 else '5 Å nested'}")
            ax.set_xlabel("x relative to source center (Å)")
            ax.set_ylabel("y relative to source center (Å)")
            ax.set_xlim(-100.0, 100.0)
            ax.set_ylim(-100.0, 100.0)
            ax.set_aspect("equal")
    fig.colorbar(image, ax=axes, label="Electron temperature rise (K)")
    fig.savefig(OUT / "nested_grid_20fs_mode_comparison.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    parameter_table = MODEL_DIR / "ttm_data" / "Ag_TTM_hybrid_figure_fit_300_20000K.csv"
    parameters = AgParameters.from_csv(parameter_table, allow_low_temperature_extension=False)
    summaries: dict[str, dict[str, dict]] = {}
    plot_data = {}
    for mode in MODE_FLAGS:
        print(f"Starting {mode} 10 Å nested run", flush=True)
        coarse_summary, coarse_model = run_case(mode, "10A_nested")
        print(f"Completed {mode} 10 Å; starting 5 Å nested run", flush=True)
        fine_summary, fine_model = run_case(mode, "5A_nested")
        summaries[mode] = {"10A_nested": coarse_summary, "5A_nested": fine_summary}
        summaries[mode]["mapped_comparison"] = [
            compare_pair(coarse_model, fine_model, parameters, time_ps)
            for time_ps in SCHEDULE_PS
        ]
        plot_data[mode] = {}
        for grid, model in (("10A_nested", coarse_model), ("5A_nested", fine_model)):
            base_shape = np.load(HERE / "qsgs_phase_fractions_10A.npz", allow_pickle=False)["ag_fraction"].shape
            center = TILE_XY // 2
            nx, ny = (base_shape[0] * (1 if grid == "10A_nested" else 2),
                      base_shape[1] * (1 if grid == "10A_nested" else 2))
            x0, y0 = center * nx, center * ny
            ag = model.lattice.ag_fraction[x0:x0 + nx, y0:y0 + ny, -1]
            te = model.field_snapshots[0.02]["electron_temperature_K"][x0:x0 + nx, y0:y0 + ny, -1]
            plot_data[mode][grid] = {
                "ag_surface_fraction": ag.copy(),
                "temperature_rise_K": (te - 300.0).copy(),
                "spacing_A_xyz": model.lattice.spacing_A.tolist(),
                "unit_size_A": (model.lattice.size_A / np.asarray([TILE_XY, TILE_XY, 1.0])).tolist(),
            }
        print(f"Completed {mode} pair", flush=True)
        del coarse_model, fine_model
        gc.collect()
    save_mode_plot(plot_data)
    output = {
        "status": "STRICTLY_NESTED_5A_10A_20FS_TRANSPORT_ABLATION",
        "grid_geometry_audit": "qsgs_nested_grid_geometry_audit.json",
        "domain_tiles_xy": TILE_XY,
        "boundary_case": "A; no bottom heat sink",
        "dt_fs": DT_PS * 1000.0,
        "heat_duration_ps": DURATION_PS,
        "snapshot_times_ps": list(SCHEDULE_PS),
        "source": {"I0_eV_ps_A2": 5.0, "R_A": 30.0, "spatial_form": "exp(-r²/R²) exp(-depth/20Å)"},
        "modes": summaries,
        "limitations": [
            "The source-only and electron-diffusion modes are diagnostic ablations, not alternate physical models.",
            "Common Ag weighting uses the minimum of the two projected cell fractions, not exact subcell phase intersections.",
            "The nested 5 Å grid has anisotropic spacing while each axis is exactly half the corresponding coarse spacing.",
            "No melting is expected in this 20 fs localization screen; this is not a melt-pool result.",
        ],
    }
    path = OUT / "nested_grid_20fs_mode_comparison.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    csv_path = OUT / "nested_grid_20fs_mode_comparison.csv"
    with csv_path.open("w", newline="") as stream:
        fields = [
            "mode", "time_ps", "raw_electron_peak_10A_K", "raw_electron_peak_5A_nested_K",
            "raw_electron_peak_abs_difference_K", "raw_electron_peak_relative_difference",
            "raw_electron_peak_rise_relative_difference", "common_ag_excess_electron_energy_relative_L1",
            "common_ag_lattice_enthalpy_relative_L1", "common_ag_coupling_energy_relative_L1",
            "ag_fraction_mapping_relative_L1", "liquid_ag_fraction_10A", "liquid_ag_fraction_5A_nested",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for mode, data in summaries.items():
            for comparison in data["mapped_comparison"]:
                writer.writerow({"mode": mode, **{key: comparison.get(key) for key in fields if key != "mode"}})
    print(path)
    for mode, data in summaries.items():
        full = data["mapped_comparison"][-1]
        print(mode, json.dumps({key: full[key] for key in (
            "raw_electron_peak_relative_difference", "raw_electron_peak_rise_relative_difference",
            "common_ag_excess_electron_energy_relative_L1",
            "common_ag_lattice_enthalpy_relative_L1", "common_ag_coupling_energy_relative_L1",
        )}))


if __name__ == "__main__":
    main()
