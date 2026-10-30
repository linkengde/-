#!/usr/bin/env python3
"""Run and audit the Ag70 1x single-temperature Case B model.

Case B uses ideal Ag/Ti3SiC2 thermal contact, heat storage in both phases,
Ag-only source deposition, periodic x/y boundaries, and insulated z faces.
It is independent of the existing Case A run and checkpoints.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


OUT_DIR = Path(__file__).resolve().parent
RUN_DIR = OUT_DIR.parent
PROJECT_DIR = OUT_DIR.parents[2]
PHASE_DIR = PROJECT_DIR / "phase_field_qsgs"
GEOMETRY = RUN_DIR / "qsgs_morphology_ag70vol.npz"
sys.path.insert(0, str(PHASE_DIR))

from qsgs_enthalpy_solver import QSGSEnthalpyGrid, write_cluster_csv  # noqa: E402
from qsgs_liquid_tracking import LiquidClusterTracker  # noqa: E402


TM_K = 1235.08
SC_TMAX_K = 1473.15
SC_WARNING_K = 1400.0
SC_POOL_HEAT_OFF_K = 1450.0
SC_NO_POOL_STOP_K = SC_TMAX_K - 1.0
AG_TMAX_K = 2300.0
LIQUID_THRESHOLD = 0.5
ENERGY_REL_TOL = 1.0e-8


def build_model() -> QSGSEnthalpyGrid:
    if not GEOMETRY.is_file():
        raise FileNotFoundError(f"Missing Ag70 QSGS geometry: {GEOMETRY}")
    return QSGSEnthalpyGrid(
        None,
        "B",
        peak_flux_eV_ps_A2=5.0,
        gaussian_radius_A=30.0,
        absorption_depth_A=20.0,
        native_geometry_npz=GEOMETRY,
        source_geometry_npz=GEOMETRY,
        mesh_refinement=1,
    )


def _volume_integral(field: np.ndarray, model: QSGSEnthalpyGrid) -> float:
    return float(np.sum(field, dtype=np.float64) * model.cell_volume_m3)


def energy_breakdown(model: QSGSEnthalpyGrid) -> dict[str, float]:
    """Return phase-resolved energies without counting latent heat twice."""
    temperature = model.state.temperature_K
    liquid = model.state.ag_liquid_fraction
    ag_mask = model.ag_fraction > 0.0
    sc_mask = model.sc_fraction > 0.0

    ag_temperature = np.where(ag_mask, temperature, 300.0)
    ag_specific_total = np.asarray(
        model.ag.specific_enthalpy_J_kg(ag_temperature, np.where(ag_mask, liquid, 0.0)),
        dtype=float,
    )
    ag_latent_density = (
        model.ag_fraction
        * model.ag.reference_density_kg_m3
        * model.ag.latent_heat_J_kg
        * liquid
    )
    ag_total_density = model.ag_fraction * model.ag.reference_density_kg_m3 * ag_specific_total
    ag_latent_J = _volume_integral(ag_latent_density, model)
    ag_total_J = _volume_integral(ag_total_density, model)
    ag_sensible_J = ag_total_J - ag_latent_J

    sc_temperature = np.where(sc_mask, temperature, 300.0)
    cp0, cp1 = (float(value) for value in model.sc.cp_mass_J_kg_K)
    t0, t1 = (float(value) for value in model.sc.temperature_K)
    slope = (cp1 - cp0) / (t1 - t0)
    delta = sc_temperature - t0
    sc_specific_h = cp0 * delta + 0.5 * slope * delta**2
    sc_specific_h = np.where(sc_temperature < t0, cp0 * delta, sc_specific_h)
    sc_density = model.sc_fraction * model.sc.density_kg_m3 * sc_specific_h
    sc_energy_J = _volume_integral(sc_density, model)
    total_J = _volume_integral(model.enthalpy_density_J_m3, model)
    return {
        "ag_sensible_enthalpy_J": ag_sensible_J,
        "ag_latent_enthalpy_J": ag_latent_J,
        "ag_total_enthalpy_J": ag_total_J,
        "ti3sic2_stored_enthalpy_J": sc_energy_J,
        "system_total_enthalpy_J": total_J,
        "cumulative_external_source_energy_J": float(model.injected_energy_J),
        "energy_balance_residual_J": total_J - float(model.injected_energy_J),
    }


def state_metrics(
    model: QSGSEnthalpyGrid,
    tracker_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    t = model.state.temperature_K
    liquid = model.state.ag_liquid_fraction
    ag_mask = model.ag_fraction > 0.0
    sc_mask = model.sc_fraction > 0.0
    ag_weight = model.ag_fraction
    sc_weight = model.sc_fraction
    ag_den = float(ag_weight.sum())
    sc_den = float(sc_weight.sum())
    metrics: dict[str, Any] = {
        "time_ps": float(model.time_s * 1.0e12),
        "ag_peak_temperature_K": float(t[ag_mask].max()) if np.any(ag_mask) else 300.0,
        "ag_mean_temperature_K": float(np.sum(ag_weight * t) / ag_den) if ag_den else 300.0,
        "ti3sic2_peak_temperature_K": float(t[sc_mask].max()) if np.any(sc_mask) else 300.0,
        "ti3sic2_mean_temperature_K": float(np.sum(sc_weight * t) / sc_den) if sc_den else 300.0,
        "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
        "ag_max_local_liquid_fraction": float(liquid.max()),
        "maximum_timestep_fs": float(model.last_run_step_limit_s * 1.0e15)
        if model.last_run_step_limit_s is not None
        else float(model.stable_timestep_s() * 1.0e15),
    }
    metrics.update(energy_breakdown(model))
    if tracker_summary:
        metrics.update(tracker_summary)
    metrics["largest_liquid_cluster_volume_nm3"] = (
        float(metrics.get("largest_liquid_cluster_volume_A3", 0.0)) / 1000.0
    )
    return metrics


def _relative_error(delta: float, reference: float) -> float:
    return abs(delta) / max(abs(reference), 1.0e-300)


def run_preflight(output_dir: Path) -> dict[str, Any]:
    """Exercise the real Ag70 1x grid, both phases, source mask, and guards."""
    model = build_model()
    ag_volume_fraction = float(model.ag_fraction.mean())
    sc_volume_fraction = float(model.sc_fraction.mean())
    source_in_sc = bool(np.any((model.source_cell_power_W > 0.0) & (model.sc_fraction > 0.0)))
    assert model.boundary_case == "B"
    assert model.temperature_shape == (100, 100, 50)
    assert abs(ag_volume_fraction - 0.699996) < 2.0e-8
    assert abs(sc_volume_fraction - 0.300004) < 2.0e-8
    assert model.source_geometry_npz == GEOMETRY
    assert model.mesh_refinement == 1 and model.tile_xy == 1
    assert model.actual_source_power_W > 0.0 and not source_in_sc
    assert np.all(model.state.temperature_K == 300.0)
    assert np.all(model.state.ag_liquid_fraction == 0.0)

    source_power = float(model.actual_source_power_W)
    stable_dt_s = float(model.stable_timestep_s())
    assert np.isfinite(stable_dt_s) and stable_dt_s > 0.0

    # At 300 K, no source means an exactly static state and zero system energy.
    initial_h = model.enthalpy_density_J_m3.copy()
    static_rows = model.run(0.001, source_on=False, sample_interval_ps=0.001)
    assert np.array_equal(model.enthalpy_density_J_m3, initial_h)
    assert np.all(model.state.temperature_K == 300.0)
    assert model.injected_energy_J == 0.0
    static_energy_J = _volume_integral(model.enthalpy_density_J_m3, model)
    assert static_energy_J == 0.0

    # A uniform two-phase temperature has no interfacial or boundary flux.
    uniform_K = 500.0
    uniform_t = np.full(model.temperature_shape, uniform_K, dtype=float)
    uniform_h = model.enthalpy_density_from_state(uniform_t, np.zeros_like(uniform_t))
    model.enthalpy_density_J_m3 = uniform_h
    model.state = model.state_from_enthalpy_density(uniform_h)
    uniform_energy_before = _volume_integral(uniform_h, model)
    uniform_rows = model.run(0.001, source_on=False, sample_interval_ps=0.001)
    uniform_energy_after = _volume_integral(model.enthalpy_density_J_m3, model)
    uniform_temp_error_K = float(np.max(np.abs(model.state.temperature_K - uniform_K)))
    assert uniform_temp_error_K < 1.0e-8
    assert _relative_error(uniform_energy_after - uniform_energy_before, uniform_energy_before) < 1.0e-12

    # Put Ag at 600 K and the adjoining skeleton at 300 K. A conservative
    # source-off step must warm the scaffold, cool Ag, and preserve total H.
    contrast_t = np.where(model.ag_fraction > 0.0, 600.0, 300.0)
    contrast_h = model.enthalpy_density_from_state(contrast_t, np.zeros_like(contrast_t))
    model.enthalpy_density_J_m3 = contrast_h
    model.state = model.state_from_enthalpy_density(contrast_h)
    energy_before = _volume_integral(contrast_h, model)
    sc_energy_before = energy_breakdown(model)["ti3sic2_stored_enthalpy_J"]
    ag_mean_before = float(np.sum(model.ag_fraction * model.state.temperature_K) / model.ag_fraction.sum())
    sc_max_before = float(model.state.temperature_K[model.sc_fraction > 0.0].max())
    dt = min(model.stable_timestep_s(), 0.01e-15)
    for _ in range(8):
        model.advance(dt, source_on=False, enforce_stability=True)
    energy_after = _volume_integral(model.enthalpy_density_J_m3, model)
    ag_mean_after = float(np.sum(model.ag_fraction * model.state.temperature_K) / model.ag_fraction.sum())
    sc_max_after = float(model.state.temperature_K[model.sc_fraction > 0.0].max())
    sc_energy_after = energy_breakdown(model)["ti3sic2_stored_enthalpy_J"]
    energy_rel_error = _relative_error(energy_after - energy_before, energy_before)
    sc_heat_gain_J = sc_energy_after - sc_energy_before
    assert sc_max_after > sc_max_before
    assert ag_mean_after < ag_mean_before
    assert sc_heat_gain_J > 0.0
    assert energy_rel_error < 1.0e-12

    # Verify the explicit stability guard and the hard Ti3SiC2 data ceiling.
    state_before_guard = model.state.temperature_K.copy()
    with np.testing.assert_raises(ValueError):
        model.advance(model.stable_timestep_s() * 1.01, source_on=False, enforce_stability=True)
    assert np.array_equal(model.state.temperature_K, state_before_guard)
    at_limit_t = np.full(model.temperature_shape, SC_TMAX_K, dtype=float)
    at_limit_h = model.enthalpy_density_from_state(at_limit_t, np.ones_like(at_limit_t))
    state_at_limit = model.state_from_enthalpy_density(at_limit_h)
    assert float(state_at_limit.temperature_K[model.sc_fraction > 0.0].max()) <= SC_TMAX_K + 1.0e-9
    over_limit_t = np.full(model.temperature_shape, SC_TMAX_K + 0.01, dtype=float)
    over_limit_h = model.enthalpy_density_from_state(over_limit_t, np.ones_like(over_limit_t))
    with np.testing.assert_raises_regex(ValueError, "1473.15 K"):
        model.state_from_enthalpy_density(over_limit_h)

    results = {
        "status": "CASE_B_PREFLIGHT_PASSED",
        "boundary_case": "B: ideal Ag/Ti3SiC2 thermal contact; both phases store and conduct heat",
        "geometry": str(GEOMETRY.relative_to(PROJECT_DIR)),
        "geometry_sha256": QSGSEnthalpyGrid._sha256_file(GEOMETRY),
        "grid_shape_xyz": list(model.temperature_shape),
        "mesh_refinement": model.mesh_refinement,
        "spacing_A": model.spacing_A.tolist(),
        "ag_volume_fraction": ag_volume_fraction,
        "ti3sic2_volume_fraction": sc_volume_fraction,
        "initial_temperature_K": 300.0,
        "source_peak_flux_W_m2": model.peak_flux_W_m2,
        "source_gaussian_radius_A": model.gaussian_radius_A,
        "source_absorption_depth_A": model.absorption_depth_A,
        "source_actual_power_W": source_power,
        "source_support_in_scaffold": source_in_sc,
        "stable_timestep_at_300K_fs": stable_dt_s * 1.0e15,
        "static_300K_no_source": {
            "status": "passed",
            "duration_ps": 0.001,
            "energy_J": static_energy_J,
            "sample_rows": len(static_rows),
        },
        "uniform_500K_no_source": {
            "status": "passed",
            "duration_ps": 0.001,
            "max_temperature_drift_K": uniform_temp_error_K,
            "energy_before_J": uniform_energy_before,
            "energy_after_J": uniform_energy_after,
            "sample_rows": len(uniform_rows),
        },
        "ag_to_scaffold_conduction": {
            "status": "passed",
            "initial_ag_temperature_K": 600.0,
            "initial_scaffold_temperature_K": 300.0,
            "steps": 8,
            "step_fs": dt * 1.0e15,
            "ag_mean_before_K": ag_mean_before,
            "ag_mean_after_K": ag_mean_after,
            "scaffold_max_before_K": sc_max_before,
            "scaffold_max_after_K": sc_max_after,
            "scaffold_stored_heat_gain_J": sc_heat_gain_J,
            "system_energy_relative_error": energy_rel_error,
        },
        "explicit_timestep_guard": "passed",
        "ti3sic2_1473_15K_hard_limit": "passed; states above the data ceiling are rejected",
        "candidate_properties": "First-pass candidate values; exact source tables are not archived and remain unverified.",
        "limitations": [
            "ideal interface contact; no finite interface thermal resistance",
            "external faces are insulated except x/y periodicity; no external heat sink",
            "single-temperature continuum model; no electron-lattice nonequilibrium",
            "no liquid flow, melt expansion, vaporization, or mechanical skeleton constraint",
        ],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "caseB_preflight.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    write_preflight_report(output_dir / "CASE_B_MODEL_VALIDATION.md", results)
    return results


def write_preflight_report(path: Path, results: dict[str, Any]) -> None:
    conduction = results["ag_to_scaffold_conduction"]
    text = f"""# Ag70 1× Case B model preflight

## Outcome

**{results['status']}** on the native Ag70 QSGS geometry. This validates the existing Case B implementation for short, source-off conduction tests; it does not establish a melt-pool result.

## Model used

- Geometry: `{results['geometry']}`; SHA-256 `{results['geometry_sha256']}`.
- Grid: `{results['grid_shape_xyz']}`; mesh refinement `{results['mesh_refinement']}×`.
- Ag/Ti₃SiC₂: `{results['ag_volume_fraction']:.6%}` / `{results['ti3sic2_volume_fraction']:.6%}` by actual grid count.
- Initial state: 300 K; x/y periodic; z insulated.
- Case B: both phases store and conduct heat with ideal thermal contact; only Ag receives the source.
- Source: peak flux `{results['source_peak_flux_W_m2']:.6g} W/m²`, Gaussian radius `{results['source_gaussian_radius_A']} Å`, absorption depth `{results['source_absorption_depth_A']} Å`; integrated power `{results['source_actual_power_W']:.9g} W`.
- Explicit diffusion limit at 300 K: `{results['stable_timestep_at_300K_fs']:.6g} fs`.

## Checks

| Check | Result |
| --- | --- |
| Native shape, Ag70 geometry identity, 1× resolution | passed |
| Heat source has no support in Ti₃SiC₂ cells | passed |
| 300 K, no source, 0.001 ps static run | passed; energy `{results['static_300K_no_source']['energy_J']:.4e} J` |
| Uniform 500 K, no source, 0.001 ps | passed; max drift `{results['uniform_500K_no_source']['max_temperature_drift_K']:.3e} K` |
| Hot Ag / cold scaffold transport | passed; scaffold max `{conduction['scaffold_max_before_K']:.6f} → {conduction['scaffold_max_after_K']:.6f} K` |
| Hot Ag cools while scaffold stores heat | passed; Ag mean `{conduction['ag_mean_before_K']:.6f} → {conduction['ag_mean_after_K']:.6f} K`; scaffold gain `{conduction['scaffold_stored_heat_gain_J']:.4e} J` |
| Closed-system energy conservation in conduction test | passed; relative change `{conduction['system_energy_relative_error']:.3e}` |
| Explicit timestep guard | passed |
| Ti₃SiC₂ upper-property limit | passed; states above 1473.15 K rejected |

## Material and physical limits

The solver uses the repository's first-pass Ag and Ti₃SiC₂ candidate properties. The Ti₃SiC₂ table ends at 1473.15 K, and the model rejects states above that limit. These candidates are not fully source-verified. Case B assumes perfect Ag/Ti₃SiC₂ contact and an insulated outer boundary; it is not a simulation of heat loss to the electrodes or ambient. A shared single temperature, no liquid flow, no vaporization, no expansion, and no mechanical constraint are included.

The next numerical gate is a 0.01–0.05 ps Ag-only source-on run from a fresh 300 K Case B state. Longer segmented heating and a source-off cooling branch should proceed only while the temperature/property limits and energy ledger remain valid.
"""
    path.write_text(text, encoding="utf-8")


def _write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    exists = path.exists() and path.stat().st_size > 0
    with path.open("a" if exists else "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        if not exists:
            writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value) if isinstance(value, (list, dict)) else value for key, value in row.items()})


def _load_or_initialize(
    model: QSGSEnthalpyGrid,
    output_dir: Path,
    restart_from: Path | None,
) -> tuple[LiquidClusterTracker, Path]:
    checkpoints = output_dir / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    tracker = LiquidClusterTracker(model.spacing_A, model.size_A, liquid_threshold=LIQUID_THRESHOLD)
    if restart_from is not None:
        info = model.load_checkpoint(restart_from.resolve(), tracker)
        if model.boundary_case != "B" or model.mesh_refinement != 1:
            raise ValueError("This driver accepts only Ag70 Case B 1x restart checkpoints")
        return tracker, restart_from.resolve()
    initial_checkpoint = checkpoints / "caseB_ag70_1x_initial_0ps.npz"
    if initial_checkpoint.exists():
        raise FileExistsError(f"Refusing to overwrite existing initial state: {initial_checkpoint}")
    tracker.update(0.0, model.state.ag_liquid_fraction, model.ag_fraction > 0.0)
    model.save_checkpoint(initial_checkpoint, tracker)
    return tracker, initial_checkpoint


def _choose_segment_ps(model: QSGSEnthalpyGrid, requested: float) -> float:
    max_sc = float(model.state.temperature_K[model.sc_fraction > 0.0].max())
    segment = requested
    if max_sc >= SC_WARNING_K:
        segment = min(segment, 0.01)
    if max_sc >= SC_POOL_HEAT_OFF_K:
        segment = min(segment, 0.005)
    if max_sc >= SC_NO_POOL_STOP_K:
        segment = min(segment, 0.001)
    return segment


def _choose_sample_ps(model: QSGSEnthalpyGrid, requested: float) -> float:
    """Increase field/cluster sampling near melt without writing a full checkpoint each sample."""
    max_ag = float(model.state.temperature_K[model.ag_fraction > 0.0].max())
    max_sc = float(model.state.temperature_K[model.sc_fraction > 0.0].max())
    sample = requested
    has_liquid = bool(np.any(model.state.ag_liquid_fraction[model.ag_fraction > 0.0] > 0.0))
    if max_ag >= TM_K - 150.0 or max_sc >= SC_WARNING_K or has_liquid:
        sample = min(sample, 0.002)
    if max_ag >= TM_K - 30.0 or max_sc >= SC_POOL_HEAT_OFF_K or has_liquid:
        sample = min(sample, 0.001)
    return sample


def dynamic_heat_stop_reason(metrics: dict[str, Any]) -> str | None:
    """Return event-driven heat-off reasons without imposing a pool lifetime."""
    sc_peak = float(metrics["ti3sic2_peak_temperature_K"])
    has_connected_pool = (
        int(metrics.get("connected_liquid_cluster_count", 0)) > 0
        and int(metrics.get("largest_liquid_cluster_voxels", 0)) >= 2
    )
    if has_connected_pool and sc_peak >= SC_POOL_HEAT_OFF_K:
        return "connected_pool_scaffold_temperature_margin_heat_off"
    if not has_connected_pool and sc_peak >= SC_NO_POOL_STOP_K:
        return "scaffold_property_ceiling_safety_stop_without_pool"
    return None


def run_transient(
    mode: str,
    duration_ps: float,
    segment_ps: float,
    sample_ps: float,
    output_dir: Path,
    restart_from: Path | None,
    pool_followup_ps: float = 0.30,
) -> dict[str, Any]:
    if mode not in {"heat", "cool"}:
        raise ValueError("mode must be heat or cool")
    if duration_ps <= 0.0 or segment_ps <= 0.0 or sample_ps <= 0.0:
        raise ValueError("duration, segment and sample intervals must be positive")
    if mode == "heat":
        validation_path = output_dir / "caseB_preflight.json"
        if not validation_path.is_file():
            raise RuntimeError("Run --mode validate successfully before starting Case B heating")
        validation = json.loads(validation_path.read_text(encoding="utf-8"))
        if validation.get("status") != "CASE_B_PREFLIGHT_PASSED":
            raise RuntimeError("Case B preflight did not pass; heating is blocked")
    model = build_model()
    tracker, last_valid_checkpoint = _load_or_initialize(model, output_dir, restart_from)
    if mode == "cool" and not np.any(model.state.ag_liquid_fraction > 0.0):
        raise ValueError("Cooling branch requires a checkpoint containing liquid Ag")

    start_time_ps = float(model.time_s * 1.0e12)
    e_off = _volume_integral(model.enthalpy_density_J_m3, model)
    e_s = _volume_integral(model._solid_at_tm_J_m3, model)
    e_l = e_s + _volume_integral(model._latent_volumetric_J_m3, model)
    energy_decision = equilibrium_energy_assessment(model, e_off) if mode == "cool" else None
    if mode == "cool" and not energy_decision["all_solid_energetically_possible"]:
        result = {
            "status": "COOLING_SKIPPED_ALL_SOLID_ENERGETICALLY_IMPOSSIBLE",
            "reason": "E_off exceeds the all-Ag-solid-at-Tm enthalpy under insulated Case B assumptions",
            "checkpoint": str(last_valid_checkpoint),
            "energy_assessment": energy_decision,
            "last_state": state_metrics(model),
        }
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "caseB_cooling_energy_assessment.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        return result

    started = time.perf_counter()
    final_time_ps = model.time_s * 1.0e12 + duration_ps
    pool_start = tracker.first_multivoxel_pool_time_ps
    rows: list[dict[str, Any]] = []
    stop_reason = "duration_completed"
    cooling_all_solid = False
    while model.time_s * 1.0e12 < final_time_ps - 1.0e-12:
        current_ps = model.time_s * 1.0e12
        segment = min(_choose_segment_ps(model, segment_ps), final_time_ps - current_ps)
        end_ps = current_ps + segment
        prefix = "caseB_ag70_1x_heat" if mode == "heat" else "caseB_ag70_1x_cool"
        try:
            sample_rows = model.run(
                segment,
                source_on=(mode == "heat"),
                sample_interval_ps=min(_choose_sample_ps(model, sample_ps), segment),
                track_liquid_connectivity=True,
                liquid_threshold=LIQUID_THRESHOLD,
                tracker=tracker,
                checkpoint_times_ps=(end_ps,),
                checkpoint_dir=output_dir / "checkpoints",
                checkpoint_prefix=prefix,
            )
        except ValueError as exc:
            message = str(exc)
            if "Ti3SiC2 candidate table" not in message and "1473.15 K" not in message:
                raise
            # The solver may have tentatively advanced the last step before
            # detecting an out-of-range scaffold state. Restore the previous
            # validated segment checkpoint; never keep or save that trial.
            model.load_checkpoint(last_valid_checkpoint, tracker)
            stop_reason = "ti3sic2_candidate_property_limit_reached"
            (output_dir / "property_limit_event.txt").write_text(
                f"Stopped before accepting an out-of-range state at {message}\n"
                f"Last complete validated checkpoint: {last_valid_checkpoint}\n",
                encoding="utf-8",
            )
            break

        if not sample_rows:
            raise RuntimeError("Transient solver returned no rows for a positive segment")
        row = state_metrics(model, sample_rows[-1])
        row["mode"] = mode
        row["segment_start_ps"] = current_ps
        row["segment_end_ps"] = end_ps
        residual = float(row["energy_balance_residual_J"])
        reference = max(abs(float(row["system_total_enthalpy_J"])), abs(float(row["cumulative_external_source_energy_J"])))
        relative_error = abs(residual) / max(reference, 1.0e-300)
        row["relative_energy_balance_error"] = relative_error
        if not np.isfinite(relative_error) or relative_error > ENERGY_REL_TOL:
            raise ArithmeticError(f"Energy balance failed at {end_ps:.9g} ps: {relative_error:.3e}")
        if not np.isfinite(row["ti3sic2_peak_temperature_K"]) or not np.isfinite(row["ag_peak_temperature_K"]):
            raise FloatingPointError("Non-finite phase temperature")
        if row["ti3sic2_peak_temperature_K"] > SC_TMAX_K + 1.0e-9:
            raise ArithmeticError("Scaffold temperature exceeded its candidate-property ceiling")

        checkpoint = output_dir / "checkpoints" / f"{prefix}_{_time_tag(end_ps)}.npz"
        if not checkpoint.is_file():
            raise FileNotFoundError(f"Expected complete segment checkpoint not found: {checkpoint}")
        row["checkpoint"] = str(checkpoint.relative_to(output_dir))
        rows.append(row)
        _write_rows(output_dir / f"{mode}_metrics.csv", [row])
        last_valid_checkpoint = checkpoint
        if tracker.history:
            write_cluster_csv(output_dir / f"{mode}_liquid_clusters.csv", tracker.history)

        if mode == "heat":
            if row["ag_peak_temperature_K"] >= AG_TMAX_K:
                stop_reason = "ag_2300K_numerical_guard_reached"
                break
            if row["ti3sic2_peak_temperature_K"] >= SC_TMAX_K - 1.0e-8:
                stop_reason = "ti3sic2_candidate_property_limit_reached"
                break
            event_stop = dynamic_heat_stop_reason(row)
            if event_stop is not None:
                stop_reason = event_stop
                break
            if tracker.first_multivoxel_pool_time_ps is not None:
                pool_start = tracker.first_multivoxel_pool_time_ps
                if end_ps >= pool_start + pool_followup_ps - 1.0e-10:
                    stop_reason = "multivoxel_pool_followup_completed"
                    break
        else:
            if model._liquid_ag_volume_fraction() == 0.0 and float(model.state.ag_liquid_fraction.max()) == 0.0:
                cooling_all_solid = True
                stop_reason = "all_ag_solid"
                break

    final_ps = model.time_s * 1.0e12
    final_cp = output_dir / "checkpoints" / f"caseB_ag70_1x_{mode}_final_{_time_tag(final_ps)}.npz"
    model.save_checkpoint(final_cp, tracker)
    if tracker.history:
        write_cluster_csv(output_dir / f"{mode}_liquid_clusters.csv", tracker.history)
    # Keep the tracked component metrics from the last accepted segment in the
    # final summary. Without this, the checkpoint summary incorrectly reported
    # a zero-size cluster even when the final validated frame contained liquid.
    final_metrics = state_metrics(model, rows[-1] if rows else None)
    result = {
        "mode": mode,
        "status": f"AG70_CASE_B_{mode.upper()}_SEGMENTED_RUN_COMPLETED",
        "boundary_case": "B: ideal contact; Ag and Ti3SiC2 store and conduct heat",
        "geometry_sha256": QSGSEnthalpyGrid._sha256_file(GEOMETRY),
        "source_on": mode == "heat",
        "restart_checkpoint": str(restart_from.resolve()) if restart_from else None,
        "final_checkpoint": str(final_cp.relative_to(output_dir)),
        "start_time_ps": start_time_ps,
        "time_end_ps": final_ps,
        "duration_requested_ps": duration_ps,
        "elapsed_wall_time_s": time.perf_counter() - started,
        "segments_completed": len(rows),
        "stop_reason": stop_reason,
        "first_liquid_time_ps": tracker.first_liquid_time_ps,
        "first_connected_liquid_time_ps": tracker.first_connected_liquid_time_ps,
        "first_multivoxel_pool_time_ps": tracker.first_multivoxel_pool_time_ps,
        "largest_cluster_sampled_continuous_duration_ps": max(
            (float(item["continuous_duration_ps"]) for item in tracker.history if item.get("is_main_cluster")),
            default=0.0,
        ),
        "cooling_reached_all_solid": cooling_all_solid,
        "heat_off_energy_assessment": equilibrium_energy_assessment(model, e_off) if mode == "cool" else None,
        "final_metrics": final_metrics,
        "candidate_property_limits": {"Ag_K": AG_TMAX_K, "Ti3SiC2_K": SC_TMAX_K},
        "physical_limitations": [
            "candidate thermophysical properties remain incompletely sourced",
            "Case B assumes ideal Ag/Ti3SiC2 contact and insulated external boundaries",
            "this single-temperature model has no melt flow, vaporization, or mechanical skeleton response",
        ],
    }
    (output_dir / f"caseB_{mode}_summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    write_transient_report(output_dir / f"CASE_B_{mode.upper()}_STATUS.md", result, rows)
    return result


def _time_tag(time_ps: float) -> str:
    return f"{time_ps:.9f}".rstrip("0").rstrip(".").replace(".", "p") + "ps"


def equilibrium_energy_assessment(model: QSGSEnthalpyGrid, e_off_J: float) -> dict[str, Any]:
    e_s = _volume_integral(model._solid_at_tm_J_m3, model)
    e_latent = _volume_integral(model._latent_volumetric_J_m3, model)
    e_l = e_s + e_latent
    tolerance = max(abs(e_off_J), abs(e_s), 1.0e-300) * 1.0e-12
    if e_off_J <= e_s + tolerance:
        lo, hi = 300.0, TM_K
        for _ in range(64):
            mid = 0.5 * (lo + hi)
            h = model.enthalpy_density_from_state(
                np.full(model.temperature_shape, mid),
                np.zeros(model.temperature_shape),
            )
            energy = _volume_integral(h, model)
            if energy < e_off_J:
                lo = mid
            else:
                hi = mid
        temp = 0.5 * (lo + hi)
        return {
            "E_off_J": e_off_J,
            "E_s_all_ag_solid_at_tm_J": e_s,
            "E_l_all_ag_liquid_at_tm_J": e_l,
            "all_solid_energetically_possible": True,
            "uniform_equilibrium_temperature_K": temp,
            "uniform_equilibrium_ag_liquid_fraction": 0.0,
            "classification": "E_off <= E_s; all-solid uniform equilibrium is energetically possible",
        }
    if e_off_J <= e_l + tolerance:
        return {
            "E_off_J": e_off_J,
            "E_s_all_ag_solid_at_tm_J": e_s,
            "E_l_all_ag_liquid_at_tm_J": e_l,
            "all_solid_energetically_possible": False,
            "uniform_equilibrium_temperature_K": TM_K,
            "uniform_equilibrium_ag_liquid_fraction": (e_off_J - e_s) / e_latent,
            "classification": "E_s < E_off <= E_l; residual liquid remains at a uniform-Tm equilibrium",
        }
    h_limit = model.enthalpy_density_from_state(
        np.full(model.temperature_shape, SC_TMAX_K),
        np.ones(model.temperature_shape),
    )
    e_limit = _volume_integral(h_limit, model)
    if e_off_J > e_limit:
        equilibrium_temperature = None
        classification = "E_off exceeds the enthalpy representable within the Ti3SiC2 candidate range"
    else:
        lo, hi = TM_K, SC_TMAX_K
        for _ in range(64):
            mid = 0.5 * (lo + hi)
            h = model.enthalpy_density_from_state(
                np.full(model.temperature_shape, mid),
                np.ones(model.temperature_shape),
            )
            if _volume_integral(h, model) < e_off_J:
                lo = mid
            else:
                hi = mid
        equilibrium_temperature = 0.5 * (lo + hi)
        classification = "E_off > E_l; all-Ag-liquid equilibrium within the Ti3SiC2 property range"
    return {
        "E_off_J": e_off_J,
        "E_s_all_ag_solid_at_tm_J": e_s,
        "E_l_all_ag_liquid_at_tm_J": e_l,
        "all_solid_energetically_possible": False,
        "uniform_equilibrium_temperature_K": equilibrium_temperature,
        "uniform_equilibrium_ag_liquid_fraction": 1.0 if equilibrium_temperature is not None else None,
        "classification": classification,
        "max_enthalpy_at_ti3sic2_property_limit_J": e_limit,
    }


def write_transient_report(path: Path, result: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    latest = result["final_metrics"]
    text = f"""# Ag70 1× Case B {result['stop_reason'].replace('_', ' ')}

## Run outcome

- Status: `{result['status']}`
- Source on during run: `{result['source_on']}`
- Stop reason: `{result['stop_reason']}`
- End time: `{result['time_end_ps']:.9g} ps`; completed segments: {result['segments_completed']}.
- Final valid checkpoint: `{result['final_checkpoint']}`
- Final Ag peak/mean temperature: {latest['ag_peak_temperature_K']:.4f} / {latest['ag_mean_temperature_K']:.4f} K.
- Final Ti₃SiC₂ peak/mean temperature: {latest['ti3sic2_peak_temperature_K']:.4f} / {latest['ti3sic2_mean_temperature_K']:.4f} K.
- Final Ag liquid volume fraction: {latest['ag_liquid_volume_fraction']:.8g}; largest sampled connected cluster: {latest['largest_liquid_cluster_volume_nm3']:.6g} nm³.
- Final system enthalpy: {latest['system_total_enthalpy_J']:.12e} J; source ledger: {latest['cumulative_external_source_energy_J']:.12e} J; residual: {latest['energy_balance_residual_J']:.4e} J.

## Sampled segment history

`{result['mode']}_metrics.csv` records phase-resolved temperature, Ag liquid fraction, connected-cluster metrics, phase enthalpy, source energy, and energy residual at completed segment boundaries. Connectivity persistence is sampled, not a proof of hydrodynamic or mechanical stability.

## Interpretation limits

The calculation uses first-pass candidate properties, ideal Ag/Ti₃SiC₂ contact, and an insulated outer boundary. It includes no finite interface resistance, external heat sink, melt flow, vaporization, thermal expansion, or mechanical skeleton response. A source-off Case B result describes internal redistribution in the composite domain, not heat rejection into an electrode or the environment.
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("validate", "heat", "cool"), required=True)
    parser.add_argument("--duration-ps", type=float, default=0.05,
                        help="Transient duration for heat/cool; heating starts from a fresh 300 K state unless restarted")
    parser.add_argument("--segment-ps", type=float, default=0.01)
    parser.add_argument("--sample-ps", type=float, default=0.005)
    parser.add_argument("--restart-from", type=Path, default=None)
    parser.add_argument("--pool-followup-ps", type=float, default=0.30)
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if args.mode == "validate":
        result = run_preflight(output_dir)
    else:
        result = run_transient(
            args.mode,
            args.duration_ps,
            args.segment_ps,
            args.sample_ps,
            output_dir,
            args.restart_from,
            args.pool_followup_ps,
        )
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
