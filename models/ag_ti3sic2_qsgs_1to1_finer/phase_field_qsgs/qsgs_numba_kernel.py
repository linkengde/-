"""Optional fused Numba kernels for large native-grid enthalpy runs."""

from __future__ import annotations

from math import log

import numpy as np
from numba import njit, prange


@njit(inline="always")
def _sc_cp(t, cp_a, cp_b, molar_mass, reference):
    effective = max(t, reference)
    return (cp_a - cp_b / effective) / molar_mass


@njit(inline="always")
def _sc_h(t, cp_a, cp_b, molar_mass, reference):
    if t < reference:
        return _sc_cp(reference, cp_a, cp_b, molar_mass, reference) * (t - reference)
    return (cp_a * (t - reference) - cp_b * log(t / reference)) / molar_mass


@njit(inline="always")
def _ag_solid_h(t, cp0, cp_slope, reference, melting):
    if t < reference:
        return cp0 * (t - reference)
    delta = min(t, melting) - reference
    return cp0 * delta + 0.5 * cp_slope * delta * delta


@njit(inline="always")
def _ag_solid_cp(t, cp0, cp_slope, reference, melting):
    return cp0 + cp_slope * (min(max(t, reference), melting) - reference)


@njit(inline="always")
def _invert_enthalpy(
    enthalpy,
    ag_fraction,
    sc_fraction,
    ag_rho,
    ag_cp0,
    ag_cp_slope,
    ag_liquid_cp,
    ag_latent,
    ag_reference,
    melting,
    sc_rho,
    sc_cp_a,
    sc_cp_b,
    sc_molar_mass,
    sc_reference,
    sc_maximum,
):
    sc_cp_ref = _sc_cp(sc_reference, sc_cp_a, sc_cp_b, sc_molar_mass, sc_reference)
    capacity_300 = ag_fraction * ag_rho * ag_cp0 + sc_fraction * sc_rho * sc_cp_ref
    if enthalpy < 0.0:
        if capacity_300 > 0.0:
            return 300.0 + enthalpy / capacity_300, 0.0, False
        return 300.0, 0.0, False

    ag_slope = ag_cp_slope
    ag_h_tm = _ag_solid_h(melting, ag_cp0, ag_slope, ag_reference, melting)
    sc_h_tm = _sc_h(melting, sc_cp_a, sc_cp_b, sc_molar_mass, sc_reference)
    solid_at_tm = ag_fraction * ag_rho * ag_h_tm + sc_fraction * sc_rho * sc_h_tm
    latent_vol = ag_fraction * ag_rho * ag_latent

    if enthalpy < solid_at_tm:
        temperature = 300.0 + enthalpy / max(capacity_300, 1.0e-300)
        temperature = min(max(temperature, 300.0), melting)
        for _ in range(8):
            h_ag = ag_fraction * ag_rho * _ag_solid_h(
                temperature, ag_cp0, ag_slope, ag_reference, melting
            )
            c_ag = ag_fraction * ag_rho * _ag_solid_cp(
                temperature, ag_cp0, ag_slope, ag_reference, melting
            )
            h_sc = sc_fraction * sc_rho * _sc_h(
                temperature, sc_cp_a, sc_cp_b, sc_molar_mass, sc_reference
            )
            c_sc = sc_fraction * sc_rho * _sc_cp(
                temperature, sc_cp_a, sc_cp_b, sc_molar_mass, sc_reference
            )
            capacity = c_ag + c_sc
            if capacity > 0.0:
                temperature -= (h_ag + h_sc - enthalpy) / capacity
            temperature = min(max(temperature, 300.0), melting)
        return temperature, 0.0, False

    liquidus = solid_at_tm + latent_vol
    if latent_vol > 0.0 and enthalpy <= liquidus:
        liquid_fraction = min(max((enthalpy - solid_at_tm) / latent_vol, 0.0), 1.0)
        return melting, liquid_fraction, False

    excess = enthalpy - liquidus
    sc_cp_max = _sc_cp(sc_maximum, sc_cp_a, sc_cp_b, sc_molar_mass, sc_reference)
    capacity_guess = ag_fraction * ag_rho * ag_liquid_cp + sc_fraction * sc_rho * sc_cp_max
    if capacity_guess <= 0.0:
        return melting, 0.0, False
    temperature = melting + excess / capacity_guess
    for _ in range(8):
        sc_temperature = min(temperature, sc_maximum)
        ag_increment = ag_fraction * ag_rho * ag_liquid_cp * (temperature - melting)
        sc_increment = sc_fraction * sc_rho * (
            _sc_h(sc_temperature, sc_cp_a, sc_cp_b, sc_molar_mass, sc_reference) - sc_h_tm
        )
        capacity = ag_fraction * ag_rho * ag_liquid_cp + sc_fraction * sc_rho * _sc_cp(
            sc_temperature, sc_cp_a, sc_cp_b, sc_molar_mass, sc_reference
        )
        if capacity > 0.0:
            temperature -= (ag_increment + sc_increment - excess) / capacity
    liquid_fraction = 1.0 if ag_fraction > 0.0 else 0.0
    over_limit = sc_fraction > 0.0 and temperature > sc_maximum + 1.0e-9
    return temperature, liquid_fraction, over_limit


@njit(parallel=True)
def _conductivity_field(
    temperature,
    liquid_fraction,
    ag_fraction,
    sc_fraction,
    ag_k0,
    ag_k_tm,
    ag_liquid_k,
    ag_reference,
    melting,
    sc_k_intercept,
    sc_k_slope,
    result,
):
    size = temperature.size
    for n in prange(size):
        t = temperature[n]
        if t < ag_reference:
            ag_k = ag_k0
        elif t < melting:
            ag_k = ag_k0 + (ag_k_tm - ag_k0) * (t - ag_reference) / (melting - ag_reference)
        elif t == melting:
            ag_k = (1.0 - liquid_fraction[n]) * ag_k_tm + liquid_fraction[n] * ag_liquid_k
        else:
            ag_k = ag_liquid_k
        sc_k = sc_k_intercept + sc_k_slope * t
        result[n] = ag_fraction[n] * ag_k + sc_fraction[n] * sc_k


@njit(parallel=True)
def _power_field(
    temperature,
    conductivity,
    enthalpy,
    ag_fraction,
    sc_fraction,
    source_power_density,
    sink_contact,
    net_power_density,
    sink_power_density,
    event_dt_cells,
    nx,
    ny,
    nz,
    dx,
    dy,
    dz,
    sink_temperature,
    source_on,
    event_enthalpy,
    event_enabled,
):
    stride_x = ny * nz
    stride_y = nz
    inverse_dx2 = 1.0 / (dx * dx)
    inverse_dy2 = 1.0 / (dy * dy)
    inverse_dz2 = 1.0 / (dz * dz)
    for n in prange(temperature.size):
        i = n // stride_x
        j = (n - i * stride_x) // stride_y
        k = n - i * stride_x - j * stride_y
        t = temperature[n]
        ki = conductivity[n]
        diffusion = 0.0

        left = n - stride_x if i > 0 else n + (nx - 1) * stride_x
        right = n + stride_x if i + 1 < nx else n - (nx - 1) * stride_x
        kl = conductivity[left]
        kr = conductivity[right]
        face_l = 2.0 * ki * kl / (ki + kl) if ki + kl > 0.0 else 0.0
        face_r = 2.0 * ki * kr / (ki + kr) if ki + kr > 0.0 else 0.0
        diffusion += face_l * (temperature[left] - t) * inverse_dx2
        diffusion += face_r * (temperature[right] - t) * inverse_dx2

        left = n - stride_y if j > 0 else n + (ny - 1) * stride_y
        right = n + stride_y if j + 1 < ny else n - (ny - 1) * stride_y
        kl = conductivity[left]
        kr = conductivity[right]
        face_l = 2.0 * ki * kl / (ki + kl) if ki + kl > 0.0 else 0.0
        face_r = 2.0 * ki * kr / (ki + kr) if ki + kr > 0.0 else 0.0
        diffusion += face_l * (temperature[left] - t) * inverse_dy2
        diffusion += face_r * (temperature[right] - t) * inverse_dy2

        if k > 0:
            left = n - 1
            kl = conductivity[left]
            face = 2.0 * ki * kl / (ki + kl) if ki + kl > 0.0 else 0.0
            diffusion += face * (temperature[left] - t) * inverse_dz2
        if k + 1 < nz:
            right = n + 1
            kr = conductivity[right]
            face = 2.0 * ki * kr / (ki + kr) if ki + kr > 0.0 else 0.0
            diffusion += face * (temperature[right] - t) * inverse_dz2

        sink_density = 0.0
        if sink_contact[n]:
            sink_density = 2.0 * ki * (t - sink_temperature) * inverse_dz2
        sink_power_density[n] = sink_density
        source_density = source_power_density[n] if source_on else 0.0
        q = diffusion + source_density - sink_density
        net_power_density[n] = q
        if event_enabled:
            event_dt_cells[n] = (
                (event_enthalpy[n] - enthalpy[n]) / q
                if sc_fraction[n] > 0.0 and q > 0.0
                else np.inf
            )


@njit(parallel=True)
def _advance_enthalpy(
    enthalpy,
    temperature,
    liquid_fraction,
    net_power_density,
    dt,
    ag_fraction,
    sc_fraction,
    ag_rho,
    ag_cp0,
    ag_cp_slope,
    ag_liquid_cp,
    ag_latent,
    ag_reference,
    melting,
    sc_rho,
    sc_cp_a,
    sc_cp_b,
    sc_molar_mass,
    sc_reference,
    sc_maximum,
):
    for n in prange(enthalpy.size):
        h_new = enthalpy[n] + dt * net_power_density[n]
        enthalpy[n] = h_new
        t, liquid, _ = _invert_enthalpy(
            h_new,
            ag_fraction[n],
            sc_fraction[n],
            ag_rho,
            ag_cp0,
            ag_cp_slope,
            ag_liquid_cp,
            ag_latent,
            ag_reference,
            melting,
            sc_rho,
            sc_cp_a,
            sc_cp_b,
            sc_molar_mass,
            sc_reference,
            sc_maximum,
        )
        temperature[n] = t
        liquid_fraction[n] = liquid


def advance_grid(
    model,
    dt_s: float,
    source_on: bool,
    stop_at_sc_temperature_K: float | None,
    collect_metrics: bool = True,
):
    """Advance a QSGSEnhthalpyGrid using fused native NumPy buffers."""
    ag = model.ag
    sc = model.sc
    state = model.state
    temperature = state.temperature_K
    liquid_fraction = state.ag_liquid_fraction
    workspace = getattr(model, "_numba_workspace", None)
    if workspace is None:
        workspace = {
            "conductivity": np.empty_like(temperature),
            "power": np.empty_like(temperature),
            "sink_power_density": np.empty_like(temperature),
            "event_dt_cells": np.empty_like(temperature),
        }
        model._numba_workspace = workspace
    conductivity = workspace["conductivity"]
    power = workspace["power"]
    sink_power_density = workspace["sink_power_density"]
    event_dt_cells = workspace["event_dt_cells"]
    temp_flat = temperature.ravel()
    liquid_flat = liquid_fraction.ravel()
    ag_fraction = model.ag_fraction.ravel()
    sc_fraction = model.sc_fraction.ravel()
    conductivity_flat = conductivity.ravel()
    _conductivity_field(
        temp_flat,
        liquid_flat,
        ag_fraction,
        sc_fraction,
        float(ag.solid_k_W_m_K[0]),
        float(ag.solid_k_W_m_K[-1]),
        float(ag.liquid_k_W_m_K[0]),
        float(ag.reference_temperature_K),
        float(ag.melting_temperature_K),
        float(sc.k_intercept_W_m_K),
        float(sc.k_slope_W_m_K2),
        conductivity_flat,
    )
    event_enabled = stop_at_sc_temperature_K is not None
    if event_enabled:
        event_enthalpy = model._temperature_event_enthalpy_cache.get(stop_at_sc_temperature_K)
        if event_enthalpy is None:
            ag_h = float(ag.specific_enthalpy_J_kg(stop_at_sc_temperature_K, 1.0))
            sc_h = float(sc.specific_enthalpy_J_kg(stop_at_sc_temperature_K))
            event_enthalpy = (
                model.ag_fraction * ag.reference_density_kg_m3 * ag_h
                + model.sc_fraction * sc.density_kg_m3 * sc_h
            )
            model._temperature_event_enthalpy_cache[stop_at_sc_temperature_K] = event_enthalpy
        event_enthalpy_flat = event_enthalpy.ravel()
    else:
        event_enthalpy_flat = np.zeros(1, dtype=float)
    _power_field(
        temp_flat,
        conductivity_flat,
        model.enthalpy_density_J_m3.ravel(),
        ag_fraction,
        sc_fraction,
        model.source_W_m3.ravel(),
        model.bottom_sink_contact_mask.ravel(),
        power.ravel(),
        sink_power_density.ravel(),
        event_dt_cells.ravel(),
        int(temperature.shape[0]),
        int(temperature.shape[1]),
        int(temperature.shape[2]),
        float(model.spacing_A[0] * 1.0e-10),
        float(model.spacing_A[1] * 1.0e-10),
        float(model.spacing_A[2] * 1.0e-10),
        float(model.bottom_sink_temperature_K or 300.0),
        bool(source_on),
        event_enthalpy_flat,
        bool(event_enabled),
    )
    minimum_event_dt = float(event_dt_cells.min()) if event_enabled else np.inf
    sink_power_w = float(sink_power_density.sum(dtype=np.float64) * model.cell_volume_m3)
    event_time_step = bool(
        event_enabled and minimum_event_dt <= dt_s * (1.0 + 1.0e-12)
    )
    if event_time_step:
        dt_s = max(0.0, min(dt_s, float(minimum_event_dt)))
    ag_cp_slope = (ag.solid_cp_J_kg_K[-1] - ag.solid_cp_J_kg_K[0]) / (
        ag.solid_temperature_K[-1] - ag.solid_temperature_K[0]
    )
    _advance_enthalpy(
        model.enthalpy_density_J_m3.ravel(),
        temp_flat,
        liquid_flat,
        power.ravel(),
        float(dt_s),
        ag_fraction,
        sc_fraction,
        float(ag.reference_density_kg_m3),
        float(ag.solid_cp_J_kg_K[0]),
        float(ag_cp_slope),
        float(ag.liquid_cp_J_kg_K[0]),
        float(ag.latent_heat_J_kg),
        float(ag.reference_temperature_K),
        float(ag.melting_temperature_K),
        float(sc.density_kg_m3),
        float(sc.cp_a_J_mol_K),
        float(sc.cp_b_J_mol),
        float(sc.molar_mass_kg_mol),
        float(sc.reference_temperature_K),
        float(sc.maximum_temperature_K),
    )
    max_t = float(temperature.max())
    max_sc_t = float(temperature[model.sc_fraction > 0.0].max())
    if max_sc_t > sc.maximum_temperature_K + 1.0e-9:
        raise ValueError(
            f"Ti3SiC2 property model upper limit is {sc.maximum_temperature_K:.2f} K"
        )
    model.time_s += dt_s
    if source_on:
        model.injected_energy_J += model.actual_source_power_W * dt_s
    model.removed_energy_J += sink_power_w * dt_s
    model.max_temperature_seen_K = max(model.max_temperature_seen_K, float(max_t))
    event_time_step = event_time_step or bool(
        event_enabled and max_sc_t >= stop_at_sc_temperature_K - 1.0e-8
    )
    row = {
        "time_ps": model.time_s * 1.0e12,
        "dt_fs": dt_s * 1.0e15,
        "max_temperature_K": float(max_t),
        "max_scaffold_temperature_K": float(max_sc_t),
        "injected_energy_J": model.injected_energy_J,
        "bottom_sink_power_W": float(sink_power_w),
        "removed_energy_J": model.removed_energy_J,
        "sc_temperature_limit_reached": 1.0 if event_time_step else 0.0,
    }
    if collect_metrics:
        total_h = float(model.enthalpy_density_J_m3.sum(dtype=np.float64))
        total_energy = total_h * model.cell_volume_m3
        row.update(
            {
                "max_ag_temperature_K": float(temperature[model.ag_fraction > 0.0].max()),
                "mean_ag_temperature_K": model._ag_weighted_mean_temperature(),
                "max_ag_liquid_fraction": float(liquid_fraction.max()),
                "ag_liquid_volume_fraction": model._liquid_ag_volume_fraction(),
                "total_energy_J": total_energy,
                "energy_balance_error_J": total_energy - model.injected_energy_J + model.removed_energy_J,
            }
        )
    return row
