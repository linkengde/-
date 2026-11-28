"""Numba CPU kernels for 3-D mixed-phase single-temperature enthalpy FVM."""

from __future__ import annotations

from math import log

import numpy as np

try:
    from numba import njit, prange
except ImportError:  # pragma: no cover - exercised in environments without optional Numba
    njit = None
    prange = range


HAVE_NUMBA = njit is not None


def require_numba() -> None:
    if not HAVE_NUMBA:
        raise RuntimeError("Numba is optional; install solver_optimization/requirements-numba.txt")


if HAVE_NUMBA:

    @njit(inline="always")
    def _fraction(t, ts, tl, width):
        if t <= ts:
            return 0.0
        if t >= tl:
            return 1.0
        return (t - ts) / width

    @njit(inline="always")
    def _ag_cp_solid(t, t0, t1, cp0, cp1):
        if t <= t0:
            return cp0
        if t >= t1:
            return cp1
        return cp0 + (cp1 - cp0) * (t - t0) / (t1 - t0)

    @njit(inline="always")
    def _ag_h_solid(t, t0, t1, cp0, cp1):
        if t < t0:
            return cp0 * (t - t0)
        if t > t1:
            t = t1
        d = t - t0
        slope = (cp1 - cp0) / (t1 - t0)
        return cp0 * d + 0.5 * slope * d * d

    @njit(inline="always")
    def _ti_cp(t, a, b, molar_mass, tref, cp_ref):
        if t < tref:
            return cp_ref
        return (a - b / t) / molar_mass

    @njit(inline="always")
    def _ti_h(t, a, b, molar_mass, tref, cp_ref):
        if t < tref:
            return cp_ref * (t - tref)
        return (a * (t - tref) - b * log(t / tref)) / molar_mass

    @njit(inline="always")
    def _cell_fields(t, ag, p):
        # p layout is assembled once from the audited CandidateMaterials object.
        t0, t1, cp0, cp1, k0, k1, rho_ag, cp_liq, k_liq, tm, ts, tl, width, cp_mid, cp_delta, h_solid_ts, h_liquid_tl, latent, rho_ti, ti_molar, ti_a, ti_b, tref, ti_cp_ref, ti_k0, ti_kslope, ti_tmax = p
        if t > ti_tmax and ag < 1.0 - 1e-14:
            return (float("nan"), float("nan"), float("nan"))
        f = _fraction(t, ts, tl, width)
        cp_s = _ag_cp_solid(t, t0, t1, cp0, cp1)
        if t < ts:
            cp_ag = cp_s
            h_ag = _ag_h_solid(t, t0, t1, cp0, cp1)
        elif t < tl:
            delta_t = t - ts
            cp_ag = cp_mid + cp_delta * f + latent / width
            h_ag = h_solid_ts + cp_mid * delta_t + 0.5 * cp_delta * delta_t * delta_t / width + f * latent
        else:
            cp_ag = cp_liq
            h_ag = h_liquid_tl + cp_liq * (t - tl)
        k_solid = k0 if t <= t0 else (k1 if t >= t1 else k0 + (k1 - k0) * (t - t0) / (t1 - t0))
        k_ag = k_liq if t >= tl else (1.0 - f) * k_solid + f * k_liq

        sc = 1.0 - ag
        h = ag * rho_ag * h_ag
        c = ag * rho_ag * cp_ag
        k = ag * k_ag
        if sc > 1.0e-14:
            h += sc * rho_ti * _ti_h(t, ti_a, ti_b, ti_molar, tref, ti_cp_ref)
            c += sc * rho_ti * _ti_cp(t, ti_a, ti_b, ti_molar, tref, ti_cp_ref)
            k += sc * (ti_k0 + ti_kslope * t)
        return h, c, k

    @njit(parallel=True, cache=False)
    def material_fields_3d(temperature, ag_fraction, parameters):
        nx, ny, nz = temperature.shape
        enthalpy = np.empty_like(temperature)
        capacity = np.empty_like(temperature)
        conductivity = np.empty_like(temperature)
        for i in prange(nx):
            for j in range(ny):
                for k in range(nz):
                    h, c, cond = _cell_fields(temperature[i, j, k], ag_fraction[i, j, k], parameters)
                    enthalpy[i, j, k] = h
                    capacity[i, j, k] = c
                    conductivity[i, j, k] = cond
        return enthalpy, capacity, conductivity

    @njit(parallel=True, cache=False)
    def conductance_faces_3d(conductivity, ag_fraction, spacing, interface_resistance, periodic):
        nx, ny, nz = conductivity.shape
        gx = conductivity.copy()
        gy = conductivity.copy()
        gz = conductivity.copy()
        for i in prange(nx):
            for j in range(ny):
                for k in range(nz):
                    ai = ag_fraction[i, j, k]
                    ki = conductivity[i, j, k]
                    if i + 1 < nx or periodic[0]:
                        ip = 0 if i + 1 == nx else i + 1
                        ap = ag_fraction[ip, j, k]
                        kp = conductivity[ip, j, k]
                        unlike = ai * (1.0 - ap) + (1.0 - ai) * ap
                        gx[i, j, k] = 1.0 / (0.5 * spacing / ki + 0.5 * spacing / kp + interface_resistance * unlike) / spacing
                    else:
                        gx[i, j, k] = 0.0
                    if j + 1 < ny or periodic[1]:
                        jp = 0 if j + 1 == ny else j + 1
                        ap = ag_fraction[i, jp, k]
                        kp = conductivity[i, jp, k]
                        unlike = ai * (1.0 - ap) + (1.0 - ai) * ap
                        gy[i, j, k] = 1.0 / (0.5 * spacing / ki + 0.5 * spacing / kp + interface_resistance * unlike) / spacing
                    else:
                        gy[i, j, k] = 0.0
                    if k + 1 < nz or periodic[2]:
                        kp_i = 0 if k + 1 == nz else k + 1
                        ap = ag_fraction[i, j, kp_i]
                        kp = conductivity[i, j, kp_i]
                        unlike = ai * (1.0 - ap) + (1.0 - ai) * ap
                        gz[i, j, k] = 1.0 / (0.5 * spacing / ki + 0.5 * spacing / kp + interface_resistance * unlike) / spacing
                    else:
                        gz[i, j, k] = 0.0
        return gx, gy, gz

    @njit(parallel=True, cache=False)
    def apply_stencil_3d(field, gx, gy, gz, periodic):
        nx, ny, nz = field.shape
        result = field.copy()
        for i in prange(nx):
            for j in range(ny):
                for k in range(nz):
                    center = field[i, j, k]
                    value = 0.0
                    if i > 0:
                        value += gx[i - 1, j, k] * (center - field[i - 1, j, k])
                    elif periodic[0]:
                        value += gx[nx - 1, j, k] * (center - field[nx - 1, j, k])
                    if i + 1 < nx:
                        value += gx[i, j, k] * (center - field[i + 1, j, k])
                    elif periodic[0]:
                        value += gx[i, j, k] * (center - field[0, j, k])
                    if j > 0:
                        value += gy[i, j - 1, k] * (center - field[i, j - 1, k])
                    elif periodic[1]:
                        value += gy[i, ny - 1, k] * (center - field[i, ny - 1, k])
                    if j + 1 < ny:
                        value += gy[i, j, k] * (center - field[i, j + 1, k])
                    elif periodic[1]:
                        value += gy[i, j, k] * (center - field[i, 0, k])
                    if k > 0:
                        value += gz[i, j, k - 1] * (center - field[i, j, k - 1])
                    elif periodic[2]:
                        value += gz[i, j, nz - 1] * (center - field[i, j, nz - 1])
                    if k + 1 < nz:
                        value += gz[i, j, k] * (center - field[i, j, k + 1])
                    elif periodic[2]:
                        value += gz[i, j, k] * (center - field[i, j, 0])
                    result[i, j, k] = value
        return result

else:

    def material_fields_3d(*args, **kwargs):  # pragma: no cover
        require_numba()

    def conductance_faces_3d(*args, **kwargs):  # pragma: no cover
        require_numba()

    def apply_stencil_3d(*args, **kwargs):  # pragma: no cover
        require_numba()
