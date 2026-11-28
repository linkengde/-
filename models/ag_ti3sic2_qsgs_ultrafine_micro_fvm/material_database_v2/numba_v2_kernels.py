#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pankratz v2物性场的Numba内核；不修改Stage 4旧内核。"""
from __future__ import annotations

import numpy as np

try:
    from numba import njit, prange
except ImportError:  # pragma: no cover
    njit, prange = None, range

HAVE_NUMBA = njit is not None


def require_numba():
    if not HAVE_NUMBA:
        raise RuntimeError("Numba依赖不可用，请安装solver_optimization/requirements-numba.txt")


if HAVE_NUMBA:
    @njit(inline="always")
    def _interp(x, xs, ys):
        if x <= xs[0]:
            return ys[0]
        if x >= xs[-1]:
            return ys[-1]
        for i in range(xs.size - 1):
            if x <= xs[i + 1]:
                f = (x - xs[i]) / (xs[i + 1] - xs[i])
                return ys[i] + f * (ys[i + 1] - ys[i])
        return ys[-1]

    @njit(inline="always")
    def _fields_one(t, ag, p):
        (m_ag, tm, ts, tl, width, latent, rho_ag, cp_liq, cp_mid,
         cp_delta, h_ts, h_tl, k_s_t, k_s, k_l_t, k_l,
         rho_ti, m_ti, ti_a, ti_b, tref, ti_cp_ref, ti_k0, ti_kslope, ti_tmax) = p
        if t < 298.15 or t > 2000.0:
            return np.nan, np.nan, np.nan
        if ag < 1.0 - 1.0e-14 and t > ti_tmax:
            return np.nan, np.nan, np.nan
        if t <= ts:
            liquid = 0.0
        elif t >= tl:
            liquid = 1.0
        else:
            liquid = (t - ts) / width
        cp_s_molar = 4.184 * (5.272 + 1.836e-3 * t + 22200.0 / (t * t))
        cp_s = cp_s_molar / m_ag
        h_solid = 4.184 * (5.272 * (t - 300.0) + 0.000918 * (t * t - 300.0 * 300.0)
                           + 22200.0 * (1.0 / 300.0 - 1.0 / t)) / m_ag
        if t < ts:
            cp_ag = cp_s
            h_ag = h_solid
        elif t < tl:
            delta_t = t - ts
            cp_ag = cp_mid + cp_delta * liquid + latent / width
            h_ag = h_ts + cp_mid * delta_t + 0.5 * cp_delta * delta_t * delta_t / width + liquid * latent
        else:
            cp_ag = cp_liq
            h_ag = h_tl + cp_liq * (t - tl)
        ks = _interp(t, k_s_t, k_s)
        kl = _interp(t, k_l_t, k_l)
        k_ag = (1.0 - liquid) * ks + liquid * kl
        sc = 1.0 - ag
        h = ag * rho_ag * h_ag
        capacity = ag * rho_ag * cp_ag
        conductivity = ag * k_ag
        if sc > 1.0e-14:
            if t < tref:
                ti_cp = ti_cp_ref
                ti_h = ti_cp_ref * (t - tref)
            else:
                ti_cp = (ti_a - ti_b / t) / m_ti
                ti_h = (ti_a * (t - tref) - ti_b * np.log(t / tref)) / m_ti
            h += sc * rho_ti * ti_h
            capacity += sc * rho_ti * ti_cp
            conductivity += sc * (ti_k0 + ti_kslope * t)
        return h, capacity, conductivity

    @njit(parallel=True, cache=False)
    def material_fields_3d_v2(temperature, ag_fraction, parameters):
        nx, ny, nz = temperature.shape
        h = np.empty_like(temperature)
        cpv = np.empty_like(temperature)
        k = np.empty_like(temperature)
        for i in prange(nx):
            for j in range(ny):
                for z in range(nz):
                    hh, cc, kk = _fields_one(temperature[i, j, z], ag_fraction[i, j, z], parameters)
                    h[i, j, z], cpv[i, j, z], k[i, j, z] = hh, cc, kk
        return h, cpv, k
else:
    def material_fields_3d_v2(*args, **kwargs):  # pragma: no cover
        require_numba()
