#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""基于Pankratz主数据的Ag/Ti₃SiC₂单温焓物性模型。"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PARAMETERS = json.loads((HERE / "material_parameters.json").read_text(encoding="utf-8"))


def _read_csv(name: str) -> list[dict[str, str]]:
    with (HERE / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


class PankratzMaterials:
    """可直接交给Stage 4隐式FVM接口的独立v2物性候选。"""

    def __init__(self, mushy_interval_width_K: float = 1.0):
        if not np.isfinite(mushy_interval_width_K) or mushy_interval_width_K <= 0.0:
            raise ValueError("相变数值温区必须为有限正数")
        ag = PARAMETERS["Ag"]
        ti = PARAMETERS["Ti3SiC2"]
        solver = PARAMETERS["solver_policy"]
        self.width_K = float(mushy_interval_width_K)
        self.Tm_K = float(ag["melting_temperature_K"])
        self.Ts_K = self.Tm_K - 0.5 * self.width_K
        self.Tl_K = self.Tm_K + 0.5 * self.width_K
        self.Lf_J_mol = float(ag["latent_heat_fusion_J_mol"])
        self.ag_M_kg_mol = float(ag["molar_mass_kg_mol"])
        self.Lf_J_kg = self.Lf_J_mol / self.ag_M_kg_mol
        self.rho_ag_solid_kg_m3 = float(ag["solid_reference_density_fixed_kg_m3"])
        self.rho_ag_liquid_candidate_kg_m3 = 9327.1
        self.ag_T_K = np.array([float(r["temperature_K"]) for r in _read_csv("Ag_solid_properties.csv")])
        self.ag_cp_solid_J_kg_K = np.array([float(r["Cp_mass_J_kg_K"]) for r in _read_csv("Ag_solid_properties.csv")])
        self.ag_k = np.array([float(r["thermal_conductivity_W_m_K"]) for r in _read_csv("Ag_solid_properties.csv")])
        self.ag_liquid_T_K = np.array([float(r["temperature_K"]) for r in _read_csv("Ag_liquid_properties.csv")])
        self.ag_k_liquid_table_W_m_K = np.array([float(r["thermal_conductivity_W_m_K"]) for r in _read_csv("Ag_liquid_properties.csv")])
        self.ag_cp_liquid_J_kg_K = float(ag["liquid_cp_J_mol_K"]) / self.ag_M_kg_mol
        self.ag_k_liquid_W_m_K = float(self.ag_k_liquid_table_W_m_K[0])
        self.ag_cp_solid_mid = float(self.ag_cp_solid_J_kg_K[-1])
        self.cp_ag_transition_delta = self.ag_cp_liquid_J_kg_K - self.ag_cp_solid_mid
        self.h_ag_solid_at_Ts_J_kg = float(self.ag_solid_enthalpy_J_kg(self.Ts_K))
        transition_sensible = 0.5 * (self.ag_cp_solid_mid + self.ag_cp_liquid_J_kg_K) * self.width_K
        self.h_ag_liquid_at_Tl_J_kg = self.h_ag_solid_at_Ts_J_kg + transition_sensible + self.Lf_J_kg

        self.ti_M_kg_mol = float(ti["molar_mass_kg_mol"])
        self.rho_ti_kg_m3 = float(ti["reference_density_kg_m3"])
        self.ti_cp_a_J_mol_K = 164.4
        self.ti_cp_b_J_mol = 16419.0
        self.ti_k_intercept_W_m_K = 38.6
        self.ti_k_slope_W_m_K2 = -0.0045
        self.ti_Tmin_K = 300.0
        self.ti_Tmax_K = float(ti["primary_solver_upper_temperature_K"])
        self.reference_temperature_K = float(solver["reference_temperature_K"])
        self._ti_cp_300_J_kg_K = self.ti_cp_J_kg_K(self.reference_temperature_K)

    @staticmethod
    def _as_float_or_array(value, original):
        return float(value) if np.asarray(original).ndim == 0 else value

    @staticmethod
    def _validate_ag_range(t):
        if np.any(~np.isfinite(t)) or np.any(t < 298.15) or np.any(t > 2000.0):
            raise ValueError("Ag热力学候选仅支持298.15–2000 K；禁止静默外推")

    def ag_cp_solid_molar_J_mol_K(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        self._validate_ag_range(t)
        cp = 4.184 * (5.272 + 1.836e-3 * t + 22200.0 / t**2)
        return self._as_float_or_array(cp, temperature_K)

    def ag_solid_enthalpy_J_kg(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        self._validate_ag_range(t)
        tref = 300.0
        h_cal_mol = (5.272 * (t - tref) + 0.000918 * (t**2 - tref**2)
                     + 22200.0 * (1.0 / tref - 1.0 / t))
        h = h_cal_mol * 4.184 / self.ag_M_kg_mol
        return self._as_float_or_array(h, temperature_K)

    def liquid_fraction(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        self._validate_ag_range(t)
        f = np.where(t <= self.Ts_K, 0.0,
                     np.where(t >= self.Tl_K, 1.0, (t - self.Ts_K) / self.width_K))
        return self._as_float_or_array(f, temperature_K)

    def ag_specific_enthalpy_J_kg(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        self._validate_ag_range(t)
        f = np.asarray(self.liquid_fraction(t))
        delta = t - self.Ts_K
        mushy = (self.h_ag_solid_at_Ts_J_kg + self.ag_cp_solid_mid * delta
                 + 0.5 * self.cp_ag_transition_delta * delta**2 / self.width_K
                 + f * self.Lf_J_kg)
        h = np.where(t <= self.Ts_K, self.ag_solid_enthalpy_J_kg(t), mushy)
        h = np.where(t >= self.Tl_K,
                     self.h_ag_liquid_at_Tl_J_kg + self.ag_cp_liquid_J_kg_K * (t - self.Tl_K), h)
        return self._as_float_or_array(h, temperature_K)

    def ag_temperature_from_enthalpy_J_kg(self, enthalpy_J_kg):
        """Invert the monotone Ag specific-enthalpy relation by bisection."""
        h = np.asarray(enthalpy_J_kg, dtype=float)
        if np.any(~np.isfinite(h)):
            raise ValueError("Ag enthalpy must be finite")
        lo = np.full_like(h, 298.15)
        hi = np.full_like(h, 2000.0)
        h_lo = self.ag_specific_enthalpy_J_kg(lo)
        h_hi = self.ag_specific_enthalpy_J_kg(hi)
        if np.any((h < h_lo) | (h > h_hi)):
            raise ValueError("Ag enthalpy is outside the supported temperature range")
        for _ in range(64):
            mid = 0.5 * (lo + hi)
            h_mid = self.ag_specific_enthalpy_J_kg(mid)
            go_right = h_mid < h
            lo = np.where(go_right, mid, lo)
            hi = np.where(go_right, hi, mid)
        t = 0.5 * (lo + hi)
        return self._as_float_or_array(t, enthalpy_J_kg)

    def ag_sensible_cp_J_kg_K(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        self._validate_ag_range(t)
        cp_s = np.asarray(self.ag_cp_solid_molar_J_mol_K(t)) / self.ag_M_kg_mol
        f = np.asarray(self.liquid_fraction(t))
        cp = np.where(t < self.Ts_K, cp_s,
                      self.ag_cp_solid_mid + self.cp_ag_transition_delta * f)
        cp = np.where(t >= self.Tl_K, self.ag_cp_liquid_J_kg_K, cp)
        return self._as_float_or_array(cp, temperature_K)

    def ag_apparent_cp_J_kg_K(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        cp = np.asarray(self.ag_sensible_cp_J_kg_K(t))
        cp += np.where((t > self.Ts_K) & (t < self.Tl_K), self.Lf_J_kg / self.width_K, 0.0)
        return self._as_float_or_array(cp, temperature_K)

    def ag_conductivity_W_m_K(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        self._validate_ag_range(t)
        # 熔点两侧分别插值；只有数值糊状温区内按液相率作明确混合。
        k_solid = np.interp(t, self.ag_T_K, self.ag_k)
        k_liquid = np.interp(t, self.ag_liquid_T_K, self.ag_k_liquid_table_W_m_K)
        f = np.asarray(self.liquid_fraction(t))
        k = (1.0 - f) * k_solid + f * k_liquid
        return self._as_float_or_array(k, temperature_K)

    def ti_cp_J_kg_K(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        if np.any(~np.isfinite(t)) or np.any(t > self.ti_Tmax_K + 1e-10):
            raise ValueError(f"Ti₃SiC₂主模型有效上限为{self.ti_Tmax_K:g} K")
        eval_t = np.maximum(t, self.ti_Tmin_K)
        cp = (self.ti_cp_a_J_mol_K - self.ti_cp_b_J_mol / eval_t) / self.ti_M_kg_mol
        cp = np.where(t < self.reference_temperature_K, self._ti_cp_300_J_kg_K if hasattr(self, "_ti_cp_300_J_kg_K") else cp, cp)
        return self._as_float_or_array(cp, temperature_K)

    def ti_specific_enthalpy_J_kg(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        if np.any(~np.isfinite(t)) or np.any(t > self.ti_Tmax_K + 1e-10):
            raise ValueError(f"Ti₃SiC₂主模型有效上限为{self.ti_Tmax_K:g} K")
        te = np.maximum(t, self.reference_temperature_K)
        h = (self.ti_cp_a_J_mol_K * (te - self.reference_temperature_K)
             - self.ti_cp_b_J_mol * np.log(te / self.reference_temperature_K)) / self.ti_M_kg_mol
        below = t < self.reference_temperature_K
        h = np.where(below, self._ti_cp_300_J_kg_K * (t - self.reference_temperature_K), h)
        return self._as_float_or_array(h, temperature_K)

    def ti_conductivity_W_m_K(self, temperature_K):
        t = np.asarray(temperature_K, dtype=float)
        if np.any(~np.isfinite(t)) or np.any(t > self.ti_Tmax_K + 1e-10):
            raise ValueError(f"Ti₃SiC₂主模型有效上限为{self.ti_Tmax_K:g} K")
        return self._as_float_or_array(self.ti_k_intercept_W_m_K + self.ti_k_slope_W_m_K2 * t, temperature_K)

    def enthalpy_density_J_m3(self, temperature_K, ag_fraction):
        t, ag = np.broadcast_arrays(np.asarray(temperature_K, dtype=float), np.asarray(ag_fraction, dtype=float))
        if np.any(~np.isfinite(ag)) or np.any((ag < 0.0) | (ag > 1.0)):
            raise ValueError("Ag volume fraction must lie in [0,1]")
        sc = 1.0 - ag
        h = ag * self.rho_ag_solid_kg_m3 * self.ag_specific_enthalpy_J_kg(t)
        if np.any(sc > 1.0e-14):
            ti_t = np.where(sc > 1.0e-14, t, self.reference_temperature_K)
            h += sc * self.rho_ti_kg_m3 * self.ti_specific_enthalpy_J_kg(ti_t)
        return h

    def volumetric_capacity_J_m3_K(self, temperature_K, ag_fraction):
        t, ag = np.broadcast_arrays(np.asarray(temperature_K, dtype=float), np.asarray(ag_fraction, dtype=float))
        sc = 1.0 - ag
        c = ag * self.rho_ag_solid_kg_m3 * self.ag_apparent_cp_J_kg_K(t)
        if np.any(sc > 1.0e-14):
            ti_t = np.where(sc > 1.0e-14, t, self.reference_temperature_K)
            c += sc * self.rho_ti_kg_m3 * self.ti_cp_J_kg_K(ti_t)
        return c

    def conductivity_W_m_K(self, temperature_K, ag_fraction):
        t, ag = np.broadcast_arrays(np.asarray(temperature_K, dtype=float), np.asarray(ag_fraction, dtype=float))
        sc = 1.0 - ag
        k = ag * self.ag_conductivity_W_m_K(t)
        if np.any(sc > 1.0e-14):
            ti_t = np.where(sc > 1.0e-14, t, self.reference_temperature_K)
            k += sc * self.ti_conductivity_W_m_K(ti_t)
        return k
