#!/usr/bin/env python3
"""候选Ag/Ti3SiC2焓、比热和热导率关系；所有单位均显式标注。"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
CONFIG = json.loads((HERE / "VALIDATION_CONFIG.json").read_text(encoding="utf-8"))
TABLE = PROJECT / CONFIG["source_candidate_table"]
SOURCE_CONFIG = PROJECT / CONFIG["source_candidate_json"]


class CandidateMaterials:
    """采用仓库候选参数，Ag相变用窄温区焓正则化。"""

    def __init__(self, mushy_interval_width_K: float | None = None):
        if mushy_interval_width_K is None:
            mushy_interval_width_K = float(CONFIG["phase_change"]["Ag_mushy_interval_width_K"])
        if not np.isfinite(mushy_interval_width_K) or mushy_interval_width_K <= 0:
            raise ValueError("Ag mushy interval width must be finite and positive")
        self.width_K = float(mushy_interval_width_K)
        source = json.loads(SOURCE_CONFIG.read_text(encoding="utf-8"))
        with TABLE.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        ag_solid = [r for r in rows if r["material"] == "Ag" and r["phase"] == "solid"]
        ag_liquid = [r for r in rows if r["material"] == "Ag" and r["phase"] == "liquid"]
        if len(ag_solid) != 2 or len(ag_liquid) != 1:
            raise ValueError("Expected two solid and one liquid Ag candidate anchors")
        ag_solid.sort(key=lambda r: float(r["temperature_K"]))
        self.ag_T_K = np.asarray([float(r["temperature_K"]) for r in ag_solid])
        self.ag_cp_molar = np.asarray([float(r["Cp_molar_J_mol_K"]) for r in ag_solid])
        self.ag_k = np.asarray([float(r["thermal_conductivity_W_m_K"]) for r in ag_solid])
        self.ag_M_kg_mol = float(source["ag"]["molar_mass_kg_mol"])
        self.rho_ag_solid_kg_m3 = float(source["ag"]["solid_reference_density_kg_m3"])
        self.rho_ag_liquid_candidate_kg_m3 = float(source["ag"]["liquid_density_at_melting_kg_m3"])
        self.Tm_K = float(source["ag"]["melting_temperature_K"])
        self.Lf_J_mol = float(source["ag"]["latent_heat_fusion_J_mol"])
        self.Lf_J_kg = float(source["ag"]["latent_heat_fusion_J_kg"])
        self.ag_cp_solid_J_kg_K = self.ag_cp_molar / self.ag_M_kg_mol
        self.ag_cp_liquid_J_kg_K = float(ag_liquid[0]["Cp_molar_J_mol_K"]) / self.ag_M_kg_mol
        self.ag_k_liquid_W_m_K = float(ag_liquid[0]["thermal_conductivity_W_m_K"])
        self.Ts_K = self.Tm_K - 0.5 * self.width_K
        self.Tl_K = self.Tm_K + 0.5 * self.width_K
        self.cp_ag_solid_mid = float(np.interp(self.Tm_K, self.ag_T_K, self.ag_cp_solid_J_kg_K))
        self.cp_ag_transition_delta = self.ag_cp_liquid_J_kg_K - self.cp_ag_solid_mid

        ti = source["ti3sic2"]
        self.ti_M_kg_mol = float(ti["molar_mass_kg_mol"])
        self.rho_ti_kg_m3 = float(ti["density_kg_m3"])
        cp = ti["cp_model"]
        self.ti_cp_a_J_mol_K = float(cp["a_J_mol_K"])
        self.ti_cp_b_J_mol = float(cp["b_J_mol"])
        k = ti["conductivity_model"]
        self.ti_k_intercept_W_m_K = float(k["intercept_W_m_K"])
        self.ti_k_slope_W_m_K2 = float(k["slope_W_m_K2"])
        self.ti_Tmin_K = float(ti["reference_temperature_K"])
        # Use the recorded near-maximum fit anchor, not the 2000 K extrapolation.
        self.ti_Tmax_K = float(CONFIG["phase_change"]["Ti3SiC2_maximum_temperature_K"])
        self.reference_temperature_K = float(CONFIG["phase_change"]["enthalpy_reference_temperature_K"])
        self._ag_cumulative_h_J_kg = np.zeros_like(self.ag_T_K)
        self._ag_cumulative_h_J_kg[1:] = np.cumsum(
            np.diff(self.ag_T_K)
            * (self.ag_cp_solid_J_kg_K[:-1] + self.ag_cp_solid_J_kg_K[1:])
            * 0.5
        )
        self.h_ag_solid_at_Ts_J_kg = float(self._ag_solid_enthalpy(self.Ts_K))
        transition_sensible = 0.5 * (self.cp_ag_solid_mid + self.ag_cp_liquid_J_kg_K) * self.width_K
        self.h_ag_liquid_at_Tl_J_kg = (
            self.h_ag_solid_at_Ts_J_kg + transition_sensible + self.Lf_J_kg
        )
        self._ti_cp_300_J_kg_K = self.ti_cp_J_kg_K(self.reference_temperature_K)

    def _ag_solid_enthalpy(self, temperature_K: np.ndarray | float) -> np.ndarray | float:
        t = np.asarray(temperature_K, dtype=float)
        clipped = np.clip(t, self.ag_T_K[0], self.ag_T_K[-1])
        index = np.clip(np.searchsorted(self.ag_T_K, clipped, side="right") - 1, 0, len(self.ag_T_K) - 2)
        dt = clipped - self.ag_T_K[index]
        slope = (self.ag_cp_solid_J_kg_K[index + 1] - self.ag_cp_solid_J_kg_K[index]) / (
            self.ag_T_K[index + 1] - self.ag_T_K[index]
        )
        h = self._ag_cumulative_h_J_kg[index] + self.ag_cp_solid_J_kg_K[index] * dt + 0.5 * slope * dt**2
        h = np.where(t < self.ag_T_K[0], self.ag_cp_solid_J_kg_K[0] * (t - self.ag_T_K[0]), h)
        return float(h) if np.asarray(temperature_K).ndim == 0 else h

    def liquid_fraction(self, temperature_K: np.ndarray) -> np.ndarray:
        t = np.asarray(temperature_K, dtype=float)
        return np.where(t <= self.Ts_K, 0.0, np.where(t >= self.Tl_K, 1.0, (t - self.Ts_K) / self.width_K))

    def ag_specific_enthalpy_J_kg(self, temperature_K: np.ndarray) -> np.ndarray:
        t = np.asarray(temperature_K, dtype=float)
        f = self.liquid_fraction(t)
        delta = t - self.Ts_K
        mushy_h = (
            self.h_ag_solid_at_Ts_J_kg
            + self.cp_ag_solid_mid * delta
            + 0.5 * self.cp_ag_transition_delta * delta**2 / self.width_K
            + f * self.Lf_J_kg
        )
        h = np.where(t <= self.Ts_K, self._ag_solid_enthalpy(t), mushy_h)
        h = np.where(t >= self.Tl_K, self.h_ag_liquid_at_Tl_J_kg + self.ag_cp_liquid_J_kg_K * (t - self.Tl_K), h)
        return float(h) if t.ndim == 0 else h

    def ag_sensible_cp_J_kg_K(self, temperature_K: np.ndarray) -> np.ndarray:
        t = np.asarray(temperature_K, dtype=float)
        f = self.liquid_fraction(t)
        cp_solid = np.interp(t, self.ag_T_K, self.ag_cp_solid_J_kg_K)
        cp = np.where(t < self.Ts_K, cp_solid, self.cp_ag_solid_mid + self.cp_ag_transition_delta * f)
        cp = np.where(t >= self.Tl_K, self.ag_cp_liquid_J_kg_K, cp)
        return float(cp) if t.ndim == 0 else cp

    def ag_apparent_cp_J_kg_K(self, temperature_K: np.ndarray) -> np.ndarray:
        t = np.asarray(temperature_K, dtype=float)
        cp = np.asarray(self.ag_sensible_cp_J_kg_K(t))
        cp = cp + np.where((t > self.Ts_K) & (t < self.Tl_K), self.Lf_J_kg / self.width_K, 0.0)
        return float(cp) if t.ndim == 0 else cp

    def ag_conductivity_W_m_K(self, temperature_K: np.ndarray) -> np.ndarray:
        t = np.asarray(temperature_K, dtype=float)
        f = self.liquid_fraction(t)
        ks = np.interp(t, self.ag_T_K, self.ag_k)
        return np.where(t >= self.Tl_K, self.ag_k_liquid_W_m_K,
                        (1.0 - f) * ks + f * self.ag_k_liquid_W_m_K)

    def ti_cp_J_kg_K(self, temperature_K: np.ndarray | float) -> np.ndarray | float:
        t = np.asarray(temperature_K, dtype=float)
        if np.any(t > self.ti_Tmax_K + 1e-10):
            raise ValueError(f"Ti3SiC2 candidate property limit exceeded: Tmax={self.ti_Tmax_K} K")
        eval_t = np.maximum(t, self.ti_Tmin_K)
        cp = (self.ti_cp_a_J_mol_K - self.ti_cp_b_J_mol / eval_t) / self.ti_M_kg_mol
        cp = np.where(t < self.ti_Tmin_K, self._ti_cp_300_J_kg_K if hasattr(self, "_ti_cp_300_J_kg_K") else cp, cp)
        return float(cp) if np.asarray(temperature_K).ndim == 0 else cp

    def ti_specific_enthalpy_J_kg(self, temperature_K: np.ndarray) -> np.ndarray:
        t = np.asarray(temperature_K, dtype=float)
        if np.any(t > self.ti_Tmax_K + 1e-10):
            raise ValueError(f"Ti3SiC2 candidate property limit exceeded: Tmax={self.ti_Tmax_K} K")
        te = np.maximum(t, self.ti_Tmin_K)
        h = (
            self.ti_cp_a_J_mol_K * (te - self.reference_temperature_K)
            - self.ti_cp_b_J_mol * np.log(te / self.reference_temperature_K)
        ) / self.ti_M_kg_mol
        below = t < self.reference_temperature_K
        h = np.where(below, self._ti_cp_300_J_kg_K * (t - self.reference_temperature_K), h)
        return float(h) if t.ndim == 0 else h

    def ti_conductivity_W_m_K(self, temperature_K: np.ndarray) -> np.ndarray:
        t = np.asarray(temperature_K, dtype=float)
        if np.any(t > self.ti_Tmax_K + 1e-10):
            raise ValueError(f"Ti3SiC2 candidate property limit exceeded: Tmax={self.ti_Tmax_K} K")
        return self.ti_k_intercept_W_m_K + self.ti_k_slope_W_m_K2 * t

    def enthalpy_density_J_m3(self, temperature_K: np.ndarray, ag_fraction: np.ndarray) -> np.ndarray:
        t, ag = np.broadcast_arrays(np.asarray(temperature_K, dtype=float), np.asarray(ag_fraction, dtype=float))
        sc = 1.0 - ag
        h = ag * self.rho_ag_solid_kg_m3 * self.ag_specific_enthalpy_J_kg(t)
        if np.any(sc > 1e-14):
            ti_t = np.where(sc > 1e-14, t, self.reference_temperature_K)
            h += sc * self.rho_ti_kg_m3 * self.ti_specific_enthalpy_J_kg(ti_t)
        return float(h) if h.ndim == 0 else h

    def volumetric_capacity_J_m3_K(self, temperature_K: np.ndarray, ag_fraction: np.ndarray) -> np.ndarray:
        t, ag = np.broadcast_arrays(np.asarray(temperature_K, dtype=float), np.asarray(ag_fraction, dtype=float))
        capacity = ag * self.rho_ag_solid_kg_m3 * self.ag_apparent_cp_J_kg_K(t)
        sc = 1.0 - ag
        if np.any(sc > 1e-14):
            ti_t = np.where(sc > 1e-14, t, self.reference_temperature_K)
            capacity += sc * self.rho_ti_kg_m3 * self.ti_cp_J_kg_K(ti_t)
        return float(capacity) if capacity.ndim == 0 else capacity

    def conductivity_W_m_K(self, temperature_K: np.ndarray, ag_fraction: np.ndarray) -> np.ndarray:
        t, ag = np.broadcast_arrays(np.asarray(temperature_K, dtype=float), np.asarray(ag_fraction, dtype=float))
        conductivity = ag * self.ag_conductivity_W_m_K(t)
        sc = 1.0 - ag
        if np.any(sc > 1e-14):
            ti_t = np.where(sc > 1e-14, t, self.reference_temperature_K)
            conductivity += sc * self.ti_conductivity_W_m_K(ti_t)
        return float(conductivity) if conductivity.ndim == 0 else conductivity
