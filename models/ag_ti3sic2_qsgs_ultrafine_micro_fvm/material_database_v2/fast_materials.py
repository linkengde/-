#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把数据库v2物性接到现有Numba/EW-safe有限体积步进器。"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).resolve().parent.parent
SOLVER = PROJECT / "solver_optimization"
PHASE = PROJECT / "phase_change_validation"
for _p in (str(PROJECT), str(PHASE), str(SOLVER)):
    if _p in sys.path:
        sys.path.remove(_p)
    sys.path.insert(0, _p)

from optimized_solver import FastCandidateMaterials
from material_database_v2.numba_v2_kernels import material_fields_3d_v2
from material_database_v2.material_model import PankratzMaterials


class FastPankratzMaterials(FastCandidateMaterials):
    """FastCandidateMaterials-compatible adapter with independent v2 tables."""

    def __init__(self, mushy_interval_width_K: float = 1.0, profile: bool = False):
        super().__init__(mushy_interval_width_K=mushy_interval_width_K, profile=profile)
        self.v2 = PankratzMaterials(mushy_interval_width_K)
        # 暴露FVM、检查点和诊断模块共同使用的物性字段接口。
        for name in (
            "width_K", "Tm_K", "Ts_K", "Tl_K", "Lf_J_mol", "Lf_J_kg",
            "ag_M_kg_mol", "rho_ag_solid_kg_m3", "rho_ag_liquid_candidate_kg_m3",
            "ag_T_K", "ag_cp_solid_J_kg_K", "ag_k", "ag_liquid_T_K",
            "ag_k_liquid_table_W_m_K", "ag_cp_liquid_J_kg_K", "ag_k_liquid_W_m_K",
            "ag_cp_solid_mid", "cp_ag_transition_delta", "h_ag_solid_at_Ts_J_kg",
            "h_ag_liquid_at_Tl_J_kg", "ti_M_kg_mol", "rho_ti_kg_m3",
            "ti_cp_a_J_mol_K", "ti_cp_b_J_mol", "ti_k_intercept_W_m_K",
            "ti_k_slope_W_m_K2", "ti_Tmin_K", "ti_Tmax_K", "reference_temperature_K",
            "_ti_cp_300_J_kg_K",
        ):
            setattr(self, name, getattr(self.v2, name))
        self._kernel_parameters = (
            self.ag_M_kg_mol, self.Tm_K, self.Ts_K, self.Tl_K, self.width_K,
            self.Lf_J_kg, self.rho_ag_solid_kg_m3, self.ag_cp_liquid_J_kg_K,
            self.ag_cp_solid_mid, self.cp_ag_transition_delta,
            self.h_ag_solid_at_Ts_J_kg, self.h_ag_liquid_at_Tl_J_kg,
            self.ag_T_K, self.ag_k, self.ag_liquid_T_K, self.ag_k_liquid_table_W_m_K,
            self.rho_ti_kg_m3, self.ti_M_kg_mol, self.ti_cp_a_J_mol_K,
            self.ti_cp_b_J_mol, self.reference_temperature_K,
            self._ti_cp_300_J_kg_K, self.ti_k_intercept_W_m_K,
            self.ti_k_slope_W_m_K2, self.ti_Tmax_K,
        )
        self._cached_temperature = None
        self._cached_ag = None
        self._cached_fields = None

    def _evaluate(self, temperature_K, ag_fraction):
        t = np.asarray(temperature_K, dtype=np.float64)
        ag = np.asarray(ag_fraction, dtype=np.float64)
        if t.ndim != 3 or t.shape != ag.shape or not t.flags.c_contiguous or not ag.flags.c_contiguous:
            return None
        if self._cached_temperature is t and self._cached_ag is ag:
            return self._cached_fields
        started = time.perf_counter() if self.profile else 0.0
        fields = material_fields_3d_v2(t, ag, self._kernel_parameters)
        if self.profile:
            self.profile_seconds += time.perf_counter() - started
            self.profile_calls += 1
        self._cached_temperature, self._cached_ag, self._cached_fields = t, ag, fields
        return fields

    def enthalpy_density_J_m3(self, temperature_K, ag_fraction):
        fields = self._evaluate(temperature_K, ag_fraction)
        return fields[0] if fields is not None else self.v2.enthalpy_density_J_m3(temperature_K, ag_fraction)

    def volumetric_capacity_J_m3_K(self, temperature_K, ag_fraction):
        fields = self._evaluate(temperature_K, ag_fraction)
        return fields[1] if fields is not None else self.v2.volumetric_capacity_J_m3_K(temperature_K, ag_fraction)

    def conductivity_W_m_K(self, temperature_K, ag_fraction):
        fields = self._evaluate(temperature_K, ag_fraction)
        return fields[2] if fields is not None else self.v2.conductivity_W_m_K(temperature_K, ag_fraction)

    def __getattr__(self, name):
        # 其他较少调用的标量函数转交给独立的v2物性模型。
        v2 = self.__dict__.get("v2")
        if v2 is not None and hasattr(v2, name):
            return getattr(v2, name)
        raise AttributeError(name)
