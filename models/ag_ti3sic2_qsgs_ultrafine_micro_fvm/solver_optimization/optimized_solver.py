"""Numba-backed implementation; base Python solver remains the reference."""

from __future__ import annotations

import time

import numpy as np

from implicit_enthalpy_fvm import BoundarySet, ImplicitEnthalpyFVM
from material_model import CandidateMaterials
from numba_kernels import (
    HAVE_NUMBA,
    apply_stencil_3d,
    conductance_faces_3d,
    material_fields_3d,
    require_numba,
)


class FastCandidateMaterials(CandidateMaterials):
    """Batch property evaluation using one fused Numba pass over the 3-D field."""

    def __init__(self, *args, profile: bool = False, **kwargs):
        require_numba()
        super().__init__(*args, **kwargs)
        self.profile = bool(profile)
        self.profile_seconds = 0.0
        self.profile_calls = 0
        self._cached_temperature = None
        self._cached_ag = None
        self._cached_fields = None
        self._kernel_parameters = (
            float(self.ag_T_K[0]), float(self.ag_T_K[-1]),
            float(self.ag_cp_solid_J_kg_K[0]), float(self.ag_cp_solid_J_kg_K[-1]),
            float(self.ag_k[0]), float(self.ag_k[-1]),
            float(self.rho_ag_solid_kg_m3), float(self.ag_cp_liquid_J_kg_K),
            float(self.ag_k_liquid_W_m_K), float(self.Tm_K), float(self.Ts_K),
            float(self.Tl_K), float(self.width_K), float(self.cp_ag_solid_mid),
            float(self.cp_ag_transition_delta), float(self.h_ag_solid_at_Ts_J_kg),
            float(self.h_ag_liquid_at_Tl_J_kg), float(self.Lf_J_kg),
            float(self.rho_ti_kg_m3), float(self.ti_M_kg_mol),
            float(self.ti_cp_a_J_mol_K), float(self.ti_cp_b_J_mol),
            float(self.reference_temperature_K), float(self._ti_cp_300_J_kg_K),
            float(self.ti_k_intercept_W_m_K), float(self.ti_k_slope_W_m_K2),
            float(self.ti_Tmax_K),
        )

    def _evaluate(self, temperature_K, ag_fraction):
        t = np.asarray(temperature_K, dtype=np.float64)
        ag = np.asarray(ag_fraction, dtype=np.float64)
        if t.ndim != 3 or t.shape != ag.shape or not t.flags.c_contiguous or not ag.flags.c_contiguous:
            return None
        if self._cached_temperature is t and self._cached_ag is ag:
            return self._cached_fields
        started = time.perf_counter() if self.profile else 0.0
        fields = material_fields_3d(t, ag, self._kernel_parameters)
        if self.profile:
            self.profile_seconds += time.perf_counter() - started
            self.profile_calls += 1
        # Strong references prevent pointer reuse from making the cache ambiguous.
        self._cached_temperature = t
        self._cached_ag = ag
        self._cached_fields = fields
        return fields

    def enthalpy_density_J_m3(self, temperature_K, ag_fraction):
        fields = self._evaluate(temperature_K, ag_fraction)
        if fields is None:
            return super().enthalpy_density_J_m3(temperature_K, ag_fraction)
        return fields[0]

    def volumetric_capacity_J_m3_K(self, temperature_K, ag_fraction):
        fields = self._evaluate(temperature_K, ag_fraction)
        if fields is None:
            return super().volumetric_capacity_J_m3_K(temperature_K, ag_fraction)
        return fields[1]

    def conductivity_W_m_K(self, temperature_K, ag_fraction):
        fields = self._evaluate(temperature_K, ag_fraction)
        if fields is None:
            return super().conductivity_W_m_K(temperature_K, ag_fraction)
        return fields[2]


class NumbaEnthalpyFVM(ImplicitEnthalpyFVM):
    """3-D uniform-grid solver with compiled, fused property and stencil kernels."""

    def __init__(self, ag_fraction, spacing_m, materials=None, **kwargs):
        require_numba()
        ag = np.ascontiguousarray(np.asarray(ag_fraction, dtype=np.float64))
        spacing = np.asarray(spacing_m, dtype=np.float64)
        if ag.ndim != 3 or not np.allclose(spacing, spacing[0], rtol=0.0, atol=0.0):
            raise ValueError("NumbaEnthalpyFVM currently supports uniform 3-D cubic cells")
        material_model = materials or FastCandidateMaterials()
        if not isinstance(material_model, FastCandidateMaterials):
            # Rebuild the accelerated adapter from the identical mushy width and
            # audited source table; no physical parameter is altered.
            material_model = FastCandidateMaterials(materials.width_K if materials else None)
        self.kernel_profile = {"conductance_seconds": 0.0, "conductance_calls": 0,
                               "stencil_seconds": 0.0, "stencil_calls": 0}
        super().__init__(ag, spacing, materials=material_model, **kwargs)
        self.periodic_flags = np.zeros(3, dtype=np.bool_)

    def _boundary_face_conductance_density(self, conductivity_W_m_K, boundary: BoundarySet):
        if self.ndim != 3:
            return super()._boundary_face_conductance_density(conductivity_W_m_K, boundary)
        faces = {}
        for key in boundary.dirichlet_K:
            axis = "xyz".index(key[0])
            index = 0 if key.endswith("low") else -1
            surface_k = np.take(conductivity_W_m_K, index, axis=axis)
            faces[key] = np.ascontiguousarray(2.0 * surface_k / self.spacing_m[axis] ** 2)
        return faces

    def _face_conductances_density(self, conductivity_W_m_K, boundary: BoundarySet):
        if self.ndim != 3 or not np.allclose(self.spacing_m, self.spacing_m[0], rtol=0.0, atol=0.0):
            return super()._face_conductances_density(conductivity_W_m_K, boundary)
        self.periodic_flags = np.asarray([axis in boundary.periodic_axes for axis in range(3)], dtype=np.bool_)
        started = time.perf_counter()
        g = conductance_faces_3d(
            np.ascontiguousarray(conductivity_W_m_K), self.ag_fraction,
            float(self.spacing_m[0]), self.interface_resistance_m2_K_W,
            self.periodic_flags,
        )
        self.kernel_profile["conductance_seconds"] += time.perf_counter() - started
        self.kernel_profile["conductance_calls"] += 1
        return list(g), self._boundary_face_conductance_density(conductivity_W_m_K, boundary)

    def _apply_laplacian(self, field, conductance_faces, boundary_faces):
        if self.ndim != 3:
            return super()._apply_laplacian(field, conductance_faces, boundary_faces)
        started = time.perf_counter()
        result = apply_stencil_3d(
            np.ascontiguousarray(field), conductance_faces[0], conductance_faces[1],
            conductance_faces[2], self.periodic_flags,
        )
        for key, gb in boundary_faces.items():
            axis = "xyz".index(key[0])
            selector = [slice(None)] * 3
            selector[axis] = 0 if key.endswith("low") else -1
            surface = result[tuple(selector)]
            surface += gb * field[tuple(selector)]
        self.kernel_profile["stencil_seconds"] += time.perf_counter() - started
        self.kernel_profile["stencil_calls"] += 1
        return result

    def _net_conduction_density(self, temperature_K, conductance_faces, boundary_faces, boundary):
        result = self._apply_laplacian(temperature_K, conductance_faces, boundary_faces)
        result *= -1.0
        for key, gb in boundary_faces.items():
            axis = "xyz".index(key[0])
            selector = [slice(None)] * 3
            selector[axis] = 0 if key.endswith("low") else -1
            result[tuple(selector)] += gb * float(boundary.dirichlet_K[key])
        return result

    def _boundary_out_rate_W(self, temperature_K, conductivity_W_m_K, boundary):
        boundary_faces = self._boundary_face_conductance_density(conductivity_W_m_K, boundary)
        total_rate = 0.0
        for key, gb in boundary_faces.items():
            axis = "xyz".index(key[0])
            selector = [slice(None)] * self.ndim
            selector[axis] = 0 if key.endswith("low") else -1
            total_rate += float(np.sum(gb * (temperature_K[tuple(selector)] - float(boundary.dirichlet_K[key])))
                                * self.volume_m3)
        return total_rate

    def advance(self, *args, **kwargs):
        record = super().advance(*args, **kwargs)
        active_ag = self.ag_fraction > 1.0e-14
        liquid = self.materials.liquid_fraction(self.temperature_K)
        record["max_ag_liquid_fraction"] = float(np.max(liquid[active_ag])) if np.any(active_ag) else 0.0
        record["ag_liquid_volume_fraction"] = float(np.sum(self.ag_fraction * liquid)
                                                    / max(float(np.sum(self.ag_fraction)), 1.0e-30))
        return record
