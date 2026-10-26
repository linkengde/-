"""Integrate the separable Gaussian/exponential source over binary QSGS Ag voxels."""

from __future__ import annotations

import numpy as np
from scipy.special import erf

ANGSTROM3_TO_M3 = 1.0e-30


def _overlap_integrals_A(
    source_count: int,
    target_count: int,
    length_A: float,
    mode: str,
    characteristic_length_A: float,
    center_A: float | None = None,
) -> np.ndarray:
    """Integrals of one spatial source factor over source/target overlaps."""
    source_edges = np.linspace(0.0, length_A, source_count + 1)
    target_edges = np.linspace(0.0, length_A, target_count + 1)
    lo = np.maximum(target_edges[:-1, None], source_edges[None, :-1])
    hi = np.minimum(target_edges[1:, None], source_edges[None, 1:])
    width = np.maximum(hi - lo, 0.0)
    length = characteristic_length_A
    if mode == "gaussian":
        assert center_A is not None
        return (
            0.5
            * length
            * np.sqrt(np.pi)
            * (erf((hi - center_A) / length) - erf((lo - center_A) / length))
        ) * (width > 0.0)
    if mode == "depth_exponential":
        # q(z)=exp[-(L-z)/length], evaluated with expm1 for narrow overlaps.
        result = length * np.exp((lo - length_A) / length) * np.expm1(width / length)
        return np.where(width > 0.0, result, 0.0)
    raise ValueError(f"Unsupported source integration mode: {mode}")


def integrate_qsgs_source_power_W(
    source_ag_mask: np.ndarray,
    target_shape_xyz: tuple[int, int, int] | np.ndarray,
    size_A: np.ndarray,
    peak_flux_W_m2: float,
    gaussian_radius_A: float,
    absorption_depth_A: float,
    source_center_xy_A: tuple[float, float] | None = None,
) -> np.ndarray:
    """Return exact cell-integrated source power for a binary QSGS Ag mask.

    The source is ``peak_flux / absorption_depth * exp(-r²/R²) * exp(-d/l)``.
    Source voxels and thermal cells are piecewise-constant rectangular boxes;
    separability lets us integrate each Gaussian/exponential factor over their
    exact overlaps and contract those weights with the original Ag mask.
    ``size_A`` and all integration lengths use Å; returned values are W/cell.
    """
    ag = np.asarray(source_ag_mask, dtype=bool)
    target_shape = np.asarray(target_shape_xyz, dtype=int)
    size = np.asarray(size_A, dtype=float)
    if ag.ndim != 3 or any(n <= 0 for n in ag.shape):
        raise ValueError("source_ag_mask must be a non-empty 3D boolean field")
    if target_shape.shape != (3,) or np.any(target_shape <= 0):
        raise ValueError("target_shape_xyz must contain three positive grid counts")
    if size.shape != (3,) or np.any(~np.isfinite(size)) or np.any(size <= 0.0):
        raise ValueError("size_A must contain three finite positive lengths")
    if not np.isfinite(peak_flux_W_m2) or peak_flux_W_m2 < 0.0:
        raise ValueError("peak_flux_W_m2 must be finite and non-negative")
    if not np.isfinite(gaussian_radius_A) or gaussian_radius_A <= 0.0:
        raise ValueError("gaussian_radius_A must be finite and positive")
    if not np.isfinite(absorption_depth_A) or absorption_depth_A <= 0.0:
        raise ValueError("absorption_depth_A must be finite and positive")

    center = source_center_xy_A or (size[0] / 2.0, size[1] / 2.0)
    x_integrals = _overlap_integrals_A(
        ag.shape[0], int(target_shape[0]), float(size[0]), "gaussian", gaussian_radius_A, float(center[0])
    )
    y_integrals = _overlap_integrals_A(
        ag.shape[1], int(target_shape[1]), float(size[1]), "gaussian", gaussian_radius_A, float(center[1])
    )
    z_integrals = _overlap_integrals_A(
        ag.shape[2], int(target_shape[2]), float(size[2]), "depth_exponential", absorption_depth_A
    )

    # Contract in stages so peak temporary memory scales with two dimensions,
    # even for the 3×3, 5 Å comparison domain.
    along_x = np.einsum("ai,ijk->ajk", x_integrals, ag, optimize=True)
    along_y = np.einsum("bj,ajk->abk", y_integrals, along_x, optimize=True)
    spatial_integral_A3 = np.einsum("ck,abk->abc", z_integrals, along_y, optimize=True)
    source_density_prefactor_W_m3 = peak_flux_W_m2 / (absorption_depth_A * 1.0e-10)
    return source_density_prefactor_W_m3 * spatial_integral_A3 * ANGSTROM3_TO_M3
