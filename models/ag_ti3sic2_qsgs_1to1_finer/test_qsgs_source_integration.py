from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.special import erf


MODEL_DIR = Path(__file__).resolve().parent
PHASE_DIR = MODEL_DIR / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_DIR))
from qsgs_enthalpy_solver import QSGSEnthalpyGrid  # noqa: E402
from qsgs_source_integration import integrate_qsgs_source_power_W  # noqa: E402


class QSGSSourceIntegrationTests(unittest.TestCase):
    def test_uniform_ag_matches_analytic_gaussian_exponential_integral(self):
        size_A = np.asarray([12.0, 9.0, 8.0])
        peak_flux_W_m2 = 8.0e13
        radius_A = 3.0
        depth_A = 2.0
        shape = (4, 3, 2)
        powers = integrate_qsgs_source_power_W(
            np.ones((6, 5, 4), dtype=bool),
            shape,
            size_A,
            peak_flux_W_m2,
            radius_A,
            depth_A,
        )
        cx, cy = size_A[:2] / 2.0
        ix = 0.5 * radius_A * np.sqrt(np.pi) * (
            erf((size_A[0] - cx) / radius_A) - erf(-cx / radius_A)
        )
        iy = 0.5 * radius_A * np.sqrt(np.pi) * (
            erf((size_A[1] - cy) / radius_A) - erf(-cy / radius_A)
        )
        iz = depth_A * (1.0 - np.exp(-size_A[2] / depth_A))
        expected = peak_flux_W_m2 / (depth_A * 1.0e-10) * ix * iy * iz * 1.0e-30
        self.assertAlmostEqual(float(powers.sum()), expected, delta=expected * 2.0e-14)

    def test_exact_source_integrals_coarsen_on_strictly_nested_grids(self):
        rng = np.random.default_rng(20261028)
        source_ag = rng.random((8, 8, 8)) < 0.5
        size_A = np.asarray([8.0, 8.0, 8.0])
        coarse = integrate_qsgs_source_power_W(
            source_ag, (2, 2, 2), size_A, 8.0e13, 3.0, 2.0
        )
        fine = integrate_qsgs_source_power_W(
            source_ag, (4, 4, 4), size_A, 8.0e13, 3.0, 2.0
        )
        fine_coarsened = fine.reshape(2, 2, 2, 2, 2, 2).sum(axis=(1, 3, 5))
        np.testing.assert_allclose(fine_coarsened, coarse, rtol=3.0e-14, atol=0.0)

    def test_exact_qsgs_source_power_is_independent_of_grid_and_lateral_tile(self):
        source_geometry = MODEL_DIR / "qsgs_morphology.npz"
        grid_10A = PHASE_DIR / "qsgs_phase_fractions_10A.npz"
        grid_5A = PHASE_DIR / "qsgs_phase_fractions_5A.npz"
        for tile_xy in (1, 3):
            ten = QSGSEnthalpyGrid(
                grid_10A,
                "A",
                tile_xy=tile_xy,
                source_geometry_npz=source_geometry,
            )
            five = QSGSEnthalpyGrid(
                grid_5A,
                "A",
                tile_xy=tile_xy,
                source_geometry_npz=source_geometry,
            )
            self.assertAlmostEqual(
                ten.actual_source_power_W,
                five.actual_source_power_W,
                delta=ten.actual_source_power_W * 2.0e-13,
            )
            self.assertAlmostEqual(
                ten.source_W_m3.sum() * ten.cell_volume_m3,
                ten.actual_source_power_W,
                delta=ten.actual_source_power_W * 2.0e-14,
            )
            self.assertEqual(
                ten.source_integration_method,
                "exact separable integral over binary QSGS Ag voxels; contributions to round-off-only projected cells are omitted",
            )


if __name__ == "__main__":
    unittest.main()
