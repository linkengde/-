from pathlib import Path
import sys
import unittest

import numpy as np


MODEL_DIR = Path(__file__).resolve().parent
PHASE_DIR = MODEL_DIR / "phase_field_qsgs"
sys.path.insert(0, str(PHASE_DIR))
from analyze_matched_grid_fields import conservative_average, weighted_relative_l1  # noqa: E402


class MatchedGridAnalysisTests(unittest.TestCase):
    def test_overlap_mapping_preserves_integrated_energy_density(self):
        rng = np.random.default_rng(20261026)
        fine_density = rng.random((8, 10, 6))
        target_shape = (4, 5, 3)
        size_A = np.asarray((8.0, 10.0, 6.0))
        coarse_density = conservative_average(fine_density, target_shape, size_A)
        fine_volume = np.prod(size_A / np.asarray(fine_density.shape))
        coarse_volume = np.prod(size_A / np.asarray(target_shape))
        self.assertAlmostEqual(
            float(fine_density.sum() * fine_volume),
            float(coarse_density.sum() * coarse_volume),
            delta=1.0e-12,
        )

    def test_overlap_mapping_preserves_uniform_fields(self):
        fine = np.full((12, 12, 8), 7.25)
        coarse = conservative_average(fine, (6, 6, 4), np.asarray((12.0, 12.0, 8.0)))
        np.testing.assert_allclose(coarse, 7.25, rtol=0.0, atol=1.0e-14)

    def test_common_ag_weighting_ignores_geometry_fraction_only_difference(self):
        coarse_phase_fraction = np.asarray([0.25, 0.80])
        mapped_fine_phase_fraction = np.asarray([0.50, 0.60])
        common_ag_volume = np.minimum(coarse_phase_fraction, mapped_fine_phase_fraction)
        coarse_intensive_field = np.asarray([3.0, 9.0])
        fine_intensive_field = coarse_intensive_field.copy()
        self.assertEqual(
            weighted_relative_l1(
                fine_intensive_field - coarse_intensive_field,
                coarse_intensive_field,
                common_ag_volume,
                volume_m3=1.0,
            ),
            0.0,
        )


if __name__ == "__main__":
    unittest.main()
