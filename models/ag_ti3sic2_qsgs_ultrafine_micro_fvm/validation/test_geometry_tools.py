"""几何处理工具的纯几何单元测试，不初始化或加热FVM。"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "geometry_audit"))
sys.path.insert(0, str(ROOT / "geometry_generation"))
from audit_geometries import chord_lengths, component_audit, interface_area_density, maximal_inscribed_sphere_diameters
from generate_qsgs_pair import correlated_exact_fraction, qsgs_growth_rank
from reference_generate_geometry import qsgs_four_parameter


class GeometryToolsTests(unittest.TestCase):
    def test_periodic_chord_joins_boundary_run_only_when_requested(self) -> None:
        phase = np.array([True, False, True], dtype=bool).reshape(3, 1, 1)
        self.assertEqual(chord_lengths(phase, axis=0, periodic=False).tolist(), [1, 1])
        self.assertEqual(chord_lengths(phase, axis=0, periodic=True).tolist(), [2])

    def test_6face_connectivity_does_not_join_diagonal_voxels(self) -> None:
        phase = np.zeros((4, 4, 4), dtype=bool)
        phase[0, 0, 0] = True
        phase[1, 1, 1] = True
        result = component_audit(phase)
        self.assertEqual(result["component_count_6face"], 2)
        self.assertAlmostEqual(result["largest_component_fraction"], 0.5)
        self.assertFalse(any(result["largest_component_spans_xyz"]))

    def test_sphere_center_diameter_and_face_interface(self) -> None:
        phase = np.zeros((9, 9, 9), dtype=bool)
        phase[2:7, 2:7, 2:7] = True
        widths, _, _ = maximal_inscribed_sphere_diameters(phase)
        self.assertAlmostEqual(float(widths.max()), 6.0)
        self.assertGreater(interface_area_density(phase, spacing_m=2.0), 0.0)

    def test_qsgs_growth_matches_repository_generator(self) -> None:
        shape = (18, 18, 18)
        rank, generations = qsgs_growth_rank(shape, 0.06, (0.22, 0.22, 0.22), 0.30, 77)
        reference, reference_generations = qsgs_four_parameter(shape, 0.06, (0.22, 0.22, 0.22), 0.30, 77)
        np.testing.assert_array_equal(rank >= 0, reference)
        self.assertEqual(generations, reference_generations)

    def test_correlated_threshold_is_deterministic_and_exact_fraction(self) -> None:
        raw, _ = qsgs_four_parameter((20, 20, 20), 0.06, (0.22, 0.22, 0.22), 0.30, 91)
        a, _ = correlated_exact_fraction(raw, 2.0)
        b, _ = correlated_exact_fraction(raw, 2.0)
        np.testing.assert_array_equal(a, b)
        self.assertEqual(int(a.sum()), int(round(raw.size * 0.30)))
        self.assertTrue(np.array_equal(~a | a, np.ones(raw.shape, dtype=bool)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
