import unittest

import numpy as np

from melt_metrics import classify_liquid_ag, largest_connected_cluster


class LiquidMetricTests(unittest.TestCase):
    def test_liquid_flag_requires_disorder_and_diffusion(self):
        cna = np.array([1, 5, 5, 5, 1])
        displacement = np.array([2.0, 2.0, 0.9, 1.0, 0.1])
        expected = np.array([False, True, False, True, False])
        np.testing.assert_array_equal(classify_liquid_ag(cna, displacement), expected)

    def test_cluster_finds_periodic_xy_connection(self):
        xyz = np.array([[0.1, 5.0, 5.0], [9.9, 5.0, 5.0], [5.0, 5.0, 5.0]])
        result = largest_connected_cluster(xyz, np.ones(3, dtype=bool), (10, 10, 10), 0.5)
        self.assertEqual(result.member_count, 3)
        self.assertEqual(result.largest_cluster_atoms, 2)
        self.assertEqual(int(result.largest_cluster_mask.sum()), 2)

    def test_cluster_does_not_wrap_across_free_z_surface(self):
        xyz = np.array([[5.0, 5.0, 0.1], [5.0, 5.0, 9.9]])
        result = largest_connected_cluster(xyz, np.ones(2, dtype=bool), (10, 10, 10), 0.5)
        self.assertEqual(result.largest_cluster_atoms, 1)

    def test_empty_liquid_set_returns_zero(self):
        result = largest_connected_cluster(np.zeros((2, 3)), np.zeros(2, dtype=bool), (10, 10, 10), 3.8)
        self.assertEqual(result.member_count, 0)
        self.assertEqual(result.largest_cluster_atoms, 0)


if __name__ == "__main__":
    unittest.main()
