import unittest

import numpy as np

from phase_field_qsgs.qsgs_liquid_connectivity import liquid_cluster_metrics


class LiquidConnectivityTests(unittest.TestCase):
    def test_connectivity_uses_shared_faces_only(self):
        mask = np.zeros((4, 4, 3), dtype=bool)
        mask[0, 0, 0] = True
        mask[1, 1, 0] = True  # Diagonal contact is not a shared face.
        mask[3, 3, 2] = True
        metrics = liquid_cluster_metrics(mask, (1.0, 1.0, 1.0), periodic_axes=())
        self.assertEqual(metrics["connected_liquid_cluster_count"], 3)
        self.assertEqual(metrics["largest_liquid_cluster_voxels"], 1)

    def test_periodic_seam_merges_cluster_and_centroid_wraps(self):
        mask = np.zeros((4, 3, 2), dtype=bool)
        mask[0, 1, 0] = True
        mask[3, 1, 0] = True
        mask[2, 0, 1] = True
        metrics = liquid_cluster_metrics(mask, (1.0, 2.0, 3.0))
        self.assertEqual(metrics["connected_liquid_cluster_count"], 2)
        self.assertEqual(metrics["liquid_voxel_count"], 3)
        self.assertEqual(metrics["largest_liquid_cluster_voxels"], 2)
        self.assertEqual(metrics["largest_liquid_cluster_volume_A3"], 12.0)
        centroid = metrics["largest_liquid_cluster_centroid_A"]
        self.assertLess(min(centroid[0], 4.0 - centroid[0]), 1.0e-12)
        self.assertAlmostEqual(centroid[1], 3.0)
        self.assertAlmostEqual(centroid[2], 1.5)

    def test_empty_mask_has_no_cluster_or_centroid(self):
        metrics = liquid_cluster_metrics(np.zeros((2, 3, 4), dtype=bool), (2.0, 3.0, 4.0))
        self.assertEqual(metrics["connected_liquid_cluster_count"], 0)
        self.assertEqual(metrics["largest_liquid_cluster_volume_A3"], 0.0)
        self.assertIsNone(metrics["largest_liquid_cluster_centroid_A"])


if __name__ == "__main__":
    unittest.main()
