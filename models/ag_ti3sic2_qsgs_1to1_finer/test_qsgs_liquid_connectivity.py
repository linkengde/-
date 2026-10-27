import unittest

import numpy as np

from phase_field_qsgs.qsgs_liquid_connectivity import liquid_cluster_labels, liquid_cluster_metrics
from phase_field_qsgs.qsgs_liquid_tracking import LiquidClusterTracker


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

    def test_component_labels_merge_periodic_seam(self):
        mask = np.zeros((4, 3, 2), dtype=bool)
        mask[0, 1, 0] = True
        mask[3, 1, 0] = True
        labels, count = liquid_cluster_labels(mask)
        self.assertEqual(count, 1)
        self.assertEqual(labels[0, 1, 0], labels[3, 1, 0])

    def test_temporal_tracker_records_split_and_merge_lineage(self):
        tracker = LiquidClusterTracker((1.0, 1.0, 1.0), (8.0, 8.0, 4.0))
        ag = np.ones((8, 8, 4), dtype=bool)
        liquid = np.zeros_like(ag, dtype=float)
        liquid[2:5, 3, 1] = 1.0
        _, birth = tracker.update(0.0, liquid, ag)
        self.assertEqual(birth[0]["event"], "birth")
        self.assertEqual(birth[0]["track_id"], 1)

        liquid[3, 3, 1] = 0.0
        _, split = tracker.update(0.1, liquid, ag)
        self.assertEqual(len(split), 2)
        self.assertEqual({item["event"] for item in split}, {"split_parent_continues", "split_child"})
        self.assertEqual({item["parent_track_ids"][0] for item in split}, {1})

        liquid[3, 3, 1] = 1.0
        _, merge = tracker.update(0.2, liquid, ag)
        self.assertEqual(len(merge), 1)
        self.assertEqual(merge[0]["event"], "merge")
        self.assertCountEqual(merge[0]["parent_track_ids"], [1, 2])
        self.assertEqual(merge[0]["track_id"], 1)
        self.assertAlmostEqual(merge[0]["continuous_duration_ps"], 0.2)
        self.assertAlmostEqual(merge[0]["lateral_size_x_A"], 3.0)
        self.assertAlmostEqual(merge[0]["depth_from_top_A"], 2.5)
        # The child first appears at one sample and has no observed interval
        # before the next frame shows it merged back into the parent.
        self.assertAlmostEqual(tracker.closed_tracks[0]["sampled_continuous_duration_ps"], 0.0)


if __name__ == "__main__":
    unittest.main()
