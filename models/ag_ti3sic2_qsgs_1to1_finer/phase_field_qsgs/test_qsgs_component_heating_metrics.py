from __future__ import annotations

from types import SimpleNamespace
import unittest

import numpy as np

from qsgs_component_heating_extend import _largest_liquid_cluster, _periodic_interval


class LiquidClusterMetricsTests(unittest.TestCase):
    def test_equal_size_disconnected_clusters_do_not_get_summed(self) -> None:
        ag = np.zeros((6, 1, 1), dtype=bool)
        ag[[0, 1, 3, 4], 0, 0] = True
        component_id = np.full(ag.shape, -1, dtype=np.int64)
        component_id[ag] = np.arange(4)
        graph = SimpleNamespace(
            component_id=component_id,
            voxel_count=np.ones(4, dtype=np.int64),
            block_index=np.argwhere(ag).astype(np.int32),
            edge_i=np.array([0, 2], dtype=np.int64),
            edge_j=np.array([1, 3], dtype=np.int64),
        )
        backend = SimpleNamespace(
            lattice_liquid_fraction=np.full(4, 0.8),
            graph=graph,
            node_volume_A3=np.ones(4),
            source_node_W=np.ones(4),
            source_power_W=4.0,
            source_ag=ag,
            spacing_A=np.array([1.0, 1.0, 1.0]),
            size_A=np.array([6.0, 1.0, 1.0]),
        )

        metrics, cluster_ids = _largest_liquid_cluster(backend)

        self.assertEqual(cluster_ids.tolist(), [0, 1])
        self.assertEqual(metrics["largest_cluster_node_count"], 2)
        self.assertEqual(metrics["largest_cluster_volume_A3"], 2.0)
        self.assertEqual(metrics["equal_size_largest_cluster_count"], 2)
        self.assertTrue(metrics["largest_cluster_is_tied"])
        self.assertEqual(metrics["bbox_voxel_xyz"]["x"]["start"], 0)
        self.assertEqual(metrics["bbox_voxel_xyz"]["x"]["end"], 1)

    def test_periodic_bbox_uses_short_interval_across_seam(self) -> None:
        interval = _periodic_interval(np.array([0, 9]), 10, 2.0)
        self.assertEqual(interval["start_index"], 9)
        self.assertEqual(interval["end_index"], 0)
        self.assertEqual(interval["extent_A"], 4.0)
        self.assertTrue(interval["wraps_boundary"])


if __name__ == "__main__":
    unittest.main()
