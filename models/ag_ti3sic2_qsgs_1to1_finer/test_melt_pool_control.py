import unittest

from melt_pool_control import HeatAction, PoolObservation, StablePoolController


def sample(t, *, cluster=800, volume=1000.0, connected=True, damage=0.0, recession=0.0, escaped=0):
    return PoolObservation(
        time_ps=t,
        largest_liquid_cluster_atoms=cluster,
        liquid_pool_volume_A3=volume,
        skeleton_connected=connected,
        skeleton_damage_fraction=damage,
        surface_recession_A=recession,
        escaped_atoms=escaped,
    )


class StablePoolControllerTests(unittest.TestCase):
    def test_heat_off_requires_persistence_and_volume_plateau(self):
        controller = StablePoolController()
        actions = [controller.update(sample(i / 10.0)).action for i in range(31)]
        self.assertTrue(all(action == HeatAction.CONTINUE for action in actions[:-1]))
        self.assertEqual(actions[-1], HeatAction.HEAT_OFF)

    def test_small_or_disconnected_liquid_cluster_does_not_turn_heat_off(self):
        controller = StablePoolController()
        for i in range(50):
            decision = controller.update(sample(i / 10.0, cluster=100))
        self.assertEqual(decision.action, HeatAction.CONTINUE)

    def test_growing_pool_does_not_turn_heat_off(self):
        controller = StablePoolController()
        for i in range(50):
            t = i / 10.0
            decision = controller.update(sample(t, volume=1000.0 + 100.0 * t))
        self.assertEqual(decision.action, HeatAction.CONTINUE)

    def test_skeleton_failure_or_escaped_atom_stops_run(self):
        controller = StablePoolController()
        decision = controller.update(sample(0.0, connected=False))
        self.assertEqual(decision.action, HeatAction.SAFE_STOP)
        self.assertIn("connectivity", decision.reason)

        other = StablePoolController()
        decision = other.update(sample(0.0, escaped=1))
        self.assertEqual(decision.action, HeatAction.SAFE_STOP)

    def test_missing_diagnostic_interval_restarts_persistence_window(self):
        controller = StablePoolController(max_sample_gap_ps=0.25)
        for t in [0.0, 0.1, 0.2, 0.3, 0.4, 1.0]:
            decision = controller.update(sample(t))
        self.assertEqual(decision.action, HeatAction.CONTINUE)
        self.assertAlmostEqual(decision.stable_duration_ps, 0.0)


if __name__ == "__main__":
    unittest.main()
