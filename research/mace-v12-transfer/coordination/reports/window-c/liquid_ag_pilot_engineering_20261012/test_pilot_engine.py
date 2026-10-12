import json
import hashlib
import tempfile
import unittest
from pathlib import Path

import numpy as np

from pilot_engine import (
    LOCK_SHA256,
    PreflightError,
    assert_segment_limits,
    make_rngs,
    pack_random_positions,
    run_authorized,
    verify_lock,
)


LOCK = Path("/tmp/LIQUID_AG_PILOT_LOCK.json")


class Rho100EngineeringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lock, cls.lock_sha = verify_lock(LOCK, LOCK_SHA256)

    def test_lock_hash_and_hold_gate_are_exact(self):
        self.assertEqual(self.lock_sha, LOCK_SHA256)
        self.assertIs(self.lock["execution_authorized"], False)
        self.assertEqual(self.lock["scope"]["master_seed"], 125010003)

    def test_three_child_pcg64_streams_and_rejection_packing_repeat(self):
        one, provenance_one = make_rngs(125010003)
        two, provenance_two = make_rngs(125010003)
        self.assertEqual(provenance_one, provenance_two)
        self.assertEqual(list(one), ["packing", "velocity_initialization", "thermostat"])
        self.assertTrue(all(type(rng.bit_generator).__name__ == "PCG64" for rng in one.values()))
        positions_one, tries_one = pack_random_positions(10.717455315295654, 64, 2.2, 200000, one["packing"])
        positions_two, tries_two = pack_random_positions(10.717455315295654, 64, 2.2, 200000, two["packing"])
        np.testing.assert_array_equal(positions_one, positions_two)
        self.assertEqual(tries_one, tries_two)
        self.assertEqual(positions_one.shape, (64, 3))
        delta = positions_one[:, None, :] - positions_one[None, :, :]
        delta -= np.rint(delta / 10.717455315295654) * 10.717455315295654
        distances = np.sqrt((delta * delta).sum(axis=2))
        distances[np.diag_indices(64)] = np.inf
        self.assertGreaterEqual(float(distances.min()), 2.2)

    def test_checkpoint_restore_is_not_fixture_release(self):
        self.assertFalse(self.lock["initial_structure_generator"]["structure_generated"])
        self.assertFalse(self.lock["execution_authorized"])

    def test_limits_accept_only_locked_boundary(self):
        assert_segment_limits(self.lock, 1000, 3600.0)
        with self.assertRaises(PreflightError):
            assert_segment_limits(self.lock, 1001, 1.0)
        with self.assertRaises(PreflightError):
            assert_segment_limits(self.lock, 1000, 3600.01)

    def test_current_hold_cannot_create_formal_output(self):
        with tempfile.TemporaryDirectory(prefix="D-rho100-hold-gate-test-") as tmp:
            root = Path(tmp)
            auth = root / "unauthorized.json"
            auth.write_text(json.dumps({
                "execution_authorized": False,
                "owner": "D",
                "task_id": "C-D-AGL-1250K-rho100-pilot01",
                "lock_sha256": LOCK_SHA256,
                "maximum_added_steps": 1000,
                "maximum_wall_seconds": 3600,
            }), encoding="utf-8")
            auth_sha = hashlib.sha256(auth.read_bytes()).hexdigest()
            output = root / "C-D-AGL-1250K-rho100-pilot01"
            with self.assertRaises(PreflightError):
                run_authorized(
                    LOCK, LOCK_SHA256,
                    Path("/workspace/WMH2006_Ag-local-validation/Ag.eam.alloy"),
                    output, auth, auth_sha,
                )
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
