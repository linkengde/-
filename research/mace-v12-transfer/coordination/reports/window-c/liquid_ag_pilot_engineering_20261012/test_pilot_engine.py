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
    walltime_reached,
    _save_exit_boundary,
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

    def test_walltime_is_a_stop_boundary_not_an_overflow_exception(self):
        assert_segment_limits(self.lock, 1000, 3600.0)
        with self.assertRaises(PreflightError):
            assert_segment_limits(self.lock, 1001, 1.0)
        self.assertFalse(walltime_reached(3599.0, 3600.0))
        self.assertTrue(walltime_reached(3600.0, 3600.0))
        self.assertTrue(walltime_reached(3600.01, 3600.0))
        # An expired wall clock must remain observable as a stop decision even
        # when elapsed time is a little over the bound.
        assert_segment_limits(self.lock, 37, 3600.01)

    def test_walltime_crossing_finishes_one_step_then_saves_nonperiodic_boundary(self):
        with tempfile.TemporaryDirectory(prefix="D-rho100-wall-boundary-test-") as tmp:
            root = Path(tmp)
            clock = {"elapsed": 3599.0}
            state = {
                "global_step": 37,
                "rng_states": {"thermostat": {"synthetic_counter": 37}},
            }
            checkpoint = root / "fixture_checkpoint.json"
            checkpoint.write_text(json.dumps(state, sort_keys=True), encoding="utf-8")
            integrate_calls = 0

            # Synthetic indivisible integrator step: it begins below the cap,
            # completes above it, and advances the state/RNG exactly once.
            while state["global_step"] < 1000:
                if walltime_reached(clock["elapsed"], 3600.0):
                    break
                integrate_calls += 1
                state["global_step"] += 1
                state["rng_states"]["thermostat"]["synthetic_counter"] += 1
                clock["elapsed"] = 3601.2
                if walltime_reached(clock["elapsed"], 3600.0):
                    break

            checkpoint.write_text(json.dumps(state, sort_keys=True), encoding="utf-8")
            saved = json.loads(checkpoint.read_text(encoding="utf-8"))
            self.assertEqual(integrate_calls, 1)
            self.assertEqual(saved["global_step"], 38)
            self.assertEqual(saved["rng_states"]["thermostat"]["synthetic_counter"], 38)
            self.assertEqual(hashlib.sha256(checkpoint.read_bytes()).hexdigest(), hashlib.sha256(
                json.dumps(saved, sort_keys=True).encode("utf-8")
            ).hexdigest())

    def test_walltime_already_expired_does_not_start_a_step_and_saves_latest_step(self):
        for elapsed in (3600.0, 3602.0):
            with self.subTest(elapsed=elapsed), tempfile.TemporaryDirectory(prefix="D-rho100-wall-expired-test-") as tmp:
                root = Path(tmp)
                state = {"global_step": 37, "rng_state": "unchanged-at-last-complete-step"}
                integrate_calls = 0
                if not walltime_reached(elapsed, 3600.0):
                    integrate_calls += 1
                    state["global_step"] += 1
                checkpoint = root / "fixture_checkpoint.json"
                checkpoint.write_text(json.dumps(state, sort_keys=True), encoding="utf-8")
                saved = json.loads(checkpoint.read_text(encoding="utf-8"))
                self.assertEqual(integrate_calls, 0)
                self.assertEqual(saved["global_step"], 37)
                self.assertEqual(saved["rng_state"], "unchanged-at-last-complete-step")

    def test_exception_after_completed_step_persists_latest_state_and_rng(self):
        with tempfile.TemporaryDirectory(prefix="D-rho100-exception-boundary-test-") as tmp:
            checkpoint = Path(tmp) / "fixture_checkpoint.json"
            state = {
                "global_step": 41,
                "rng_states": {"thermostat": {"synthetic_counter": 41}},
            }
            # Simulate a completed step followed by an output exception before
            # the periodic checkpoint; exercise the same production exit-save
            # helper used by run_authorized().
            state["global_step"] += 1
            state["rng_states"]["thermostat"]["synthetic_counter"] += 1

            def save_boundary(status):
                state["boundary_status"] = status
                checkpoint.write_text(json.dumps(state, sort_keys=True), encoding="utf-8")
                return hashlib.sha256(checkpoint.read_bytes()).hexdigest(), status

            status, failure, checkpoint_sha = _save_exit_boundary(
                save_boundary, "RUNNING", None, OSError("synthetic thermal-log write failure"),
            )

            saved = json.loads(checkpoint.read_text(encoding="utf-8"))
            self.assertEqual(status, "STOPPED_AFTER_EXCEPTION_AT_LAST_COMPLETED_STEP")
            self.assertIn("synthetic thermal-log write failure", failure)
            self.assertEqual(checkpoint_sha, hashlib.sha256(checkpoint.read_bytes()).hexdigest())
            self.assertEqual(saved["boundary_status"], status)
            self.assertEqual(saved["global_step"], 42)
            self.assertEqual(saved["rng_states"]["thermostat"]["synthetic_counter"], 42)

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
