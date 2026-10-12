#!/usr/bin/env python3
"""Non-integrating checks for Stage A continuation from archived step 21000."""
import copy
import json
from pathlib import Path
import sys
import unittest

import numpy as np

PACKAGE_DIR = Path(__file__).resolve().parent
REPO = PACKAGE_DIR.parents[5]
sys.path.insert(0, str(PACKAGE_DIR))
import followup_runner_step21000_extension as runner

CONFIG_PATH = PACKAGE_DIR / "FOLLOWUP_EXECUTION_CONFIG_STEP21000_60000.json"


class Step21000ExtensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        cls.archive = REPO / cls.config["resume_checkpoint"]["archive_path"]
        cls.checkpoint = json.loads((cls.archive / "checkpoint.json").read_text(encoding="utf-8"))

    def test_exact_resume_and_global_bounds(self):
        phase = self.config["schedule"]["phase_1_hot_liquidification"]
        policy = self.config["step_policy"]
        self.assertEqual(policy["resume_global_step"], 21000)
        self.assertEqual(policy["maximum_global_step"], 60000)
        self.assertEqual(policy["maximum_added_steps"], 39000)
        self.assertEqual(phase["maximum_global_step"], 60000)
        self.assertEqual(phase["maximum_additional_steps_from_current_checkpoint"], 39000)
        self.assertEqual(phase["maximum_total_duration_from_step0_ps"], 30.0)
        self.assertEqual(policy["maximum_chunk_steps"], 1000)
        self.assertEqual(policy["maximum_chunk_wall_seconds"], 3600)
        self.assertFalse(self.config["execution_authorized"])

    def test_stage_a_pending_continues_and_pass_transitions(self):
        state = {"stage": "A_HOT_1600K"}
        pending = {"gate_status": "PENDING"}
        passed = {"gate_status": "PASS"}
        held = {"gate_status": "QA_HOLD"}
        self.assertEqual(runner.apply_monitor_transition(copy.deepcopy(state), 22000, pending), "CHUNK_COMPLETE")
        transitioned = copy.deepcopy(state)
        self.assertEqual(runner.apply_monitor_transition(transitioned, 22000, passed), "CHUNK_COMPLETE")
        self.assertEqual(transitioned["stage"], "B_COOL_1600_TO_1250K")
        self.assertEqual(transitioned["stage_b_start_step"], 22000)
        self.assertEqual(runner.apply_monitor_transition(copy.deepcopy(state), 22000, held), "HOLD_QA_REQUIRED_INPUT_OR_METRIC")

    def test_global_safety_cap_is_fail_closed(self):
        state = {"stage": "A_HOT_1600K"}
        self.assertEqual(runner.apply_completed_chunk_limits(copy.deepcopy(state), 59999, "CHUNK_COMPLETE", 60000, 60000), "CHUNK_COMPLETE")
        self.assertEqual(runner.apply_completed_chunk_limits(copy.deepcopy(state), 60000, "CHUNK_COMPLETE", 60000, 60000), "HOLD_STAGE_A_MAXIMUM_WITHOUT_LIQUID_GATE")
        self.assertEqual(runner.apply_completed_chunk_limits(copy.deepcopy(state), 60001, "CHUNK_COMPLETE", 60000, 60000), "HOLD_STAGE_A_MAXIMUM_WITHOUT_LIQUID_GATE")

    def test_step21000_checkpoint_positions_momenta_and_rng_restore_exactly(self):
        restored = runner._load_resume_checkpoint(REPO, self.config)
        cp = restored["checkpoint"]
        self.assertEqual(cp["global_step"], 21000)
        self.assertEqual(restored["checkpoint_sha256"], self.config["resume_checkpoint"]["checkpoint_sha256"])
        self.assertEqual(restored["state"]["stage"], "A_HOT_1600K")
        self.assertTrue(np.array_equal(restored["atoms"].get_positions(), np.asarray(cp["positions_A"])))
        self.assertTrue(np.array_equal(restored["atoms"].get_momenta() * runner.fs, np.asarray(cp["momenta_amu_A_per_fs"])))
        self.assertEqual(set(restored["rngs"]), set(runner.RNG_NAMES))
        for name in runner.RNG_NAMES:
            self.assertEqual(runner.canonical_sha(restored["rngs"][name].bit_generator.state), runner.canonical_sha(cp["rng_states"][name]))

    def test_restore_does_not_touch_output_or_create_integrator(self):
        output = Path(self.config["output_directory"])
        self.assertFalse(output.exists())
        restored = runner._load_resume_checkpoint(REPO, self.config)
        self.assertEqual(restored["checkpoint"]["global_step"], 21000)
        self.assertEqual(restored["qa_monitor"]["gate_status"], "PENDING")
        self.assertFalse(output.exists())

    def test_graceful_signal_after_archive_does_not_duplicate_boundary(self):
        # A signal at an already archived segment boundary must terminate
        # without asking the immutable archive backend to write that same path.
        passing_archive = {k: True for k in (
            "archive_upload_pass", "archive_fresh_fetch_pass", "archive_restore_pass",
            "checkpoint_restore_pass", "sha256_match_pass",
        )}
        self.assertTrue(runner._should_continue_after_archive("CHUNK_COMPLETE", passing_archive, False))
        self.assertFalse(runner._should_continue_after_archive("CHUNK_COMPLETE", passing_archive, True))


if __name__ == "__main__":
    unittest.main(verbosity=2)
