#!/usr/bin/env python3
"""Non-integrating tests for the Ag64 step-20000 continuation entry."""
import copy
import json
from pathlib import Path
import sys
import unittest

import numpy as np

PACKAGE_DIR = Path(__file__).resolve().parent
REPO = PACKAGE_DIR.parents[5]
sys.path.insert(0, str(PACKAGE_DIR))
import followup_runner_step20000_extension as runner

CONFIG_PATH = PACKAGE_DIR / "FOLLOWUP_EXECUTION_CONFIG_STEP20000_40000.json"


class Step20000ExtensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        cls.checkpoint_path = REPO / cls.config["resume_checkpoint"]["archive_path"] / "checkpoint.json"
        cls.checkpoint = json.loads(cls.checkpoint_path.read_text(encoding="utf-8"))

    def test_exact_extension_bounds(self):
        phase = self.config["schedule"]["phase_1_hot_liquidification"]
        policy = self.config["step_policy"]
        self.assertEqual(phase["maximum_global_step"], 40000)
        self.assertEqual(phase["maximum_additional_steps_from_current_checkpoint"], 20000)
        self.assertEqual(phase["maximum_total_duration_from_step0_ps"], 20.0)
        self.assertEqual(phase["target_K"], 1600.0)
        self.assertEqual(policy["resume_global_step"], 20000)
        self.assertEqual(policy["maximum_global_step"], 40000)
        self.assertEqual(policy["maximum_added_steps"], 20000)
        self.assertEqual(policy["maximum_chunk_steps"], 1000)
        self.assertEqual(policy["maximum_chunk_wall_seconds"], 3600)
        self.assertFalse(self.config["execution_authorized"])

    def test_stage_a_boundary_and_global_cap(self):
        state = {"stage": "A_HOT_1600K"}
        self.assertEqual(runner.apply_completed_chunk_limits(copy.deepcopy(state), 39999, "CHUNK_COMPLETE", 40000, 40000), "CHUNK_COMPLETE")
        self.assertEqual(runner.apply_completed_chunk_limits(copy.deepcopy(state), 40000, "CHUNK_COMPLETE", 40000, 40000), "HOLD_STAGE_A_MAXIMUM_WITHOUT_LIQUID_GATE")
        self.assertEqual(runner.apply_completed_chunk_limits(copy.deepcopy(state), 40001, "CHUNK_COMPLETE", 40000, 40000), "HOLD_STAGE_A_MAXIMUM_WITHOUT_LIQUID_GATE")
        other_stage = {"stage": "B_COOL_1600_TO_1250K"}
        self.assertEqual(runner.apply_completed_chunk_limits(other_stage, 40000, "CHUNK_COMPLETE", 40000, 40000), "HOLD_GLOBAL_STEP_LIMIT")

    def test_step20000_restore_is_bitwise_state_preserving(self):
        restored = runner._load_resume_checkpoint(REPO, self.config)
        checkpoint = restored["checkpoint"]
        self.assertEqual(checkpoint["global_step"], 20000)
        self.assertEqual(restored["checkpoint_sha256"], self.config["resume_checkpoint"]["checkpoint_sha256"])
        self.assertEqual(restored["state"]["global_step"], 20000)
        self.assertEqual(restored["state"]["stage"], "A_HOT_1600K")
        self.assertEqual(restored["state"]["recovered_from_checkpoint_sha256"], self.config["resume_checkpoint"]["checkpoint_sha256"])
        self.assertTrue(np.array_equal(restored["atoms"].get_positions(), np.asarray(checkpoint["positions_A"])))
        self.assertTrue(np.array_equal(restored["atoms"].get_momenta() * runner.fs, np.asarray(checkpoint["momenta_amu_A_per_fs"])))
        self.assertEqual(set(restored["rngs"]), set(runner.RNG_NAMES))
        for name in runner.RNG_NAMES:
            self.assertEqual(
                runner.canonical_sha(restored["rngs"][name].bit_generator.state),
                runner.canonical_sha(checkpoint["rng_states"][name]),
            )
        self.assertEqual(checkpoint["potential_sha256"], runner.EXPECTED_POTENTIAL_SHA256)
        self.assertEqual(checkpoint["lock_sha256"], runner.SCIENTIFIC_LOCK_SHA256)
        self.assertEqual(checkpoint["parameter_sha256"], self.config["source_checkpoint"]["parameter_sha256"])

    def test_no_new_structure_or_integrator_state_is_created_by_restore(self):
        restored = runner._load_resume_checkpoint(REPO, self.config)
        self.assertEqual(restored["checkpoint"]["global_step"], 20000)
        self.assertEqual(restored["qa_monitor"]["gate_status"], "PENDING")
        self.assertIsNone(restored["qa_monitor"].get("MSD_A2"))
        self.assertFalse((Path(self.config["output_directory"])).exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
