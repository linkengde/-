import copy
import json
import unittest
from pathlib import Path

from segment_gate_proposal import SegmentHold, next_step_guard, validate_segment


HERE = Path(__file__).parent
PROPOSED = json.loads((HERE / "PROPOSED_SEGMENT_1000_2000.json").read_text())
STATE = {
    "task_id": "C-D-AGL-1250K-rho097-seed01", "step": 1000,
    "scientific_parameters_sha256": "8f792d6206c1b7bae7b900e7ea81445d6cd3e837de2b5979c88bb977650ab7b2",
    "source_manifest_sha256": "9a694ea8f93655efa7a1c8695669ab2f17520c2d1b90f6df28115b2b1f9c2ad6",
    "rng_state": {name: {"bit_generator": "PCG64"} for name in ("packing", "velocity", "thermostat")},
}
STATE_SHA = "e08574678b3f8231c0a52ab25febb13c9d06e815b28149eb306d4b451c83ef6f"
ARCHIVE_SHA = "b62a818e58a2cf423e4cb1b9ed0468ae73ad578a75c4b2cb9fb3f28baf75e7b1"
CONFIG_SHA = "a" * 64


class SegmentGateProposalTests(unittest.TestCase):
    def signed(self):
        cfg = copy.deepcopy(PROPOSED)
        cfg["release_status"] = "SIGNED_BY_C"
        return cfg

    def validate(self, config, state=STATE, requested_stop_step=2000,
                 checkpoint_state_sha256=STATE_SHA, checkpoint_archive_sha256=ARCHIVE_SHA,
                 actual_config_sha256=CONFIG_SHA, externally_pinned_config_sha256=CONFIG_SHA,
                 owner="D", task_id=STATE["task_id"]):
        return validate_segment(config, state, owner=owner, task_id=task_id,
                                requested_stop_step=requested_stop_step,
                                checkpoint_state_sha256=checkpoint_state_sha256,
                                checkpoint_archive_sha256=checkpoint_archive_sha256,
                                actual_config_sha256=actual_config_sha256,
                                externally_pinned_config_sha256=externally_pinned_config_sha256)

    def test_unsigned_proposal_holds(self):
        with self.assertRaisesRegex(SegmentHold, "NOT_SIGNED"):
            self.validate(PROPOSED)

    def test_explicitly_signed_synthetic_exact_segment_validates_without_running_md(self):
        result = self.validate(self.signed())
        self.assertEqual(result["status"], "VALIDATOR_PASS_ONLY_NOT_COMPUTE_AUTHORIZATION")
        self.assertEqual((result["start_step"], result["end_step"], result["added_steps"]), (1000, 2000, 1000))

    def test_implicit_continue_endpoint_holds(self):
        with self.assertRaisesRegex(SegmentHold, "STOP_STEP"):
            self.validate(self.signed(), requested_stop_step=None)

    def test_default_120000_endpoint_holds(self):
        with self.assertRaisesRegex(SegmentHold, "STOP_STEP"):
            self.validate(self.signed(), requested_stop_step=120000)

    def test_wrong_owner_or_task_holds(self):
        with self.assertRaisesRegex(SegmentHold, "OWNER_OR_TASK"):
            self.validate(self.signed(), owner="E")

    def test_checkpoint_hash_mismatch_holds(self):
        with self.assertRaisesRegex(SegmentHold, "CHECKPOINT_STATE_SHA"):
            self.validate(self.signed(), checkpoint_state_sha256="0" * 64)

    def test_wrong_start_step_holds(self):
        bad = dict(STATE, step=1001)
        with self.assertRaisesRegex(SegmentHold, "START_STEP"):
            self.validate(self.signed(), state=bad)

    def test_missing_rng_state_holds(self):
        bad = dict(STATE, rng_state={"thermostat": {"bit_generator": "PCG64"}})
        with self.assertRaisesRegex(SegmentHold, "RNG_STATE"):
            self.validate(self.signed(), state=bad)

    def test_changed_scientific_identity_holds(self):
        bad = dict(STATE, scientific_parameters_sha256="f" * 64)
        with self.assertRaisesRegex(SegmentHold, "SCIENTIFIC_PARAMETER"):
            self.validate(self.signed(), state=bad)

    def test_source_manifest_identity_mismatch_holds(self):
        bad = dict(STATE, source_manifest_sha256="0" * 64)
        with self.assertRaisesRegex(SegmentHold, "SOURCE_MANIFEST"):
            self.validate(self.signed(), state=bad)

    def test_checkpoint_archive_and_release_config_hash_mismatch_hold(self):
        with self.assertRaisesRegex(SegmentHold, "ARCHIVE_SHA"):
            self.validate(self.signed(), checkpoint_archive_sha256="0" * 64)
        with self.assertRaisesRegex(SegmentHold, "CONFIG_SHA"):
            self.validate(self.signed(), actual_config_sha256="b" * 64)

    def test_wrong_end_or_more_than_1000_steps_holds(self):
        bad = self.signed()
        bad["authorized_end_step_if_later_signed"] = 3000
        bad["maximum_added_steps_if_later_signed"] = 2000
        with self.assertRaisesRegex(SegmentHold, "EXACTLY_1000"):
            self.validate(bad, requested_stop_step=3000)

    def test_per_step_guard_allows_only_through_2000(self):
        self.assertEqual(next_step_guard(1999, 2000), 2000)
        with self.assertRaisesRegex(SegmentHold, "END_REACHED"):
            next_step_guard(2000, 2000)
        with self.assertRaisesRegex(SegmentHold, "ONLY_STEP_2000"):
            next_step_guard(1999, 120000)


if __name__ == "__main__":
    unittest.main(verbosity=2)
