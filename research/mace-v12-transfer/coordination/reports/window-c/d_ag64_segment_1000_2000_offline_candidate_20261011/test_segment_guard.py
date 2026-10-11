"""Synthetic-only tests for segment_guard.py; never imports an atomistic engine."""

from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

import segment_guard as sg


def rng_state(seed: int) -> dict:
    return copy.deepcopy(np.random.Generator(np.random.PCG64(seed)).bit_generator.state)


def initial_state() -> dict:
    return {
        "owner": "D", "task_id": sg.EXPECTED_TASK, "step": 1000,
        "positions_A": [[float(i), 0.0, 0.0] for i in range(64)],
        "momenta": [[0.0, float(i) / 10.0, 0.0] for i in range(64)],
        "velocities_A_per_fs": [[0.0, float(i) / 100.0, 0.0] for i in range(64)],
        "cell_A": [[10.826824658074, 0, 0], [0, 10.826824658074, 0], [0, 0, 10.826824658074]],
        "pbc": [True, True, True], "atom_ids": list(range(1, 65)),
        "rng_state": {name: rng_state(110 + i) for i, name in enumerate(("packing", "velocity", "thermostat"))},
        "scientific_parameters_sha256": sg.EXPECTED_SCIENCE_SHA256,
        "source_manifest_sha256": sg.EXPECTED_SOURCE_MANIFEST_SHA256,
    }


def test_payload() -> dict[str, bytes]:
    state = initial_state()
    return {
        "checkpoint/state-step-0002000.json": json.dumps(state, sort_keys=True).encode(),
        "trajectory/frames.jsonl": b'{"step":2000,"synthetic":true}\n',
        "logs/thermal.csv": b"step,T_actual_K\n2000,1600.0\n",
        "source/configuration.json": b'{"synthetic_test":true}\n',
    }


class SegmentBoundaryTests(unittest.TestCase):
    def test_exact_1000_to_2000_boundary(self):
        state, saved, advanced = initial_state(), [], []

        def advance(candidate):
            advanced.append(candidate["step"])
            candidate["positions_A"][0][0] += 0.001

        result = sg.run_bounded_segment(state, advance, saved.append, stop_at_step=2000,
                                        max_added_steps=1000, max_wall_seconds=600)
        self.assertEqual(result["status"], "SYNTHETIC_SEGMENT_BOUNDARY_PASS_ONLY")
        self.assertEqual((result["start_step"], result["end_step"], len(advanced)), (1000, 2000, 1000))
        self.assertEqual(advanced[0], 1000)
        self.assertEqual(saved[-1]["step"], 2000)

    def test_step_2001_is_forbidden(self):
        with self.assertRaisesRegex(sg.GateError, "STOP_AT_STEP_REACHED"):
            sg.assert_step_allowed(2000, start_step=1000, stop_at_step=2000, max_added_steps=1000)
        with self.assertRaisesRegex(sg.GateError, "ONLY_EXACT"):
            sg.assert_step_allowed(2000, start_step=1000, stop_at_step=2001, max_added_steps=1001)

    def test_bad_start_or_step_budget_is_rejected(self):
        state = initial_state()
        state["step"] = 999
        with self.assertRaisesRegex(sg.GateError, "RESUME_STATE_OWNER_TASK_OR_STEP"):
            sg.run_bounded_segment(state, lambda s: None, lambda s: None, stop_at_step=2000,
                                  max_added_steps=1000, max_wall_seconds=600)
        with self.assertRaisesRegex(sg.GateError, "SEGMENT_LIMIT"):
            sg.run_bounded_segment(initial_state(), lambda s: None, lambda s: None, stop_at_step=120000,
                                  max_added_steps=119000, max_wall_seconds=600)
        bad = initial_state()
        bad["task_id"] = "wrong-task"
        with self.assertRaisesRegex(sg.GateError, "OWNER_TASK_OR_STEP"):
            sg.run_bounded_segment(bad, lambda s: None, lambda s: None, stop_at_step=2000,
                                  max_added_steps=1000, max_wall_seconds=600)
        bad_rng = initial_state()
        del bad_rng["rng_state"]["thermostat"]
        with self.assertRaisesRegex(sg.GateError, "RNG_STREAMS_INCOMPLETE"):
            sg.run_bounded_segment(bad_rng, lambda s: None, lambda s: None, stop_at_step=2000,
                                  max_added_steps=1000, max_wall_seconds=600)

    def test_rng_state_restoration_repeats_sequence(self):
        original = np.random.Generator(np.random.PCG64(999))
        state = copy.deepcopy(original.bit_generator.state)
        expected = original.random(32)
        restored_bitgen = np.random.PCG64()
        restored_bitgen.state = state
        restored = np.random.Generator(restored_bitgen)
        np.testing.assert_array_equal(expected, restored.random(32))

    def test_timeout_and_interrupt_leave_last_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "new-segment"
            out.mkdir()
            callback = sg.write_checkpoint_callback(out)
            callback(initial_state())
            before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()}
            ticks = iter((0.0, 1.0, 1.0))
            with self.assertRaisesRegex(sg.GateError, "WALL_LIMIT"):
                sg.run_bounded_segment(initial_state(), lambda s: None, callback, stop_at_step=2000,
                                       max_added_steps=1000, max_wall_seconds=1, clock=lambda: next(ticks))
            self.assertEqual(before, {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()})
            calls = {"n": 0}

            def cancel_after_one_checkpoint():
                calls["n"] += 1
                return calls["n"] > 1

            with self.assertRaisesRegex(sg.GateError, "CANCELLED"):
                sg.run_bounded_segment(initial_state(), lambda s: None, callback, stop_at_step=2000,
                                       max_added_steps=1000, max_wall_seconds=600,
                                       cancel_requested=cancel_after_one_checkpoint, save_every=1)
            self.assertTrue((out / "state-step-0001000.json").exists())

    def test_step_finishing_past_wall_budget_is_discarded(self):
        state, saved, calls = initial_state(), [], {"clock": 0.0, "steps": 0}

        def clock():
            return calls["clock"]

        def slow_synthetic_step(candidate):
            calls["steps"] += 1
            candidate["positions_A"][0][0] += 1.0
            calls["clock"] = 2.0

        with self.assertRaisesRegex(sg.GateError, "WALL_LIMIT_AFTER_STEP"):
            sg.run_bounded_segment(state, slow_synthetic_step, saved.append, stop_at_step=2000,
                                   max_added_steps=1000, max_wall_seconds=1, clock=clock)
        self.assertEqual(calls["steps"], 1)
        self.assertEqual(state["step"], 1000)
        self.assertEqual(state["positions_A"][0][0], 0.0)
        self.assertEqual(saved, [])

    def test_nonfinite_or_raising_step_does_not_replace_prior_state(self):
        state, saved = initial_state(), []

        def nan_step(candidate):
            candidate["positions_A"][0][0] = math.nan

        with self.assertRaisesRegex(sg.GateError, "NONFINITE_STATE"):
            sg.run_bounded_segment(state, nan_step, saved.append, stop_at_step=2000,
                                   max_added_steps=1000, max_wall_seconds=600)
        self.assertEqual(state["step"], 1000)
        self.assertEqual(saved, [])

        def interrupted_step(candidate):
            candidate["positions_A"][0][0] = 12345
            raise KeyboardInterrupt()

        with self.assertRaises(KeyboardInterrupt):
            sg.run_bounded_segment(state, interrupted_step, saved.append, stop_at_step=2000,
                                   max_added_steps=1000, max_wall_seconds=600)
        self.assertEqual(state["positions_A"][0][0], 0.0)
        self.assertEqual(saved, [])


class ArchiveTests(unittest.TestCase):
    def test_archive_manifest_sha_and_local_git_fresh_readback(self):
        archive = sg.make_archive_bytes(test_payload(), task_id=sg.EXPECTED_TASK, owner="D", final_step=2000)
        internal = sg.verify_successor_archive(archive)
        self.assertEqual(internal["final_step"], 2000)
        with tempfile.TemporaryDirectory() as tmp:
            result = sg.local_git_roundtrip(archive, tmp)
            self.assertEqual(result["status"], "LOCAL_GIT_FRESH_CLONE_READBACK_PASS_NOT_GITHUB")
            self.assertEqual(result["archive_sha256"], hashlib.sha256(archive).hexdigest())

    def test_archive_refuses_potential_and_oversize_member(self):
        with self.assertRaisesRegex(sg.GateError, "POTENTIAL_FILE"):
            sg.make_archive_bytes({"Ag.eam.alloy": b"synthetic"}, task_id=sg.EXPECTED_TASK, owner="D", final_step=2000)
        with self.assertRaisesRegex(sg.GateError, "OVER_50MIB"):
            sg.make_archive_bytes({"large.bin": b"x" * (sg.MAX_ARCHIVE_MEMBER_BYTES + 1)},
                                  task_id=sg.EXPECTED_TASK, owner="D", final_step=2000)

    def test_archive_failure_retains_only_valid_local_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            checkpoint = Path(tmp) / "step-001000.json"
            checkpoint.write_bytes(b"immutable-original-checkpoint")
            digest = sg.sha256_file(checkpoint)
            archive = sg.make_archive_bytes(test_payload(), task_id=sg.EXPECTED_TASK, owner="D", final_step=2000)
            corrupt = b"not-a-valid-gzip-or-tar-archive"
            with self.assertRaises(Exception):
                sg.verify_successor_archive(bytes(corrupt))
            self.assertTrue(checkpoint.exists())
            self.assertEqual(sg.sha256_file(checkpoint), digest)

    def test_atomic_single_file_write_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "first.bin"
            sg.atomic_write_bytes(path, b"first")
            with self.assertRaisesRegex(sg.GateError, "OVERWRITE"):
                sg.atomic_write_bytes(path, b"second")
            self.assertEqual(path.read_bytes(), b"first")

    def test_atomic_checkpoint_bundle_and_no_overwrite(self):
        members = {
            "sampler_state.json": b'{"step":2000,"synthetic":true}\n',
            "state.json": b'{"step":2000}\n',
            "state_arrays.npz": b"synthetic-array-bytes",
            "state_atoms.traj": b"synthetic-atoms-bytes",
            "thermal.csv": b"step,T\n2000,1600\n",
            "trajectory.traj": b"synthetic-trajectory",
        }
        with tempfile.TemporaryDirectory() as tmp:
            result = sg.atomic_publish_checkpoint_bundle(tmp, "segment-1000-2000", members)
            self.assertEqual(result["status"], "ATOMIC_CHECKPOINT_BUNDLE_PUBLISHED")
            final = Path(tmp) / "segment-1000-2000"
            self.assertEqual((final / "state_arrays.npz").read_bytes(), members["state_arrays.npz"])
            manifest = json.loads((final / "MANIFEST.json").read_text())
            self.assertEqual(len(manifest["members"]), len(members))
            with self.assertRaisesRegex(sg.GateError, "OVERWRITE"):
                sg.atomic_publish_checkpoint_bundle(tmp, "segment-1000-2000", members)
            self.assertEqual((final / "state.json").read_bytes(), members["state.json"])

    def test_atomic_checkpoint_bundle_failure_publishes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(sg.GateError, "OVER_50MIB"):
                sg.atomic_publish_checkpoint_bundle(tmp, "failed-stage", {"x.bin": b"x" * (sg.MAX_ARCHIVE_MEMBER_BYTES + 1)})
            self.assertFalse((Path(tmp) / "failed-stage").exists())
            self.assertFalse(any("staging" in child.name for child in Path(tmp).iterdir()))

    def test_remote_roundtrip_failure_keeps_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            checkpoint = Path(tmp) / "immutable-step1000.json"
            checkpoint.write_bytes(b"do-not-delete")
            checkpoint_sha = sg.sha256_file(checkpoint)
            payload = test_payload()
            archive = sg.make_archive_bytes(payload, task_id=sg.EXPECTED_TASK, owner="D", final_step=2000)
            root = Path(tmp) / "git-roundtrip"
            original_run = sg.subprocess.run

            def fail_fresh_clone(args, *a, **kw):
                if isinstance(args, list) and "clone" in args and any("fresh-reader" in str(x) for x in args):
                    raise sg.subprocess.CalledProcessError(1, args, stderr=b"injected readback failure")
                return original_run(args, *a, **kw)

            with patch.object(sg.subprocess, "run", side_effect=fail_fresh_clone):
                with self.assertRaises(sg.subprocess.CalledProcessError):
                    sg.local_git_roundtrip(archive, root)
            self.assertTrue(checkpoint.exists())
            self.assertEqual(sg.sha256_file(checkpoint), checkpoint_sha)


class ReleaseGateTests(unittest.TestCase):
    def _fixture(self, root: Path, status: str = "SIGNED_BY_C"):
        code = root / "segment_guard.py"
        code.write_text("synthetic execution code identity\n")
        config = {
            "release_status": status, "owner": "D", "task_id": sg.EXPECTED_TASK,
            "start_step": 1000, "stop_at_step_required": 2000, "max_added_steps": 1000,
            "stage_continue_default_endpoint_warning": 120000,
            "execution_code_sha256": {"segment_guard.py": sg.sha256_file(code)},
            "checkpoint_archive_sha256": sg.EXPECTED_ARCHIVE_SHA256,
            "checkpoint_archive_manifest_sha256": sg.EXPECTED_ARCHIVE_MANIFEST_SHA256,
            "potential_sha256_must_remain": sg.EXPECTED_POTENTIAL_SHA256,
            "scientific_parameters_sha256_must_remain": sg.EXPECTED_SCIENCE_SHA256,
        }
        config_path = root / "config.json"
        config_path.write_text(json.dumps(config, sort_keys=True) + "\n")
        lock = {
            "approval_status": status,
            "signed_commit": "synthetic-test-commit" if status == "SIGNED_BY_C" else None,
            "proposal_sha256": sg.sha256_file(config_path),
            "code_sha256": {"segment_guard.py": sg.sha256_file(code)},
        }
        lock_path = root / "release-lock.json"
        lock_path.write_text(json.dumps(lock, sort_keys=True) + "\n")
        potential = root / "potential.synthetic"
        potential.write_bytes(b"test fixture only; not a potential")
        return config_path, lock_path, code, potential

    def _run(self, config, lock, code, potential, *, owner="D", task=sg.EXPECTED_TASK, stop=2000, lock_sha=None):
        return sg.verify_release_inputs(
            proposal_path=config, release_lock_path=lock,
            archive_path="synthetic-archive", archive_manifest_path="synthetic-manifest",
            code_paths={"segment_guard.py": code}, owner=owner, task_id=task,
            requested_stop_step=stop, expected_release_lock_sha256=lock_sha or sg.sha256_file(lock),
            potential_path=potential,
        )

    def test_unsigned_proposal_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._fixture(Path(tmp), status="PROPOSED_UNSIGNED_NO_MD")
            with self.assertRaisesRegex(sg.GateError, "NOT_SIGNED"):
                self._run(*paths)

    def test_code_config_lock_task_owner_stop_sha_mismatches_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            config, lock, code, potential = self._fixture(Path(tmp))
            with self.assertRaisesRegex(sg.GateError, "RELEASE_LOCK_SHA"):
                self._run(config, lock, code, potential, lock_sha="0" * 64)
            with self.assertRaisesRegex(sg.GateError, "CALLER_OWNER_OR_TASK"):
                self._run(config, lock, code, potential, owner="E")
            with self.assertRaisesRegex(sg.GateError, "CALLER_OWNER_OR_TASK"):
                self._run(config, lock, code, potential, task="C-E-AGL-1250K-rho100-seed02")
            with self.assertRaisesRegex(sg.GateError, "EXPLICIT_STOP"):
                self._run(config, lock, code, potential, stop=None)
            original_config = config.read_bytes()
            config.write_text(config.read_text() + " ")
            with self.assertRaisesRegex(sg.GateError, "ACTUAL_CONFIG_SHA"):
                self._run(config, lock, code, potential)
            # Restore the locked fixture, then prove code is read from bytes too.
            config.write_bytes(original_config)
            code.write_text("tampered after lock\n")
            with self.assertRaisesRegex(sg.GateError, "EXECUTION_CODE_SHA"):
                self._run(config, lock, code, potential)

    def test_preflight_refuses_archive_or_potential_identity_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            config, lock, code, potential = self._fixture(Path(tmp))
            fake_archive = Path(tmp) / "archive.tar.gz"
            fake_manifest = Path(tmp) / "archive-manifest.json"
            fake_archive.write_bytes(b"different archive")
            fake_manifest.write_text("{}")
            with patch.object(sg, "verify_original_archive", side_effect=sg.GateError("ORIGINAL_ARCHIVE_SHA256_MISMATCH")):
                with self.assertRaisesRegex(sg.GateError, "ORIGINAL_ARCHIVE_SHA"):
                    sg.verify_release_inputs(proposal_path=config, release_lock_path=lock,
                        archive_path=fake_archive, archive_manifest_path=fake_manifest,
                        code_paths={"segment_guard.py": code}, owner="D", task_id=sg.EXPECTED_TASK,
                        requested_stop_step=2000, expected_release_lock_sha256=sg.sha256_file(lock),
                        potential_path=potential)

    def test_local_potential_file_is_rehashed_and_mismatch_holds(self):
        with tempfile.TemporaryDirectory() as tmp:
            config, lock, code, potential = self._fixture(Path(tmp))
            result = {"archive_sha256": sg.EXPECTED_ARCHIVE_SHA256,
                      "manifest_sha256": sg.EXPECTED_ARCHIVE_MANIFEST_SHA256}
            with patch.object(sg, "verify_original_archive", return_value=result):
                with self.assertRaisesRegex(sg.GateError, "LOCAL_POTENTIAL_SHA256_MISMATCH"):
                    self._run(config, lock, code, potential)

    def test_complete_signed_fixture_passes_only_with_actual_local_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            config, lock, code, potential = self._fixture(Path(tmp))
            fake_archive_result = {
                "status": "SYNTHETIC_ARCHIVE_BYTES_VERIFIED",
                "archive_sha256": sg.EXPECTED_ARCHIVE_SHA256,
                "manifest_sha256": sg.EXPECTED_ARCHIVE_MANIFEST_SHA256,
                "checkpoint_step": 1000, "task_id": sg.EXPECTED_TASK,
            }
            real_sha = sg.sha256_file

            def sha_with_unavailable_source_fixture(path):
                if Path(path) == potential:
                    return sg.EXPECTED_POTENTIAL_SHA256
                return real_sha(path)

            with patch.object(sg, "verify_original_archive", return_value=fake_archive_result), \
                 patch.object(sg, "sha256_file", side_effect=sha_with_unavailable_source_fixture):
                result = self._run(config, lock, code, potential)
            self.assertEqual(result["status"], "PREFLIGHT_PASS_NOT_EXECUTION_AUTHORIZATION")
            self.assertEqual(result["stop_step"], 2000)


if __name__ == "__main__":
    unittest.main(verbosity=2)
