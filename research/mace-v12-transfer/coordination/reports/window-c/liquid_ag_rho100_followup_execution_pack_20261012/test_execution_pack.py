#!/usr/bin/env python3
"""Non-scientific tests; no structure generation, calculator, or MD is used."""
from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

import archive_backend
import followup_runner
import liquid_criteria


class AnalysisTests(unittest.TestCase):
    def test_autocorrelation_and_bootstrap_are_finite_and_repeatable(self):
        values = np.sin(np.arange(64) * 0.3) + np.arange(64) * 0.01
        self.assertGreaterEqual(liquid_criteria.integrated_autocorrelation_time(values), 0.5)
        self.assertEqual(
            liquid_criteria.block_bootstrap_mean_ci(values),
            liquid_criteria.block_bootstrap_mean_ci(values),
        )

    def test_multiple_origin_msd_uses_no_periodic_wrapping_or_randomness(self):
        t = np.arange(21, dtype=float) * 0.05
        x = np.zeros((len(t), 2, 3), dtype=float)
        x[:, 0, 0] = t
        x[:, 1, 0] = 2.0 * t
        lag, msd = liquid_criteria.multiple_origin_msd(x, 0.05)
        self.assertEqual(len(lag), 10)
        self.assertAlmostEqual(msd[0], (0.05**2 + 0.1**2) / 2)
        self.assertTrue(np.isfinite(msd).all())

    def test_rdf_shape_policy_is_explicit_and_honored(self):
        r = np.arange(0.025, 5.4, 0.05)
        g = np.exp(-((r - 2.8) / 0.30) ** 2) + 0.05
        frames = [
            liquid_criteria.FrameMetrics(
                step=i * 100, temperature_K=1600.0,
                potential_energy_eV_per_atom=-2.0, global_q6=0.1,
                crystalline_fraction=0.0, largest_crystalline_cluster_atoms=0,
                coordination_number=11.0, rdf_first_peak_A=2.8,
                rdf_first_minimum_A=3.75, minimum_pair_distance_A=2.0,
                rdf_r_A=tuple(r), rdf_g=tuple(g), stress_eV_A3=(0.0,) * 6,
            ) for i in range(8)
        ]
        policy = {
            "minimum_first_peak_fwhm_A": 0.5,
            "secondary_search_start_A": 4.2,
            "maximum_secondary_peak_prominence_fraction": 0.2,
        }
        out = liquid_criteria.rdf_shape_assess(frames, policy)
        self.assertTrue(out["rdf_shape_pass"])
        self.assertGreater(out["first_peak_fwhm_A"], 0.5)


class ScheduleTests(unittest.TestCase):
    def test_cooling_boundaries_and_segment_cap(self):
        self.assertEqual(followup_runner.cooling_target_K(0), 1600.0)
        self.assertEqual(followup_runner.cooling_target_K(100), 1596.5)
        self.assertEqual(followup_runner.cooling_target_K(9900), 1253.5)
        self.assertEqual(followup_runner.cooling_target_K(10000), 1250.0)
        self.assertEqual(followup_runner.next_chunk_end(1000, 22000), 2000)
        self.assertEqual(followup_runner.next_chunk_end(19500, 20000), 20000)
        self.assertEqual(followup_runner.next_chunk_end(59900, 70000), 60000)
        with self.assertRaises(followup_runner.GateError):
            followup_runner.cooling_target_K(10001)

    def test_automatic_chunk_sequence_never_crosses_phase_or_global_boundary(self):
        step = 1000
        ends = []
        while step < 4500:
            step = followup_runner.next_chunk_end(step, 4500)
            ends.append(step)
        self.assertEqual(ends, [2000, 3000, 4000, 4500])
        self.assertLessEqual(max(b - a for a, b in zip([1000] + ends[:-1], ends)), 1000)
        self.assertEqual(followup_runner.next_chunk_end(59000, 70000), 60000)

    def test_frozen_phase_transitions_are_automatic_and_bounded(self):
        state = {"stage": "A_HOT_1600K"}
        self.assertEqual(
            followup_runner.apply_monitor_transition(state, 9000, {"gate_status": "PASS"}),
            "CHUNK_COMPLETE",
        )
        self.assertEqual(state["stage"], "B_COOL_1600_TO_1250K")
        self.assertEqual(state["liquid_onset_step"], 9000)
        followup_runner.advance_completed_phase(state, 19000)
        self.assertEqual(state["stage"], "C_HOLD_1250K")
        self.assertEqual(state["stage_c_start_step"], 19000)
        followup_runner.apply_monitor_transition(
            state, 27000, {"gate_status": "PASS"}, {"snapshot_spacing_steps": 2000}
        )
        self.assertEqual(state["stage"], "D_CANDIDATE_1250K")
        self.assertEqual(state["candidate_snapshot_spacing_steps"], 2000)
        status = followup_runner.apply_monitor_transition(
            state, 28000, {"gate_status": "FAIL_WINDOW"}
        )
        self.assertEqual(status, "HOLD_CANDIDATE_WINDOW_QA_FAILED")
        self.assertEqual(state["stage"], status)

    def test_phase_or_global_limit_stops_without_extending(self):
        state = {"stage": "A_HOT_1600K"}
        self.assertEqual(
            followup_runner.apply_completed_chunk_limits(state, 20000, "CHUNK_COMPLETE", 60000),
            "HOLD_STAGE_A_MAXIMUM_WITHOUT_LIQUID_GATE",
        )
        state = {"stage": "C_HOLD_1250K", "stage_c_start_step": 40000}
        self.assertEqual(
            followup_runner.apply_completed_chunk_limits(state, 60000, "CHUNK_COMPLETE", 60000),
            "HOLD_STAGE_C_MAXIMUM_WITHOUT_LIQUID_GATE",
        )
        state = {"stage": "D_CANDIDATE_1250K", "stage_d_start_step": 50000,
                 "last_monitor": {"gate_status": "PASS"}}
        self.assertEqual(
            followup_runner.apply_completed_chunk_limits(state, 60000, "CHUNK_COMPLETE", 60000),
            "CANDIDATE_WINDOW_COMPLETE_FINAL_QA_PASS",
        )

    def test_archive_failure_and_abnormal_chunk_exception_stop_continuation(self):
        passed = {
            "archive_upload_pass": True,
            "archive_fresh_fetch_pass": True,
            "archive_restore_pass": True,
            "checkpoint_restore_pass": True,
            "sha256_match_pass": True,
        }
        failed = {**passed, "archive_fresh_fetch_pass": False}
        starts = []
        current = 1000
        result = failed
        while current < 5000:
            starts.append(current)
            current = followup_runner.next_chunk_end(current, 5000)
            if not archive_backend.segment_archive_passed(result):
                break
            result = passed
        self.assertEqual(starts, [1000])

        starts = []
        current = 1000
        try:
            while current < 5000:
                starts.append(current)
                current = followup_runner.next_chunk_end(current, 5000)
                if current == 3000:
                    raise OSError("synthetic interrupted segment")
                self.assertTrue(archive_backend.segment_archive_passed(passed))
        except OSError:
            pass
        self.assertEqual(starts, [1000, 2000])
        self.assertEqual(current, 3000)

    def test_archived_step1000_checkpoint_restores_without_calculation(self):
        repo = followup_runner._git_root(Path(__file__).resolve().parent)
        checkpoint_path = repo / "research/mace-v12-transfer/coordination/reports/window-c/liquid_ag_rho100_pilot_results_20261012/C-D-AGL-1250K-rho100-pilot01/checkpoint.json"
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        atoms, streams = followup_runner._restore_atoms(checkpoint)
        self.assertEqual(len(atoms), 64)
        self.assertEqual(atoms.get_chemical_formula(), "Ag64")
        self.assertEqual(set(streams), {"packing", "velocity_initialization", "thermostat"})
        for name, stream in streams.items():
            self.assertEqual(stream.bit_generator.state, checkpoint["rng_states"][name])


class ArchiveTests(unittest.TestCase):
    def test_missing_or_failed_remote_gate_never_passes(self):
        self.assertFalse(archive_backend.segment_archive_passed(None))
        self.assertFalse(archive_backend.segment_archive_passed({"archive_upload_pass": True}))
        self.assertFalse(archive_backend.segment_archive_passed({
            "archive_upload_pass": True,
            "archive_fresh_fetch_pass": True,
            "archive_restore_pass": True,
            "checkpoint_restore_pass": False,
            "sha256_match_pass": True,
        }))

    def test_local_bare_remote_push_fresh_fetch_and_restore(self):
        with tempfile.TemporaryDirectory(prefix="ag64-archive-test-") as tmp_name:
            root = Path(tmp_name)
            repo = root / "repo"
            remote = root / "remote.git"
            source = root / "bundle"
            restore = root / "restored"
            repo.mkdir()
            source.mkdir()
            subprocess.run(["git", "init", "-q", str(repo)], check=True)
            subprocess.run(["git", "-C", str(repo), "config", "user.name", "Test"], check=True)
            subprocess.run(["git", "-C", str(repo), "config", "user.email", "test@example.invalid"], check=True)
            (repo / "base.txt").write_text("base\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(repo), "add", "base.txt"], check=True)
            subprocess.run(["git", "-C", str(repo), "commit", "-qm", "base"], check=True)
            subprocess.run(["git", "init", "--bare", "-q", str(remote)], check=True)
            subprocess.run(["git", "-C", str(repo), "branch", "-M", "test-owner"], check=True)
            subprocess.run(["git", "-C", str(repo), "remote", "add", "origin", str(remote)], check=True)
            subprocess.run(["git", "-C", str(repo), "push", "-q", "-u", "origin", "test-owner"], check=True)
            (source / "README.md").write_text("synthetic archive fixture\n", encoding="utf-8")
            (source / "manifest.json").write_text(json.dumps({"fixture_only": True}) + "\n", encoding="utf-8")
            (source / "checkpoint.json").write_text(json.dumps({
                "global_step": 7,
                "positions_A": [[0.0, 0.0, 0.0]],
                "momenta_amu_A_per_fs": [[0.0, 0.0, 0.0]],
                "rng_states": {"fixture": {"state": 1}},
            }) + "\n", encoding="utf-8")
            (source / "trajectory_sim.txt").write_bytes(b"not a trajectory\n")
            head = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip()
            result = archive_backend.archive_push_fresh_restore(
                repo_root=repo,
                source_dir=source,
                archive_relpath="reports/test-bundle",
                branch="test-owner",
                remote="origin",
                expected_remote_head=head,
                restore_dir=restore,
                commit_message="test synthetic archive",
            )
            self.assertTrue(result["archive_upload_pass"])
            self.assertTrue(result["archive_fresh_fetch_pass"])
            self.assertTrue(result["archive_restore_pass"])
            self.assertTrue(result["checkpoint_restore_pass"])
            self.assertTrue(result["sha256_match_pass"])
            self.assertEqual(archive_backend.verify_sha256sums(restore), archive_backend.verify_sha256sums(repo / "reports/test-bundle"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
