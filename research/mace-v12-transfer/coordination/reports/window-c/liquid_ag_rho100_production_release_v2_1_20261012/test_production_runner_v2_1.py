#!/usr/bin/env python3
from __future__ import annotations

import importlib
import csv
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import Trajectory

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[5]
sys.path.insert(0, str(HERE))
import liquid_criteria_v2 as qa
import archive_backend_v2 as archive_backend
import followup_runner_v2_1 as runner

HISTORY = REPO / "research/mace-v12-transfer/coordination/reports/window-c/liquid_ag_rho100_pilot_results_20261012/C-D-AGL-1250K-rho100-pilot01/trajectory.traj"
OLD_CRITERIA_DIR = REPO / "research/mace-v12-transfer/coordination/reports/window-c/liquid_ag_rho100_followup_execution_pack_20261012"


def synthetic_counts(frame_count: int = 160) -> np.ndarray:
    edges = qa._locked_edges()
    r = (edges[:-1] + edges[1:]) / 2.0
    g = 0.65 + 1.6 * np.exp(-0.5 * ((r - 2.85) / 0.18) ** 2) - 0.45 * np.exp(-0.5 * ((r - 3.7) / 0.20) ** 2)
    g = np.maximum(g, 0.02)
    shell = 4.0 * np.pi / 3.0 * (edges[1:] ** 3 - edges[:-1] ** 3)
    per_frame_expectation = g * 64 * 63 * shell / float(qa.LOCK["finite_cell"]["volume_angstrom3"])
    return np.tile(np.rint(per_frame_expectation).astype(np.int64), (frame_count, 1))


def synthetic_metrics(frame_count: int = 160) -> list[qa.FrameMetrics]:
    counts = synthetic_counts(frame_count)
    result = []
    for i in range(frame_count):
        result.append(qa.FrameMetrics(
            step=20000 + 100 * i, time_ps=10.0 + 0.05 * i,
            temperature_K=1250.0, potential_energy_eV_per_atom=-2.0,
            global_q6=0.08,
            mean_local_q6=0.10, ptm_status="PASS", ptm_crystalline_fraction=0.0,
            ptm_largest_cluster_atoms=0, ptm_fcc_fraction=0.0, ptm_hcp_fraction=0.0,
            ptm_mean_rmsd=0.0, ptm_version="synthetic", ptm_rmsd_cutoff=0.10,
            pbc_verified=True, coordination_number=None, minimum_pair_distance_A=2.5,
            rdf_counts=tuple(int(x) for x in counts[i]), stress_eV_A3=(0.0,) * 6,
        ))
    return result


def synthetic_stress(frame_count: int = 400, start_ps: float = 14.0,
                     linear_drift: float = 0.0) -> dict:
    times = start_ps + np.arange(frame_count, dtype=float) * 0.01
    values = np.zeros((frame_count, 6), dtype=float)
    # ASE tensile-positive convention: negative diagonal stress gives positive
    # hydrostatic compression under P = -trace(sigma)/3.
    values[:, :3] = -0.001
    values[:, 0] += linear_drift * (times - times[0])
    return {"times_ps": times, "stress_eV_A3": values}


class ProductionRunnerV21Tests(unittest.TestCase):
    def test_history_q6_regression_reproduces_legacy_global_and_reports_local(self):
        self.assertTrue(HISTORY.is_file())
        sys.path.insert(0, str(OLD_CRITERIA_DIR))
        old = importlib.import_module("liquid_criteria")
        with Trajectory(str(HISTORY)) as reader:
            frames = list(reader)
        self.assertEqual(len(frames), 11)
        local = [qa.local_q6(atoms)[0] for atoms in frames]
        old_global = [old.global_q6(atoms) for atoms in frames]
        new_global = [qa.global_q6(atoms) for atoms in frames]
        np.testing.assert_allclose(new_global, old_global, rtol=1e-12, atol=1e-12)
        report = qa.q6_regression_report(new_global, local)
        self.assertAlmostEqual(float(np.mean(old_global)), 0.08871920638262232, places=10)
        self.assertAlmostEqual(float(np.mean(new_global)), 0.08871920638262232, places=10)
        self.assertAlmostEqual(report["legacy_global_q6_mean"], 0.08871920638262232, places=10)
        self.assertAlmostEqual(report["local_q6_mean_report_only"], 0.3516369060394911, places=10)
        self.assertTrue(report["regression_pass"])
        self.assertEqual(report["status"], "PASS")

    def test_q6_roles_and_strict_global_gate(self):
        frames = synthetic_metrics()
        frames[0] = qa.FrameMetrics(**{**qa.monitor_record_to_json(frames[0]), "global_q6": 0.21})
        passed = qa._window_order_gate(frames, [10.0, 12.0])
        self.assertEqual(passed["status"], "PASS")
        self.assertEqual(passed["global_q6_frame_pass_fraction"], 0.975)
        self.assertIn("mean_local_q6_report_only", passed)
        for i in range(33):
            frames[i] = qa.FrameMetrics(**{**qa.monitor_record_to_json(frames[i]), "global_q6": 0.21})
        failed = qa._window_order_gate(frames, [10.0, 12.0])
        self.assertEqual(failed["status"], "FAIL")
        # A high local diagnostic must never independently fail the global-Q6 screen.
        local_only = [qa.FrameMetrics(**{**qa.monitor_record_to_json(f), "mean_local_q6": 0.99}) for f in synthetic_metrics()]
        self.assertEqual(qa._window_order_gate(local_only, [10.0, 12.0])["status"], "PASS")

    def test_real_ovito_ptm_fcc_and_hcp_and_periodicity(self):
        fcc = bulk("Ag", "fcc", a=4.08, cubic=True).repeat((2, 2, 2))
        hcp = bulk("Ag", "hcp", a=2.89, c=4.71).repeat((2, 2, 2))
        result_fcc = qa.ptm_classify(fcc)
        result_hcp = qa.ptm_classify(hcp)
        self.assertEqual(result_fcc["status"], "PASS", result_fcc)
        self.assertEqual(result_hcp["status"], "PASS", result_hcp)
        self.assertEqual(result_fcc["fcc_fraction"], 1.0)
        self.assertEqual(result_hcp["hcp_fraction"], 1.0)
        nonperiodic = fcc.copy(); nonperiodic.pbc = (True, True, False)
        self.assertEqual(qa.ptm_classify(nonperiodic)["status"], "PERIODICITY_FAIL")
        self.assertAlmostEqual(qa.LOCK["finite_cell"]["rdf_rmax_angstrom"], 5.35)
        self.assertLess(5.35, 10.717455315295654 / 2.0)

    def test_rdf_normalization_and_bootstrap_gate_on_synthetic_data(self):
        frames = synthetic_metrics()
        result = qa.analyze_rdf_phase(frames, [[14.0, 16.0], [16.0, 18.0]])
        self.assertIn(result["status"], {"PASS", "FAIL"})
        self.assertTrue(result.get("rdf_gate_pass"), result)
        self.assertGreater(result["peak_valley_delta_one_sided_95pct_lower"], 0.0)
        self.assertTrue(9.0 <= result["coordination_number"] <= 13.0)
        self.assertEqual(result["bootstrap_replicates"], 2000)

    def test_complete_liquid_gate_and_q6_hold_path(self):
        frames = synthetic_metrics()
        rng = np.random.default_rng(4)
        increments = rng.normal(0.0, 0.12, size=(80, 64, 3))
        unwrapped = np.concatenate((np.zeros((1, 64, 3)), np.cumsum(increments, axis=0)))
        thresholds = {
            "temperature_mean_tolerance_percent": 5.0,
            "adjacent_window_mean_difference_tolerance_percent": 3.0,
            "absolute_temperature_trend_limit_K_per_ps": 25.0,
            "potential_energy_per_atom_adjacent_window_difference_limit_eV": 0.02,
        }
        passed_mock = qa.assess_phase(frames, [[14.0, 16.0], [16.0, 18.0]], unwrapped,
                                      1250.0, {"regression_pass": True}, thresholds,
                                      pressure_stress=synthetic_stress())
        self.assertEqual(passed_mock["gate_status"], "PASS", passed_mock)
        drift_mock = qa.assess_phase(frames, [[14.0, 16.0], [16.0, 18.0]], unwrapped,
                                     1250.0, {"regression_pass": True}, thresholds,
                                     pressure_stress=synthetic_stress(linear_drift=0.005))
        self.assertNotEqual(drift_mock["gate_status"], "PASS")
        self.assertEqual(drift_mock["reason"], "PRESSURE_STRESS_DRIFT")
        regression_hold = qa.assess_phase(frames, [[14.0, 16.0], [16.0, 18.0]], unwrapped,
                                          1250.0, {"regression_pass": False}, thresholds)
        self.assertEqual(regression_hold["gate_status"], "QA_HOLD")
        self.assertEqual(regression_hold["reason"], "LEGACY_GLOBAL_Q6_REGRESSION_HOLD")

    def test_insufficient_and_unresolved_rdf_fail_closed(self):
        short = synthetic_metrics(40)
        insufficient = qa.analyze_rdf_phase(short, [[14.0, 16.0], [16.0, 18.0]])
        self.assertEqual(insufficient["status"], "RDF_UNCERTAINTY_INSUFFICIENT")
        self.assertEqual(qa._rdf_point(np.zeros_like(synthetic_counts(1)))["status"], "RDF_UNRESOLVED")
        pending = qa.assess_phase(short, [[14.0, 16.0], [16.0, 18.0]], [], 1250.0,
                                  {"regression_pass": True}, {"temperature_mean_tolerance_percent": 5,
                                  "adjacent_window_mean_difference_tolerance_percent": 3,
                                  "absolute_temperature_trend_limit_K_per_ps": 25,
                                  "potential_energy_per_atom_adjacent_window_difference_limit_eV": 0.02})
        self.assertIn(pending["gate_status"], {"PENDING", "QA_HOLD"})
        self.assertNotEqual(pending["gate_status"], "PASS")

    def test_ptm_unavailable_returns_hold_status(self):
        atoms = bulk("Ag", "fcc", a=4.08, cubic=True).repeat((2, 2, 2))
        real_import = __import__
        def blocked(name, *args, **kwargs):
            if name.startswith("ovito"):
                raise ModuleNotFoundError("synthetic missing OVITO")
            return real_import(name, *args, **kwargs)
        with patch("builtins.__import__", side_effect=blocked):
            self.assertEqual(qa.ptm_classify(atoms)["status"], "PTM_UNAVAILABLE")

    def test_control_caps_and_archive_failure_prevent_advance(self):
        self.assertEqual(runner.next_chunk_end(1000, 20000), 2000)
        self.assertEqual(runner.next_chunk_end(59000, 60000), 60000)
        self.assertEqual(runner.cooling_target_K(0), 1600.0)
        self.assertEqual(runner.cooling_target_K(100), 1596.5)
        self.assertEqual(runner.cooling_target_K(10000), 1250.0)
        self.assertEqual(runner.next_chunk_end(59900, 60000), 60000)
        self.assertTrue(runner.can_start_step(1000, 2000, 3599.0))
        self.assertFalse(runner.can_start_step(1000, 2000, 3600.0))
        self.assertFalse(runner.can_start_step(1000, 2000, 3600.1))
        self.assertFalse(runner.can_start_step(2000, 2000, 3599.0))
        self.assertFalse(runner.can_start_step(60000, 60001, 10.0))
        with self.assertRaises(runner.GateError):
            runner.cooling_target_K(10001)
        good = {"archive_upload_pass": True, "archive_fresh_fetch_pass": True,
                "archive_restore_pass": True, "checkpoint_restore_pass": True,
                "sha256_match_pass": True}
        bad = dict(good, archive_restore_pass=False)
        self.assertFalse(runner._archive_allows_continuation("CHUNK_COMPLETE", bad))
        self.assertFalse(runner._archive_allows_continuation("CHUNK_COMPLETE", None))
        self.assertTrue(runner._archive_allows_continuation("CHUNK_COMPLETE", good))
        self.assertFalse(runner._archive_allows_continuation("HOLD_QA", good))
        state = {"stage": "A_HOT_1600K"}
        self.assertEqual(runner.apply_monitor_transition(state, 5000, {"gate_status": "PENDING"}), "CHUNK_COMPLETE")
        self.assertEqual(state["stage"], "A_HOT_1600K")
        self.assertEqual(runner.apply_monitor_transition(state, 5000, {"gate_status": "PASS"}), "CHUNK_COMPLETE")
        self.assertEqual(state["stage"], "B_COOL_1600_TO_1250K")
        state = {"stage": "B_COOL_1600_TO_1250K", "stage_b_start_step": 5000}
        runner.advance_completed_phase(state, 15000)
        self.assertEqual(state["stage"], "C_HOLD_1250K")
        state = {"stage": "C_HOLD_1250K", "stage_c_start_step": 15000}
        self.assertEqual(runner.apply_monitor_transition(state, 25000, {"gate_status": "PASS"},
                                                        {"snapshot_spacing_steps": 2000}), "CHUNK_COMPLETE")
        self.assertEqual(state["stage"], "D_CANDIDATE_1250K")
        self.assertEqual(state["stage_d_start_step"], 25000)
        state = {"stage": "A_HOT_1600K"}
        self.assertEqual(runner.apply_monitor_transition(state, 6000, {"gate_status": "QA_HOLD"}), "HOLD_QA_REQUIRED_INPUT_OR_METRIC")
        state = {"stage": "A_HOT_1600K"}
        self.assertEqual(runner.apply_completed_chunk_limits(state, 20000, "CHUNK_COMPLETE", 60000), "HOLD_STAGE_A_MAXIMUM_WITHOUT_LIQUID_GATE")
        state = {"stage": "C_HOLD_1250K", "stage_c_start_step": 20000}
        self.assertEqual(runner.apply_completed_chunk_limits(state, 40000, "CHUNK_COMPLETE", 60000), "HOLD_STAGE_C_MAXIMUM_WITHOUT_LIQUID_GATE")
        state = {"stage": "D_CANDIDATE_1250K", "stage_d_start_step": 50000, "last_monitor": {"gate_status": "PASS"}}
        self.assertEqual(runner.apply_completed_chunk_limits(state, 60000, "CHUNK_COMPLETE", 60000), "CANDIDATE_WINDOW_COMPLETE_FINAL_QA_PASS")
        state = {"stage": "D_CANDIDATE_1250K", "stage_d_start_step": 50000, "last_monitor": {"gate_status": "PENDING"}}
        self.assertEqual(runner.apply_completed_chunk_limits(state, 60000, "CHUNK_COMPLETE", 60000), "HOLD_CANDIDATE_FINAL_QA_FAILED")

    def test_source_checkpoint_positions_velocities_and_all_rng_streams_restore_exactly(self):
        checkpoint_path = REPO / "research/mace-v12-transfer/coordination/reports/window-c/liquid_ag_rho100_pilot_results_20261012/C-D-AGL-1250K-rho100-pilot01/checkpoint.json"
        before_sha = runner.sha256_file(checkpoint_path)
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        atoms, streams = runner._restore_atoms(checkpoint)
        np.testing.assert_array_equal(atoms.get_positions(), np.asarray(checkpoint["positions_A"]))
        np.testing.assert_allclose(atoms.get_velocities(), np.asarray(checkpoint["velocities_A_per_fs"]), rtol=1e-14, atol=0.0)
        self.assertEqual(set(streams), {"packing", "velocity_initialization", "thermostat"})
        for name, stream in streams.items():
            self.assertEqual(stream.bit_generator.state, checkpoint["rng_states"][name])
            expected = np.random.Generator(np.random.PCG64())
            expected.bit_generator.state = checkpoint["rng_states"][name]
            self.assertEqual(stream.random(), expected.random())
        self.assertEqual(before_sha, runner.sha256_file(checkpoint_path))

    def test_production_remote_boundary_restore_preserves_rng_and_state(self):
        checkpoint_path = REPO / "research/mace-v12-transfer/coordination/reports/window-c/liquid_ag_rho100_pilot_results_20261012/C-D-AGL-1250K-rho100-pilot01/checkpoint.json"
        source = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        live_atoms, live_rngs = runner._restore_atoms(source, expected_global_step=1000)
        state = {"stage": "A_HOT_1600K", "global_step": 2000, "liquid_onset_step": None,
                 "last_monitor": {"gate_status": "PENDING"}}
        remote_copy = dict(source)
        remote_copy["global_step"] = 2000
        remote_copy["stage_state"] = dict(state)
        restored_atoms, restored_rngs, dynamics = runner._restore_verified_segment(
            remote_copy, live_atoms, live_rngs, state, calculator=None
        )
        np.testing.assert_array_equal(restored_atoms.get_positions(), live_atoms.get_positions())
        np.testing.assert_allclose(restored_atoms.get_momenta(), live_atoms.get_momenta(), rtol=1e-14, atol=0.0)
        for name in live_rngs:
            self.assertEqual(restored_rngs[name].bit_generator.state, live_rngs[name].bit_generator.state)
        self.assertEqual(dynamics.atoms, restored_atoms)
        corrupt = dict(remote_copy)
        corrupt["rng_states"] = dict(remote_copy["rng_states"])
        corrupt["rng_states"]["thermostat"] = {"bit_generator": "PCG64", "state": {"state": 1, "inc": 1}, "has_uint32": 0, "uinteger": 0}
        with self.assertRaises(runner.GateError):
            runner._restore_verified_segment(corrupt, live_atoms, live_rngs, state, calculator=None)

    def test_packaged_lock_is_preparation_only(self):
        self.assertFalse(qa.LOCK["execution_authorized"])
        self.assertEqual(qa.SCIENCE_LOCK_SHA256, "848f8db72d83cea91dcf4cd42ee9fbbf49a80e7b5e3759cb4f535be8cc5873d2")

    def test_pressure_stress_constant_pass_pressure_sign_units_and_summary(self):
        result = qa.pressure_stress_gate(**synthetic_stress())
        self.assertEqual(result["status"], "PASS", result)
        self.assertTrue(result["pressure_stress_gate_pass"])
        self.assertAlmostEqual(result["series"]["pressure"]["mean"], 0.001, places=14)
        self.assertAlmostEqual(result["pressure_mean_GPa"], 0.16021766208, places=10)
        self.assertEqual(result["stress_unit"], "eV/Angstrom^3")
        self.assertEqual(result["series"]["pressure"]["sample_count"], 400)
        self.assertAlmostEqual(result["series"]["pressure"]["slope_per_ps"], 0.0, places=14)
        self.assertEqual(result["series"]["pressure"]["slope_ci95_per_ps"], [0.0, 0.0])

    def test_pressure_stress_known_linear_drift_fails_closed(self):
        result = qa.pressure_stress_gate(**synthetic_stress(linear_drift=0.005))
        self.assertEqual(result["status"], "FAIL_DRIFT", result)
        self.assertFalse(result["pressure_stress_gate_pass"])
        self.assertTrue(result["series"]["pressure"]["statistically_resolved_drift"])
        self.assertGreater(result["series"]["stress_xx"]["slope_ci95_per_ps"][0], 0.0)

    def test_pressure_stress_correlated_noise_uses_block_uncertainty(self):
        rng = np.random.default_rng(20261012)
        noise = np.zeros(2000)
        innovations = rng.normal(0.0, 1.0, size=len(noise))
        for i in range(1, len(noise)):
            noise[i] = 0.9 * noise[i - 1] + innovations[i]
        times = np.arange(len(noise), dtype=float) * 0.01
        summary = qa._stress_observable_summary(times, noise, "correlated_fixture", 0.01, 2000)
        iid_half_width = 1.96 * float(np.std(noise, ddof=1)) / math.sqrt(len(noise))
        self.assertGreater(summary["block_length_samples"], 1)
        self.assertGreater(summary["complete_effective_blocks"], 8)
        self.assertGreater(summary["mean_block_uncertainty_95pct_half_width"], iid_half_width)

    def test_pressure_stress_insufficient_missing_nan_and_bad_components_hold(self):
        short = synthetic_stress(frame_count=300)
        self.assertEqual(qa.pressure_stress_gate(**short)["status"], "PENDING")
        innovations = np.random.default_rng(912).normal(size=400)
        correlated = np.zeros(400)
        for i in range(1, len(correlated)):
            correlated[i] = 0.99 * correlated[i - 1] + innovations[i]
        times = np.arange(400, dtype=float) * 0.01
        underblocked = qa.pressure_stress_gate(times, np.repeat(correlated[:, None], 6, axis=1))
        self.assertEqual(underblocked["status"], "PENDING", underblocked)
        self.assertLess(underblocked["series"]["pressure"]["complete_effective_blocks"], 8)
        wrong_cadence = qa.pressure_stress_gate(np.arange(400, dtype=float) * 0.05,
                                                np.zeros((400, 6), dtype=float))
        self.assertEqual(wrong_cadence["status"], "QA_HOLD")
        self.assertEqual(qa.pressure_stress_gate(None, None)["status"], "QA_HOLD")
        missing = np.zeros((400, 5))
        self.assertEqual(qa.pressure_stress_gate(np.arange(400) * 0.01, missing)["status"], "QA_HOLD")
        nan_values = synthetic_stress()
        nan_values["stress_eV_A3"][5, 2] = np.nan
        self.assertEqual(qa.pressure_stress_gate(**nan_values)["status"], "QA_HOLD")

    def test_pressure_stress_restored_archive_matches_continuous_result(self):
        with tempfile.TemporaryDirectory(prefix="ag64-stress-archive-") as td:
            base = Path(td)
            bare, repo = base / "remote.git", base / "repo"
            source, restored = base / "source", base / "restored"
            subprocess.run(["git", "init", "--bare", str(bare)], check=True, capture_output=True)
            subprocess.run(["git", "init", "-b", "cloud-d-b1-dft", str(repo)], check=True, capture_output=True)
            subprocess.run(["git", "-C", str(repo), "config", "user.name", "QA Test"], check=True)
            subprocess.run(["git", "-C", str(repo), "config", "user.email", "qa-test@example.invalid"], check=True)
            (repo / "README.md").write_text("synthetic stress archive test\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(repo), "add", "README.md"], check=True)
            subprocess.run(["git", "-C", str(repo), "commit", "-m", "fixture"], check=True, capture_output=True)
            subprocess.run(["git", "-C", str(repo), "remote", "add", "origin", str(bare)], check=True)
            subprocess.run(["git", "-C", str(repo), "push", "-u", "origin", "cloud-d-b1-dft"], check=True, capture_output=True)
            head = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
            source.mkdir()
            (source / "checkpoint.json").write_text(json.dumps({
                "global_step": 36000, "positions_A": [[0.0, 0.0, 0.0]],
                "momenta_amu_A_per_fs": [[0.0, 0.0, 0.0]], "rng_states": {"test": {"state": 3}},
            }, sort_keys=True) + "\n", encoding="utf-8")
            fields = runner.STRESS_CSV_COLUMNS
            thermal_fields = runner.THERMAL_CSV_COLUMNS
            data = synthetic_stress()
            with (source / "stress.csv").open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream); writer.writerow(fields)
                for i, (time_ps, row) in enumerate(zip(data["times_ps"], data["stress_eV_A3"])):
                    step = 28000 + 20 * i
                    writer.writerow([step, step * 0.5, *row])
            with (source / "thermal.csv").open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream); writer.writerow(thermal_fields)
                for i in range(len(data["times_ps"])):
                    step = 28000 + 20 * i
                    writer.writerow([step, step * 0.5, 1250.0, 1250.0, -2.0, 0.1, -1.9])
            live_input = runner._eligible_pressure_stress_samples(source, 1250.0, 27999, 36000)
            live_result = qa.pressure_stress_gate(live_input["times_ps"], live_input["stress_eV_A3"])
            archived = archive_backend.archive_push_fresh_restore(
                repo_root=repo, source_dir=source, archive_relpath="synthetic/stress/step_36000",
                branch="cloud-d-b1-dft", remote="origin", expected_remote_head=head,
                restore_dir=restored, commit_message="synthetic stress archive QA",
            )
            self.assertTrue(archive_backend.segment_archive_passed(archived), archived)
            restored_input = runner._eligible_pressure_stress_samples(restored, 1250.0, 27999, 36000)
            restored_result = qa.pressure_stress_gate(restored_input["times_ps"], restored_input["stress_eV_A3"])
            self.assertEqual(live_input, restored_input)
            self.assertEqual(live_result, restored_result)

    def test_synthetic_git_archive_fresh_fetch_restore(self):
        with tempfile.TemporaryDirectory(prefix="ag64-qa-v2-archive-") as td:
            base = Path(td)
            bare = base / "remote.git"
            repo = base / "repo"
            source = base / "bundle"
            restore = base / "restored"
            subprocess.run(["git", "init", "--bare", str(bare)], check=True, capture_output=True)
            subprocess.run(["git", "init", "-b", "cloud-d-b1-dft", str(repo)], check=True, capture_output=True)
            subprocess.run(["git", "-C", str(repo), "config", "user.name", "QA Test"], check=True)
            subprocess.run(["git", "-C", str(repo), "config", "user.email", "qa-test@example.invalid"], check=True)
            (repo / "README.md").write_text("synthetic archive test only\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(repo), "add", "README.md"], check=True)
            subprocess.run(["git", "-C", str(repo), "commit", "-m", "fixture"], check=True, capture_output=True)
            subprocess.run(["git", "-C", str(repo), "remote", "add", "origin", str(bare)], check=True)
            subprocess.run(["git", "-C", str(repo), "push", "-u", "origin", "cloud-d-b1-dft"], check=True, capture_output=True)
            base_head = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
            source.mkdir()
            checkpoint = {
                "global_step": 1000, "positions_A": [[0.0, 0.0, 0.0]],
                "momenta_amu_A_per_fs": [[0.0, 0.0, 0.0]], "rng_states": {"synthetic": {"state": 7}},
            }
            (source / "checkpoint.json").write_text(json.dumps(checkpoint, sort_keys=True), encoding="utf-8")
            (source / "trajectory.mock").write_bytes(b"synthetic trajectory member\n")
            (source / "README.md").write_text("No scientific simulation was run.\n", encoding="utf-8")
            result = archive_backend.archive_push_fresh_restore(
                repo_root=repo, source_dir=source,
                archive_relpath="archive/test/step_01000", branch="cloud-d-b1-dft",
                remote="origin", expected_remote_head=base_head, restore_dir=restore,
                commit_message="synthetic QA archive test",
            )
            self.assertTrue(archive_backend.segment_archive_passed(result))
            self.assertEqual(result["member_count"], 3)
            self.assertEqual(json.loads((restore / "checkpoint.json").read_text()), checkpoint)
            forbidden = base / "forbidden"
            forbidden.mkdir()
            (forbidden / "Ag.eam.alloy").write_text("never archive a potential", encoding="ascii")
            with self.assertRaises(archive_backend.ArchiveError):
                archive_backend.write_sha256sums(forbidden)


if __name__ == "__main__":
    unittest.main(verbosity=2)
