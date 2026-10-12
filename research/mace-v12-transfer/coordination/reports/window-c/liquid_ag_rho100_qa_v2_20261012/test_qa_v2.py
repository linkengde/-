#!/usr/bin/env python3
from __future__ import annotations

import importlib
import json
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
from adaptive_control_v2 import advance_after_archive, archive_gate, can_start_next_step, next_chunk_end

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
            mean_local_q6=0.10, ptm_status="PASS", ptm_crystalline_fraction=0.0,
            ptm_largest_cluster_atoms=0, ptm_fcc_fraction=0.0, ptm_hcp_fraction=0.0,
            ptm_mean_rmsd=0.0, ptm_version="synthetic", ptm_rmsd_cutoff=0.10,
            pbc_verified=True, coordination_number=None, minimum_pair_distance_A=2.5,
            rdf_counts=tuple(int(x) for x in counts[i]), stress_eV_A3=(0.0,) * 6,
        ))
    return result


class QAV2Tests(unittest.TestCase):
    def test_history_q6_regression_holds_without_tuning(self):
        self.assertTrue(HISTORY.is_file())
        sys.path.insert(0, str(OLD_CRITERIA_DIR))
        old = importlib.import_module("liquid_criteria")
        with Trajectory(str(HISTORY)) as reader:
            frames = list(reader)
        self.assertEqual(len(frames), 11)
        local = [qa.local_q6(atoms)[0] for atoms in frames]
        old_global = [old.global_q6(atoms) for atoms in frames]
        report = qa.q6_regression_report(local, old_global)
        self.assertAlmostEqual(report["v2_mean_frame_Q6"], 0.3516369060394911, places=10)
        self.assertAlmostEqual(report["old_global_Q6_mean_for_diagnosis_only"], 0.08871920638262232, places=10)
        self.assertFalse(report["regression_pass"])
        self.assertEqual(report["status"], "Q6_THRESHOLD_COMPARABILITY_HOLD")

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
                                      1250.0, {"regression_pass": True}, thresholds)
        self.assertEqual(passed_mock["gate_status"], "PASS", passed_mock)
        regression_hold = qa.assess_phase(frames, [[14.0, 16.0], [16.0, 18.0]], unwrapped,
                                          1250.0, {"regression_pass": False}, thresholds)
        self.assertEqual(regression_hold["gate_status"], "QA_HOLD")
        self.assertEqual(regression_hold["reason"], "Q6_THRESHOLD_COMPARABILITY_HOLD")

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
        self.assertEqual(next_chunk_end(1000, 20000), 2000)
        self.assertEqual(next_chunk_end(59000, 60000), 60000)
        self.assertTrue(can_start_next_step(1000, 2000, 3599.0))
        self.assertFalse(can_start_next_step(1000, 2000, 3600.0))
        self.assertFalse(can_start_next_step(1000, 2000, 3600.1))
        self.assertFalse(can_start_next_step(2000, 2000, 3599.0))
        good = {"archive_upload_pass": True, "archive_fresh_fetch_pass": True,
                "archive_restore_pass": True, "checkpoint_restore_pass": True,
                "sha256_match_pass": True}
        bad = dict(good, archive_restore_pass=False)
        self.assertFalse(archive_gate(bad))
        state = {"stage": "A_HOT_1600K"}
        status = advance_after_archive(state, 5000, "PASS", bad)
        self.assertEqual(status, "HOLD_ARCHIVE_OR_RESTORE_FAILURE")
        self.assertIsNone(state.get("stage_b_start_step"))
        state = {"stage": "A_HOT_1600K"}
        self.assertEqual(advance_after_archive(state, 5000, "PENDING", good), "CHUNK_COMPLETE")
        self.assertEqual(state["stage"], "A_HOT_1600K")
        self.assertEqual(advance_after_archive(state, 5000, "PASS", good), "CHUNK_COMPLETE")
        self.assertEqual(state["stage"], "B_COOL_1600_TO_1250K")
        state = {"stage": "B_COOL_1600_TO_1250K", "stage_b_start_step": 5000}
        advance_after_archive(state, 15000, "PENDING", good)
        self.assertEqual(state["stage"], "C_HOLD_1250K")
        state = {"stage": "A_HOT_1600K"}
        self.assertEqual(advance_after_archive(state, 6000, "FAIL_WINDOW", good), "HOLD_QA_FAIL_WINDOW")
        state = {"stage": "A_HOT_1600K"}
        self.assertEqual(advance_after_archive(state, 20000, "PENDING", good), "HOLD_STAGE_A_CAP_WITHOUT_LIQUID_GATE")
        state = {"stage": "C_HOLD_1250K", "stage_c_start_step": 20000}
        self.assertEqual(advance_after_archive(state, 40000, "PENDING", good), "HOLD_STAGE_C_CAP_WITHOUT_LIQUID_GATE")
        state = {"stage": "D_CANDIDATE_1250K", "stage_d_start_step": 50000}
        self.assertEqual(advance_after_archive(state, 60000, "PENDING", good), "CANDIDATE_WINDOW_COMPLETE_FINAL_QA_REQUIRED")

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
