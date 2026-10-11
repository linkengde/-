#!/usr/bin/env python3
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import codex_release_gate as gate
import archive_checkpoint
import numpy as np

BASE = Path(__file__).parent
CONFIG_PATH = BASE / "CODEX_WARMUP_RELEASE_CONFIG_PRE_RELEASE_HOLD.json"
SIGNED_CONFIG_PATH = BASE / "CODEX_WARMUP_RELEASE_CONFIG.json"
EXPECTED_POTENTIAL = gate.EXPECTED_POTENTIAL_SHA256


def ready_config(tmp: Path) -> dict:
    config = json.loads(CONFIG_PATH.read_text())
    config.update({"config_status": "C_SIGNED_FOR_WARMUP",
                   "release_decision": "APPROVE_WARMUP_ONLY",
                   "source_release": True, "density_reference_pass": True,
                   "storage_release": True, "warmup_release": True,
                   "continuation_release": False,
                   "warmup_task_ids": ["C-D-AGL-1250K-rho097-seed01"]})
    # Bind all release roots into the writable temporary test area, outside the checkout.
    for owner in ("D", "E", "F"):
        config["hosts"][owner]["sampler_work_root"] = str(tmp / ("sampler-" + owner))
        config["hosts"][owner]["archive_staging_root"] = str(tmp / ("archive-" + owner))
        config["hosts"][owner]["C_host_release"] = True
        config["hosts"][owner]["machine_fingerprint_sha256"] = "a" * 64
    config["hosts"]["D"]["remote_archive"]["roundtrip_status"] = "PASS"
    config["hosts"]["D"]["remote_archive"]["synthetic_fixture_verified"] = True
    config["hosts"]["D"]["remote_archive"]["expected_remote_head"] = "1" * 40
    config["hosts"]["D"]["potential_path"] = "/nonexistent/Ag.eam.alloy"
    config["hosts"]["D"]["potential_local_sha_verified"] = True
    config["hosts"]["D"]["potential_local_sha_receipt"] = {"sha256": EXPECTED_POTENTIAL}
    return config


def invoke(config, *, owner="D", task="C-D-AGL-1250K-rho097-seed01", head="1" * 40, potential=None):
    return gate.validate_ephemeral_release(
        config, owner=owner, task_id=task, stage="warmup", repository_root=Path("/workspace/-"),
        current_fingerprint="a" * 64, potential_path=potential, remote_branch_head=head)


class ArchiveIntegrationTests(unittest.TestCase):
    def test_non_scientific_local_git_archive_fresh_fetch_restore(self):
        # This exercises archive code only with synthetic bytes and a local bare Git remote.
        with tempfile.TemporaryDirectory(prefix="codex-archive-smoke-", dir="/workspace") as tmp_value:
            tmp = Path(tmp_value)
            bare = tmp / "remote.git"
            subprocess.run(["git", "init", "--bare", "--initial-branch=cloud-d-b1-dft", str(bare)],
                           check=True, capture_output=True)
            seed = tmp / "seed"
            seed.mkdir()
            subprocess.run(["git", "init", "--initial-branch=cloud-d-b1-dft"], cwd=seed,
                           check=True, capture_output=True)
            subprocess.run(["git", "config", "user.name", "test"], cwd=seed, check=True, capture_output=True)
            subprocess.run(["git", "config", "user.email", "test@example.invalid"], cwd=seed,
                           check=True, capture_output=True)
            (seed / "seed.txt").write_text("non-science branch seed\n")
            subprocess.run(["git", "add", "seed.txt"], cwd=seed, check=True, capture_output=True)
            subprocess.run(["git", "commit", "-m", "seed"], cwd=seed, check=True, capture_output=True)
            subprocess.run(["git", "remote", "add", "origin", str(bare)], cwd=seed, check=True, capture_output=True)
            subprocess.run(["git", "push", "-u", "origin", "cloud-d-b1-dft"], cwd=seed,
                           check=True, capture_output=True)
            expected_head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=seed,
                                           check=True, capture_output=True, text=True).stdout.strip()
            task_id = "C-D-AGL-1250K-rho097-seed01"
            config_sha = "a" * 64
            task = tmp / "task" / task_id
            bundle = task / "restart" / "step-000001000-test"
            (task / "logs").mkdir(parents=True)
            bundle.mkdir(parents=True)
            (task / "source_manifest.json").write_text(json.dumps({"task_id": task_id, "origin": "synthetic-test"}) + "\n")
            (task / "trajectory.traj").write_bytes(b"synthetic trajectory payload\n")
            (task / "thermal.csv").write_text("step,time_fs,T_K\n1000,500,1250\n")
            (task / "logs" / "warmup-performance-test.json").write_text("{}\n")
            (task / "stage_wall_budget.json").write_text("{}\n")
            (bundle / "state_atoms.traj").write_bytes(b"synthetic atom state\n")
            (bundle / "trajectory.traj").write_bytes((task / "trajectory.traj").read_bytes())
            (bundle / "thermal.csv").write_bytes((task / "thermal.csv").read_bytes())
            arrays_path = bundle / "state_arrays.npz"
            np.savez_compressed(arrays_path, numbers=np.ones(64, dtype=np.int32),
                                positions=np.zeros((64, 3)), velocities_A_per_fs=np.zeros((64, 3)),
                                momenta=np.zeros((64, 3)), masses=np.ones(64),
                                cell=np.eye(3), pbc=np.ones(3, dtype=bool),
                                atom_ids=np.arange(1, 65))
            import hashlib
            file_hashes = {name: hashlib.sha256((bundle / name).read_bytes()).hexdigest()
                           for name in ("state_atoms.traj", "state_arrays.npz", "trajectory.traj", "thermal.csv")}
            state = {"task_id": task_id, "owner": "D", "step": 1000,
                     "status": "WARMUP_COMPLETE_C_REVIEW_REQUIRED", "time_fs": 500.0,
                     "restart_bundle": "restart/step-000001000-test",
                     "rng_state": {"packing": {"seed": 1}, "velocity": {"seed": 2}, "thermostat": {"seed": 3}},
                     "bundle_sha256": file_hashes, "C_release_config_sha256_at_checkpoint": config_sha,
                     "task_matrix_sha256": "b" * 64,
                     "source_manifest_sha256": hashlib.sha256((task / "source_manifest.json").read_bytes()).hexdigest()}
            (bundle / "state.json").write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")
            (task / "sampler_state.json").write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")
            config = {"repository": {"url": str(bare)}, "tasks_sha256": "b" * 64,
                      "potential": {"sha256": EXPECTED_POTENTIAL},
                      "sampling": {"min_free_bytes": 1, "checkpoint_reserve_bytes": 1},
                      "codex_archive_policy": {"max_single_file_bytes": 52428800,
                                               "max_archive_bytes": 52428800}}
            host = {"owner": "D", "remote_archive": {"branch": "cloud-d-b1-dft",
                    "expected_remote_head": expected_head, "synthetic_fixture_verified": True}}
            receipt = archive_checkpoint.archive_and_restore(
                task, tmp / "staging", config, host, task_id, config_sha,
                {"run_sampling_codex.py": "c" * 64})
            self.assertTrue(receipt["STORAGE_RELEASE_PROOF"].startswith("PASS"))
            self.assertTrue(receipt["outer_sha_pass"] and receipt["inner_sha_pass"])
            self.assertGreaterEqual(receipt["final_restore_member_count"], 10)


class GateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="codex-release-gate-test-", dir="/workspace")
        self.root = Path(self.tmp.name)
        self.config = ready_config(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def assertRefused(self, expected, fn):
        with self.assertRaises(gate.GateRefusal) as caught:
            fn()
        self.assertEqual(str(caught.exception), expected)

    def test_actual_runner_refuses_unsigned_config_before_science_or_network(self):
        config_sha = hashlib.sha256(CONFIG_PATH.read_bytes()).hexdigest()
        sampler_root = Path("/workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot")
        before = sampler_root.exists()
        script = BASE / "run_sampling_codex.py"
        result = subprocess.run([
            sys.executable, str(script), "--config", str(CONFIG_PATH),
            "--config-sha256", config_sha, "--owner", "D",
            "--task-id", "C-D-AGL-1250K-rho097-seed01", "--stage", "warmup",
        ], capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 2)
        self.assertIn("CODEX_CONFIG_NOT_C_SIGNED", result.stderr)
        bad_sha = subprocess.run([
            sys.executable, str(script), "--config", str(CONFIG_PATH),
            "--config-sha256", "0" * 64, "--owner", "D",
            "--task-id", "C-D-AGL-1250K-rho097-seed01", "--stage", "warmup",
        ], capture_output=True, text=True, timeout=20)
        self.assertEqual(bad_sha.returncode, 2)
        self.assertIn("C_FROZEN_CONFIG_SHA_MISMATCH", bad_sha.stderr)
        self.assertEqual(sampler_root.exists(), before)

    def test_actual_runner_refuses_execution_code_sha_mismatch(self):
        config = json.loads(CONFIG_PATH.read_text())
        config["execution_code_sha256"]["inspect_eam.py"] = "0" * 64
        with tempfile.TemporaryDirectory(prefix="wrong-code-pin-", dir="/workspace") as tmp_value:
            config_path = Path(tmp_value) / "config.json"
            config_path.write_text(json.dumps(config))
            result = subprocess.run([
                sys.executable, str(BASE / "run_sampling_codex.py"),
                "--config", str(config_path),
                "--config-sha256", hashlib.sha256(config_path.read_bytes()).hexdigest(),
                "--owner", "D", "--task-id", "C-D-AGL-1250K-rho097-seed01", "--stage", "warmup",
            ], capture_output=True, text=True, timeout=20)
            self.assertEqual(result.returncode, 2)
            self.assertIn("EXECUTION_CODE_SHA_MISMATCH", result.stderr)

    def test_runner_wrong_owner_and_wrong_task_refuse_early(self):
        script = BASE / "run_sampling_codex.py"
        config_sha = hashlib.sha256(CONFIG_PATH.read_bytes()).hexdigest()
        result = subprocess.run([
            sys.executable, str(script), "--config", str(CONFIG_PATH),
            "--config-sha256", config_sha, "--owner", "X",
            "--task-id", "C-D-AGL-1250K-rho097-seed01", "--stage", "warmup",
        ], capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 2)
        self.assertIn("invalid choice", result.stderr)
        with tempfile.TemporaryDirectory(prefix="wrong-task-config-", dir="/workspace") as tmp_value:
            tmp = Path(tmp_value)
            cfg = json.loads(CONFIG_PATH.read_text())
            cfg.update({"config_status": "C_SIGNED_FOR_WARMUP",
                        "release_decision": "APPROVE_WARMUP_ONLY",
                        "source_release": True, "density_reference_pass": True,
                        "storage_release": True, "warmup_release": True})
            config_path = tmp / "config.json"
            config_path.write_text(json.dumps(cfg))
            digest = hashlib.sha256(config_path.read_bytes()).hexdigest()
            result = subprocess.run([
                sys.executable, str(script), "--config", str(config_path),
                "--config-sha256", digest, "--owner", "E",
                "--task-id", "C-D-AGL-1250K-rho097-seed01", "--stage", "warmup",
            ], capture_output=True, text=True, timeout=20)
            self.assertEqual(result.returncode, 2)
            self.assertIn("TASK_NOT_ALLOWED_FOR_OWNER", result.stderr)

    def test_new_config_does_not_change_frozen_science_or_task_matrix(self):
        old = json.loads((BASE / "PERSISTENT_REFERENCE_CONFIG.json").read_text())
        new = json.loads(SIGNED_CONFIG_PATH.read_text())
        for key, value in old["sampling"].items():
            if key == "max_wall_definition":
                continue  # explanatory resource-policy text, not a numerical sampler parameter
            self.assertEqual(new["sampling"][key], value, key)
        self.assertIn("warmup <=1 hour", new["sampling"]["max_wall_definition"])
        self.assertIn("Continuation is not authorized", new["sampling"]["max_wall_definition"])
        self.assertEqual(new["potential"]["sha256"], old["potential"]["sha256"])
        self.assertEqual(new["potential"]["cutoff_A"], old["potential"]["cutoff_A"])
        self.assertEqual((BASE / "AUTHORIZED_SIX_TASKS.csv").read_bytes(),
                         subprocess.check_output([
                             "git", "show",
                             "2959fd3e2396a7319bd1f5c7f4fec27ce442abe:research/mace-v12-transfer/coordination/reports/window-c/liquid_ag_six_pilot_authorized_dispatch_20261011/AUTHORIZED_SIX_TASKS.csv",
                         ], cwd=BASE))

    def test_signed_config_releases_only_first_D_warmup(self):
        config = json.loads(SIGNED_CONFIG_PATH.read_text())
        self.assertEqual(config["config_status"], "C_SIGNED_FOR_WARMUP")
        self.assertEqual(config["release_decision"], "APPROVE_WARMUP_ONLY")
        self.assertTrue(config["storage_release"] and config["warmup_release"])
        self.assertFalse(config["continuation_release"])
        self.assertEqual(config["warmup_task_ids"], ["C-D-AGL-1250K-rho097-seed01"])
        self.assertTrue(config["hosts"]["D"]["C_host_release"])
        self.assertFalse(config["hosts"]["E"]["C_host_release"])
        self.assertFalse(config["hosts"]["F"]["C_host_release"])
        self.assertEqual(config["sampling"]["warmup_steps"], 1000)
        self.assertEqual(config["sampling"]["dt_fs"], 0.5)
        self.assertFalse(config["hosts"]["D"]["persistence_verified"])
        self.assertEqual(config["potential"]["radial_tail_policy"],
                         "C_APPROVED_PUBLISHED_SETFL_ONE_BIN_ZERO_TAIL")
        self.assertEqual(config["potential"]["radial_tail_review_evidence"]["potential_sha256"],
                         EXPECTED_POTENTIAL)
        self.assertEqual(config["potential"]["radial_tail_review_evidence"]["decision"],
                         "APPROVED_FOR_THIS_EXACT_FILE_AND_ASE_3_29_INTERPOLATION_ONLY")
        self.assertEqual(config["hosts"]["D"]["identity_evidence_commit"],
                         "2c8b6f41ec8a7605ddc06f7218014c55957166ac")

    def test_runner_legacy_persistent_path_still_requires_true(self):
        old = json.loads((BASE / "PERSISTENT_REFERENCE_CONFIG.json").read_text())
        old["source_release"] = True
        old["storage_release"] = True
        old["warmup_release"] = True
        old["density_reference_pass"] = True
        old.pop("storage_mode", None)
        old["execution_code_sha256"] = json.loads(CONFIG_PATH.read_text())["execution_code_sha256"]
        with tempfile.TemporaryDirectory(prefix="legacy-persist-check-", dir="/workspace") as tmp_value:
            tmp = Path(tmp_value)
            config_path = tmp / "legacy.json"
            config_path.write_text(json.dumps(old))
            (tmp / "AUTHORIZED_SIX_TASKS.csv").write_bytes((BASE / "AUTHORIZED_SIX_TASKS.csv").read_bytes())
            config_sha = hashlib.sha256(config_path.read_bytes()).hexdigest()
            result = subprocess.run([
                sys.executable, str(BASE / "run_sampling_codex.py"), "--config", str(config_path),
                "--config-sha256", config_sha, "--owner", "D",
                "--task-id", "C-D-AGL-1250K-rho097-seed01", "--stage", "warmup",
            ], capture_output=True, text=True, timeout=20)
            self.assertEqual(result.returncode, 2)
            self.assertIn("PERSISTENT_WORK_ROOT_UNVERIFIED", result.stderr)

    def test_unsigned_release_rejected(self):
        cfg = copy.deepcopy(self.config)
        cfg["config_status"] = "C_HOLD_NOT_SIGNED_FOR_SCIENTIFIC_EXECUTION"
        self.assertRefused("CODEX_CONFIG_NOT_C_SIGNED", lambda: invoke(cfg, potential=Path("/missing")))

    def test_config_sha_mismatch_rejected(self):
        self.assertRefused("CODEX_CONFIG_SHA_MISMATCH",
                           lambda: gate.load_signed_config(CONFIG_PATH, "0" * 64))

    def test_wrong_owner_rejected(self):
        self.assertRefused("OWNER_NOT_REGISTERED", lambda: invoke(self.config, owner="X"))

    def test_task_from_other_owner_rejected(self):
        self.assertRefused("TASK_NOT_ALLOWED_FOR_OWNER",
                           lambda: invoke(self.config, owner="E", task="C-D-AGL-1250K-rho097-seed01"))

    def test_missing_remote_archive_rejected(self):
        cfg = copy.deepcopy(self.config)
        cfg["hosts"]["D"].pop("remote_archive")
        self.assertRefused("REMOTE_ARCHIVE_CHANNEL_MISSING", lambda: invoke(cfg))

    def test_e_report_file_roundtrip_does_not_release_synthetic_gate(self):
        cfg = copy.deepcopy(self.config)
        cfg["warmup_task_ids"] = ["C-E-AGL-1250K-rho097-seed02"]
        self.assertRefused("REMOTE_SYNTHETIC_ROUNDTRIP_NOT_VERIFIED",
                           lambda: invoke(cfg, owner="E", task="C-E-AGL-1250K-rho097-seed02",
                                          head=cfg["hosts"]["E"]["remote_archive"]["expected_remote_head"]))

    def test_wrong_potential_sha_rejected_without_substitution(self):
        candidate = self.root / "Ag.eam.alloy"
        candidate.write_bytes(b"synthetic test data only")
        self.assertRefused("POTENTIAL_SHA_MISMATCH_NO_ALTERNATIVE",
                           lambda: invoke(self.config, potential=candidate))

    def test_exact_potential_without_node_receipt_rejected(self):
        candidate = self.root / "Ag.eam.alloy"
        candidate.write_bytes(b"not a potential; synthetic gate fixture")
        cfg = copy.deepcopy(self.config)
        cfg["hosts"]["D"]["potential_local_sha_verified"] = False
        with patch.object(gate, "sha256_file", return_value=EXPECTED_POTENTIAL):
            self.assertRefused("NODE_LOCAL_POTENTIAL_RECEIPT_MISSING",
                               lambda: invoke(cfg, potential=candidate))

    def test_persistent_mode_still_requires_true(self):
        self.assertRefused("PERSISTENT_WORK_ROOT_UNVERIFIED",
                           lambda: gate.storage_mode_precheck(gate.MODE_PERSISTENT,
                                                             {"persistence_verified": False}))
        self.assertEqual(gate.MODE_PERSISTENT,
                         gate.storage_mode_precheck(gate.MODE_PERSISTENT,
                                                   {"persistence_verified": True}))

    def test_ephemeral_mode_keeps_persistence_false(self):
        self.assertEqual(gate.MODE_EPHEMERAL,
                         gate.storage_mode_precheck(gate.MODE_EPHEMERAL,
                                                   {"persistence_verified": False}))
        self.assertRefused("EPHEMERAL_MODE_MUST_KEEP_PERSISTENCE_FALSE",
                           lambda: gate.storage_mode_precheck(gate.MODE_EPHEMERAL,
                                                              {"persistence_verified": True}))

    def test_failed_preflight_creates_no_sampler_root(self):
        root = Path(self.config["hosts"]["D"]["sampler_work_root"])
        self.assertFalse(root.exists())
        self.assertRefused("POTENTIAL_NOT_AVAILABLE_NO_ALTERNATIVE",
                           lambda: invoke(self.config, potential=Path("/nonexistent/Ag.eam.alloy")))
        self.assertFalse(root.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
