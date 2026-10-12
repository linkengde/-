#!/usr/bin/env python3
"""Synthetic-only checks for the R2 recovery runner and collector."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

import collector
import runner


class FakeWorld:
    rank = 0
    size = 1

    def barrier(self):
        return None

    def sum_scalar(self, value):
        return value

    def broadcast(self, value, root):
        return value


class FakeCalculator:
    def __init__(self, fail=False):
        self.fail = fail

    def write(self, filename, mode=""):
        if self.fail:
            raise OSError("synthetic checkpoint write failure")
        if mode != "all":
            raise AssertionError("wavefunction-inclusive mode was not requested")
        Path(filename).write_bytes(b"SYNTHETIC_GPW_WITH_WAVEFUNCTIONS")


def fake_reader(path):
    match = re.search(r"iter_(\d+)", Path(path).name)
    return {"wavefunctions": True, "iteration": int(match.group(1))}


def synthetic_atoms():
    from ase import Atoms
    return Atoms("H", positions=[[0.0, 0.0, 0.0]],
                 cell=[[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]],
                 pbc=[True, True, True])


class RecoveryPackageTests(unittest.TestCase):
    def setUp(self):
        self.config, self.config_sha = runner.load_config()
        self.lock, self.lock_path = runner.verify_parameter_lock(self.config)
        self.world = FakeWorld()

    def test_phase2_lock_is_the_collector_lock_and_recovery_deltas_are_bounded(self):
        self.assertEqual(runner.sha256_file(self.lock_path), runner.EXPECTED_LOCK_SHA)
        self.assertEqual(self.config["method"]["max_scf_iterations"], 200)
        self.assertEqual(self.config["method"]["baseline_max_scf_iterations"], 160)
        self.assertEqual(self.config["method"]["eigensolver"], {"name": "PPCG", "niter": 4})
        self.assertFalse(self.config["execution_authorized"])

    def test_ppcg_niter_four_is_passed_to_gpaw_without_calculation(self):
        import gpaw_data
        from gpaw.new.ase_interface import GPAW

        old_setup = os.environ.get("GPAW_SETUP_PATH")
        old_omp = os.environ.get("OMP_NUM_THREADS")
        os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
        os.environ["OMP_NUM_THREADS"] = "1"
        try:
            with tempfile.TemporaryDirectory() as tmp:
                kwargs = runner.build_calculator_kwargs(self.config, Path(tmp) / "synthetic.log")
                self.assertEqual(kwargs["eigensolver"], {"name": "ppcg", "niter": 4})
                calc = GPAW(**kwargs)
                self.assertEqual(calc.params.eigensolver.niter, 4)
                self.assertIsNone(calc._dft)
                del calc
        finally:
            if old_setup is None:
                os.environ.pop("GPAW_SETUP_PATH", None)
            else:
                os.environ["GPAW_SETUP_PATH"] = old_setup
            if old_omp is None:
                os.environ.pop("OMP_NUM_THREADS", None)
            else:
                os.environ["OMP_NUM_THREADS"] = old_omp

    def test_real_gpaw_gpw_write_and_read_roundtrip_uses_wavefunctions_no_scf(self):
        """Build a synthetic one-H initialization state; do not run a single SCF step."""
        import gpaw_data
        from gpaw import FermiDirac, PW
        from gpaw.dft import PPCG
        from gpaw.mixer import Mixer
        from gpaw.new.ase_interface import GPAW

        old_setup = os.environ.get("GPAW_SETUP_PATH")
        old_omp = os.environ.get("OMP_NUM_THREADS")
        os.environ["GPAW_SETUP_PATH"] = str(gpaw_data.datapath())
        os.environ["OMP_NUM_THREADS"] = "1"
        try:
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp) / "synthetic-r2-checkpoint-test"
                root.mkdir()
                calc = GPAW(
                    mode=PW(100), xc="PBE", kpts=(1, 1, 1),
                    occupations=FermiDirac(0.1),
                    convergence={"energy": 1e-5, "density": 1e-5, "eigenstates": 1e-5},
                    maxiter=200, mixer=Mixer(0.05, 8, 100),
                    eigensolver={"name": "ppcg", "niter": 4},
                    txt=str(root / "synthetic.log"),
                )
                atoms = synthetic_atoms()
                calc.initialize(atoms)
                self.assertEqual(calc.dft.scf_loop.niter, 0)
                manager = runner.CheckpointManager(
                    calc, atoms, root, self.world, self.config, self.config_sha,
                    runner.EXPECTED_LOCK_SHA,
                )
                manager.on_scf_step(type("Context", (), {"niter": 40})())
                checkpoint = manager.latest_checkpoint
                self.assertIsNotNone(checkpoint)
                sidecar = json.loads(checkpoint.with_suffix(".json").read_text())
                self.assertEqual(sidecar["completed_scf_iteration"], 40)
                self.assertEqual(sidecar["checkpoint_sha256"], runner.sha256_file(checkpoint))
                self.assertTrue(sidecar["wavefunctions_included"])
                self.assertTrue(sidecar["readback_verified"])
                self.assertEqual(calc.dft.scf_loop.niter, 0)
                self.assertFalse((root / "RESULTS_RAW.json").exists())
                del manager, calc, atoms
        finally:
            if old_setup is None:
                os.environ.pop("GPAW_SETUP_PATH", None)
            else:
                os.environ["GPAW_SETUP_PATH"] = old_setup
            if old_omp is None:
                os.environ.pop("OMP_NUM_THREADS", None)
            else:
                os.environ["OMP_NUM_THREADS"] = old_omp

    def test_checkpoint_write_failure_stops_without_retry_or_valid_receipt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "fail-case"
            root.mkdir()
            manager = runner.CheckpointManager(
                FakeCalculator(fail=True), synthetic_atoms(), root, self.world,
                self.config, self.config_sha, runner.EXPECTED_LOCK_SHA,
                gpw_reader=fake_reader,
            )
            with self.assertRaises(runner.CheckpointFailure):
                manager.save(40)
            self.assertEqual(manager.saved_steps, set())
            self.assertIsNone(manager.latest_checkpoint)
            self.assertFalse((root / "checkpoint" / "checkpoint_iter_000040.gpw").exists())

    def test_wall_clock_guard_saves_at_safe_iteration_then_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "wall-case"
            root.mkdir()
            manager = runner.CheckpointManager(
                FakeCalculator(), synthetic_atoms(), root, self.world,
                self.config, self.config_sha, runner.EXPECTED_LOCK_SHA,
                clock=lambda: 12601.0, started_at=0.0, gpw_reader=fake_reader,
            )
            with self.assertRaises(runner.SafeStop) as caught:
                manager.on_scf_step(type("Context", (), {"niter": 41})())
            self.assertEqual(caught.exception.reason, "WALL_CLOCK_SAFE_STOP")
            self.assertEqual(caught.exception.step, 41)
            self.assertEqual(manager.saved_steps, {41})
            self.assertTrue(manager.latest_checkpoint.is_file())

    def test_unconverged_or_exception_run_never_exports_partial_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "RUN_STATUS.json").write_text(json.dumps({
                "run_status": "FAILED_RUNTIME",
                "scf_converged": False,
                "scf_steps": 40,
                "failure_reason": "synthetic exception",
            }))
            (root / "RESULTS_RAW.json").write_text(json.dumps({
                "status": "FAILED_SCF", "scf_converged": False,
                "total_energy_eV": -100.0,
                "force_consistent_free_energy_eV": -100.1,
                "forces_eV_A": [[0.0, 0.0, 0.0]] * 20,
                "stress_eV_A3": [[0.0] * 3 for _ in range(3)],
            }))
            (root / "gpaw.log").write_text("SCF step 40; no convergence marker\n")
            result = collector.collect(root)
            self.assertEqual(result["status"], "FAILED_SCF")
            self.assertIsNone(result["total_energy_eV"])
            self.assertIsNone(result["force_consistent_free_energy_eV"])
            self.assertIsNone(result["forces_eV_A"])
            self.assertIsNone(result["stress_eV_A3"])

    def test_collector_accepts_only_complete_converged_synthetic_record(self):
        source = self.config["source_bindings"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            status = {"run_status": "SUCCESS", "scf_converged": True,
                      "scf_steps": 42, "authorization_sha256": "a" * 64,
                      "authorization_id": "C-TEST-ONLY"}
            raw = {
                "status": "SUCCESS", "registry_id": "R2",
                "role": "REGISTRY_SCREENING_ONLY_DIAGNOSTIC_COPY",
                "scf_converged": True, "scf_steps": 42,
                "total_energy_eV": -10.0,
                "force_consistent_free_energy_eV": -10.1,
                "forces_eV_A": [[0.0, 0.0, 0.0] for _ in range(20)],
                "stress_eV_A3": [[0.0, 0.0, 0.0] for _ in range(3)],
                "input_sha256": source["input_sha256"],
                "structure_sha256": source["structure_sha256"],
                "parameter_lock_sha256": source["parameter_lock_sha256"],
                "config_sha256": self.config_sha,
                "authorization_sha256": status["authorization_sha256"],
                "authorization_id": status["authorization_id"],
            }
            (root / "RUN_STATUS.json").write_text(json.dumps(status))
            (root / "RESULTS_RAW.json").write_text(json.dumps(raw))
            (root / "gpaw.log").write_text("Converged in 42 steps\n")
            result = collector.collect(root)
            self.assertEqual(result["status"], "SUCCESS")
            self.assertEqual(result["scf_steps"], 42)
            self.assertEqual(result["total_energy_eV"], -10.0)

    def test_run_gate_rejects_missing_c_authorization_before_creating_work_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            work_root = base / "must-not-be-created"
            auth_path = base / "missing-c-authorization.json"
            proc = subprocess.run(
                [sys.executable, str(runner.HERE / "runner.py"), "--run",
                 "--authorization", str(auth_path), "--authorization-sha256", "0" * 64,
                 "--work-root", str(work_root)],
                cwd=runner.REPO_ROOT,
                text=True, capture_output=True, timeout=30,
                env={**os.environ,
                     "GPAW_SETUP_PATH": os.environ.get("GPAW_SETUP_PATH", "/missing"),
                     "OMP_NUM_THREADS": "1"},
            )
            self.assertEqual(proc.returncode, 2)
            self.assertIn("C execution authorization file is missing", proc.stderr)
            self.assertFalse(work_root.exists())

    def test_launcher_has_outer_four_hour_deadline_and_orderly_stop_buffer(self):
        launcher = (runner.HERE / "launch_r2.sh").read_text(encoding="utf-8")
        self.assertIn("timeout --signal=TERM --kill-after=1s 14399s", launcher)
        self.assertIn('"$MPIEXEC" -n 4 "$GPAW" python -u', launcher)
        self.assertEqual(self.config["method"]["max_wall_clock_seconds"], 14400)
        self.assertEqual(self.config["method"]["wall_clock_exit_margin_seconds"], 1800)
        self.assertEqual(self.config["checkpoint_policy"]["interval_completed_scf_iterations"], 40)


if __name__ == "__main__":
    unittest.main(verbosity=2)
