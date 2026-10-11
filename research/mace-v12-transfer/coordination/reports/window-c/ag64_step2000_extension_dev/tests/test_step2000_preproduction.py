"""Offline pre-production tests for the step-2000 continuation.

These tests deliberately stub ASE because the current validation runtime does
not have ASE 3.29.0 installed. They never attach a calculator or evaluate real
forces or energies. Timeout recovery uses deterministic array-only virtual
dynamics. Set AG64_STEP2000_ARCHIVE to the verified step-2000 archive to run
the source-checkpoint cases.
"""
from __future__ import annotations

import copy
import json
import os
import signal
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

DEV_ROOT = Path(__file__).resolve().parents[1]
CORE_ROOT = DEV_ROOT / "core"
sys.path.insert(0, str(DEV_ROOT))
sys.path.insert(0, str(CORE_ROOT))


class _MockCell:
    def __init__(self, value):
        self.array = np.asarray(value).copy()


class _MockAtoms:
    """Array-only ASE stand-in for checkpoint restoration tests."""

    def __init__(self, *, numbers, positions, cell, pbc):
        self.numbers = np.asarray(numbers).copy()
        self.positions = np.asarray(positions).copy()
        self.cell = _MockCell(cell)
        self.pbc = np.asarray(pbc).copy()
        self._masses = None
        self._momenta = None
        self._arrays = {}
        self.calc = None

    def set_masses(self, values):
        self._masses = np.asarray(values).copy()

    def get_masses(self):
        return self._masses.copy()

    def set_momenta(self, values, apply_constraint=True):
        self._momenta = np.asarray(values).copy()

    def get_momenta(self):
        return self._momenta.copy()

    def get_velocities(self):
        return self._momenta / self._masses[:, None]

    def set_positions(self, values, apply_constraint=True):
        self.positions = np.asarray(values).copy()

    def get_positions(self):
        return self.positions.copy()

    def new_array(self, name, values):
        self._arrays[name] = np.asarray(values).copy()

    def get_array(self, name):
        return self._arrays[name].copy()

    def set_constraint(self, constraint):
        self.constraint = constraint

    def __len__(self):
        return len(self.numbers)

    def get_potential_energy(self):
        return 0.0

    def get_kinetic_energy(self):
        return float(np.sum(self._momenta ** 2) / 2.0)

    def get_stress(self, voigt=False):
        return np.zeros((3, 3), dtype=float)

    def get_volume(self):
        return float(np.linalg.det(self.cell.array))

    def get_temperature(self):
        return 1.0


class _MockLangevin:
    step_calls = 0

    def __init__(self, *args, **kwargs):
        self.nsteps = 0

    def step(self, forces=None):
        type(self).step_calls += 1
        raise AssertionError("mock integrator must not be advanced in pre-production tests")


def _install_ase_test_stubs():
    ase = types.ModuleType("ase")
    ase.__version__ = "MOCK_FOR_OFFLINE_TESTS"
    ase.Atoms = _MockAtoms
    ase.units = types.SimpleNamespace(fs=1.0)
    constraints = types.ModuleType("ase.constraints")
    constraints.FixCom = type("FixCom", (), {})
    md = types.ModuleType("ase.md")
    langevin = types.ModuleType("ase.md.langevin")
    langevin.Langevin = _MockLangevin
    ase.md = md
    md.langevin = langevin
    sys.modules.update({
        "ase": ase,
        "ase.constraints": constraints,
        "ase.md": md,
        "ase.md.langevin": langevin,
    })


_install_ase_test_stubs()
import segment_core as core  # noqa: E402
import continuation_adapter as adapter  # noqa: E402

ARCHIVE_PATH = Path(os.environ.get("AG64_STEP2000_ARCHIVE", "/nonexistent"))
EXPECTED_EXECUTOR_SHA = "3033c3234c4074358f8ef54a6783c3c4c9495549c12c93641f4f777f3451ed76"
EXPECTED_STATE_SHA = "5d03864f07a98f177c3c0cbdcb850158ed89a64fa2a6ab7f879e8c7ab003be2a"
EXPECTED_ARRAYS_SHA = "0f7830ebeba5f9a83da77b15e153546b24fb3efaf0299de767d9b4fbe4639f78"
EXPECTED_ARCHIVE_SHA = "c3a88a759e0bd258e74632aba0a68e4f1a765fcbe8324034d5bf127b5aa69965"


@unittest.skipUnless(ARCHIVE_PATH.is_file(), "set AG64_STEP2000_ARCHIVE to the verified archive")
class Step2000CheckpointTests(unittest.TestCase):
    def _restore(self):
        with tempfile.TemporaryDirectory(prefix="ag64-step2000-copy-") as tmp_name:
            archive_copy = Path(tmp_name) / "step2000.tar.gz"
            archive_copy.write_bytes(ARCHIVE_PATH.read_bytes())
            self.assertEqual(adapter.sha(archive_copy), EXPECTED_ARCHIVE_SHA)
            with tempfile.TemporaryDirectory(prefix="ag64-step2000-restore-") as restore_name:
                bundle, result = adapter._extract_bundle_verified(
                    archive_copy,
                    Path(restore_name),
                    {
                        "source_archive_sha256": EXPECTED_ARCHIVE_SHA,
                        "source_checkpoint_state_sha256": EXPECTED_STATE_SHA,
                        "source_checkpoint_arrays_sha256": EXPECTED_ARRAYS_SHA,
                        "source_executor_code_sha256": EXPECTED_EXECUTOR_SHA,
                        "start_step": 2000,
                    },
                )
                # Copy the result out of the temporary extraction directory
                # while preserving only immutable test input for the caller.
                self._bundle = bundle
                self._restored = result["restored"]
                self._manifest = result["manifest"]
                self._source_files = {
                    name: (bundle / name).read_bytes()
                    for name in ("state.json", "rng_state.json")
                }

    def test_restored_global_step_is_2000(self):
        self._restore()
        self.assertEqual(self._restored.global_step, 2000)
        self.assertEqual(self._manifest["global_step"], 2000)

    def test_rng_streams_restore_without_state_change(self):
        self._restore()
        expected = json.loads(self._source_files["rng_state.json"])
        self.assertEqual(self._restored.rng_state_all, expected)
        self.assertEqual(set(expected), {"packing", "velocity", "thermostat"})
        for stream_state in expected.values():
            rng = np.random.Generator(np.random.PCG64())
            rng.bit_generator.state = copy.deepcopy(stream_state)


class BoundaryAndHashTests(unittest.TestCase):
    def test_candidate_config_parses_and_2000_to_3000_boundary_passes(self):
        config_path = DEV_ROOT / "CANDIDATE_CONTINUATION_CONFIG_STEP2000.json"
        cfg = json.loads(config_path.read_text(encoding="utf-8"))
        self.assertEqual(
            core.validate_segment_bounds(cfg["start_step"], cfg["end_step"],
                                         cfg["max_added_steps"]),
            (2000, 3000, 1000),
        )
        self.assertIs(cfg["execution_authorized"], False)


class _VirtualSegmentDynamics:
    """Deterministic integer-step harness with no ASE integrator or calculator."""

    def __init__(self, atoms, rng, *, nsteps=2000, end_step=3000):
        self.atoms = atoms
        self.rng = rng
        self.nsteps = nsteps
        self.segment_start = 2000
        self.stop_at_step = end_step
        self.on_step = None
        self.visited_steps = []
        self.signal_step = None
        self.interrupt_step = None

    def step(self):
        if self.nsteps >= self.stop_at_step:
            raise AssertionError("virtual runner exceeded its configured boundary")
        noise = self.rng.normal(size=self.atoms.positions.shape)
        self.atoms.set_positions(self.atoms.positions + noise * 1.0e-5)
        self.atoms.set_momenta(self.atoms.get_momenta() + noise * 1.0e-3,
                               apply_constraint=False)
        self.nsteps += 1
        self.visited_steps.append(self.nsteps)
        if self.on_step is not None:
            self.on_step(self, "hot_melt_candidate", 1250.0, np.zeros((64, 3)))
        if self.signal_step == self.nsteps:
            handler = signal.getsignal(signal.SIGTERM)
            handler(signal.SIGTERM, None)
        if self.interrupt_step == self.nsteps:
            raise KeyboardInterrupt()


def _virtual_segment_state(seed=2000):
    numbers = np.full(64, 47, dtype=np.int64)
    positions = np.arange(192, dtype=float).reshape(64, 3) / 100.0
    atoms = _MockAtoms(numbers=numbers, positions=positions,
                       cell=np.eye(3) * core.EXPECTED_CELL_EDGE_A,
                       pbc=np.ones(3, dtype=bool))
    atoms.set_masses(np.full(64, 107.8682))
    atoms.set_momenta(np.full((64, 3), 0.25), apply_constraint=False)
    atoms.new_array("atom_ids", np.arange(1, 65, dtype=np.int64))
    rng = np.random.Generator(np.random.PCG64(seed))
    all_rng = {}
    for name, stream_seed in (("packing", seed + 1), ("velocity", seed + 2),
                              ("thermostat", seed)):
        stream = np.random.Generator(np.random.PCG64(stream_seed))
        all_rng[name] = copy.deepcopy(stream.bit_generator.state)
    return atoms, rng, all_rng


class TimeoutRecoveryTests(unittest.TestCase):
    def _run(self, dyn, rng_states, output_root, *, clock=None,
             stop_requested=None, thermal_rows=None, trajectory_frames=None):
        kwargs = {}
        if clock is not None:
            kwargs["clock"] = clock
        if stop_requested is not None:
            kwargs["stop_requested"] = stop_requested
        return core.advance_segment_with_outputs(
            dyn, rng_state_all=rng_states,
            source_provenance={"synthetic": True, "execution_authorized": False},
            output_root=output_root,
            remote_roundtrip=lambda archive, digest: self.assertEqual(
                core.sha256_file(archive), digest),
            existing_thermal_rows=thermal_rows,
            existing_trajectory_frames=trajectory_frames,
            max_wall_s=3600.0, **kwargs)

    def test_walltime_checkpoint_restore_resume_matches_continuous_virtual_run(self):
        atoms, rng, rng_states = _virtual_segment_state()
        timed = _VirtualSegmentDynamics(atoms, rng)
        with tempfile.TemporaryDirectory(prefix="ag64-timeout-split-") as temp:
            root = Path(temp)
            with self.assertRaises(core.SegmentCheckpointSaved) as caught:
                self._run(timed, rng_states, root / "split",
                          clock=lambda: 3601.0 if timed.nsteps >= 2500 else 0.0)
            saved = caught.exception.checkpoint
            self.assertEqual(caught.exception.reason, "SEGMENT_WALLTIME_BUDGET_EXCEEDED")
            self.assertEqual(saved["step"], 2500)
            checkpoint = Path(saved["bundle"])
            self.assertTrue(checkpoint.is_dir())
            self.assertTrue(checkpoint.name.startswith("checkpoint_step2500"))

            positions_2500 = timed.atoms.get_positions().copy()
            velocities_2500 = timed.atoms.get_velocities().copy()
            rng_2500 = copy.deepcopy(timed.rng.bit_generator.state)
            restored = core.restore_segment_bundle_with_executor(
                checkpoint, core.sha256_file(core.__file__), expected_global_step=2500)
            np.testing.assert_array_equal(restored.atoms.get_positions(), positions_2500)
            np.testing.assert_array_equal(restored.atoms.get_velocities(), velocities_2500)
            self.assertEqual(restored.rng_state_all["thermostat"], rng_2500)
            self.assertEqual(restored.rng_state_all, json.loads(
                (checkpoint / "rng_state.json").read_text(encoding="utf-8")))

            resumed = _VirtualSegmentDynamics(restored.atoms, restored.thermostat_rng,
                                              nsteps=restored.global_step)
            resumed_states = copy.deepcopy(restored.rng_state_all)
            outputs = self._run(resumed, resumed_states, root / "resumed",
                                thermal_rows=restored.thermal_rows,
                                trajectory_frames=restored.trajectory_frames)
            self.assertEqual(resumed.nsteps, 3000)
            self.assertEqual([item["step"] for item in outputs], [3000])
            all_steps = timed.visited_steps + resumed.visited_steps
            self.assertEqual(all_steps, list(range(2001, 3001)))
            self.assertEqual(len(all_steps), len(set(all_steps)))
            self.assertNotIn(3001, all_steps)

            control_atoms, control_rng, control_states = _virtual_segment_state()
            control = _VirtualSegmentDynamics(control_atoms, control_rng)
            self._run(control, control_states, root / "control")
            np.testing.assert_array_equal(resumed.atoms.get_positions(),
                                          control.atoms.get_positions())
            np.testing.assert_array_equal(resumed.atoms.get_velocities(),
                                          control.atoms.get_velocities())
            self.assertEqual(resumed.rng.bit_generator.state,
                             control.rng.bit_generator.state)
            self.assertEqual(resumed.visited_steps, list(range(2501, 3001)))

    def test_signal_interrupt_and_graceful_stop_each_save_current_checkpoint(self):
        scenarios = ("signal", "interrupt", "graceful")
        for scenario in scenarios:
            with self.subTest(scenario=scenario), tempfile.TemporaryDirectory(
                    prefix=f"ag64-stop-{scenario}-") as temp:
                atoms, rng, rng_states = _virtual_segment_state()
                dyn = _VirtualSegmentDynamics(atoms, rng)
                output = Path(temp) / "out"
                stop = None
                if scenario == "signal":
                    dyn.signal_step = 2250
                elif scenario == "interrupt":
                    dyn.interrupt_step = 2250
                else:
                    stop = lambda: dyn.nsteps >= 2250
                with self.assertRaises(core.SegmentCheckpointSaved) as caught:
                    self._run(dyn, rng_states, output, stop_requested=stop)
                checkpoint = caught.exception.checkpoint
                self.assertEqual(checkpoint["step"], 2250)
                self.assertTrue(Path(checkpoint["bundle"]).is_dir())
                if scenario == "signal":
                    self.assertEqual(caught.exception.reason,
                                     "SEGMENT_SIGNAL_STOP_REQUESTED_SIGTERM")
                elif scenario == "interrupt":
                    self.assertEqual(caught.exception.reason, "SEGMENT_INTERRUPTED")
                else:
                    self.assertEqual(caught.exception.reason,
                                     "SEGMENT_GRACEFUL_STOP_REQUESTED")
    def test_step_3001_is_refused_without_integrator_call(self):
        _MockLangevin.step_calls = 0
        dyn = object.__new__(core.BoundedLangevin)
        dyn.nsteps = 3000
        dyn.segment_start = 2000
        dyn.stop_at_step = 3000
        with self.assertRaisesRegex(core.SegmentRefused,
                                    "SIGNED_STOP_AT_STEP_REACHED_OR_STEP_OUT_OF_RANGE"):
            dyn.step()
        self.assertEqual(_MockLangevin.step_calls, 0)

    def test_wrong_potential_sha_is_rejected_before_source_access(self):
        with tempfile.TemporaryDirectory(prefix="ag64-wrong-potential-") as tmp_name:
            tmp = Path(tmp_name)
            config_path = tmp / "candidate.json"
            config_path.write_text("{}", encoding="utf-8")
            potential_path = tmp / "wrong.eam.alloy"
            potential_path.write_text("not the pinned potential", encoding="utf-8")
            output_root = tmp / "new-output"
            cfg = {"new_output_root": str(output_root)}
            with mock.patch.object(adapter, "validate_release_config", return_value=cfg):
                with self.assertRaisesRegex(adapter.ContinuationRefused,
                                            "POTENTIAL_FILE_MISSING_OR_SHA_MISMATCH"):
                    adapter.preflight(config_path, "0" * 64, potential_path, output_root)
            self.assertFalse(output_root.exists())

    @unittest.skipUnless(ARCHIVE_PATH.is_file(), "set AG64_STEP2000_ARCHIVE to the verified archive")
    def test_wrong_checkpoint_sha_is_rejected(self):
        with tempfile.TemporaryDirectory(prefix="ag64-bad-checkpoint-sha-") as tmp_name:
            root = Path(tmp_name)
            archive_copy = root / "step2000.tar.gz"
            archive_copy.write_bytes(ARCHIVE_PATH.read_bytes())
            with self.assertRaisesRegex(adapter.ContinuationRefused,
                                        "SOURCE_STATE_JSON_SHA_MISMATCH"):
                adapter._extract_bundle_verified(
                    archive_copy,
                    root,
                    {
                        "source_archive_sha256": EXPECTED_ARCHIVE_SHA,
                        "source_checkpoint_state_sha256": "0" * 64,
                        "source_checkpoint_arrays_sha256": EXPECTED_ARRAYS_SHA,
                        "source_executor_code_sha256": EXPECTED_EXECUTOR_SHA,
                        "start_step": 2000,
                    },
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
