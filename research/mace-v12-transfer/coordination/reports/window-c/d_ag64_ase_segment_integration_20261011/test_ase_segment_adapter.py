from __future__ import annotations

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixCom

import ase_segment_adapter as adapter

SYNTHETIC_PROVENANCE = {"fixture": "synthetic-virtual-calculator-only"}


class VirtualHarmonic(Calculator):
    """Synthetic-only deterministic calculator; no Ag model is loaded."""

    implemented_properties = ["energy", "forces", "stress"]
    codex_synthetic_test_calculator = True

    def __init__(self, reference, *, nan_on_call=None):
        super().__init__()
        self.reference = np.array(reference, dtype=float, copy=True)
        self.k = 0.025
        self.calls = 0
        self.nan_on_call = nan_on_call

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.calls += 1
        delta = atoms.positions - self.reference
        self.results["energy"] = 0.5 * self.k * float(np.sum(delta * delta))
        self.results["stress"] = np.zeros(6, dtype=float)
        forces = -self.k * delta
        if self.nan_on_call == self.calls:
            forces = np.full_like(forces, np.nan)
        self.results["forces"] = forces


def initial_state():
    rng = np.random.Generator(np.random.PCG64(20261011))
    edge = adapter.EXPECTED_CELL_EDGE_A
    positions = rng.uniform(0.3, edge - 0.3, (64, 3))
    masses = np.full(64, 107.8682)
    velocity_rng = np.random.Generator(np.random.PCG64(111))
    thermostat_rng = np.random.Generator(np.random.PCG64(222))
    packing_rng = np.random.Generator(np.random.PCG64(333))
    momenta = velocity_rng.normal(0.0, 0.02, (64, 3)) * masses[:, None]
    atoms = Atoms(numbers=np.full(64, 47), positions=positions,
                  cell=np.eye(3) * edge, pbc=True)
    atoms.set_masses(masses)
    atoms.set_momenta(momenta, apply_constraint=False)
    atoms.new_array("atom_ids", np.arange(1, 65, dtype=np.int64))
    atoms.set_constraint(FixCom())
    atoms.set_momenta(momenta, apply_constraint=False)
    atoms.calc = VirtualHarmonic(positions)
    streams = {"packing": packing_rng.bit_generator.state,
               "velocity": velocity_rng.bit_generator.state,
               "thermostat": thermostat_rng.bit_generator.state}
    return atoms, thermostat_rng, streams


def make_dyn(atoms, rng, start=1000, on_step=None):
    return adapter.BoundedLangevin(atoms, rng, initial_step=start, on_step=on_step,
                                   synthetic_test_mode=True)


def state_hash(atoms, step, rng, streams):
    payload = {
        "step": step,
        "positions": hashlib.sha256(np.ascontiguousarray(atoms.positions).tobytes()).hexdigest(),
        "momenta": hashlib.sha256(np.ascontiguousarray(atoms.get_momenta()).tobytes()).hexdigest(),
        "rng": copy.deepcopy(rng.bit_generator.state),
        "streams": copy.deepcopy(streams),
    }
    return payload


def recorded_step(dyn, streams, trace):
    streams["thermostat"] = copy.deepcopy(dyn.rng.bit_generator.state)
    trace[dyn.nsteps] = state_hash(dyn.atoms, dyn.nsteps, dyn.rng, streams)


def restore_synthetic(bundle, calculator_reference):
    restored = adapter.restore_segment_bundle(bundle)
    atoms = restored.atoms
    # The calculator is a test fixture and is deliberately reconstructed;
    # calculator state is not a physical restart variable.
    atoms.calc = VirtualHarmonic(calculator_reference)
    dyn = make_dyn(atoms, restored.thermostat_rng, start=restored.global_step)
    return dyn, restored.rng_state_all


class SegmentAdapterTests(unittest.TestCase):
    def test_ase_version_and_schedule_are_frozen_at_global_step(self):
        self.assertTrue(adapter.ase_version_is_329())
        self.assertEqual(adapter.schedule_target_K(999), ("hot_melt_candidate", 1600.0))
        self.assertEqual(adapter.schedule_target_K(1000), ("hot_melt_candidate", 1600.0))
        self.assertEqual(adapter.schedule_target_K(20000), ("cooling", 1600.0))

    def test_continuous_and_bundle_interrupt_resume_match_every_step(self):
        full_trace = {}
        atoms_a, rng_a, streams_a = initial_state()
        dyn_a = make_dyn(atoms_a, rng_a, on_step=lambda d, _p, _t, _f: recorded_step(d, streams_a, full_trace))
        self.assertEqual(adapter.advance_to_stop(dyn_a), 2000)
        self.assertEqual(len(full_trace), 1000)

        split_trace = {}
        atoms_b, rng_b, streams_b = initial_state()
        dyn_b = make_dyn(atoms_b, rng_b)
        dyn_b.on_step = lambda d, _p, _t, _f: recorded_step(d, streams_b, split_trace)
        while dyn_b.nsteps < 1400:
            dyn_b.step()
        snapshot = adapter.snapshot_state(dyn_b.atoms, dyn_b.nsteps, streams_b, dyn_b.rng,
                                          SYNTHETIC_PROVENANCE)
        with tempfile.TemporaryDirectory() as tmp:
            bundle = adapter.write_checkpoint_bundle(tmp, "synthetic-interruption", snapshot)
            # Bundle restoration uses actual SHA checks before resuming.
            dyn_c, streams_c = restore_synthetic(bundle, atoms_b.calc.reference.copy())
            self.assertEqual(dyn_c.nsteps, 1400)
            dyn_c.on_step = lambda d, _p, _t, _f: recorded_step(d, streams_c, split_trace)
            self.assertEqual(adapter.advance_to_stop(dyn_c), 2000)
        self.assertEqual(set(full_trace), set(split_trace))
        for step in full_trace:
            self.assertEqual(full_trace[step], split_trace[step], f"different RNG/state at step {step}")
        np.testing.assert_array_equal(dyn_a.atoms.positions, dyn_c.atoms.positions)
        np.testing.assert_array_equal(dyn_a.atoms.get_momenta(), dyn_c.atoms.get_momenta())
        self.assertEqual(rng_a.bit_generator.state, dyn_c.rng.bit_generator.state)

    def test_step_2001_is_refused_before_force_or_rng_consumption(self):
        atoms, rng, _streams = initial_state()
        dyn = make_dyn(atoms, rng)
        adapter.advance_to_stop(dyn)
        before_positions = atoms.positions.copy()
        before_momenta = atoms.get_momenta().copy()
        before_rng = copy.deepcopy(rng.bit_generator.state)
        before_calls = atoms.calc.calls
        with self.assertRaisesRegex(adapter.SegmentRefused, "STOP_AT_STEP_REACHED"):
            dyn.step()
        np.testing.assert_array_equal(atoms.positions, before_positions)
        np.testing.assert_array_equal(atoms.get_momenta(), before_momenta)
        self.assertEqual(rng.bit_generator.state, before_rng)
        self.assertEqual(atoms.calc.calls, before_calls)

    def test_run_and_irun_cannot_bypass_owned_step_gate(self):
        atoms, rng, _ = initial_state()
        dyn = make_dyn(atoms, rng)
        with self.assertRaisesRegex(adapter.SegmentRefused, "ASE_RUN_API_DISABLED"):
            dyn.run(120000)
        with self.assertRaisesRegex(adapter.SegmentRefused, "ASE_IRUN_API_DISABLED"):
            dyn.irun(120000)

    def test_non_test_integrator_requires_separately_hashed_signed_lock(self):
        atoms, rng, _ = initial_state()
        atoms.calc = None
        with self.assertRaisesRegex(adapter.SegmentRefused, "SIGNED_RELEASE_LOCK_REQUIRED"):
            adapter.BoundedLangevin(atoms, rng)

    def test_wrong_initial_or_stop_boundary_refused(self):
        atoms, rng, _ = initial_state()
        with self.assertRaisesRegex(adapter.SegmentRefused, "ONLY_SIGNED_1000_TO_2000"):
            adapter.BoundedLangevin(atoms, rng, stop_at_step=120000)
        with self.assertRaisesRegex(adapter.SegmentRefused, "SEGMENT_MUST_START_WITHIN"):
            adapter.BoundedLangevin(atoms, rng, initial_step=2001)

    def test_nonfinite_force_rolls_back_coordinates_momenta_and_rng(self):
        atoms, rng, _ = initial_state()
        atoms.calc = VirtualHarmonic(atoms.positions.copy(), nan_on_call=2)
        dyn = make_dyn(atoms, rng)
        before = (atoms.positions.copy(), atoms.get_momenta().copy(), copy.deepcopy(rng.bit_generator.state))
        with self.assertRaisesRegex(adapter.SegmentRefused, "INTEGRATOR_NONFINITE_STATE"):
            dyn.step()
        np.testing.assert_array_equal(atoms.positions, before[0])
        np.testing.assert_array_equal(atoms.get_momenta(), before[1])
        self.assertEqual(rng.bit_generator.state, before[2])
        self.assertEqual(dyn.nsteps, 1000)

    def test_timeout_stops_before_next_step(self):
        atoms, rng, _ = initial_state()
        dyn = make_dyn(atoms, rng)
        before = state_hash(atoms, dyn.nsteps, rng, {})
        ticks = iter([0.0, 1.0])
        with self.assertRaisesRegex(adapter.SegmentRefused, "WALLTIME_BUDGET"):
            adapter.advance_to_stop(dyn, max_wall_s=0.5, clock=lambda: next(ticks))
        after = state_hash(atoms, dyn.nsteps, rng, {})
        self.assertEqual(before, after)

    def test_checkpoint_manifest_and_no_overwrite(self):
        atoms, rng, streams = initial_state()
        dyn = make_dyn(atoms, rng)
        dyn.step()
        snapshot = adapter.snapshot_state(atoms, dyn.nsteps, streams, rng, SYNTHETIC_PROVENANCE)
        with tempfile.TemporaryDirectory() as tmp:
            bundle = adapter.write_checkpoint_bundle(tmp, "immutable", snapshot)
            restored = adapter.restore_segment_bundle(bundle)
            self.assertEqual(restored.global_step, 1001)
            with self.assertRaisesRegex(adapter.SegmentRefused, "ALREADY_EXISTS"):
                adapter.write_checkpoint_bundle(tmp, "immutable", snapshot)

    def test_remote_archive_failure_keeps_new_bundle_and_original_step1000(self):
        atoms, rng, streams = initial_state()
        dyn = make_dyn(atoms, rng)
        dyn.step()
        snapshot = adapter.snapshot_state(atoms, dyn.nsteps, streams, rng, SYNTHETIC_PROVENANCE)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = root / "immutable-step1000.bin"
            original.write_bytes(b"original-step1000-mock")
            original_sha = adapter.sha256_file(original)
            bundle = adapter.write_checkpoint_bundle(root / "new", "segment", snapshot)

            def remote_failure(_archive, _sha):
                raise OSError("injected fresh-fetch failure")

            with self.assertRaisesRegex(OSError, "injected fresh-fetch failure"):
                adapter.publish_and_verify(bundle, root / "archives" / "step1001.tar.gz", remote_failure)
            self.assertTrue(bundle.is_dir())
            self.assertTrue((root / "archives" / "step1001.tar.gz").is_file())
            self.assertEqual(adapter.sha256_file(original), original_sha)
            self.assertEqual(original.read_bytes(), b"original-step1000-mock")

    def test_local_archive_fresh_extract_and_internal_restore(self):
        atoms, rng, streams = initial_state()
        dyn = make_dyn(atoms, rng)
        dyn.step()
        snapshot = adapter.snapshot_state(atoms, dyn.nsteps, streams, rng, SYNTHETIC_PROVENANCE)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bundle = adapter.write_checkpoint_bundle(root / "origin", "fresh-readback", snapshot)
            archive_path = root / "archives" / "step1001.tar.gz"

            def local_fresh_extract(source_archive, expected_sha):
                self.assertEqual(adapter.sha256_file(source_archive), expected_sha)
                extract_root = root / "fresh-clone"
                extract_root.mkdir()
                import tarfile
                with tarfile.open(source_archive, "r:gz") as tf:
                    tf.extractall(extract_root, filter="data")
                recovered = adapter.restore_segment_bundle(extract_root / bundle.name)
                self.assertEqual(recovered.global_step, 1001)
                np.testing.assert_array_equal(recovered.atoms.positions, snapshot["positions"])
                np.testing.assert_array_equal(recovered.atoms.get_momenta(), snapshot["momenta"])

            sha = adapter.publish_and_verify(bundle, archive_path, local_fresh_extract)
            self.assertEqual(sha, adapter.sha256_file(archive_path))

    def test_output_workflow_stops_on_remote_failure_then_resumes_from_fresh_checkpoint(self):
        atoms_cont, rng_cont, streams_cont = initial_state()
        dyn_cont = make_dyn(atoms_cont, rng_cont)
        adapter.advance_to_stop(dyn_cont)

        atoms, rng, streams = initial_state()
        calculator_reference = atoms.calc.reference.copy()
        dyn = make_dyn(atoms, rng)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_checkpoint = root / "step1000-source.bin"
            source_checkpoint.write_bytes(b"immutable-source-step-1000")
            source_hash = adapter.sha256_file(source_checkpoint)

            def fail_once(_archive, _sha):
                raise OSError("synthetic remote archive interruption")

            with self.assertRaisesRegex(OSError, "synthetic remote archive interruption"):
                adapter.advance_segment_with_outputs(
                    dyn, rng_state_all=streams, source_provenance=SYNTHETIC_PROVENANCE,
                    output_root=root / "continuation", remote_roundtrip=fail_once)
            self.assertEqual(dyn.nsteps, 1100)
            saved = root / "continuation" / "segment-step-000001100"
            self.assertTrue(saved.is_dir())
            self.assertEqual(source_checkpoint.read_bytes(), b"immutable-source-step-1000")
            self.assertEqual(adapter.sha256_file(source_checkpoint), source_hash)

            restored = adapter.restore_segment_bundle(saved)
            restored.atoms.calc = VirtualHarmonic(calculator_reference)
            resumed = make_dyn(restored.atoms, restored.thermostat_rng, start=restored.global_step)

            def fresh_readback(archive_path, expected_sha):
                self.assertEqual(adapter.sha256_file(archive_path), expected_sha)
                fresh = root / "fresh" / f"{archive_path.stem}-{archive_path.stat().st_size}"
                fresh.mkdir(parents=True)
                import tarfile
                with tarfile.open(archive_path, "r:gz") as tf:
                    tf.extractall(fresh, filter="data")
                step_text = archive_path.name.removeprefix("step-").removesuffix(".tar.gz")
                recovered = adapter.restore_segment_bundle(fresh / f"segment-step-{int(step_text):09d}")
                self.assertGreaterEqual(recovered.global_step, 1200)

            outputs = adapter.advance_segment_with_outputs(
                resumed, rng_state_all=restored.rng_state_all,
                source_provenance=SYNTHETIC_PROVENANCE,
                output_root=root / "continuation", remote_roundtrip=fresh_readback,
                existing_thermal_rows=restored.thermal_rows,
                existing_trajectory_frames=restored.trajectory_frames)
            self.assertEqual(outputs[-1]["step"], 2000)
            self.assertEqual(resumed.nsteps, 2000)
            np.testing.assert_array_equal(resumed.atoms.positions, dyn_cont.atoms.positions)
            np.testing.assert_array_equal(resumed.atoms.get_momenta(), dyn_cont.atoms.get_momenta())
            self.assertEqual(resumed.rng.bit_generator.state, dyn_cont.rng.bit_generator.state)
            final_bundle = root / "continuation" / "segment-step-000002000"
            final_state = adapter.restore_segment_bundle(final_bundle)
            self.assertEqual(len(final_state.thermal_rows), 50)
            self.assertEqual(len(final_state.trajectory_frames), 11)
            self.assertEqual(final_state.thermal_rows[-1]["step"], 2000)
            self.assertEqual(final_state.trajectory_frames[-1]["step"], 2000)

    def test_release_code_hash_and_unsigned_lock_gates(self):
        with tempfile.TemporaryDirectory() as tmp:
            lock = Path(tmp) / "lock.json"
            document = {"release_status": "PREPARED_UNSIGNED_NOT_EXECUTABLE",
                        "execution_authorized": False,
                        "execution_code_sha256": adapter.sha256_file(Path(adapter.__file__)),
                        "task_id": adapter.TASK_ID, "owner": "D",
                        "start_step": 1000, "stop_at_step": 2000, "max_added_steps": 1000,
                        "source_scientific_parameters_sha256": adapter.EXPECTED_SCIENTIFIC_PARAMETERS_SHA256}
            lock.write_bytes(adapter.canonical_json_bytes(document))
            with self.assertRaisesRegex(adapter.SegmentRefused, "UNSIGNED_OR_NOT_AUTHORIZED"):
                adapter.validate_release_lock(lock, adapter.sha256_file(lock), document["execution_code_sha256"])
            with self.assertRaisesRegex(adapter.SegmentRefused, "LOCK_SHA_MISMATCH"):
                adapter.validate_release_lock(lock, "0" * 64, document["execution_code_sha256"])
            with self.assertRaisesRegex(adapter.SegmentRefused, "EXECUTION_CODE_SHA_MISMATCH"):
                adapter.validate_release_lock(lock, adapter.sha256_file(lock), "0" * 64)
            document["release_status"] = "C_SIGNED_FOR_SEGMENT_1000_2000"
            document["execution_authorized"] = True
            document["owner"] = "E"
            lock.write_bytes(adapter.canonical_json_bytes(document))
            with self.assertRaisesRegex(adapter.SegmentRefused, "TASK_OR_OWNER_MISMATCH"):
                adapter.validate_release_lock(lock, adapter.sha256_file(lock), document["execution_code_sha256"])
            document["owner"] = "D"
            document["stop_at_step"] = 120000
            lock.write_bytes(adapter.canonical_json_bytes(document))
            with self.assertRaisesRegex(adapter.SegmentRefused, "BOUNDARY_MISMATCH"):
                adapter.validate_release_lock(lock, adapter.sha256_file(lock), document["execution_code_sha256"])

    def test_remote_fresh_fetch_gate_is_mandatory_even_for_candidate_signed_shape(self):
        pins = {
            "source_archive_sha256": adapter.EXPECTED_ARCHIVE_SHA256,
            "source_archive_manifest_sha256": adapter.EXPECTED_ARCHIVE_MANIFEST_SHA256,
            "source_checkpoint_state_sha256": adapter.EXPECTED_STATE_SHA256,
            "source_state_arrays_sha256": adapter.EXPECTED_ARRAYS_SHA256,
            "source_manifest_sha256": adapter.EXPECTED_SOURCE_MANIFEST_SHA256,
            "source_executor_code_sha256": adapter.EXPECTED_OLD_RUNNER_SHA256,
            "potential_sha256": adapter.EXPECTED_POTENTIAL_SHA256,
        }
        with tempfile.TemporaryDirectory() as tmp:
            lock = Path(tmp) / "lock.json"
            document = {"release_status": "C_SIGNED_FOR_SEGMENT_1000_2000",
                        "execution_authorized": True,
                        "execution_code_sha256": adapter.sha256_file(Path(adapter.__file__)),
                        "task_id": adapter.TASK_ID, "owner": "D",
                        "start_step": 1000, "stop_at_step": 2000, "max_added_steps": 1000,
                        "source_scientific_parameters_sha256": adapter.EXPECTED_SCIENTIFIC_PARAMETERS_SHA256,
                        "sampling": adapter.SIGNED_SAMPLING_PARAMETERS,
                        "remote_archive_verified": False,
                        "remote_archive_method": "GITHUB_PUSH_FRESH_FETCH_SHA256",
                        **pins}
            lock.write_bytes(adapter.canonical_json_bytes(document))
            with self.assertRaisesRegex(adapter.SegmentRefused, "REMOTE_GITHUB_FRESH_FETCH"):
                adapter.validate_release_lock(lock, adapter.sha256_file(lock), document["execution_code_sha256"])
            document["remote_archive_verified"] = True
            lock.write_bytes(adapter.canonical_json_bytes(document))
            verified = adapter.validate_release_lock(lock, adapter.sha256_file(lock), document["execution_code_sha256"])
            self.assertEqual(verified["stop_at_step"], 2000)


if __name__ == "__main__":
    unittest.main(verbosity=2)
