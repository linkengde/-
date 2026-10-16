"""No GPAW import and no scientific calculation."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import time
import unittest
from types import SimpleNamespace

P = Path(__file__).with_name('tsc48_segmented.py')
spec = importlib.util.spec_from_file_location('tsc48_segmented', P)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class OfflineTests(unittest.TestCase):
    def test_dual_slot_fallback_and_tamper(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            lock = {'scope':'TSC48_STATIC_SINGLEPOINT','input_sha256':'X',
                    'h_config_sha256':'Y','runner_sha256':'Z'}
            for sequence, slot in ((1,'B'),(2,'A')):
                data = f'fake checkpoint {slot}'.encode()
                (root/f'slot_{slot}.gpw').write_bytes(data)
                mod.write_atomic(root/f'slot_{slot}.json', {
                    **lock,'mode':'all','slot':slot,'sequence':sequence,
                    'checkpoint_finalized':True,
                    'cumulative_elapsed_seconds':12.0,
                    'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
            found = mod.verified_slots(root,lock)
            self.assertEqual([x['sequence'] for x,_ in found],[2,1])
            (root/'slot_A.gpw').write_bytes(b'bad')
            found = mod.verified_slots(root,lock)
            self.assertEqual([x['sequence'] for x,_ in found],[1])
            (root/'slot_B.writing.gpw').write_bytes(b'incomplete')
            self.assertEqual(len(mod.verified_slots(root,lock)),1)
            # A byte-valid, not-yet-finalized slot is not recoverable.
            meta = mod.read_json(root/'slot_B.json')
            meta['checkpoint_finalized'] = False
            mod.write_atomic(root/'slot_B.json',meta)
            self.assertEqual(mod.verified_slots(root,lock),[])

    def test_mpi_shared_decision_root_only(self):
        import numpy as np
        class Bus:
            packet = None
        class FakeComm:
            def __init__(self, rank, bus):
                self.rank, self.bus = rank, bus
            def broadcast(self, buf, root):
                self_ref = self
                if self_ref.rank == root:
                    self_ref.bus.packet = buf.copy()
                else:
                    buf[:] = self_ref.bus.packet
        bus = Bus()
        root = FakeComm(0, bus)
        follower = FakeComm(1, bus)
        # Root crossed wallclock threshold; rank 1 did not.
        a = mod.shared_decision(root, 3, 10, 20, 97, 90,
                                False, 300, 100)
        b = mod.shared_decision(follower, 3, 10, 20, 40, 90,
                                False, 999999, 100)
        self.assertEqual(a, (True, True, False))
        self.assertEqual(b, a)
        # A low-space event blocks all ranks even when non-root sees free disk.
        c = mod.shared_decision(root, 10, 10, 20, 10, 90,
                                False, 30, 100)
        d = mod.shared_decision(follower, 10, 10, 20, 10, 90,
                                False, 1000000, 100)
        self.assertEqual(c, (True, False, True))
        self.assertEqual(d, c)
        self.assertEqual(bus.packet.dtype, np.int32)

    def test_append_only_ledger_source_linkage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mod.append_ledger(root, {'event':'SEGMENT_STARTED',
                                     'resume_basis_checkpoint_sha256':'oldsha'})
            mod.append_ledger(root, {'event':'CHECKPOINT_COMMITTED',
                                     'source_checkpoint_sha256':'oldsha',
                                     'checkpoint_sha256':'newsha',
                                     'cumulative_elapsed_seconds':33.5,
                                     'cumulative_logged_scf_steps':20})
            events = [json.loads(x) for x in
                      (root/'segment_ledger.jsonl').read_text().splitlines()]
            self.assertEqual(len(events),2)
            self.assertEqual(events[1]['source_checkpoint_sha256'],
                             events[0]['resume_basis_checkpoint_sha256'])
            self.assertEqual(events[1]['cumulative_logged_scf_steps'],20)
            self.assertEqual(mod.latest_verified_counters({
                'cumulative_elapsed_seconds':33.5,
                'cumulative_logged_scf_steps':20}), (33.5,20))

    def test_gpaw_26_7_context_mapping(self):
        from ase.units import Ha
        ctx = SimpleNamespace(
            niter=7,
            ham=SimpleNamespace(e_total_extrapolated=2.0),
            dens=SimpleNamespace(error=0.004),
            wfs=SimpleNamespace(
                nvalence=8.0,
                eigensolver=SimpleNamespace(error=0.003)),
            eig_error=0.005)
        # Deliberately diverge the SCF counter and local hook-call number,
        # preventing the regression missed by the V3 tests.
        sample = mod.scf_context_sample(ctx, 3, 12.5)
        self.assertEqual(sample['local_iteration'], 3)
        self.assertEqual(sample['scf_iteration_from_context'], 7)
        self.assertEqual(sample['elapsed_seconds'], 12.5)
        self.assertAlmostEqual(sample['total_extrapolated_energy_eV'],
                               2.0 * Ha)
        self.assertAlmostEqual(sample['density_error_per_valence'], 0.004 / 8.0)
        self.assertAlmostEqual(sample['eigenstates_error_eV2_per_valence'],
                               0.003 * Ha**2 / 8.0)
        self.assertAlmostEqual(sample['eigenvalue_error_eV'], 0.005 * Ha)
        self.assertNotIn('e_error', sample)
        self.assertNotIn('density_error', sample)
        self.assertNotIn('eigenstates_error', sample)

    def test_checkpoint_cumulative_includes_write_and_publication(self):
        class FakeWorld:
            rank = 0
            def __init__(self):
                self.value = None
            def barrier(self):
                # Explicitly delay the MPI barrier so the test can detect
                # whether it is incorrectly counted as calc.write duration.
                time.sleep(0.035)
            def broadcast(self, buf, root):
                if root == 0:
                    self.value = buf.copy()
                else:
                    buf[:] = self.value

        class FakeCalc:
            def write(self, filename, mode=''):
                self.mode = mode
                time.sleep(0.02)
                Path(filename).write_bytes(b'fake-full-wavefunction-gpw')

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            world = FakeWorld()
            lock = {'scope':mod.SCOPE,'input_sha256':'I',
                    'h_config_sha256':'C','runner_sha256':'R'}
            started = time.monotonic()
            cumulative = mod.save_checkpoint(
                FakeCalc(), world, root, lock, 1, 20, 'SOURCE',
                1, started, 100.0, 20, 'gpaw.start.log')
            meta = mod.read_json(root/'slot_B.json')
            self.assertEqual(meta['mode'], 'all')
            self.assertEqual(meta['source_checkpoint_sha256'], 'SOURCE')
            self.assertTrue(meta['checkpoint_finalized'])
            self.assertGreater(meta['checkpoint_write_elapsed_seconds'], 0.015)
            self.assertLess(meta['checkpoint_write_elapsed_seconds'], 0.04)
            self.assertGreater(meta['cumulative_elapsed_seconds'], 100.07)
            self.assertGreater(meta['cumulative_elapsed_seconds'],
                               meta['checkpoint_write_elapsed_seconds']+100.05)
            self.assertGreaterEqual(
                cumulative, meta['cumulative_elapsed_seconds'] - 0.01)
            self.assertGreater(meta['cumulative_elapsed_seconds'], 100.015)
            self.assertTrue((root/'slot_B.gpw').is_file())
            self.assertEqual(len(mod.verified_slots(root,lock)),1)
            ledger=(root/'segment_ledger.jsonl').read_text()
            self.assertIn('CHECKPOINT_PUBLISHING',ledger)
            self.assertFalse((root/'slot_B.writing.gpw').exists())

    def test_resource_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir=Path(tmp)/'task'
            clearance=Path(tmp)/'clearance.json'
            clearance.write_text(json.dumps({'scope':mod.SCOPE,'resource_ready':False}))
            with self.assertRaises(mod.Blocked):
                mod.approved_resources(clearance,workdir,3600,600)

    def test_atomic_json_no_pending(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'meta.json'
            mod.write_atomic(p,{'n':1})
            self.assertEqual(mod.read_json(p)['n'],1)
            self.assertFalse((Path(tmp)/'meta.json.pending').exists())

if __name__=='__main__':
    unittest.main()
