"""No GPAW import and no scientific calculation."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

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
                    'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
            found = mod.verified_slots(root,lock)
            self.assertEqual([x['sequence'] for x,_ in found],[2,1])
            (root/'slot_A.gpw').write_bytes(b'bad')
            found = mod.verified_slots(root,lock)
            self.assertEqual([x['sequence'] for x,_ in found],[1])
            (root/'slot_B.writing.gpw').write_bytes(b'incomplete')
            self.assertEqual(len(mod.verified_slots(root,lock)),1)

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
