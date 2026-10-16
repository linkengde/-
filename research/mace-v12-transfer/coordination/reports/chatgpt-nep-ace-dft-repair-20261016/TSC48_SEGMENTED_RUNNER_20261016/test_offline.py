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
