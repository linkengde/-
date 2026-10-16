import hashlib
import importlib
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[7]
MINIMAL=REPO/'research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016/TSC48_MINIMAL_RELAY_20261016'
sys.path.insert(0,str(MINIMAL))
rr=importlib.import_module('relay_runner')

class LocalRotationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.folder=Path(self.tmp.name)
        self.lock={'scope':'TSC48_STATIC_SINGLEPOINT','input_sha256':'input',
                   'h_config_sha256':'config','runner_sha256':'v4',
                   'relay_adapter_sha256':rr.PREVIOUS_ADAPTER_SHA,
                   'local_rotation_adapter_sha256':'local'}
        self.old_bytes=b'previous-complete-gpw-fixture'
        self.new_bytes=b'new-complete-gpw-fixture'
        self.old_path=self.folder/'checkpoint.00000002.gpw'
        self.new_path=self.folder/'checkpoint.00000003.gpw'
        self.old_path.write_bytes(self.old_bytes)
        self.new_path.write_bytes(self.new_bytes)
        self.old={'file':self.old_path.name,'sequence':2,'bytes':len(self.old_bytes),
                  'sha256':hashlib.sha256(self.old_bytes).hexdigest(),
                  'checkpoint_finalized':True,'mode':'all',**self.lock}
        self.new={'file':self.new_path.name,'sequence':3,'bytes':len(self.new_bytes),
                  'sha256':hashlib.sha256(self.new_bytes).hexdigest(),
                  'checkpoint_finalized':True,'mode':'all',
                  'previous_file':self.old_path.name,'previous_sha256':self.old['sha256'],
                  'cumulative_logged_scf_steps':3,**self.lock}
        (self.folder/'current.json').write_text(rr.base.json.dumps(self.new))
        (self.folder/'lock.json').write_text(rr.base.json.dumps(self.lock))
        (self.folder/'ACTIVE_RUN.json').write_text(rr.base.json.dumps(
            {'source_checkpoint_sha256':self.old['sha256'],
             'latest_checkpoint_sequence':2}))
        (self.folder/'checkpoint.writing.gpw').unlink(missing_ok=True)
        self.env=patch.dict(os.environ,{'TSC48_LOCAL_ONLY_ROTATION':'YES'})
        self.env.start()

    def tearDown(self):
        self.env.stop()
        self.tmp.cleanup()

    def test_bounded_hash_matches_standard_sha256(self):
        p=self.folder/'sample.bin'
        p.write_bytes(bytes(range(256))*8193)
        self.assertEqual(rr.bounded_sha256(p),hashlib.sha256(p.read_bytes()).hexdigest())

    def test_retire_after_new_current_and_archive(self):
        records=[]
        with patch.object(rr,'archive_record',side_effect=lambda r: records.append(r) or {'commit':'fixture'}):
            rr.retire_previous_local(self.folder,self.lock,self.old,self.new)
        self.assertFalse(self.old_path.exists())
        self.assertTrue(self.new_path.exists())
        self.assertEqual(len(records),2)
        current=rr.base.read_json(self.folder/'current.json')
        self.assertIsNone(current['previous_file'])
        self.assertEqual(current['sha256'],self.new['sha256'])
        self.assertEqual(records[1]['retired_checkpoint_sha256'],self.old['sha256'])

    def test_open_predecessor_is_never_deleted(self):
        records=[]
        with self.old_path.open('rb') as held:
            with patch.object(rr,'archive_record',side_effect=lambda r: records.append(r) or {'commit':'fixture'}):
                with self.assertRaises(rr.base.Blocked):
                    rr.retire_previous_local(self.folder,self.lock,self.old,self.new)
        self.assertTrue(self.old_path.exists())
        self.assertTrue(self.new_path.exists())

    def test_wrong_active_pointer_preserves_predecessor(self):
        (self.folder/'ACTIVE_RUN.json').write_text(rr.base.json.dumps(
            {'latest_checkpoint_sha256':'different'}))
        with self.assertRaises(rr.base.Blocked):
            rr.retire_previous_local(self.folder,self.lock,self.old,self.new)
        self.assertTrue(self.old_path.exists())

if __name__=='__main__':
    unittest.main()
