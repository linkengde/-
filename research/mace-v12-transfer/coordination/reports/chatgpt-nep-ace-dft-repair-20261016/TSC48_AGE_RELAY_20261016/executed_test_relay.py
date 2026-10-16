import tempfile,unittest,json
from pathlib import Path
from unittest.mock import patch
import relay_runner as r

class World:
    rank=0
    def barrier(self): pass
    def broadcast(self,a,root): pass
class Calc:
    def __init__(self,fail=False): self.fail=fail
    def write(self,p,mode):
        assert mode=='all'
        Path(p).write_bytes(b'wavefunction-test'*100)
        if self.fail: raise OSError('injected disk failure')
class RelayTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.p=Path(self.temp.name);self.lock={'scope':'test'}
    def tearDown(self): self.temp.cleanup()
    def save(self,n=1,fail=False):
        r.mpi_save(Calc(fail),World(),self.p,self.lock,n,3,None,n,0,0,3,'test.log')
        return r.read_current(self.p,self.lock)[0]
    def test_mpi_writer_publishes(self):
        m,p=self.save();self.assertTrue(m['checkpoint_finalized']);self.assertEqual(m['bytes'],p.stat().st_size)
        self.assertGreater(m['cumulative_elapsed_seconds'],m['checkpoint_write_elapsed_seconds'])
    def test_failed_write_preserves_current(self):
        m,p=self.save()
        with self.assertRaises(OSError):self.save(2,True)
        self.assertEqual(r.read_current(self.p,self.lock)[0][0]['sha256'],m['sha256']);self.assertTrue(p.exists())
    def test_pointer_failure_preserves_current(self):
        m,p=self.save();original=r.base.write_atomic
        def fail(path,data):
            if Path(path).name=='current.json':raise OSError('pointer publication failure')
            return original(path,data)
        with patch.object(r.base,'write_atomic',fail):
            with self.assertRaises(r.base.Blocked): self.save(2)
        self.assertEqual(r.read_current(self.p,self.lock)[0][0]['sequence'],1)
        with self.assertRaises(r.base.Blocked):r.ensure_two_file_budget(self.p,self.lock)
    def test_two_files_then_block_third(self):
        self.save();self.save(2)
        self.assertEqual(len(list(self.p.glob('*.gpw'))),2)
        with self.assertRaises(r.base.Blocked):self.save(3)
    def test_missing_receipt_no_retirement(self):
        self.save();self.save(2);receipt=self.p/'receipt.json';receipt.write_text('{}')
        with self.assertRaises(r.base.Blocked):r.retire_previous(self.p,self.lock,receipt)
        self.assertEqual(len(list(self.p.glob('*.gpw'))),2)
    def test_verified_attestation_only_removes_previous(self):
        self.save();m,p=self.save(2);receipt=self.p/'receipt.json'
        receipt.write_text(json.dumps({'independent_reader_verified':True,'sha256':m['sha256'],
            'bytes':m['bytes'],'reader_identity':'synthetic-test-reader','remote_immutable_id':'synthetic-only'}))
        r.retire_previous(self.p,self.lock,receipt);self.assertTrue(p.exists());self.assertEqual(len(list(self.p.glob('*.gpw'))),1)
    def test_corrupt_current_rejected(self):
        m,p=self.save();p.write_bytes(b'bad')
        with self.assertRaises(r.base.Blocked):r.read_current(self.p,self.lock)
    def test_unfinalized_rejected(self):
        self.save();q=self.p/'current.json';m=json.loads(q.read_text());m['checkpoint_finalized']=False;q.write_text(json.dumps(m))
        with self.assertRaises(r.base.Blocked):r.read_current(self.p,self.lock)
if __name__=='__main__':unittest.main()

class RuntimePolicyTests(unittest.TestCase):
    def test_ppcg_named_dict_preserves_parameters(self):
        from gpaw.dft import PPCG
        from types import SimpleNamespace
        got={}
        def factory(**kw):got.update(kw);return SimpleNamespace(hooks={})
        solver=PPCG(niter=2)
        r.compatible_factory(factory,eigensolver=solver,maxiter=160)
        self.assertEqual(got['eigensolver'],{'name':'ppcg',**solver.todict()})
        self.assertEqual(got['maxiter'],160)
    def test_actual_new_factory_has_hook(self):
        from gpaw import GPAW,PW
        from gpaw.dft import PPCG
        c=r.compatible_factory(GPAW,mode=PW(500),eigensolver=PPCG(niter=2),txt=None)
        self.assertTrue(hasattr(c,'hooks'));self.assertEqual(c.params.eigensolver.niter,2)
    def test_regular_save_then_pause_with_two_file_peak(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            # every=5, limit huge, elapsed=0, wallclock far away, plenty disk.
            args=(5,5,2**30,0,10000,False,1000,1)
            self.assertEqual(r.periodic_decision(folder,World(),*args),(True,False,False))
            (folder/'checkpoint.00000001.gpw').write_bytes(b'synthetic')
            self.assertEqual(r.periodic_decision(folder,World(),*args),(True,True,False))
            args=(4,5,2**30,0,10000,False,1000,1)
            self.assertEqual(r.periodic_decision(folder,World(),*args),(False,False,False))
