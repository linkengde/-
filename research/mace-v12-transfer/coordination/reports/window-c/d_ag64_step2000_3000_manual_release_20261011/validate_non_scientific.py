"""Counter-only state/transport tests; no ASE integration or EAM calculation."""
import copy,json,sys,tempfile,platform,hashlib
from pathlib import Path
import numpy as np
import ase
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'core'))
import continuation_adapter as adapter
import segment_core as core
import github_archive_backend as backend
SOURCE_CFG=ROOT.parent/'d_ag64_step2000_3000_config_review_20261011/resource_budget_revision/UNSIGNED_STEP2000_3000_CONFIG_BUDGET3600.json'
def same(a,b):
    return (np.array_equal(a.atoms.positions,b.atoms.positions)
            and np.array_equal(a.atoms.get_momenta(),b.atoms.get_momenta())
            and a.rng_state_all==b.rng_state_all)
class VirtualAtoms:
    def __init__(self,atoms): self.a=atoms
    def __getattr__(self,n): return getattr(self.a,n)
    def get_potential_energy(self):return 0.0
    def get_kinetic_energy(self):return 0.0
    def get_temperature(self):return 1600.0
    def get_stress(self,voigt=False):return np.zeros((3,3))
class Counter:
    def __init__(self,restored,start=2000,end=3000):
        self.atoms=VirtualAtoms(restored.atoms.copy());self.rng=np.random.Generator(np.random.PCG64())
        self.rng.bit_generator.state=copy.deepcopy(restored.thermostat_rng.bit_generator.state)
        self.nsteps=self.segment_start=start;self.stop_at_step=end;self.on_step=None
    def step(self):
        if self.nsteps>=self.stop_at_step:raise core.SegmentRefused('VIRTUAL_STOP_BOUNDARY')
        self.rng.random();self.nsteps+=1
        if self.on_step:self.on_step(self,'hot_melt_candidate',1600.,None)
def main():
    cfg=json.loads(SOURCE_CFG.read_text());report={'scientific_MD_performed':False,'real_ASE_integrator_steps':0,'EAM_calls':0}
    with tempfile.TemporaryDirectory(prefix='c-manual-counter-review-') as td:
        tmp=Path(td)
        _,e=adapter._extract_bundle_verified(Path('/tmp/c-step2000-review.tar.gz'),tmp/'source',cfg)
        source=e['restored'];report['real_ASE_source_restore']='PASS'
        clock_values=iter([0.]+[0.]*351+[3601.])
        interrupted=Counter(source)
        stopped=tmp/'interrupted'
        local_calls=[]
        def local_verify(p,h):
            m,f=backend._read_archive_members(p,h);local_calls.append(m['global_step'])
        try:
            core.advance_segment_with_outputs(interrupted,rng_state_all=source.rng_state_all,
                source_provenance={'synthetic':True,'execution_authorized':False},output_root=stopped,
                remote_roundtrip=local_verify,clock=lambda:next(clock_values),max_wall_s=3600.)
        except core.SegmentRefused as exc:
            assert str(exc)=='SEGMENT_WALLTIME_BUDGET_EXCEEDED'
        assert interrupted.nsteps==2351 and local_calls[-1]==2351
        recovered=core.restore_segment_bundle(stopped/'segment-step-000002351',expected_global_step=2351)
        assert recovered.sampler_state['source_provenance']['segment_result']=='INCOMPLETE'
        report['timeout_atomic_checkpoint_and_restore']='PASS'
        continued=Counter(recovered,start=2351)
        core.advance_segment_with_outputs(continued,rng_state_all=recovered.rng_state_all,
            source_provenance={'synthetic':True,'execution_authorized':False},output_root=tmp/'continued',
            remote_roundtrip=local_verify,checkpoint_every=100)
        continuous=Counter(source)
        outputs=core.advance_segment_with_outputs(continuous,rng_state_all=source.rng_state_all,
            source_provenance={'synthetic':True,'execution_authorized':False},output_root=tmp/'continuous',
            remote_roundtrip=local_verify,checkpoint_every=100)
        r1=core.restore_segment_bundle(tmp/'continuous/segment-step-000003000',expected_global_step=3000)
        r2=core.restore_segment_bundle(tmp/'continued/segment-step-000003000',expected_global_step=3000)
        assert same(r1,r2)
        assert len(outputs)==10
        report['continuous_vs_interrupted_state_and_RNG']='PASS'
        probe=object.__new__(core.BoundedLangevin);probe.nsteps=3000;probe.segment_start=2000;probe.stop_at_step=3000
        try:probe.step()
        except core.SegmentRefused:report['step3001_rejected_before_ASE']='PASS'
        else:raise AssertionError('boundary not enforced')
        b=backend.GitHubArchiveBackend(remote_url='https://github.com/linkengde/-.git',
            branch='cloud-c-ag64-manual-v2-test-20261011',expected_head='bc3ff0631db3b7d78dcc17556befa98c2f454480',
            archive_prefix='C_SYNTHETIC_TRANSPORT_ONLY',author_name='Codex Cloud C',author_email='codex-cloud-c@users.noreply.github.com')
        arc=Path(outputs[-1]['archive'])
        receipt=b.upload(arc,outputs[-1]['archive_sha256'],start_step=2000,end_step=3000,max_added_steps=1000,source_checkpoint_sha256=cfg['source_checkpoint_state_sha256'])
        verified=b.verify(receipt);restored=b.restore(receipt,tmp/'fresh-remote-restored')
        remote=core.restore_segment_bundle(tmp/'fresh-remote-restored',expected_global_step=3000)
        assert same(remote,r1)
        report['backend_real_upload_verify_fresh_fetch_restore']='PASS'
        report['backend_receipt']=receipt;report['backend_verify']=verified;report['backend_restore']=restored
    report['code_sha256']={n:backend.sha256_file(ROOT/n) for n in ['continuation_adapter.py','core/segment_core.py','github_archive_backend.py']}
    report['runtime']={'python':platform.python_version(),'ase':ase.__version__,'numpy':np.__version__}
    (ROOT/'NON_SCIENTIFIC_VALIDATION.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['backend_receipt','backend_verify','backend_restore']},indent=2))
if __name__=='__main__':main()
