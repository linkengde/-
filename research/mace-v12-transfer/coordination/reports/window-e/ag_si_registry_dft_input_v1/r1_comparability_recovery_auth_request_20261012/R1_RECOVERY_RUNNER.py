#!/usr/bin/env python3
"""Run one frozen R1 recovery point only after C issues matching authorization."""
from __future__ import annotations
import hashlib, json, math, os, shutil, signal, sys, time, traceback, threading
from collections import Counter
from importlib.metadata import version
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent
REPO = Path('/workspace/-')
ROOT = Path(os.environ.get('R1_RECOVERY_OUTPUT_ROOT', '/workspace/e-r1-energy-comparability-recovery-authorized-20261012')).resolve()
CONFIG_PATH = PACKAGE / 'R1_RECOVERY_CONFIG.json'
AUTH_PATH = PACKAGE / 'R1_RECOVERY_EXEC_AUTH.json'
CONFIG = json.loads(CONFIG_PATH.read_text())
CONFIG_SHA = hashlib.sha256(CONFIG_PATH.read_bytes()).hexdigest()
INPUT = REPO / CONFIG['source_bindings']['input_path']
STRUCTURE = REPO / CONFIG['source_bindings']['structure_path']
LOCK = REPO / CONFIG['source_bindings']['parameter_lock_path']
R2_CONFIG = REPO / CONFIG['source_bindings']['matched_recovery_reference_config_path']
PAW_SHA = CONFIG['paw_sha256_by_element']
EXPECTED_CELL = [[3.071,0.0,0.0],[-1.5355,2.6595640150220112,0.0],[0.0,0.0,71.26826600017428]]
EXPECTED_COUNTS = {'Ag':7,'C':4,'Si':3,'Ti':6}

class SafeStop(Exception): pass
class CheckpointFailure(Exception): pass

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def write_json(path: Path, data: dict) -> None:
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    os.replace(tmp,path)

def cgroup_cpu_quota() -> float:
    q,p=Path('/sys/fs/cgroup/cpu.max').read_text().split()
    return int(q)/int(p)

def resources() -> dict:
    lim=int(Path('/sys/fs/cgroup/memory.max').read_text())
    cur=int(Path('/sys/fs/cgroup/memory.current').read_text())
    du=shutil.disk_usage(ROOT)
    st=os.statvfs(ROOT)
    free_inodes=st.f_favail
    estimate=CONFIG['checkpoint_policy']['per_checkpoint_estimate_bytes']
    reserve=CONFIG['checkpoint_policy']['free_space_reserve_bytes']
    minimum=CONFIG['checkpoint_policy']['initial_required_free_bytes']
    return {
      'cpu_quota_cores':cgroup_cpu_quota(),
      'memory_limit_bytes':lim,'memory_current_bytes':cur,'memory_headroom_bytes':lim-cur,
      'estimated_working_set_bytes':estimate,
      'disk_total_bytes':du.total,'disk_free_bytes':du.free,'free_inodes':free_inodes,
      'checkpoint_estimate_bytes':estimate,'free_space_reserve_bytes':reserve,
      'minimum_initial_free_bytes':minimum,
      'cpu_pass':cgroup_cpu_quota() >= 4.0,
      'memory_pass':lim-cur >= estimate,
      'disk_pass':du.free >= minimum,
      'inode_pass':free_inodes > 10000,
    }

def validate_inputs(check_mpi=False) -> tuple[dict,object,dict]:
    import ase, gpaw
    from ase.io import read
    from gpaw.setup_data import search_for_file
    from gpaw_data import datapath
    require(sha256(INPUT)==CONFIG['source_bindings']['input_sha256'],'R1 input SHA mismatch')
    require(sha256(STRUCTURE)==CONFIG['source_bindings']['structure_sha256'],'R1 structure SHA mismatch')
    require(sha256(LOCK)==CONFIG['source_bindings']['parameter_lock_sha256'],'Phase2 parameter-lock SHA mismatch')
    require(sha256(R2_CONFIG)==CONFIG['source_bindings']['matched_recovery_reference_config_sha256'],'R2 recovery reference config SHA mismatch')
    require(CONFIG['scope']['registry_id']=='R1' and CONFIG['method']['eigensolver']=={'name':'PPCG','niter':4},'scope/PPCG mismatch')
    require(CONFIG['method']['max_scf_iterations']==200,'maxiter must be 200')
    require(gpaw.__version__=='26.7.0' and ase.__version__=='3.29.0' and version('gpaw-data')=='1.2.1','engine version mismatch')
    expected_path=Path(datapath()).resolve()
    require(os.environ.get('GPAW_SETUP_PATH') is not None and Path(os.environ['GPAW_SETUP_PATH']).resolve()==expected_path,'GPAW_SETUP_PATH mismatch')
    require(os.environ.get('OMP_NUM_THREADS')=='1','OMP_NUM_THREADS must equal 1')
    atoms=read(STRUCTURE,format='extxyz')
    require(len(atoms)==20 and dict(sorted(Counter(atoms.get_chemical_symbols()).items()))==EXPECTED_COUNTS,'R1 composition mismatch')
    require(atoms.pbc.tolist()==[True,True,False],'R1 PBC mismatch')
    require(atoms.cell.array.tolist()==EXPECTED_CELL,'R1 cell mismatch')
    require(all(math.isfinite(float(x)) for x in atoms.positions.ravel()),'nonfinite R1 coordinates')
    observed={}
    for el,expected in PAW_SHA.items():
        filename,_=search_for_file(f'{el}.PBE')
        got=sha256(Path(filename)); require(got==expected,f'{el} PAW SHA mismatch')
        observed[el]={'path':str(filename),'sha256':got}
    mpi_size=None
    if check_mpi:
        from gpaw.mpi import world
        mpi_size=world.size
        require(mpi_size==4,'MPI world size must equal 4')
    root_sha=subprocess_git_head()
    import subprocess
    ancestor=subprocess.run(['git','-C',str(REPO),'merge-base','--is-ancestor',CONFIG['source_bindings']['fresh_fetched_e_commit'],root_sha])
    require(ancestor.returncode==0,'frozen E source commit is not an ancestor of the checked-out branch')
    return atoms,observed,{'python':sys.version.split()[0],'gpaw':gpaw.__version__,'ase':ase.__version__,'gpaw_data':version('gpaw-data'),'gpaw_setup_path':str(expected_path),'omp_threads':os.environ['OMP_NUM_THREADS'],'mpi_world_size':mpi_size,'paw_files':observed,'source_commit':root_sha}

def subprocess_git_head() -> str:
    import subprocess
    return subprocess.run(['git','-C',str(REPO),'rev-parse','HEAD'],check=True,capture_output=True,text=True).stdout.strip()

def require(ok: bool,msg: str) -> None:
    if not ok: raise RuntimeError(msg)

def resource_preflight() -> dict:
    r=resources()
    require(r['cpu_pass'] and r['memory_pass'] and r['disk_pass'] and r['inode_pass'],'resource preflight failed: '+json.dumps(r))
    return r

def authorization_gate() -> dict:
    require(AUTH_PATH.is_file(),'C execution authorization is missing; HOLD')
    auth=json.loads(AUTH_PATH.read_text())
    require(auth.get('execution_authorized') is True,'C authorization does not set execution_authorized=true')
    require(auth.get('registry_id')=='R1','C authorization registry_id must be R1')
    require(auth.get('scope')=='R1_ENERGY_COMPARABILITY_RECOVERY_SINGLE_POINT','C authorization scope mismatch')
    bindings=auth.get('bindings',{})
    expected={
      'config_sha256':CONFIG_SHA,
      'runner_sha256':sha256(Path(__file__)),
      'input_sha256':CONFIG['source_bindings']['input_sha256'],
      'structure_sha256':CONFIG['source_bindings']['structure_sha256'],
      'parameter_lock_sha256':CONFIG['source_bindings']['parameter_lock_sha256'],
    }
    for key,value in expected.items(): require(bindings.get(key)==value,f'C authorization {key} mismatch')
    require(auth.get('method')==CONFIG['method'],'C authorization method mismatch')
    return {'authorization_file':str(AUTH_PATH),'authorization_sha256':sha256(AUTH_PATH),'authorization_commit':auth.get('authorization_commit'),'bindings':bindings}

class Checkpoints:
    def __init__(self,calc,atoms,world,start_time,stop_event):
        self.calc=calc; self.atoms=atoms; self.world=world; self.start_time=start_time
        self.stop_event=stop_event; self.dir=ROOT/'checkpoints'; self.saved={}; self.latest=None
    def verify_readback(self,path):
        from gpaw.new.gpw import read_gpw
        from gpaw.mpi import serial_comm
        a,dft,_=read_gpw(path,comm=serial_comm)
        require((a.numbers==self.atoms.numbers).all(),'checkpoint atom numbers differ')
        require((a.positions==self.atoms.positions).all(),'checkpoint coordinates differ')
        require((a.cell.array==self.atoms.cell.array).all(),'checkpoint cell differs')
        require((a.pbc==self.atoms.pbc).all(),'checkpoint PBC differs')
        require(dft.ibzwfs.has_wave_functions(),'checkpoint lacks wavefunctions')
        require(dft.scf_loop.eigensolver.niter==4,'checkpoint PPCG niter mismatch')
        del a,dft
    def save(self,step,final=False,converged=False):
        self.dir.mkdir(parents=True,exist_ok=True)
        stem='checkpoint_final' if final else f'checkpoint_iter_{step:06d}'
        final_path=self.dir/(stem+'.gpw'); tmp=self.dir/('.'+stem+'.gpw.partial')
        local_error=None
        try:
            r=resources()
            require(r['disk_free_bytes']>=CONFIG['checkpoint_policy']['per_checkpoint_estimate_bytes']+CONFIG['checkpoint_policy']['free_space_reserve_bytes'],'checkpoint disk reserve no longer fits')
            require(not final_path.exists() and not tmp.exists(),'checkpoint path collision')
            self.calc.write(str(tmp),mode='all')
        except Exception as e:
            local_error=f'{type(e).__name__}: {e}'
            self.world.barrier()
        self.world.barrier()
        errors=self.world.sum_scalar(int(local_error is not None))
        payload=None
        if self.world.rank==0:
            if errors:
                payload={'ok':False,'error':local_error or 'MPI rank checkpoint write failed'}
            else:
                try:
                    require(tmp.is_file() and tmp.stat().st_size>0,'checkpoint file missing or empty')
                    size=tmp.stat().st_size
                    require(size<=CONFIG['checkpoint_policy']['per_checkpoint_estimate_bytes'],'checkpoint exceeded size estimate')
                    digest=sha256(tmp)
                    self.verify_readback(tmp)
                    os.replace(tmp,final_path)
                    sidecar=final_path.with_suffix('.json')
                    metadata={'status':'CHECKPOINT_VALID','registry_id':'R1','completed_scf_iteration':int(step),'scf_converged':bool(converged),'checkpoint_sha256':digest,'checkpoint_size_bytes':size,'checkpoint_path':final_path.name,'wavefunctions_included':True,'readback_verified':True,'config_sha256':CONFIG_SHA,'parameter_lock_sha256':CONFIG['source_bindings']['parameter_lock_sha256'],'input_sha256':CONFIG['source_bindings']['input_sha256'],'structure_sha256':CONFIG['source_bindings']['structure_sha256']}
                    write_json(sidecar,metadata)
                    payload={'ok':True,'path':str(final_path),'sidecar':str(sidecar),'sha256':digest,'size_bytes':size}
                except Exception as e:
                    payload={'ok':False,'error':f'{type(e).__name__}: {e}'}
        from gpaw.mpi import broadcast
        payload=broadcast(payload,0,comm=self.world)
        if not payload['ok']: raise CheckpointFailure(payload['error'])
        self.saved[int(step)]=payload; self.latest=payload
        return payload
    def on_step(self,ctx):
        step=int(ctx.niter)
        if step%40==0 and step not in self.saved: self.save(step)
        elapsed=time.monotonic()-self.start_time
        if self.stop_event.is_set() or elapsed>=CONFIG['method']['orderly_stop_seconds']:
            if step not in self.saved: self.save(step)
            raise SafeStop(f'orderly stop after checkpoint at SCF step {step}')
    def save_final(self,step):
        if int(step) in self.saved:
            return self.saved[int(step)]
        return self.save(int(step),final=True,converged=True)


def main():
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('--preflight',action='store_true')
    ap.add_argument('--run',action='store_true')
    args=ap.parse_args()
    require(args.preflight != args.run,'choose exactly one mode')
    ROOT.mkdir(parents=True,exist_ok=True)
    auth=authorization_gate() if args.run else None
    atoms,paws,env=validate_inputs(check_mpi=args.run)
    res=resource_preflight()
    pre={'status':'PREFLIGHT_PASS','authorization_status':'C_AUTHORIZED' if auth else 'PENDING_C_AUTH','execution_authorized':bool(auth),'r1_dft_started':False,'config_sha256':CONFIG_SHA,'source_bindings':CONFIG['source_bindings'],'atom_count':len(atoms),'formula':atoms.get_chemical_formula(),'cell_A':atoms.cell.array.tolist(),'pbc':atoms.pbc.tolist(),'runtime':env,'resources':res,'method':CONFIG['method'],'authorization':auth}
    if args.preflight:
        write_json(ROOT/'PRECHECK.json',pre)
        print(json.dumps(pre,indent=2)); return 0
    from gpaw.mpi import world
    require(not (ROOT/'RUN_STATUS.json').exists(),'refusing to overwrite a prior R1 run')
    world.barrier()
    if world.rank==0: write_json(ROOT/'RUN_STATUS.json',{'run_status':'RUNNING','scf_converged':False,'authorization_status':'C_AUTHORIZED','authorization_sha256':auth['authorization_sha256'],'config_sha256':CONFIG_SHA})
    world.barrier()
    from gpaw.new.ase_interface import GPAW
    from gpaw import FermiDirac,PW
    from gpaw.mixer import Mixer
    pos0=atoms.positions.copy(); cell0=atoms.cell.array.copy(); pbc0=atoms.pbc.copy()
    log=ROOT/'gpaw.log'
    kwargs={'mode':PW(500),'xc':'PBE','kpts':(8,8,1),'occupations':FermiDirac(0.1),'convergence':{'energy':1e-5,'density':1e-5,'eigenstates':1e-5},'maxiter':200,'mixer':Mixer(0.05,8,100),'eigensolver':{'name':'ppcg','niter':4},'symmetry':{'point_group':True,'time_reversal':True},'poissonsolver':{},'parallel':{'sl_auto':True},'txt':str(log)}
    calc=GPAW(**kwargs)
    start=time.monotonic(); stop_event=threading.Event(); checkpoints=Checkpoints(calc,atoms,world,start,stop_event)
    signal.signal(signal.SIGTERM,lambda *_: stop_event.set())
    signal.signal(signal.SIGINT,lambda *_: stop_event.set())
    calc.hooks['scf_step']=checkpoints.on_step
    atoms.calc=calc
    result=None; failure=None
    try:
        total=float(atoms.get_potential_energy(force_consistent=False))
        free=float(atoms.get_potential_energy(force_consistent=True))
        forces=atoms.get_forces(apply_constraint=False)
        stress=atoms.get_stress(voigt=False,apply_constraint=False)
        require((atoms.positions==pos0).all(),'R1 coordinates changed')
        require((atoms.cell.array==cell0).all(),'R1 cell changed')
        require((atoms.pbc==pbc0).all(),'R1 PBC changed')
        require(math.isfinite(total) and math.isfinite(free),'nonfinite energy')
        require(all(math.isfinite(float(x)) for x in forces.ravel()),'nonfinite forces')
        require(all(math.isfinite(float(x)) for x in stress.ravel()),'nonfinite stress')
        steps=int(calc.dft.scf_loop.niter)
        cp=checkpoints.save_final(steps)
        result={'status':'SUCCESS','registry_id':'R1','role':'REGISTRY_SCREENING_ONLY_COMPARABILITY_RECOVERY','scf_converged':True,'scf_steps':steps,'total_energy_eV':total,'force_consistent_free_energy_eV':free,'forces_eV_A':forces.tolist(),'forces_status':'AVAILABLE','stress_eV_A3':stress.tolist(),'stress_status':'AVAILABLE','checkpoint':cp,'input_sha256':CONFIG['source_bindings']['input_sha256'],'structure_sha256':CONFIG['source_bindings']['structure_sha256'],'parameter_lock_sha256':CONFIG['source_bindings']['parameter_lock_sha256'],'config_sha256':CONFIG_SHA,'authorization_sha256':auth['authorization_sha256'],'method':CONFIG['method'],'runtime':env,'paw_files':paws,'geometry_unchanged':True,'elapsed_seconds':time.monotonic()-start,'failure_reason':None}
    except SafeStop as e:
        failure={'run_status':'SAFE_STOP_UNCONVERGED','scf_converged':False,'scf_steps':int(getattr(calc.dft.scf_loop,'niter',0)),'failure_reason':str(e),'latest_checkpoint':checkpoints.latest}
    except BaseException as e:
        kind='FAILED_SCF' if type(e).__name__=='KohnShamConvergenceError' else ('FAILED_CHECKPOINT' if isinstance(e,CheckpointFailure) else 'FAILED_RUNTIME')
        failure={'run_status':kind,'scf_converged':False,'scf_steps':int(getattr(calc.dft.scf_loop,'niter',0)),'failure_reason':f'{type(e).__name__}: {e}','traceback':traceback.format_exc(),'latest_checkpoint':checkpoints.latest}
    if world.rank==0:
        if result is not None:
            write_json(ROOT/'R1_RECOVERY_RESULT.json',result)
            status={'run_status':'SUCCESS','scf_converged':True,'scf_steps':result['scf_steps'],'elapsed_seconds':result['elapsed_seconds'],'final_checkpoint_sha256':result['checkpoint']['sha256']}
        else:
            write_json(ROOT/'R1_FAILURE.json',failure)
            status={**failure,'authorization_status':'C_AUTHORIZED','authorization_sha256':auth['authorization_sha256'],'config_sha256':CONFIG_SHA}
        write_json(ROOT/'RUN_STATUS.json',status)
    world.barrier()
    if world.rank==0: print(json.dumps(status,indent=2))
    return 0 if result is not None else 5

if __name__=='__main__':
    try:
        code=main()
    except Exception as e:
        print(f'HOLD: {type(e).__name__}: {e}',file=sys.stderr)
        raise SystemExit(2)
    raise SystemExit(code)
