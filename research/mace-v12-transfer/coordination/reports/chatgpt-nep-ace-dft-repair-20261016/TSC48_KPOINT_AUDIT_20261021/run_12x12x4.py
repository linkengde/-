from __future__ import annotations
import hashlib, json, os, sys, time, traceback
from pathlib import Path
from ase.io import read
import gpaw, ase, relay_runner
from gpaw import PW
from gpaw.dft import PPCG
from gpaw.mixer import Mixer
from gpaw.mpi import world
from gpaw.occupations import FermiDirac

root = Path('/workspace/h-dft-runs/tsc48-kpoint-audit-20261021')
out = root / 'grid_12x12x4_retry'
source = root / 'input.extxyz'
expected_input = 'da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b'
expected_paw = {
    'C.PBE.gz': 'fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983',
    'Si.PBE.gz': '2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0',
    'Ti.PBE.gz': 'bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f',
}
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()

assert world.size == 4, f'expected 4 MPI ranks, got {world.size}'
assert gpaw.__version__ == '26.7.0', gpaw.__version__
assert ase.__version__ == '3.29.0', ase.__version__
assert sha(source) == expected_input
setup_dir=Path(os.environ['GPAW_SETUP_PATH'])
for name, expected in expected_paw.items():
    assert sha(setup_dir/name) == expected, f'PAW mismatch: {name}'
if world.rank == 0:
    lock=out/'ACTIVE_RUN.json'
    fd=os.open(lock, os.O_CREAT|os.O_EXCL|os.O_WRONLY, 0o600)
    with os.fdopen(fd,'w') as f: json.dump({'pid':os.getpid(),'mesh':[12,12,4],'started_unix':time.time()},f)
world.barrier()
atoms=read(source)
assert len(atoms)==48 and all(atoms.pbc)
calc=relay_runner.compatible_factory(gpaw.GPAW,
    mode=PW(500), xc='PBE', kpts=(12,12,4),
    occupations=FermiDirac(0.10), eigensolver=PPCG(niter=2),
    mixer=Mixer(0.05,8,100),
    convergence={'energy':1e-5,'density':1e-5,'eigenstates':1e-5},
    maxiter=160, spinpol=False, txt=str(out/'gpaw_12x12x4.log'))
atoms.calc=calc
t0=time.monotonic()
try:
    e_native=float(atoms.get_potential_energy())
    e_free=float(atoms.get_potential_energy(force_consistent=True))
    forces=atoms.get_forces().tolist()
    stress=atoms.get_stress().tolist()
    elapsed=time.monotonic()-t0
    record={
      'status':'SCF_RESULT_RETURNED',
      'mesh':[12,12,4],
      'mesh_convention':'standard Monkhorst-Pack; gamma omitted; no offset',
      'atom_count':len(atoms),
      'input_sha256':sha(source),
      'gpaW_version':gpaw.__version__,
      'ase_version':ase.__version__,
      'mpi_ranks':world.size,
      'xc':'PBE','cutoff_eV':500,
      'occupations':'FermiDirac(0.10 eV)',
      'eigensolver':'PPCG(niter=2)',
      'mixer':'Mixer(0.05,8,100)',
      'convergence':{'energy':1e-5,'density':1e-5,'eigenstates':1e-5},
      'maxiter':160,'spinpol':False,'fixed_geometry':True,
      'energy_native_eV':e_native,'free_energy_eV':e_free,
      'forces_eV_A':forces,'stress_voigt_eV_A3':stress,
      'elapsed_seconds':elapsed,
      'paw_sha256':{k:sha(setup_dir/k) for k in expected_paw},
    }
    if world.rank==0:
        (out/'result.json').write_text(json.dumps(record,indent=2)+'\n')
except BaseException as e:
    if world.rank==0:
        (out/'failure.json').write_text(json.dumps({'status':'FAILED','error':repr(e),'traceback':traceback.format_exc()},indent=2)+'\n')
    raise
finally:
    world.barrier()
    if world.rank==0:
        try: (out/'ACTIVE_RUN.json').unlink()
        except FileNotFoundError: pass
