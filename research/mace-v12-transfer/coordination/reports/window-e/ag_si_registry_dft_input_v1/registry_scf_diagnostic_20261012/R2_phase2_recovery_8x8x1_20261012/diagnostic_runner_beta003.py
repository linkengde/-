from __future__ import annotations
import hashlib, importlib.metadata, json, os, re, sys, time, traceback
from datetime import datetime, timezone
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from ase import __version__ as ase_version
from gpaw import __version__ as gpaw_version
from gpaw import FermiDirac, GPAW, PW
from gpaw.mixer import Mixer
from gpaw.mpi import world
from gpaw.setup_data import search_for_file
from gpaw_data import datapath
from input import load_frozen_structure, PAW_SHA256, STRUCTURE_SHA256

HERE = Path(__file__).resolve().parent
RUN = HERE / 'run_maxiter300_mixer_beta005_nmaxold003'
LOCK_SHA = 'd97de3fe5a26ab6529d5e6e8ceef15adc1f3a2f548f03b59f9f092361c2c0c78'
BASELINE_COMMIT = '0dd8628e4a91e88f31aba3fd7bbce1c8a40a2aa2'
EXPECTED_INPUT_SHA = '7c1327f403c94240e6181da990fe7b078f420cc32f92a68bee3282792b0bcc66'
EXPECTED_MESH = [8, 8, 1]
MAXITER = 300
MIXER = {'name': 'Mixer', 'beta': 0.05, 'nmaxold': 3, 'weight': 100}

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def parse_scf_log(path: Path):
    text = path.read_text(encoding='utf-8', errors='replace') if path.exists() else ''
    converged = list(re.finditer(r'Converged in\s+(\d+)\s+steps', text, re.I))
    if converged:
        return True, int(converged[-1].group(1))
    iterations = [int(x) for x in re.findall(r'\|iter:\s*(\d+)\|', text)]
    return (False, max(iterations)) if iterations else (None, None)

if os.environ.get('OMP_NUM_THREADS') != '1':
    raise RuntimeError('OMP_NUM_THREADS must equal 1')
if world.size != 4:
    raise RuntimeError(f'expected MPI world size 4, got {world.size}')
if gpaw_version != '26.7.0' or ase_version != '3.29.0':
    raise RuntimeError(f'version mismatch GPAW={gpaw_version} ASE={ase_version}')
if importlib.metadata.version('gpaw-data') != '1.2.1':
    raise RuntimeError('gpaw-data version mismatch')
setup_path = os.environ.get('GPAW_SETUP_PATH')
if not setup_path or Path(setup_path).resolve() != Path(datapath()).resolve():
    raise RuntimeError('GPAW_SETUP_PATH mismatch')

atoms, raw, input_checks = load_frozen_structure()
input_sha = sha(HERE / 'input.py')
structure_sha = sha(HERE / 'structure.xyz')
if input_sha != EXPECTED_INPUT_SHA:
    raise RuntimeError(f'input SHA mismatch: {input_sha}')
if structure_sha != STRUCTURE_SHA256:
    raise RuntimeError(f'structure SHA mismatch: {structure_sha}')
if atoms.pbc.tolist() != [True, True, False] or len(atoms) != 20:
    raise RuntimeError('frozen geometry count/PBC mismatch')
if world.rank == 0:
    RUN.mkdir(parents=False, exist_ok=False)
world.barrier()

paw_files = {}
for element, expected in PAW_SHA256.items():
    paw_path, _ = search_for_file(f'{element}.PBE')
    actual = sha(Path(paw_path))
    if actual != expected:
        raise RuntimeError(f'PAW SHA mismatch for {element}: {actual}')
    paw_files[element] = {'path': str(paw_path), 'sha256': actual}

if world.rank == 0:
    preflight = {
        'status': 'PASS', 'baseline_phase2_commit': BASELINE_COMMIT,
        'input_sha256': input_sha, 'structure_sha256': structure_sha,
        'phase2_parameter_lock_sha256': LOCK_SHA, 'role': 'REGISTRY_SCREENING_ONLY_DIAGNOSTIC_COPY',
        'kpoint_mesh': EXPECTED_MESH, 'fixed_geometry': True,
        'mpi_ranks': world.size, 'omp_threads': 1,
        'python': sys.version, 'gpaw': gpaw_version, 'ase': ase_version,
        'gpaw_data': importlib.metadata.version('gpaw-data'),
        'paw_files': paw_files, 'geometry_check': input_checks,
        'diagnostic_parameters': {'maxiter': MAXITER, 'mixer': MIXER},
    }
    (RUN / 'DIAGNOSTIC_PREFLIGHT.json').write_text(json.dumps(preflight, indent=2) + '\n')
world.barrier()

before_positions = atoms.positions.copy()
before_cell = atoms.cell.array.copy()
before_pbc = atoms.pbc.copy()
log_path = RUN / 'gpaw.log'
result_path = RUN / 'result.json'
started = datetime.now(timezone.utc).isoformat()
t0 = time.monotonic()
energy = free_energy = forces = stress = None
error = None
calc = GPAW(
    mode=PW(500), xc='PBE', kpts=tuple(EXPECTED_MESH),
    occupations=FermiDirac(0.1),
    convergence={'energy': 1e-5, 'density': 1e-5, 'eigenstates': 1e-5},
    maxiter=MAXITER, mixer=Mixer(MIXER['beta'], MIXER['nmaxold'], MIXER['weight']),
    symmetry={'point_group': True, 'time_reversal': True},
    poissonsolver={}, parallel={'sl_auto': True}, txt=str(log_path),
)
atoms.calc = calc
try:
    energy = float(atoms.get_potential_energy(force_consistent=False))
    free_energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = np.asarray(atoms.get_forces(apply_constraint=False), dtype=float)
    stress = np.asarray(atoms.get_stress(voigt=False, apply_constraint=False), dtype=float)
except Exception as exc:
    error = f'{type(exc).__name__}: {exc}'
    if world.rank == 0:
        (RUN / 'stderr.log').write_text(traceback.format_exc())
world.barrier()
geometry_unchanged = (np.array_equal(atoms.positions, before_positions)
    and np.array_equal(atoms.cell.array, before_cell) and np.array_equal(atoms.pbc, before_pbc))
if not geometry_unchanged:
    raise RuntimeError('geometry changed during fixed-geometry diagnostic')
converged, steps = parse_scf_log(log_path)
elapsed = time.monotonic() - t0
world.barrier()
if world.rank == 0:
    status = 'SUCCESS' if converged is True and error is None else ('FAILED_SCF' if converged is False else 'FAILED_RUNTIME')
    result = {
        'structure_id': 'R2', 'diagnostic_only': True, 'status': status,
        'baseline_phase2_commit': BASELINE_COMMIT,
        'input_sha256': input_sha, 'structure_sha256': structure_sha,
        'phase2_parameter_lock_sha256': LOCK_SHA,
        'method': 'PBE/PAW fixed-geometry single point',
        'kpoint_mesh': EXPECTED_MESH, 'maxiter': MAXITER, 'mixer': MIXER,
        'total_energy_eV': energy, 'force_consistent_free_energy_eV': free_energy,
        'forces_eV_A': None if forces is None else forces.tolist(),
        'stress_eV_A3': None if stress is None else stress.tolist(),
        'forces_status': 'AVAILABLE' if forces is not None else 'NOT_AVAILABLE',
        'stress_status': 'AVAILABLE' if stress is not None else 'NOT_AVAILABLE',
        'scf_converged': converged, 'scf_steps': steps, 'failure_reason': error,
        'geometry_unchanged': geometry_unchanged, 'role': 'REGISTRY_SCREENING_ONLY_DIAGNOSTIC_COPY',
        'paw_files': paw_files, 'started_utc': started,
        'finished_utc': datetime.now(timezone.utc).isoformat(),
        'elapsed_seconds': elapsed,
    }
    result_path.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': status, 'scf_steps': steps, 'elapsed_seconds': elapsed,
        'energy_eV': energy, 'free_energy_eV': free_energy, 'failure_reason': error}, indent=2))
world.barrier()
