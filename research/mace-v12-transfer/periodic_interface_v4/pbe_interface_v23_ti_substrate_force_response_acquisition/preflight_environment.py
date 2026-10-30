#!/usr/bin/env python3
"""Read-only geometry, provenance, MPI and resource preflight for v23."""
from __future__ import annotations
import fcntl, hashlib, importlib.metadata as metadata, io, json, os, re, shutil, subprocess
from pathlib import Path
import numpy as np
from ase.io import read
from geometry_utils import compare_geometry, extxyz_geometry_frames

ROOT = Path(__file__).resolve().parent
REPO = Path('/workspace/jeitschko-main-run')
CGROUP = Path('/sys/fs/cgroup')
PYTHON = '/workspace/.venvs/gpaw-mpi/bin/python'
LIMIT = 34359738368
START_LIMIT = 3 * 1024**3
MIN_DISK = 5 * 1024**3
LOCK = Path('/workspace/.mace-v19-dft.lock')
DEP = '/workspace/.local/gpaw-mpi-deps/usr/lib/x86_64-linux-gnu'
V22_COMMIT = '5d1e734c16d17dc9b460e3184eab4670ccc2739c'
TRAIN49_COMMIT = 'e1a6c134b240a8afa9b44febc0bb646a4015ca3b'
V22_TASK = 'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v22_force_response_acquisition'
TRAIN49_TASK = 'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v21_jeitschko_interface_acquisition/training_candidate_v19_train49_20261030_A'

def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def git_blob(commit: str, path: str) -> bytes:
    return subprocess.check_output(['git', '-C', str(REPO), 'show', f'{commit}:{path}'])

def main() -> None:
    manifest = json.loads((ROOT / 'input_manifest.json').read_text())
    audit = json.loads((ROOT / 'geometry_audit.json').read_text())
    if manifest.get('status') != 'V23_GEOMETRY_SCREENED_ACQUISITION_CANDIDATES' or audit.get('status') != 'PASS':
        raise SystemExit('Unexpected v23 manifest/audit status')
    jobs = manifest['preflight_jobs']
    records = manifest['records']
    if len(jobs) != 6 or len(records) != 6 or len({x['label'] for x in jobs}) != 6:
        raise SystemExit('Expected six unique v23 jobs')
    train_dir = REPO / TRAIN49_TASK / 'dataset'
    train_manifest = json.loads((train_dir / 'manifest.json').read_text())
    train_raw = (train_dir / 'train.extxyz').read_bytes()
    dev_raw = (train_dir / 'valid.extxyz').read_bytes()
    if digest(train_raw) != train_manifest['train_sha256'] or digest(dev_raw) != train_manifest['dev_sha256']:
        raise SystemExit('train49 or frozen dev3 hash mismatch')
    ref = manifest['reference_datasets']
    if (digest((train_dir / 'manifest.json').read_bytes()) != ref['train49_manifest_sha256']
            or digest(train_raw) != ref['train49_train_sha256']
            or digest(dev_raw) != ref['frozen_dev3_sha256']):
        raise SystemExit('Pinned train49/development references changed')
    if subprocess.run(['git', '-C', str(REPO), 'merge-base', '--is-ancestor', V22_COMMIT, 'HEAD']).returncode != 0:
        raise SystemExit('Pinned v22 input commit is not an ancestor of current HEAD')
    if subprocess.run(['git', '-C', str(REPO), 'merge-base', '--is-ancestor', TRAIN49_COMMIT, 'HEAD']).returncode != 0:
        raise SystemExit('Pinned train49 geometry-screen commit is not an ancestor of current HEAD')
    v22 = json.loads(git_blob(V22_COMMIT, f'{V22_TASK}/input_manifest.json'))
    train_geos = extxyz_geometry_frames(train_dir / 'train.extxyz')
    dev_geos = extxyz_geometry_frames(train_dir / 'valid.extxyz')
    if len(train_geos) != 49 or len(dev_geos) != 3:
        raise SystemExit('Geometry-only training/development frame count mismatch')
    by_label = {row['label']: row for row in records}
    if set(by_label) != {x['label'] for x in jobs}:
        raise SystemExit('Record/job label mismatch')
    for job in jobs:
        record = by_label[job['label']]
        input_path = ROOT / job['input']
        raw = input_path.read_bytes()
        if digest(raw) != job['input_sha256'] or record['input_sha256'] != job['input_sha256']:
            raise SystemExit(f'Input hash mismatch: {job["label"]}')
        if not job['launch_enabled'] or job['kpts'] != [6, 6, 1] or job['smearing_eV'] != 0.10:
            raise SystemExit(f'Unexpected method/launch flag: {job["label"]}')
        parent_raw = git_blob(V22_COMMIT, record['parent_input'])
        if digest(parent_raw) != record['parent_input_sha256']:
            raise SystemExit(f'Parent input hash mismatch: {job["label"]}')
        parent = read(__import__('io').StringIO(parent_raw.decode()), format='extxyz')
        candidate = read(input_path)
        if len(candidate) != 52 or candidate.get_chemical_formula() != 'C16Ag4Si8Ti24':
            raise SystemExit(f'Unexpected candidate composition: {job["label"]}')
        if (candidate.get_chemical_symbols() != parent.get_chemical_symbols()
                or not np.array_equal(candidate.cell.array, parent.cell.array)
                or not np.array_equal(candidate.pbc, parent.pbc)):
            raise SystemExit(f'Parent structure/cell/PBC changed: {job["label"]}')
        symbols = np.asarray(parent.get_chemical_symbols())
        ti = np.flatnonzero(symbols == 'Ti')
        zmax = float(parent.positions[ti, 2].max())
        target = ti[np.isclose(parent.positions[ti, 2], zmax, rtol=0, atol=1e-6)]
        ids = np.asarray(parent.arrays.get('lammps_id', np.arange(len(parent))))
        target_ids = sorted(map(int, ids[target]))
        if len(target) != 4 or target_ids != sorted(record['target_atom_ids']):
            raise SystemExit(f'Target Ti orbit differs from registered selector: {job["label"]}')
        delta = candidate.positions - parent.positions
        expected = np.zeros_like(delta)
        expected[target, 2] = float(record['target_displacement_z_A'])
        if not np.allclose(delta, expected, rtol=0, atol=1e-8):
            raise SystemExit(f'Unexpected coordinate displacement: {job["label"]}')
        if ('REF_energy' in candidate.info or 'PW_PBE_energy_eV' in candidate.info
                or 'REF_forces' in candidate.arrays or 'PW_PBE_forces' in candidate.arrays):
            raise SystemExit(f'Candidate input already contains labels: {job["label"]}')
        parent_matches = []
        for idx, frame in enumerate(train_geos):
            comp = compare_geometry(parent, frame)
            if comp is not None and comp['exact_duplicate']:
                parent_matches.append(idx)
        if parent_matches != [record['intentional_parent_train49_geometry_frame_index']]:
            raise SystemExit(f'Parent-frame geometry mapping changed: {job["label"]}: {parent_matches}')
        geom = extxyz_geometry_frames(input_path)
        if len(geom) != 1:
            raise SystemExit(f'Input geometry frame count mismatch: {job["label"]}')
        train_comparisons = [(idx, compare_geometry(geom[0], frame)) for idx, frame in enumerate(train_geos)]
        train_comparisons = [(idx, x) for idx, x in train_comparisons if x is not None]
        train_exact = [idx for idx, x in train_comparisons if x['exact_duplicate']]
        train_near = [idx for idx, x in train_comparisons if x['near_duplicate']]
        parent_idx = record['intentional_parent_train49_geometry_frame_index']
        if train_exact or train_near != [parent_idx]:
            raise SystemExit(f'Unexpected train49 geometry hits for {job["label"]}: exact={train_exact}, near={train_near}')
        dev_comparisons = [compare_geometry(geom[0], frame) for frame in dev_geos]
        dev_comparisons = [x for x in dev_comparisons if x is not None]
        if any(x['exact_duplicate'] or x['near_duplicate'] for x in dev_comparisons):
            raise SystemExit(f'Unexpected frozen dev3 geometry hit: {job["label"]}')
    if int((CGROUP / 'memory.max').read_text()) != LIMIT:
        raise SystemExit('Unexpected cgroup memory ceiling')
    memory = int((CGROUP / 'memory.current').read_text())
    if memory >= START_LIMIT:
        raise SystemExit(f'Starting memory exceeds 3 GiB: {memory}')
    free = shutil.disk_usage(ROOT).free
    if free < MIN_DISK:
        raise SystemExit(f'Free disk below 5 GiB: {free}')
    active = []
    for proc in Path('/proc').iterdir():
        if proc.name.isdigit():
            try:
                cmdline = (proc / 'cmdline').read_bytes().decode(errors='ignore')
            except OSError:
                continue
            if 'run_preflight_label.py' in cmdline and int(proc.name) != os.getpid():
                active.append(proc.name)
    if active:
        raise SystemExit(f'Another DFT runner is active: {active}')
    env = os.environ.copy()
    env.update({'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1'})
    env['LD_LIBRARY_PATH'] = ':'.join([DEP + '/openmpi/lib', DEP,
        '/workspace/.local/gpaw-libs/usr/lib/x86_64-linux-gnu',
        '/workspace/.local/openblas/usr/lib/x86_64-linux-gnu/openblas-pthread', env.get('LD_LIBRARY_PATH', '')])
    info = subprocess.check_output([str(Path(PYTHON).with_name('gpaw')), 'info'], env=env, text=True)
    if not re.search(r'MPI enabled\s+yes', info) or not re.search(r'scalapack\s+yes', info):
        raise SystemExit('GPAW MPI/ScaLAPACK preflight failed')
    method = manifest['method_candidate']
    probe = "import hashlib,importlib.metadata as m,gpaw_data,pathlib; p=pathlib.Path(gpaw_data.datapath()); print(m.version('gpaw'),m.version('ase'),m.version('gpaw-data')); [print(n,hashlib.sha256((p/n).read_bytes()).hexdigest()) for n in " + repr(list(method['PAW_setup_sha256'])) + "]"
    pkg = subprocess.check_output([PYTHON, '-c', probe], env=env, text=True).splitlines()
    if pkg[0] != '26.7.0 3.29.0 1.2.1':
        raise SystemExit(f'GPAW/ASE/data version mismatch: {pkg[0]}')
    paw_hashes = dict(x.split() for x in pkg[1:])
    if paw_hashes != method['PAW_setup_sha256']:
        raise SystemExit('Installed PAW hashes differ from registered method')
    with LOCK.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        fcntl.flock(lock, fcntl.LOCK_UN)
    for job in jobs:
        if (ROOT / 'calculations' / job['label']).exists() or (ROOT / 'calculations' / ('.running_' + job['label'])).exists():
            raise SystemExit(f'Existing calculation/scratch found; refusing to overwrite: {job["label"]}')
    result = {'status': 'V23_GEOMETRY_AND_ENVIRONMENT_PREFLIGHT_PASS',
              'candidate_count': len(jobs), 'all_input_hashes_pass': True,
              'all_train49_dev3_geometry_screens_clear': True,
              'intentional_train49_parent_near_hits': 6,
              'train49_sha256': ref['train49_train_sha256'],
              'frozen_dev3_sha256': ref['frozen_dev3_sha256'],
              'cgroup_memory_limit_bytes': LIMIT, 'memory_current_bytes': memory,
              'free_disk_bytes': free, 'mpi_enabled': True, 'scalapack_enabled': True,
              'GPAW_version': '26.7.0', 'ASE_version': '3.29.0',
              'PAW_setup_hashes_pass': True, 'shared_DFT_lock_available': True,
              'maximum_concurrent_jobs': 1, 'automatic_retries': 0,
              'official_split_mutated': False, 'test_or_sealed_data_opened': False}
    (ROOT / 'preflight_environment.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
