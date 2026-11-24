#!/workspace/venvs/gpaw/bin/python
"""External force-only check of V14 on archived Ag cooling labels."""
import json, os, subprocess
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
OLD = Path('/workspace/potential-fit/round7_rebuild/labels')
FIT = Path(os.environ.get('V14_FIT_DIR', str(ROOT / 'fit_v14_targeted_ag')))
COEFF = Path(os.environ.get('V14_COEFF_FILE', str(FIT / 'v14_targeted_ag_coefficients.pod')))
OUT = Path(os.environ.get('V14_LEGACY_OUT', str(ROOT / 'v14_stage2_p0_ag_dft_diagnostics' / 'legacy_ag_cooling_v14_force_checks')))
OUT.mkdir(parents=True, exist_ok=True)
LMP = os.environ.get('LAMMPS_BIN', '/workspace/venvs/podfit/bin/lmp')
rows = []
for src in sorted(OLD.glob('ag_cooling_*_step*_pw350.extxyz')):
    label_json = src.with_suffix('.json')
    if not label_json.exists():
        continue
    ref = read(src, format='extxyz')
    f_dft = np.asarray(ref.get_forces(), dtype=float)
    tag = src.stem
    data, inp, log, dump = (OUT / f'{tag}.{ext}' for ext in ('data', 'in', 'log', 'dump'))
    write(data, ref, format='lammps-data', atom_style='atomic', units='metal',
          specorder=['Ag', 'C', 'Si', 'Ti'])
    inp.write_text(f'''clear
units metal
atom_style atomic
boundary p p p
atom_modify sort 0 0
read_data {data}
mass 1 107.8682
mass 2 12.011
mass 3 28.0855
mass 4 47.867
pair_style pod
pair_coeff * * {FIT/'AgTiSiC_POD_v14.pod'} {COEFF} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
thermo_style custom step atoms pe press
thermo 1
run 0
write_dump all custom {dump} id type x y z fx fy fz modify sort id
''')
    env = os.environ.copy()
    env.update({'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
    p = subprocess.run([LMP, '-in', str(inp), '-log', str(log), '-screen', 'none'],
                       cwd=OUT, env=env, capture_output=True, text=True)
    if p.returncode:
        raise SystemExit(f'LAMMPS failed for {src.name}: {p.stdout[-1000:]} {p.stderr[-1000:]}')
    lines = dump.read_text().splitlines(); k = lines.index('ITEM: ATOMS id type x y z fx fy fz')
    n = int(lines[3]); values = np.array([[float(x) for x in line.split()] for line in lines[k+1:k+1+n]])
    forces = values[np.argsort(values[:, 0].astype(int))][:, 5:8]
    delta = forces - f_dft; norm = np.linalg.norm(delta, axis=1)
    meta = json.loads(label_json.read_text())
    rows.append({'label_file': src.name,
                 'historical_role': 'training' if 'train' in src.name else 'validation',
                 'scf_iterations_recorded': meta.get('scf_iterations'),
                 'historical_protocol': meta.get('method'),
                 'historical_energy_shifted': meta.get('energy_aligned_eV') is not None,
                 'force_rmse_component_eV_A': float(np.sqrt(np.mean(delta**2))),
                 'force_rmse_vector_norm_eV_A': float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
                 'force_error_p95_eV_A': float(np.percentile(norm, 95)),
                 'force_error_max_eV_A': float(norm.max()),
                 'minimum_pair_distance_A': float(np.min(ref.get_all_distances(mic=True) + np.eye(len(ref))*1e6)),
                 'natoms': len(ref)})
summary = {'purpose': 'external force-only diagnostic; labels are old Gamma-point data with historical energy alignment, not V14 training data',
           'potential': 'V14 targeted Ag correction POD', 'labels_directory': str(OLD), 'results': rows}
(OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
