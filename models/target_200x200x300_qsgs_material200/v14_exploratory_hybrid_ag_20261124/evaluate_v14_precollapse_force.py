#!/workspace/venvs/gpaw/bin/python
"""Recheck the pre-collapse Ag DFT candidate with fitted V14 forces."""
import json, os, subprocess
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
FIT = Path(os.environ.get('V14_FIT_DIR', str(ROOT / 'fit_v14_targeted_ag')))
COEFF = Path(os.environ.get('V14_COEFF_FILE', str(FIT / 'v14_targeted_ag_coefficients.pod')))
SRC = ROOT / 'v14_stage2_p0_ag_dft_diagnostics' / 'ag_liquid_precollapse' / 'label.extxyz'
OUT = Path(os.environ.get('V14_LIQUID_OUT', str(ROOT / 'v14_stage2_p0_ag_dft_diagnostics' / 'ag_liquid_precollapse_v14')))
OUT.mkdir(parents=True, exist_ok=True)
LMP = os.environ.get('LAMMPS_BIN', '/workspace/venvs/podfit/bin/lmp')
ref = read(SRC, format='extxyz')
f_dft = np.asarray(ref.get_forces(), dtype=float)
data, inp, log, dump = [OUT / f'candidate.{x}' for x in ('data', 'in', 'log', 'dump')]
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
    raise SystemExit(p.stdout[-1000:] + p.stderr[-1000:])
lines = dump.read_text().splitlines()
k = lines.index('ITEM: ATOMS id type x y z fx fy fz')
n = int(lines[3])
val = np.array([[float(x) for x in row.split()] for row in lines[k+1:k+1+n]])
f_v14 = val[np.argsort(val[:, 0].astype(int))][:, 5:8]
delta = f_v14 - f_dft
norm = np.linalg.norm(delta, axis=1)
d = ref.get_all_distances(mic=True)
np.fill_diagonal(d, np.inf)
i, j = np.unravel_index(np.argmin(d), d.shape)
rij = ref.positions[j] - ref.positions[i]
rij -= np.round(rij @ np.linalg.inv(ref.cell.array)) @ ref.cell.array
rhat = rij / np.linalg.norm(rij)
summary = {
    'candidate_id': 'V14_AG_LIQUID_PRECOLLAPSE_0060',
    'natoms': len(ref),
    'scf_converged': True,
    'force_rmse_component_eV_A': float(np.sqrt(np.mean(delta**2))),
    'force_rmse_vector_norm_eV_A': float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
    'force_error_p95_eV_A': float(np.percentile(norm, 95)),
    'force_error_max_eV_A': float(norm.max()),
    'minimum_pair': {
        'indices_0based': [int(i), int(j)], 'distance_A': float(d[i, j]),
        'dft_outward_force_on_j_eV_A': float(np.dot(f_dft[j], rhat)),
        'v14_outward_force_on_j_eV_A': float(np.dot(f_v14[j], rhat)),
        'dft_outward_force_on_i_eV_A': float(np.dot(f_dft[i], -rhat)),
        'v14_outward_force_on_i_eV_A': float(np.dot(f_v14[i], -rhat))
    },
    'note': 'This is a training candidate. Passing it does not constitute an independent liquid validation.'
}
(OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
