#!/workspace/venvs/gpaw/bin/python
"""Static Ag FCC EOS check for the targeted V14 POD fit."""
import json, os, subprocess
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.io import write

ROOT = Path(__file__).resolve().parent
FIT = Path(os.environ.get('V14_FIT_DIR', str(ROOT / 'fit_v14_targeted_ag')))
COEFF = Path(os.environ.get('V14_COEFF_FILE', str(FIT / 'v14_targeted_ag_coefficients.pod')))
OUT = Path(os.environ.get('V14_EOS_OUT', str(ROOT / 'v14_stage2_p0_ag_dft_diagnostics' / 'ag_fcc_v14_static')))
LMP = os.environ.get('LAMMPS_BIN', '/workspace/venvs/podfit/bin/lmp')
OUT.mkdir(parents=True, exist_ok=True)
rows = []
for a in [4.09, 3.90, 3.75, 3.60, 3.30]:
    tag = f'a_{a:.2f}A'
    data, inp, log = (OUT / f'{tag}.{ext}' for ext in ('data', 'in', 'log'))
    atoms = bulk('Ag', 'fcc', a=a, cubic=True).repeat((4, 4, 4))
    write(data, atoms, format='lammps-data', atom_style='atomic', units='metal',
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
pair_coeff * * {FIT / 'AgTiSiC_POD_v14.pod'} {COEFF} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo_style custom step atoms pe press c_fmax lx
thermo 1
run 0
''')
    env = os.environ.copy()
    env.update({'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
    p = subprocess.run([LMP, '-in', str(inp), '-log', str(log), '-screen', 'none'],
                       cwd=OUT, env=env, text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
    if p.returncode:
        raise SystemExit(p.stdout[-2000:])
    header = datum = None
    for line in log.read_text(errors='replace').splitlines():
        fields = line.split()
        if fields and fields[0] == 'Step' and 'PotEng' in fields:
            header = fields
        elif header and len(fields) == len(header):
            try:
                datum = dict(zip(header, map(float, fields)))
            except ValueError:
                pass
    if datum is None:
        raise RuntimeError(f'No thermo row for {tag}')
    rows.append({'a_A': a, 'relative_energy_to_a4.09_eV_per_atom': None,
                 'energy_eV_per_atom': datum['PotEng'] / 256,
                 'pressure_bar': datum['Press'], 'fmax_eV_A': datum['c_fmax'],
                 'nearest_neighbor_A': a / np.sqrt(2), 'volume_per_atom_A3': a**3 / 4,
                 'natoms': 256, 'role': 'diagnostic_only' if a == 3.30 else 'train_or_heldout_check'})
base = next(r['energy_eV_per_atom'] for r in rows if r['a_A'] == 4.09)
for r in rows:
    r['relative_energy_to_a4.09_eV_per_atom'] = r['energy_eV_per_atom'] - base
summary = {'potential': 'V14 targeted Ag correction, POD coefficients fit on 31 configs',
           'status': 'static diagnostic; not an MD-use approval',
           'rows': rows}
(OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
