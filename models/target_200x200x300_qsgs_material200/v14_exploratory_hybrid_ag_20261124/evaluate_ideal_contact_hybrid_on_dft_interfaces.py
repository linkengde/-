#!/workspace/venvs/gpaw/bin/python
"""Force-only check of the literature Ag EAM + MAX TSCBOP + Morse v0 interface route."""
import json, os, subprocess
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
V13 = ROOT / 'fit_v13'
SOURCES = [V13 / 'training' / f'V13_IFACE_{i:02d}.xyz' for i in range(1, 6)]
SOURCES += [V13 / 'validation' / 'V13_VAL_INTERFACE_01.xyz']
OUT = ROOT / 'v14_stage2_p0_ag_dft_diagnostics' / 'hybrid_interface_force_checks'
OUT.mkdir(parents=True, exist_ok=True)
POT = Path(os.environ.get('V14_POT_DIR', '/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/potentials'))
LMP = os.environ.get('LAMMPS_BIN', '/workspace/venvs/podfit/bin/lmp')
rows = []
for src in SOURCES:
    ref = read(src, format='extxyz'); f_dft = np.asarray(ref.get_forces(), dtype=float)
    tag = src.stem
    data, inp, log, dump = [OUT / f'{tag}.{x}' for x in ('data', 'in', 'log', 'dump')]
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
pair_style hybrid eam tersoff morse 8.0
pair_coeff 1 1 eam {POT/'Ag_u3.eam'}
pair_coeff * * tersoff {POT/'TSCBOP.tersoff'} NULL C Si Ti
pair_coeff 1 2 morse 1.35077230 2.0 2.15255924 8.0
pair_coeff 1 3 morse 0.93944948 2.0 2.40407837 8.0
pair_coeff 1 4 morse 0.94814194 2.0 2.33708864 8.0
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
thermo_style custom step atoms pe press
thermo 1
run 0
write_dump all custom {dump} id type x y z fx fy fz modify sort id
''')
    env = os.environ.copy(); env.update({'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
    p = subprocess.run([LMP, '-in', str(inp), '-log', str(log), '-screen', 'none'],
                       cwd=OUT, env=env, capture_output=True, text=True)
    if p.returncode:
        raise SystemExit(f'{tag}: {p.stdout[-1000:]} {p.stderr[-1000:]}')
    lines = dump.read_text().splitlines(); k = lines.index('ITEM: ATOMS id type x y z fx fy fz')
    n = int(lines[3]); a = np.array([[float(x) for x in row.split()] for row in lines[k+1:k+1+n]])
    f_model = a[np.argsort(a[:, 0].astype(int))][:, 5:8]
    diff = f_model - f_dft; err = np.linalg.norm(diff, axis=1)
    by_element = {}
    symbols = np.array(ref.get_chemical_symbols())
    for s in sorted(set(symbols)):
        idx = symbols == s
        by_element[s] = {'count': int(idx.sum()),
                         'force_RMSE_vector_eV_A': float(np.sqrt(np.mean(np.sum(diff[idx]**2, axis=1)))),
                         'force_error_P95_eV_A': float(np.percentile(err[idx], 95)),
                         'force_error_max_eV_A': float(err[idx].max())}
    rows.append({'configuration': tag, 'natoms': len(ref), 'formula': ref.get_chemical_formula(),
                 'min_Ag_skeleton_A': float(ref.get_all_distances(mic=True)[np.ix_(np.where(symbols=='Ag')[0], np.where(symbols!='Ag')[0])].min()),
                 'force_RMSE_component_eV_A': float(np.sqrt(np.mean(diff**2))),
                 'force_RMSE_vector_eV_A': float(np.sqrt(np.mean(np.sum(diff**2, axis=1)))),
                 'force_error_P95_eV_A': float(np.percentile(err, 95)),
                 'force_error_max_eV_A': float(err.max()), 'by_element': by_element})
summary = {'purpose': 'exploratory force-only check; verifies the idealized interface approximation against existing DFT interface labels',
           'model': 'Foiles Ag_u3 EAM + TSCBOP Ti3SiC2 + calibration-v0 cross Morse',
           'cross_parameters_are_independently_validated': False,
           'results': rows}
(OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
