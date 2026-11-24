#!/workspace/venvs/gpaw/bin/python
"""Same-geometry 32-atom NVE check: V14 POD force5 versus literature Ag EAM."""
import json, os, subprocess
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
SRC = ROOT / 'v14_stage2_p0_ag_dft_diagnostics' / 'ag_liquid_precollapse' / 'label.extxyz'
BASE = ROOT / 'v14_stage2_p0_ag_dft_diagnostics' / 'ag_precollapse_nve_compare'
FIT = ROOT / 'fit_v14_targeted_ag_force5'
COEFF = FIT / 'v14_targeted_ag_f5_coefficients.pod'
POT = Path(os.environ.get('V14_AG_EAM', '/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/potentials/Ag_u3.eam'))
LMP = os.environ.get('LAMMPS_BIN', '/workspace/venvs/podfit/bin/lmp')
atoms = read(SRC, format='extxyz')
methods = {
    'V14_force5': f'''pair_style pod
pair_coeff * * {FIT/'AgTiSiC_POD_v14.pod'} {COEFF} Ag C Si Ti''',
    'Ag_u3_EAM': f'''pair_style eam
pair_coeff * * {POT}'''
}
summaries = {}
for name, pair in methods.items():
    out = BASE / name
    out.mkdir(parents=True, exist_ok=True)
    data, inp, log, dump = [out / f'candidate.{x}' for x in ('data', 'in', 'log', 'dump')]
    order = ['Ag'] if name == 'Ag_u3_EAM' else ['Ag', 'C', 'Si', 'Ti']
    write(data, atoms, format='lammps-data', atom_style='atomic', units='metal',
          specorder=order)
    masses = 'mass 1 107.8682' if name == 'Ag_u3_EAM' else '''mass 1 107.8682
mass 2 12.011
mass 3 28.0855
mass 4 47.867'''
    inp.write_text(f'''clear
units metal
atom_style atomic
boundary p p p
atom_modify sort 0 0
read_data {data}
{masses}
{pair}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo 10
thermo_style custom step time atoms temp pe ke etotal press c_fmax
thermo_modify lost error flush yes format float %.14g
timestep 0.001
run 0
velocity all create 1810.056942965326 291107 mom yes rot no dist gaussian
write_dump all custom {out / 'start.dump'} id type x y z modify sort id
fix integ all nve
dump trj all custom 10 {dump} id type x y z
dump_modify trj sort id
run 100
undump trj
unfix integ
''')
    env = os.environ.copy()
    env.update({'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
    p = subprocess.run([LMP, '-in', str(inp), '-log', str(log), '-screen', 'none'],
                       cwd=out, env=env, capture_output=True, text=True)
    if p.returncode:
        raise SystemExit(f'{name} failed: {p.stdout[-1000:]} {p.stderr[-1000:]}')
    rows = []; header = None
    for line in log.read_text(errors='replace').splitlines():
        fields = line.split()
        if fields and fields[0] == 'Step' and 'Temp' in fields:
            header = fields
        elif header and len(fields) == len(header):
            try:
                rows.append(dict(zip(header, map(float, fields))))
            except ValueError:
                pass
    frame_paths = [out / 'start.dump', dump]
    frames = []
    for path in frame_paths:
        frames.extend(read(path, index=':', format='lammps-dump-text'))
    dmin = []
    for fr in frames:
        ds = fr.get_all_distances(mic=True); np.fill_diagonal(ds, np.inf)
        dmin.append(float(ds.min()))
    rows = [r for r in rows if r['Temp'] > 1.0e-8]
    temps = [r['Temp'] for r in rows]
    energies = [r['TotEng'] for r in rows]
    summaries[name] = {'natoms_initial': int(rows[0]['Atoms']), 'natoms_final': int(rows[-1]['Atoms']),
                       'temperature_initial_K': temps[0], 'temperature_final_K': temps[-1],
                       'temperature_max_K': max(temps), 'min_Ag_Ag_initial_A': dmin[0],
                       'min_Ag_Ag_over_0p1ps_A': min(dmin), 'min_Ag_Ag_final_A': dmin[-1],
                       'Fmax_initial_eV_A': rows[0]['c_fmax'],
                       'Fmax_max_eV_A': max(r['c_fmax'] for r in rows),
                       'total_energy_drift_eV_per_atom': (energies[-1]-energies[0])/32,
                       'dump_frames': len(frames)}
result = {'initial_structure': str(SRC), 'initial_temperature_K': 1810.056942965326,
          'timestep_fs': 1.0, 'duration_ps': 0.1,
          'role': 'small-system exploratory stability comparison; no phase classification or long MD approval',
          'results': summaries}
BASE.mkdir(parents=True, exist_ok=True)
(BASE / 'summary.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
