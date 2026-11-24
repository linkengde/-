#!/workspace/venvs/gpaw/bin/python
"""Stage-1 V14 baseline: relax and run an Ag FCC crystal without a heat source."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.io import read, write

WORK = Path(__file__).resolve().parent
FIT = WORK / 'fit_v13'
if not FIT.exists():
    FIT = WORK / 'fit_v13_weight100'
OUT = WORK / 'v14_stage1_pure_ag_fcc'
LMP = os.environ.get('LAMMPS_EXECUTABLE', '/workspace/venvs/podfit/bin/lmp')
SPECIES = ['Ag', 'C', 'Si', 'Ti']


def sha(p: Path):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run_lmp(run: Path, input_name: str):
    env = os.environ.copy()
    env.update({'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
    p = subprocess.run([LMP, '-in', input_name, '-log', 'log.lammps', '-screen', 'none'],
                       cwd=run, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (run / 'screen.out').write_text(p.stdout)
    return p.returncode, p.stdout


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    crystal = bulk('Ag', 'fcc', a=4.09, cubic=True).repeat((4, 4, 4))
    crystal.pbc = True
    start = OUT / 'fcc_256_initial.data'
    write(start, crystal, format='lammps-data', atom_style='atomic', units='metal', specorder=SPECIES)
    # The POD pair_coeff maps all four model species, so LAMMPS data must
    # declare four atom types even though this baseline contains Ag only.
    data_text = start.read_text()
    if '\n4 atom types\n' not in data_text:
        if '\n1 atom types\n' not in data_text:
            raise RuntimeError('Unexpected ASE LAMMPS data header; expected one or four atom types')
        start.write_text(data_text.replace('\n1 atom types\n', '\n4 atom types\n', 1))
    relax = OUT / 'zero_pressure_relax'
    relax.mkdir(exist_ok=True)
    opt_data = relax / 'fcc_256_relaxed.data'
    (relax / 'in.relax.lmp').write_text(f'''clear
units metal
atom_style atomic
boundary p p p
atom_modify sort 0 0
read_data {start}
mass 1 107.8682
mass 2 12.011
mass 3 28.0855
mass 4 47.867
pair_style pod
pair_coeff * * {FIT / 'AgTiSiC_POD_v13.pod'} {FIT / 'v13_coefficients.pod'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo 10
thermo_style custom step atoms pe press pxx pyy pzz lx ly lz c_fmax
thermo_modify flush yes lost error
fix relax all box/relax iso 0.0 vmax 0.0005
min_style cg
minimize 1.0e-8 1.0e-6 1000 10000
unfix relax
run 0
write_data {opt_data}
''')
    code, stdout = run_lmp(relax, 'in.relax.lmp')
    if code:
        raise SystemExit(f'FCC zero-pressure relaxation failed: {stdout[-2000:]}')
    optimized = read(opt_data, format='lammps-data', atom_style='atomic', units='metal',
                     Z_of_type={1: 47, 2: 6, 3: 14, 4: 22})
    if not set(optimized.get_chemical_symbols()) == {'Ag'}:
        raise RuntimeError('FCC relaxed data type mapping is not pure Ag')
    test = OUT / 'nve_300K_1ps'
    test.mkdir(exist_ok=True)
    (test / 'in.nve.lmp').write_text(f'''clear
units metal
atom_style atomic
boundary p p p
atom_modify sort 0 0
read_data {opt_data}
mass 1 107.8682
mass 2 12.011
mass 3 28.0855
mass 4 47.867
pair_style pod
pair_coeff * * {FIT / 'AgTiSiC_POD_v13.pod'} {FIT / 'v13_coefficients.pod'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
compute T all temp
timestep 0.001
thermo 100
thermo_style custom step time atoms c_T pe ke etotal press c_fmax
thermo_modify temp T flush yes lost error
run 0
velocity all create 300.0 291201 mom yes rot no dist gaussian
run 0
fix int all nve
dump traj all custom 100 {test / 'trajectory.dump'} id type x y z vx vy vz fx fy fz
dump_modify traj sort id format float %.12g flush yes
run 1000
unfix int
undump traj
write_data {test / 'end.data'}
''')
    code, stdout = run_lmp(test, 'in.nve.lmp')
    log = test / 'log.lammps'
    text = log.read_text(errors='replace') if log.exists() else stdout
    rows, head = [], None
    for line in text.splitlines():
        f = line.split()
        if f and f[0] == 'Step' and 'TotEng' in f:
            head = f
        elif head and len(f) == len(head):
            try:
                rows.append(dict(zip(head, map(float, f))))
            except ValueError:
                pass
    zero = [r for r in rows if int(r['Step']) == 0]
    final = rows[-1] if rows else {}
    min_d = None
    if (test / 'trajectory.dump').exists():
        last = read(test / 'trajectory.dump', format='lammps-dump-text', index=-1)
        d = last.get_all_distances(mic=True)
        np.fill_diagonal(d, np.inf)
        min_d = float(d.min())
    summary = {'status': 'completed' if code == 0 and rows and int(final.get('Step', -1)) == 1000 else 'failed',
               'natoms': len(crystal), 'initial_fcc_lattice_A': 4.09, 'relaxed_cell_A': np.asarray(optimized.cell).tolist(),
               'relaxed_fcc_lattice_parameter_A': float(optimized.cell.lengths()[0] / 4.0),
               'relaxed_volume_per_atom_A3': float(optimized.get_volume() / len(optimized)),
               'initial_volume_per_atom_A3': float(crystal.get_volume() / len(crystal)),
               'relaxed_minimum_pair_distance_A': float(np.min(np.where(np.eye(len(optimized), dtype=bool), np.inf, optimized.get_all_distances(mic=True)))),
               'zero_pressure_minimize_log': str(relax / 'log.lammps'), 'nve_temperature_K': 300.0,
               'nve_dt_ps': 0.001, 'nve_steps': 1000, 'nve_duration_ps': 1.0,
               'initial_total_energy_eV': zero[-1].get('TotEng') if zero else None,
               'final_total_energy_eV': final.get('TotEng'),
               'energy_drift_eV_per_atom': ((final['TotEng'] - zero[-1]['TotEng']) / len(crystal)) if zero and final else None,
               'initial_fmax_eV_A': zero[-1].get('c_fmax') if zero else None,
               'final_fmax_eV_A': final.get('c_fmax'), 'final_temperature_K': final.get('c_T'),
               'final_minimum_pair_distance_A': min_d, 'return_code': code,
               'physical_assessment': 'FAIL: zero-pressure relaxation contracts FCC Ag to a nonphysical cell; see Stage-1 report',
               'input_sha256': sha(test / 'in.nve.lmp'), 'log_sha256': sha(log) if log.exists() else None,
               'fcc_relax_input_sha256': sha(relax / 'in.relax.lmp'), 'fcc_relaxed_data_sha256': sha(opt_data)}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    if code:
        raise SystemExit(code)


if __name__ == '__main__':
    main()
