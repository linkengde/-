#!/workspace/venvs/gpaw/bin/python
"""Short no-heat 1 fs NVE screen on independent V13 validation structures."""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import numpy as np
from ase.io import read, write

WORK = Path('/workspace/potential-fit/round7_pod_snap/models/target_200x200x300_qsgs_material200/v13_source_and_fit_20261123')
FIT = WORK / 'fit_v13'
OUT = WORK / 'small_cold_screen_v13'
LMP = Path('/workspace/venvs/podfit/bin/lmp')
SPECIES = ['Ag', 'C', 'Si', 'Ti']
CASES = [
    ('V13_VAL_AG_SOLID_01', 'validation', 300.0),
    ('V13_VAL_AG_LIQUID_01', 'validation', 1300.0),
    ('V13_TSC_VIB_01_ANISO', 'training', 300.0),
    ('V13_VAL_INTERFACE_01', 'validation', 900.0),
    ('V13_VAL_AG_RICH_SI_01', 'validation', 1100.0),
]
STEPS = 100
DT_PS = 0.001


def thermo_rows(path: Path):
    header = None
    rows = []
    for line in path.read_text(errors='replace').splitlines():
        fields = line.split()
        if fields and fields[0] == 'Step' and 'TotEng' in fields:
            header = fields
            continue
        if header and len(fields) == len(header):
            try:
                rows.append(dict(zip(header, map(float, fields))))
            except ValueError:
                continue
    return rows


def minimum_distance(atoms):
    d = atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, np.inf)
    return float(np.min(d))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results = []
    for index, (cid, role, target_temp) in enumerate(CASES, start=1):
        source = FIT / role / f'{cid}.xyz'
        atoms = read(source)
        run = OUT / cid
        run.mkdir(parents=True, exist_ok=True)
        data = run / 'start.data'
        write(data, atoms, format='lammps-data', atom_style='atomic', units='metal', specorder=SPECIES)
        dump = run / 'end.dump'
        input_file = run / 'in.nve_0p1ps_1fs.lmp'
        seed = 291104 + index
        input_file.write_text(f'''clear
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
pair_coeff * * {FIT / 'AgTiSiC_POD_v13.pod'} {FIT / 'v13_coefficients.pod'} Ag C Si Ti
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
compute T all temp
timestep {DT_PS:.6f}
thermo 1
thermo_style custom step time atoms c_T pe ke etotal c_fmax
thermo_modify temp T lost error flush yes format float %.16g
run 0
velocity all create {target_temp:.1f} {seed} mom yes rot no dist gaussian
run 0
fix integ all nve
dump traj all custom {STEPS} {dump} id type x y z vx vy vz fx fy fz
dump_modify traj sort id format float %.12g
run {STEPS}
unfix integ
undump traj
write_data {run / 'end.data'}
''')
        env = os.environ.copy()
        env.update({'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
        proc = subprocess.run([str(LMP), '-in', str(input_file), '-log', 'log.lammps', '-screen', 'none'],
                              cwd=run, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True)
        (run / 'screen.out').write_text(proc.stdout)
        if proc.returncode != 0:
            results.append({'candidate_id': cid, 'status': 'failed', 'target_temperature_K': target_temp,
                            'return_code': proc.returncode,
                            'output_tail': proc.stdout[-3000:]})
            (OUT / 'summary.json').write_text(json.dumps({'results': results}, indent=2) + '\n')
            print(json.dumps(results[-1], ensure_ascii=False), flush=True)
            continue
        rows = thermo_rows(run / 'log.lammps')
        zero = [r for r in rows if int(r['Step']) == 0]
        if len(zero) < 2 or int(rows[-1]['Step']) != STEPS:
            raise RuntimeError(f'{cid}: unexpected thermo rows or incomplete run')
        start, end = zero[-1], rows[-1]
        final = read(dump, format='lammps-dump-text', index=-1)
        types = np.asarray(final.arrays['type'], dtype=int)
        type_map = {i + 1: symbol for i, symbol in enumerate(SPECIES)}
        final.set_chemical_symbols([type_map[int(t)] for t in types])
        for key in list(final.arrays):
            if key.lower() in {'id', 'type', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz'}:
                final.arrays.pop(key, None)
        results.append({
            'candidate_id': cid, 'status': 'completed', 'target_temperature_K': target_temp,
            'natoms': len(final),
            'steps': STEPS, 'dt_ps': DT_PS, 'duration_ps': STEPS * DT_PS,
            'initial_temperature_K': start['c_T'], 'final_temperature_K': end['c_T'],
            'initial_total_energy_eV': start['TotEng'], 'final_total_energy_eV': end['TotEng'],
            'energy_drift_eV_per_atom': (end['TotEng'] - start['TotEng']) / len(final),
            'initial_fmax_eV_A': start['c_fmax'], 'final_fmax_eV_A': end['c_fmax'],
            'initial_min_distance_A': minimum_distance(atoms),
            'final_min_distance_A': minimum_distance(final),
            'atom_count_conserved': len(final) == len(atoms),
            'source_file': str(source), 'input_file': str(input_file),
        })
        (OUT / 'summary.json').write_text(json.dumps({'results': results}, indent=2) + '\n')
        print(json.dumps(results[-1], ensure_ascii=False), flush=True)
    print('SUMMARY', OUT / 'summary.json', flush=True)


if __name__ == '__main__':
    main()
