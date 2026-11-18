#!/workspace/venvs/dft-labels/bin/python
"""Controlled static minimization of a copy of the 347k model with POD v10."""
from __future__ import annotations

import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path('/workspace/potential-fit/round7_pod_snap')
MODEL = ROOT / 'models/final_v10'
DATA = ROOT / 'models/target_200x200x300_qsgs/md_runs/attempt_008_POD_TTM_mpi_big/start_wrapped.data'
OUT = MODEL / 'pod_minimize_347k'
LMP = Path('/workspace/software/lammps-pod-snap-clean/build/lmp')


def read_initial(path: Path):
    lines = path.read_text().splitlines()
    begin = next(i for i, line in enumerate(lines) if line.startswith('Atoms')) + 1
    rows = []
    for line in lines[begin:]:
        if not line.strip():
            if rows:
                break
            continue
        fields = line.split()
        if len(fields) < 5 or not fields[0].isdigit():
            continue
        rows.append((int(fields[0]), int(fields[1]), float(fields[2]), float(fields[3]), float(fields[4])))
    rows.sort()
    return np.asarray(rows, dtype=float)


def read_dump(path: Path):
    lines = path.read_text().splitlines()
    head = next(i for i, line in enumerate(lines) if line.startswith('ITEM: ATOMS '))
    cols = lines[head].split()[2:]
    rows = [dict(zip(cols, line.split())) for line in lines[head+1:] if line.strip() and not line.startswith('ITEM:')]
    rows.sort(key=lambda r: int(r['id']))
    return (np.asarray([int(r['id']) for r in rows]),
            np.asarray([int(r['type']) for r in rows]),
            np.asarray([[float(r[k]) for k in ('x', 'y', 'z')] for r in rows]),
            np.asarray([[float(r[k]) for k in ('fx', 'fy', 'fz')] for r in rows]))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    inp = OUT / 'in.minimize'
    dump = OUT / 'minimized.dump'
    input_lines = [
        'clear', 'units metal', 'atom_style atomic', 'boundary p p m',
        'atom_modify sort 0 0', f'read_data {DATA}',
        'group source_Ti type 1', 'group source_Si type 2',
        'group source_C type 3', 'group source_Ag type 4',
        'set group source_Ti type 4', 'set group source_Si type 3',
        'set group source_C type 2', 'set group source_Ag type 1',
        'mass 1 107.8682', 'mass 2 12.011', 'mass 3 28.085', 'mass 4 47.867',
        'pair_style pod',
        f"pair_coeff * * {MODEL/'POD/AgTiSiC_POD.pod'} {MODEL/'POD/pod_coefficients.pod'} Ag C Si Ti",
        'neighbor 2.0 bin', 'neigh_modify delay 0 every 1 check yes',
        'region fixed_z block INF INF INF INF 0 20 units box',
        'group bottom_z region fixed_z',
        'group skeleton union source_Ti source_Si source_C',
        'group fixed intersect bottom_z skeleton',
        'group mobile subtract all fixed',
        'variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)',
        'compute fmax all reduce max v_fmag',
        'compute fmax_mobile mobile reduce max v_fmag',
        'compute fmax_fixed fixed reduce max v_fmag',
        'thermo 10', 'thermo_style custom step atoms pe etotal c_fmax c_fmax_mobile c_fmax_fixed',
        'thermo_modify lost error flush yes format float %.16g',
        'run 0', 'fix hold fixed setforce 0.0 0.0 0.0',
        'min_style fire', 'min_modify dmax 0.1',
        'minimize 0.0 0.3 1000 5000',
        'unfix hold', 'run 0',
        f"print \"FINAL $(pe:%.16g) $(c_fmax:%.16g) $(c_fmax_mobile:%.16g) $(c_fmax_fixed:%.16g)\" file {OUT/'final_metrics.txt'} screen no",
        f'dump d all custom 1 {dump} id type x y z fx fy fz',
        'dump_modify d sort id', f'write_data {OUT/"minimized.data"}', '']
    inp.write_text('\n'.join(input_lines))
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '1'
    started = time.time()
    result = subprocess.run([str(LMP), '-in', str(inp), '-log', 'log.lammps'],
                            cwd=OUT, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    (OUT / 'stdout.log').write_text(result.stdout)
    if result.returncode:
        raise SystemExit(result.stdout[-5000:])
    ids, types, pos, force = read_dump(dump)
    initial = read_initial(DATA)
    if len(ids) != len(initial):
        raise ValueError('Atom count changed during minimize')
    delta = pos - initial[:, 2:5]
    delta[:, 0] -= 200 * np.rint(delta[:, 0] / 200)
    delta[:, 1] -= 200 * np.rint(delta[:, 1] / 200)
    displacement = np.linalg.norm(delta, axis=1)
    vals = [float(x) for x in (OUT / 'final_metrics.txt').read_text().split()[1:]]
    log_lines = (OUT / 'log.lammps').read_text(errors='replace').splitlines()
    summary = {
        'status': 'MINIMIZATION_COMPLETED',
        'method': 'POD v10; FIRE minimization; fixed ceramic atoms with 0<=z<20 A; no MD integration',
        'natoms_before_after': [int(len(initial)), int(len(ids))],
        'iterations_and_evaluations': next((line.strip() for line in reversed(log_lines) if 'iterations,' in line and 'force evaluations' in line), None),
        'final_pe_eV': vals[0], 'final_fmax_all_eV_A': vals[1],
        'final_fmax_mobile_eV_A': vals[2], 'final_fmax_fixed_skeleton_eV_A': vals[3],
        'max_displacement_A': float(displacement.max()),
        'rms_displacement_A': float(np.sqrt(np.mean(displacement ** 2))),
        'wall_time_s': time.time() - started,
        'output_data': str(OUT / 'minimized.data'),
        'note': 'Connectivity and shortest-neighbor distances still require explicit post-minimization analysis.'}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
