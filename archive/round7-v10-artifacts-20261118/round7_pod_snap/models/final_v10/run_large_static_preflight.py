#!/workspace/venvs/dft-labels/bin/python
"""Run LAMMPS run-0 force screens for POD/SNAP v10 on the 347k model."""
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
OUT = MODEL / 'static_347k_preflight'
LAMMPS = Path('/workspace/software/lammps-pod-snap-clean/build/lmp')
ELEMENTS = {1: 'Ag', 2: 'C', 3: 'Si', 4: 'Ti'}
MASSES = {1: 107.8682, 2: 12.011, 3: 28.085, 4: 47.867}


def input_text(name: str, folder: Path) -> str:
    if name == 'POD':
        pair = 'pair_style pod\npair_coeff * * {} {} Ag C Si Ti'.format(
            MODEL / 'POD/AgTiSiC_POD.pod', MODEL / 'POD/pod_coefficients.pod')
    elif name == 'SNAP':
        pair = 'pair_style snap\npair_coeff * * {} {} Ag C Si Ti'.format(
            MODEL / 'SNAP/AgTiSiC_SNAP_final_v10.snapcoeff',
            MODEL / 'SNAP/AgTiSiC_SNAP_final_v10.snapparam')
    else:
        raise ValueError(name)
    return f'''clear
units metal
atom_style atomic
boundary p p m
atom_modify sort 0 0
read_data {DATA}
# Source data types Ti/Si/C/Ag -> model types Ag/C/Si/Ti.
group source_Ti type 1
group source_Si type 2
group source_C type 3
group source_Ag type 4
set group source_Ti type 4
set group source_Si type 3
set group source_C type 2
set group source_Ag type 1
mass 1 {MASSES[1]}
mass 2 {MASSES[2]}
mass 3 {MASSES[3]}
mass 4 {MASSES[4]}
{pair}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo 1
thermo_style custom step atoms temp pe etotal c_fmax
thermo_modify lost error flush yes format float %.16g
dump d all custom 1 {folder / 'forces.dump'} id type x y z fx fy fz
dump_modify d sort id
run 0
undump d
'''


def read_dump(path: Path):
    lines = path.read_text().splitlines()
    head = next(i for i, line in enumerate(lines) if line.startswith('ITEM: ATOMS '))
    cols = lines[head].split()[2:]
    rows = [dict(zip(cols, line.split())) for line in lines[head+1:] if line.strip() and not line.startswith('ITEM:')]
    rows.sort(key=lambda r: int(r['id']))
    f = np.asarray([[float(r[k]) for k in ('fx', 'fy', 'fz')] for r in rows])
    pos = np.asarray([[float(r[k]) for k in ('x', 'y', 'z')] for r in rows])
    types = np.asarray([int(r['type']) for r in rows])
    mag = np.linalg.norm(f, axis=1)
    if len(rows) != 347464:
        raise ValueError(f'Expected 347464 atoms, got {len(rows)}')
    def stats(mask):
        a = mag[mask]
        return {'atoms': int(mask.sum()), 'Fmax_eV_A': float(a.max()) if len(a) else None,
                'F99_eV_A': float(np.quantile(a, .99)) if len(a) else None,
                'F95_eV_A': float(np.quantile(a, .95)) if len(a) else None,
                'counts_gt_eV_A': {str(t): int((a > t).sum()) for t in (10, 50, 100, 500)}}
    z = pos[:, 2]
    regions = {'fixed_0_20': (z >= 0) & (z < 20), 'thermostat_20_40': (z >= 20) & (z < 40),
               'Newton_40_200': (z >= 40) & (z < 200), 'sputter_200_250': (z >= 200) & (z < 250),
               'evaporation_250_300': (z >= 250) & (z < 300), 'outside_0_300': (z < 0) | (z >= 300)}
    top = int(np.argmax(mag))
    return {'natoms': len(rows), 'force_overall': stats(np.ones(len(rows), bool)),
            'by_element': {e: stats(types == t) for t, e in ELEMENTS.items()},
            'by_user_z_region': {n: stats(m) for n, m in regions.items()},
            'max_force_atom': {'rank': 1, 'id': int(top+1), 'type': int(types[top]),
                               'element': ELEMENTS[int(types[top])], 'xyz_A': pos[top].tolist(),
                               'force_eV_A': f[top].tolist(), 'F_eV_A': float(mag[top])}}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    ld = '/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:' + str(LAMMPS.parent)
    env['LD_LIBRARY_PATH'] = ld + ':' + env.get('LD_LIBRARY_PATH', '')
    env['OMP_NUM_THREADS'] = '1'
    summaries = {}
    for model in ('POD', 'SNAP'):
        folder = OUT / model
        folder.mkdir(parents=True, exist_ok=True)
        (folder / 'in.static').write_text(input_text(model, folder))
        started = time.time()
        result = subprocess.run([str(LAMMPS), '-in', str(folder / 'in.static'), '-log', 'log.lammps'],
                                cwd=folder, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True)
        (folder / 'stdout.log').write_text(result.stdout)
        if result.returncode:
            summaries[model] = {'status': 'FAILED', 'returncode': result.returncode,
                                'tail': result.stdout[-5000:]}
            continue
        parsed = read_dump(folder / 'forces.dump')
        parsed.update({'status': 'STATIC_RUN_0_COMPLETED', 'method': 'LAMMPS run 0; no MD integration',
                       'wall_time_s': time.time() - started, 'lammps': str(LAMMPS),
                       'start_data': str(DATA), 'potential': model})
        (folder / 'summary.json').write_text(json.dumps(parsed, indent=2) + '\n')
        summaries[model] = parsed
    (OUT / 'summary.json').write_text(json.dumps(summaries, indent=2) + '\n')
    print(json.dumps({k: {'status': v.get('status'), 'force_overall': v.get('force_overall'), 'wall_time_s': v.get('wall_time_s'), 'tail': v.get('tail')} for k,v in summaries.items()}, indent=2))


if __name__ == '__main__':
    main()
