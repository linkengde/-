#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, json, os, re, shutil, subprocess, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
LMP = Path('/workspace/.local/lammps-static-30Sep2026/lammps-static/bin/lmp')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/agti-liquid-control-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/agti-liquid-control-cache')
sys.path.insert(0, str(ROOT / 'audit'))
import audit_crater as ac
import audit_connectivity as conn

MASS = {1: 47.867, 2: 28.0855, 3: 12.0107, 4: 107.87}
KB = 8.617333262e-5
M_VV2E = 1.0364269e-4


def write_input(window):
    template = (ROOT / 'inputs/control_template.lmp').read_text()
    template = re.sub(r'^create_atoms 4 single .* units box\n', '', template, flags=re.M)
    template = template.replace('group projectiles id 4306:4750', 'group projectiles id 4306:4751')
    template = template.replace('group new_projectile id 4751', 'group new_projectile id 4752')
    template = re.sub(r'^velocity new_projectile set .*\n', '', template, flags=re.M)
    template = template.replace('624146', str(623700 + window))
    template = template.replace(
        'dump tr all custom 4000 dumps/impact_446_62eV.dump id type x y z vx vy vz fx fy fz',
        f'dump tr all custom 100 dumps/control_{window}_62eV.dump id type x y z xu yu zu vx vy vz fx fy fz')
    template = template.replace('write_data structures/state_current.data',
                                f'write_data structures/control_{window}_62eV.data')
    template = template.replace('write_restart restart/state_current.restart',
                                f'write_restart restart/control_{window}_62eV.restart')
    if re.search(r'^create_atoms\b|^velocity new_projectile\b', template, flags=re.M):
        raise RuntimeError('control input still creates or launches an atom')
    path = ROOT / 'inputs' / f'control_{window}_62eV.lmp'
    path.write_text(template)
    return path


def read_last_dump(path):
    lines = path.read_text(errors='replace').splitlines()
    starts = [i for i, s in enumerate(lines) if s == 'ITEM: TIMESTEP']
    i = starts[-1]
    step, n = int(lines[i+1]), int(lines[i+3])
    j = i + 4
    while not lines[j].startswith('ITEM: ATOMS'):
        j += 1
    cols = lines[j].split()[2:]
    rows = np.array([[float(x) for x in line.split()] for line in lines[j+1:j+1+n]])
    return step, cols, rows


def one_window(window, threads=5):
    inp = write_input(window)
    name = f'control_{window}_62eV'
    log = ROOT / 'logs' / f'{name}.lammps'
    console = ROOT / 'logs' / f'{name}.console.txt'
    dump = ROOT / 'dumps' / f'{name}.dump'
    data = ROOT / 'structures' / f'{name}.data'
    restart = ROOT / 'restart' / f'{name}.restart'
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(threads)
    cmd = [str(LMP), '-sf', 'omp', '-pk', 'omp', str(threads), '-in', str(inp.relative_to(ROOT)), '-log', str(log.relative_to(ROOT))]
    with console.open('w') as f:
        p = subprocess.run(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, env=env)
    if p.returncode != 0:
        raise RuntimeError(f'LAMMPS failed in control window {window}; see {console}')
    if not all(pth.is_file() for pth in (log, dump, data, restart)):
        raise RuntimeError(f'missing output in control window {window}')
    logtext = log.read_text(errors='replace')
    if ('ERROR:' in logtext or 'Lost atoms' in logtext or re.search(r'\b(?:nan|inf)\b', logtext)
            or 'Total wall time:' not in logtext):
        raise RuntimeError(f'control log check failed in window {window}')

    atomdata = ac.read_data(data)
    if len(atomdata['ids']) != 4751:
        raise RuntimeError(f'control atom count changed in window {window}: {len(atomdata["ids"])}')
    step, cols, rows = read_last_dump(dump)
    if step != 4000 or len(rows) != 4751 or not np.isfinite(rows).all():
        raise RuntimeError(f'bad final dump in control window {window}')
    ids = rows[:, cols.index('id')].astype(int)
    typ = rows[:, cols.index('type')].astype(int)
    vel = rows[:, [cols.index(x) for x in ('vx', 'vy', 'vz')]]
    frc = rows[:, [cols.index(x) for x in ('fx', 'fy', 'fz')]]
    fn = np.linalg.norm(frc, axis=1)
    k = int(np.argmax(fn))
    xy = np.array([atomdata['lx'], 0.0, 0.0])
    xyb = np.array([atomdata['tilt'][0], atomdata['ly'], 0.0])
    connectivity = {}
    for key, mask, cutoff in (('TSC', atomdata['typ'] != 4, 3.1), ('Ag', atomdata['typ'] == 4, 3.35)):
        largest, clusters, frac = conn.largest_component(atomdata['xyz'][mask], atomdata['ids'][mask], cutoff, xy, xyb)
        connectivity[key] = {'largest': int(largest), 'total': int(mask.sum()), 'clusters': int(clusters),
                             'fraction': float(frac), 'cutoff_A': cutoff}
    xyz = rows[:, [cols.index(x) for x in ('x', 'y', 'z')]]
    # Local temperature uses the same fixed impact center and z window as the impact branch.
    dxy = xyz[:, :2] - np.array([9.9573898076, 24.2578814220])
    shifts = np.array([i * np.array([atomdata['lx'], 0.0]) + j * np.array([atomdata['tilt'][0], atomdata['ly']])
                       for i in (-1, 0, 1) for j in (-1, 0, 1)])
    img = dxy[:, None, :] + shifts[None, :, :]
    r = np.sqrt(np.min(np.sum(img * img, axis=2), axis=1))
    local = (r < 8.0) & (xyz[:, 2] >= 12.0) & (xyz[:, 2] < 55.0) & (ids <= 4305)
    def temp(mask):
        ix = np.flatnonzero(mask)
        if len(ix) < 2:
            return None, len(ix)
        m = np.array([MASS[int(typ[i])] for i in ix])
        v = vel[ix] - np.average(vel[ix], axis=0, weights=m)
        ke = 0.5 * M_VV2E * np.sum(m[:, None] * v * v)
        return float(2 * ke / ((3 * len(ix) - 3) * KB)), len(ix)
    t_ag, n_ag = temp(local & (typ == 4))
    t_tsc, n_tsc = temp(local & (typ != 4))
    result = {'window': window, 'time_ps': 100.0 + 0.2 * (window - 447 + 1), 'step': step,
              'atoms': len(ids), 'type_counts': {str(int(t)): int(np.sum(typ == t)) for t in np.unique(typ)},
              'max_force_eV_per_A': float(fn[k]), 'max_force_atom_id': int(ids[k]),
              'local_kinetic_temperature_K': {'Ag': {'N': int(n_ag), 'T_K': t_ag},
                                              'TSC': {'N': int(n_tsc), 'T_K': t_tsc}},
              'connectivity': connectivity}
    (ROOT / 'audit' / f'audit_control_{window}.json').write_text(json.dumps(result, indent=2) + '\n')
    shutil.copy2(data, ROOT / 'structures/state_current.data')
    shutil.copy2(restart, ROOT / 'restart/state_current.restart')
    shutil.copy2(data, ROOT / 'checkpoints' / f'after_control_window_{window:03d}.data')
    shutil.copy2(restart, ROOT / 'checkpoints' / f'after_control_window_{window:03d}.restart')
    gz = dump.with_suffix(dump.suffix + '.gz')
    with dump.open('rb') as src, gzip.open(gz, 'wb', compresslevel=5) as dst:
        shutil.copyfileobj(src, dst)
    dump.unlink()
    if fn[k] > 150.0 or connectivity['TSC']['fraction'] < 0.99 or connectivity['Ag']['fraction'] < 0.95:
        raise RuntimeError(f'control threshold crossed at window {window}; inspect saved checkpoint')
    row = [window, f'{result["time_ps"]:.3f}', len(ids), f'{fn[k]:.6f}',
           '' if t_ag is None else f'{t_ag:.3f}', '' if t_tsc is None else f'{t_tsc:.3f}',
           f'{connectivity["TSC"]["fraction"]:.6f}', f'{connectivity["Ag"]["fraction"]:.6f}']
    with (ROOT / 'run_progress.csv').open('a', newline='') as f:
        csv.writer(f).writerow(row)
    print(f'PASS control={window} time={row[1]}ps N={len(ids)} maxF={fn[k]:.1f} '
          f'TAg={t_ag or 0:.0f}K TTSC={t_tsc or 0:.0f}K connectivity='
          f'{connectivity["TSC"]["fraction"]:.4f}/{connectivity["Ag"]["fraction"]:.4f}', flush=True)


def main():
    if not LMP.is_file():
        raise FileNotFoundError(LMP)
    with (ROOT / 'run_progress.csv').open('w', newline='') as f:
        csv.writer(f).writerow(['window', 'end_time_ps', 'atoms', 'max_force_eV_per_A',
                                'Ag_local_T_r8_K', 'TSC_local_T_r8_K', 'TSC_largest_cluster_fraction',
                                'Ag_largest_cluster_fraction'])
    for window in range(447, 472):
        one_window(window)
    print('NO_IMPACT_CONTROL_COMPLETE windows=25 duration=5ps', flush=True)


if __name__ == '__main__':
    main()
