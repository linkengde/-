#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, json, math, os, re, shutil, subprocess, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
LMP = Path('/workspace/.local/lammps-static-30Sep2026/lammps-static/bin/lmp')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/agti-liquid-diagnostic-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/agti-liquid-diagnostic-cache')
CENTER = np.array([9.9573898076, 24.2578814220])
KB = 8.617333262e-5
MIN_CLEARANCE_A = 2.5
sys.path.insert(0, str(ROOT / 'audit'))
import audit_crater as ac


def candidate_rows(first=472, last=476, seed=20261012):
    rng = np.random.default_rng(seed)
    R = 7.5
    lmax = math.log1p((R / 6.0) ** 2)
    for _ in range(first - 1):
        rng.random()
        rng.uniform(0.0, 2.0 * math.pi)
    for hit in range(first, last + 1):
        u = float(rng.random())
        radius = 6.0 * math.sqrt(math.exp(u * lmax) - 1.0)
        theta = float(rng.uniform(0.0, 2.0 * math.pi))
        x, y = CENTER + radius * np.array([math.cos(theta), math.sin(theta)])
        yield (hit, radius, theta, float(x), float(y), 60.0, 62.0, -105.3164, 200.0)


def launch_clearance(row):
    d = ac.read_data(ROOT / 'structures/state_current.data')
    _, _, _, x, y, z, *_ = row
    target = np.linalg.solve(d['H'], np.array([x - d['bounds']['x'][0], y - d['bounds']['y'][0]]))
    df = d['frac'] - target
    a = np.array([d['lx'], 0.0])
    b = np.array([d['tilt'][0], d['ly']])
    images = np.stack([df + np.array([i, j]) for i in (-1, 0, 1) for j in (-1, 0, 1)])
    dxy = images @ d['H'].T
    dxy2 = np.min(np.sum(dxy * dxy, axis=2), axis=0)
    dist = np.sqrt(dxy2 + (d['xyz'][:, 2] - z) ** 2)
    k = int(np.argmin(dist))
    return float(dist[k]), int(d['ids'][k]), int(d['typ'][k])


def safe_launch(row):
    hit, radius, theta, x, y, z, energy, vz, spacing = row
    clearance, near_id, near_type = launch_clearance(row)
    if clearance >= MIN_CLEARANCE_A:
        return row, clearance, 0, near_id, near_type
    rng = np.random.default_rng(np.random.SeedSequence([20261012, hit, 3232026]))
    R = 7.5
    lmax = math.log1p((R / 6.0) ** 2)
    for attempt in range(1, 10001):
        u = float(rng.random())
        radius = 6.0 * math.sqrt(math.exp(u * lmax) - 1.0)
        theta = float(rng.uniform(0.0, 2.0 * math.pi))
        x, y = CENTER + radius * np.array([math.cos(theta), math.sin(theta)])
        trial = (hit, radius, theta, float(x), float(y), z, energy, vz, spacing)
        clearance, near_id, near_type = launch_clearance(trial)
        if clearance >= MIN_CLEARANCE_A:
            return trial, clearance, attempt, near_id, near_type
    raise RuntimeError(f'could not find safe launch position for hit {hit}')


def write_input(hit, row):
    _, _, _, x, y, z, _, _, _ = row
    path = ROOT / 'inputs' / f'impact_{hit}_62eV.lmp'
    s = (ROOT / 'inputs/impact_446_template.lmp').read_text()
    repl = {
        'group projectiles id 4306:4750': f'group projectiles id 4306:{4304 + hit}',
        'group new_projectile id 4751': f'group new_projectile id {4305 + hit}',
        '624146': str(623700 + hit),
        'dump tr all custom 4000 dumps/impact_446_62eV.dump id type x y z vx vy vz fx fy fz':
            f'dump tr all custom 100 dumps/impact_{hit}_62eV.dump id type x y z xu yu zu vx vy vz fx fy fz',
        'write_data structures/state_current.data': f'write_data structures/impact_{hit}_62eV.data',
        'write_restart restart/state_current.restart': f'write_restart restart/state_after_hit_{hit:03d}.restart',
    }
    for old, new in repl.items():
        if old not in s:
            raise RuntimeError(f'template token missing: {old}')
        s = s.replace(old, new)
    s = re.sub(r'^create_atoms 4 single .* units box$',
               f'create_atoms 4 single {x:.10f} {y:.10f} {z:.1f} units box', s, flags=re.M)
    if 'impact_446_62eV.dump' in s or 'write_data structures/state_current.data' in s:
        raise RuntimeError('old output name remains in generated input')
    path.write_text(s)
    return path


def sampled_dump_peak(path):
    peak = {'force_eV_per_A': -1.0, 'atom_id': None, 'type': None, 'step': None, 'frames': 0, 'atoms_per_frame': None}
    with path.open() as f:
        while True:
            line = f.readline()
            if not line: break
            if line.strip() != 'ITEM: TIMESTEP': continue
            step = int(f.readline()); f.readline(); n = int(f.readline()); f.readline()
            for _ in range(3): f.readline()
            cols = f.readline().split()[2:]
            ix = {name: cols.index(name) for name in ('id','type','fx','fy','fz')}
            peak['frames'] += 1; peak['atoms_per_frame'] = n
            for _ in range(n):
                row = f.readline().split()
                fx,fy,fz = (float(row[ix[k]]) for k in ('fx','fy','fz'))
                fn = math.sqrt(fx*fx+fy*fy+fz*fz)
                if not math.isfinite(fn): raise RuntimeError(f'non-finite force in {path} at step {step}')
                if fn > peak['force_eV_per_A']:
                    peak.update(force_eV_per_A=fn, atom_id=int(row[ix['id']]), type=int(row[ix['type']]), step=step)
    return peak


def run_hit(hit, row, threads=5):
    row, clearance, resamples, near_id, near_type = safe_launch(row)
    _, radius, theta, x, y, z, energy, vz, spacing = row
    with (ROOT / 'impact_positions_472_476.csv').open('a', newline='') as f:
        csv.writer(f).writerow([hit, f'{radius:.10f}', f'{theta:.12f}', f'{x:.10f}', f'{y:.10f}',
                                f'{z:.3f}', energy, f'{vz:.4f}', spacing, f'{clearance:.6f}',
                                resamples, near_id, near_type])
        f.flush()
        os.fsync(f.fileno())

    inp = write_input(hit, row)
    name = f'impact_{hit}_62eV'
    log = ROOT / 'logs' / f'{name}.lammps'
    console = ROOT / 'logs' / f'{name}.console.txt'
    dump = ROOT / 'dumps' / f'{name}.dump'
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(threads)
    env['MPLCONFIGDIR'] = '/tmp/agti-liquid-diagnostic-mpl'
    env['XDG_CACHE_HOME'] = '/tmp/agti-liquid-diagnostic-cache'
    data = ROOT / 'structures' / f'{name}.data'
    restart = ROOT / 'restart' / f'state_after_hit_{hit:03d}.restart'
    existing = all(p.is_file() for p in (log, dump, data, restart))
    existing_log = log.read_text(errors='replace') if existing else ''
    complete_existing = ('Total wall time:' in existing_log and 'ERROR:' not in existing_log
                         and 'Lost atoms' not in existing_log
                         and not re.search(r'\b(?:nan|inf)\b', existing_log))
    if not complete_existing:
        cmd = [str(LMP), '-sf', 'omp', '-pk', 'omp', str(threads), '-in', str(inp.relative_to(ROOT)), '-log', str(log.relative_to(ROOT))]
        with console.open('w') as f:
            proc = subprocess.run(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, env=env)
        if proc.returncode != 0:
            raise RuntimeError(f'LAMMPS returned {proc.returncode} at hit {hit}; see {console}')
    if not all(p.is_file() for p in (log, dump, data, restart)):
        raise RuntimeError(f'missing output at hit {hit}')
    logtext = log.read_text(errors='replace')
    if ('ERROR:' in logtext or 'Lost atoms' in logtext
            or re.search(r'\b(?:nan|inf)\b', logtext)
            or 'Total wall time:' not in logtext):
        raise RuntimeError(f'LAMMPS log check failed at hit {hit}; inspect log before continuing')

    audit_txt = ROOT / 'audit' / f'audit_hit_{hit}.txt'
    with audit_txt.open('w') as f:
        proc = subprocess.run([sys.executable, 'audit/audit_stage.py', str(hit)], cwd=ROOT,
                              stdout=f, stderr=subprocess.STDOUT, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f'audit failed at hit {hit}; see {audit_txt}')
    audit_path = ROOT / 'audit' / f'audit_hit_{hit:02d}.json'
    result = json.loads(audit_path.read_text())
    sampled_peak = sampled_dump_peak(dump)
    if sampled_peak['frames'] != 41 or sampled_peak['atoms_per_frame'] != result['atoms']:
        raise RuntimeError(f'incomplete high-frequency dump at hit {hit}: {sampled_peak}')
    sampled_peak['cumulative_time_ps'] = 100.0 + 0.2 * (hit - 446) + sampled_peak['step'] * 0.00005
    result['sampled_dump_max_force'] = sampled_peak
    audit_path.write_text(json.dumps(result, indent=2) + '\n')
    if result['atoms'] != 4305 + hit:
        raise RuntimeError(f'atom count mismatch at hit {hit}: {result["atoms"]}')

    # Preserve each completed state before allowing the next projectile.
    shutil.copy2(data, ROOT / 'structures/state_current.data')
    shutil.copy2(restart, ROOT / 'restart/state_current.restart')
    shutil.copy2(data, ROOT / 'checkpoints' / f'after_hit_{hit:03d}.data')
    shutil.copy2(restart, ROOT / 'checkpoints' / f'after_hit_{hit:03d}.restart')

    fmax = result['max_force_eV_per_A']
    tsc = result['connectivity']['TSC']['fraction']
    ag = result['connectivity']['Ag']['fraction']
    rowout = [hit, f'{100.0 + 0.2 * (hit - 446):.3f}', result['atoms'], f'{fmax:.6f}',
              f'{sampled_peak["force_eV_per_A"]:.6f}', f'{sampled_peak["cumulative_time_ps"]:.5f}', sampled_peak['atom_id'],
              f'{result["local_kinetic_temperature_K"]["r_lt_8"]["Ag_matrix"]["T_K"] or 0:.3f}',
              f'{result["local_kinetic_temperature_K"]["r_lt_8"]["TSC"]["T_K"] or 0:.3f}',
              f'{tsc:.6f}', f'{ag:.6f}']
    with (ROOT / 'run_progress.csv').open('a', newline='') as f:
        csv.writer(f).writerow(rowout)
        f.flush()
        os.fsync(f.fileno())

    # Compression keeps the 5 ps, 5 fs-resolution trajectory practical to archive.
    gz = dump.with_suffix(dump.suffix + '.gz')
    with dump.open('rb') as src, gzip.open(gz, 'wb', compresslevel=5) as dst:
        shutil.copyfileobj(src, dst)
    dump.unlink()
    # High sampled forces are recorded for review; in this model, an Ag projectile–target collision activates the ZBL short-range term.
    if tsc < 0.99 or ag < 0.95:
        raise RuntimeError(f'connectivity threshold crossed at hit {hit}: TSC={tsc:.5f}, Ag={ag:.5f}')
    print(f'PASS hit={hit} cumulative={rowout[1]}ps N={result["atoms"]} endpointF={fmax:.1f} sampledF={sampled_peak["force_eV_per_A"]:.1f} '
          f'TAg={rowout[7]}K TTSC={rowout[8]}K connectivity={tsc:.4f}/{ag:.4f} '
          f'launch_clearance={clearance:.2f}A frames={sum(1 for _ in gzip.open(gz,"rt") if _.startswith("ITEM: TIMESTEP"))}',
          flush=True)


def main():
    if not LMP.is_file():
        raise FileNotFoundError(LMP)
    if not (ROOT / 'structures/state_current.data').is_file():
        raise FileNotFoundError('missing state_current.data')
    posfile = ROOT / 'impact_positions_472_476.csv'
    with posfile.open('w', newline='') as f:
        csv.writer(f).writerow(['hit', 'radius_A', 'theta_rad', 'x_A', 'y_A', 'launch_z_A', 'energy_eV',
                                'Ag_velocity_z_A_per_ps', 'spacing_fs', 'min_initial_distance_A',
                                'resamples', 'nearest_initial_atom_id', 'nearest_initial_atom_type'])
    with (ROOT / 'run_progress.csv').open('w', newline='') as f:
        csv.writer(f).writerow(['hit', 'cumulative_time_ps', 'atoms', 'endpoint_max_force_eV_per_A',
                                'sampled_max_force_eV_per_A', 'sampled_max_force_time_ps', 'sampled_max_force_atom_id',
                                'Ag_local_T_r8_K', 'TSC_local_T_r8_K', 'TSC_largest_cluster_fraction',
                                'Ag_largest_cluster_fraction'])
    for row in candidate_rows(first=472, last=476):
        run_hit(row[0], row)
    print('IMPACT_EXTENSION_COMPLETE hits=472-476 added=1.0ps cumulative=106.0ps', flush=True)


if __name__ == '__main__':
    main()
