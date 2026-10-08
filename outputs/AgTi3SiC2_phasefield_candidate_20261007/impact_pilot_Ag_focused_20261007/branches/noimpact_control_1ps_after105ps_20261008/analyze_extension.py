#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, json, os
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.special import sph_harm_y
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DUMP_PREFIX = os.environ.get('LIQUID_DUMP_PREFIX', 'impact')
CENTER = np.array([9.9573898076, 24.2578814220])
KB = 8.617333262e-5
MASS = {1: 47.867, 2: 28.0855, 3: 12.0107, 4: 107.87}
M_VV2E = 1.0364269e-4
SAMPLE_EVERY_FRAMES = 10  # 50 fs between local-order samples
LAG_STEPS = [1, 2, 4, 10, 20, 40, 100, 200, 400, 600, 800, 1000]


def read_dump(path):
    frames = []
    with gzip.open(path, 'rt', errors='replace') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.strip() != 'ITEM: TIMESTEP':
                continue
            step = int(f.readline())
            if f.readline().strip() != 'ITEM: NUMBER OF ATOMS':
                raise RuntimeError(f'bad dump header in {path}')
            n = int(f.readline())
            boxline = f.readline().split()
            if boxline[:2] != ['ITEM:', 'BOX']:
                raise RuntimeError(f'missing box in {path}')
            box_style = boxline[3:]
            b1 = list(map(float, f.readline().split()))
            b2 = list(map(float, f.readline().split()))
            b3 = list(map(float, f.readline().split()))
            atomline = f.readline().split()
            if atomline[:2] != ['ITEM:', 'ATOMS']:
                raise RuntimeError(f'missing ATOMS in {path}')
            cols = atomline[2:]
            rows = np.empty((n, len(cols)), dtype=np.float64)
            for i in range(n):
                vals = np.fromstring(f.readline(), sep=' ')
                if len(vals) != len(cols):
                    raise RuntimeError(f'bad atom row {i} in {path}')
                rows[i] = vals
            # Restricted-triclinic dump bounds are bounding-box values; recover true bounds.
            if box_style[:3] == ['xy', 'xz', 'yz']:
                xy, xz, yz = b1[2], b2[2], b3[2]
                xlo = b1[0] - min(0.0, xy, xz, xy + xz)
                xhi = b1[1] - max(0.0, xy, xz, xy + xz)
                ylo = b2[0] - min(0.0, yz)
                yhi = b2[1] - max(0.0, yz)
                lx, ly = xhi - xlo, yhi - ylo
            else:
                raise RuntimeError(f'expected triclinic box in {path}')
            H = np.array([[lx, xy], [0.0, ly]])
            frames.append({'step': step, 'cols': cols, 'rows': rows,
                           'bounds': (xlo, ylo, b3[0], xhi, yhi, b3[1]),
                           'H': H, 'a': np.array([lx, 0.0]), 'b': np.array([xy, ly])})
    return frames


def column(frame, name):
    return frame['rows'][:, frame['cols'].index(name)]


def min_image_xy(xy, frame):
    delta = xy - CENTER
    shifts = np.array([i * frame['a'] + j * frame['b']
                       for i in (-1, 0, 1) for j in (-1, 0, 1)])
    cand = delta[:, None, :] + shifts[None, :, :]
    d2 = np.sum(cand * cand, axis=2)
    return cand[np.arange(len(xy)), np.argmin(d2, axis=1)]


def local_temperature(types, ids, velocity, mask):
    ix = np.flatnonzero(mask)
    if len(ix) < 2:
        return float('nan'), int(len(ix))
    masses = np.array([MASS[int(t)] for t in types[ix]])
    v = velocity[ix] - np.average(velocity[ix], axis=0, weights=masses)
    ke = 0.5 * M_VV2E * np.sum(masses[:, None] * v * v)
    return float(2.0 * ke / ((3 * len(ix) - 3) * KB)), int(len(ix))


def local_ag_order(frame, sample_time):
    rows = frame['rows']
    ids = column(frame, 'id').astype(np.int64)
    typ = column(frame, 'type').astype(np.int64)
    xyz = rows[:, [frame['cols'].index(c) for c in ('x', 'y', 'z')]]
    ag = np.flatnonzero((typ == 4) & (xyz[:, 2] >= 12.0) & (xyz[:, 2] < 55.0))
    dxy = min_image_xy(xyz[ag, :2], frame)
    rho = np.linalg.norm(dxy, axis=1)
    centers = ag[rho < 8.0]
    if len(centers) == 0:
        return {'time_ps': sample_time, 'n_center': 0, 'mean_cn35': np.nan, 'mean_q6': np.nan,
                'median_q6': np.nan, 'q6_n': 0, 'nn_distances': np.array([])}

    agxyz = xyz[ag]
    shifts = np.array([i * np.array([frame['a'][0], frame['a'][1], 0.0]) +
                        j * np.array([frame['b'][0], frame['b'][1], 0.0])
                        for i in (-1, 0, 1) for j in (-1, 0, 1)])
    points = np.concatenate([agxyz + shift for shift in shifts], axis=0)
    point_ids = np.tile(ids[ag], len(shifts))
    tree = cKDTree(points)
    cn = []
    q6s = []
    nn = []
    for ci in centers:
        candidates = tree.query_ball_point(xyz[ci], 3.6)
        if not candidates:
            continue
        cand = np.asarray(candidates, dtype=np.int64)
        keep = point_ids[cand] != ids[ci]
        cand = cand[keep]
        vec = points[cand] - xyz[ci]
        dist = np.linalg.norm(vec, axis=1)
        valid = (dist > 1.0e-8) & (dist <= 3.6)
        vec, dist = vec[valid], dist[valid]
        cn.append(int(np.sum(dist <= 3.5)))
        nn.extend(dist[dist <= 4.2].tolist())
        if len(dist) >= 5:
            theta = np.arccos(np.clip(vec[:, 2] / dist, -1.0, 1.0))
            phi = np.arctan2(vec[:, 1], vec[:, 0])
            ylm = np.array([sph_harm_y(6, m, theta, phi) for m in range(-6, 7)])
            q6s.append(float(np.sqrt(4.0 * np.pi / 13.0 * np.sum(np.abs(np.mean(ylm, axis=1)) ** 2))))
    return {'time_ps': sample_time, 'n_center': int(len(centers)),
            'mean_cn35': float(np.mean(cn)) if cn else np.nan,
            'mean_q6': float(np.mean(q6s)) if q6s else np.nan,
            'median_q6': float(np.median(q6s)) if q6s else np.nan,
            'q6_n': int(len(q6s)), 'nn_distances': np.asarray(nn, dtype=float)}


def unwrap_increment(prev, current, frame):
    raw = current[:, :2] - prev[:, :2]
    shifts = np.array([i * frame['a'] + j * frame['b']
                       for i in (-1, 0, 1) for j in (-1, 0, 1)])
    cand = raw[:, None, :] + shifts[None, :, :]
    d2 = np.sum(cand * cand, axis=2)
    xy = cand[np.arange(len(raw)), np.argmin(d2, axis=1)]
    return np.column_stack((xy, current[:, 2] - prev[:, 2]))


def main():
    frame_times = []
    frame_hits = []
    order_rows = []
    nn_all = []
    local_temps = []
    frame_force_peaks = []
    tracer_ids = None
    tracer_wrapped = []
    tracer_unwrapped = []
    previous_tracer_wrapped = None
    cumulative = None
    boundary_jumps = []
    expected_index = 0
    previous_segment_last = None

    for hit in range(472, 477):
        path = ROOT / 'dumps' / f'{DUMP_PREFIX}_{hit}_62eV.dump.gz'
        frames = read_dump(path)
        if len(frames) != 41 or frames[0]['step'] != 0 or frames[-1]['step'] != 4000:
            raise RuntimeError(f'{path.name}: expected 41 frames from step 0 to 4000, got {len(frames)}')
        for j, frame in enumerate(frames):
            if hit > 472 and j == 0:
                # Boundary duplicate at t=end(previous hit)=start(next hit); retain previous endpoint.
                old_last = previous_segment_last
                ids0 = column(old_last, 'id').astype(np.int64)
                ids1 = column(frame, 'id').astype(np.int64)
                common, ix0, ix1 = np.intersect1d(ids0, ids1, assume_unique=True, return_indices=True)
                p0 = old_last['rows'][ix0][:, [old_last['cols'].index(c) for c in ('xu', 'yu', 'zu')]]
                p1 = frame['rows'][ix1][:, [frame['cols'].index(c) for c in ('xu', 'yu', 'zu')]]
                jump = np.linalg.norm(p1 - p0, axis=1)
                boundary_jumps.append({'between_hits': f'{hit-1}-{hit}', 'common_atoms': int(len(common)),
                                       'median_A': float(np.median(jump)), 'max_A': float(np.max(jump))})
                continue
            time_ps = 105.0 + 0.2 * (hit - 472) + frame['step'] * 0.00005
            frame_times.append(time_ps)
            frame_hits.append(hit)

            ids = column(frame, 'id').astype(np.int64)
            typ = column(frame, 'type').astype(np.int64)
            xyz = frame['rows'][:, [frame['cols'].index(c) for c in ('x', 'y', 'z')]]
            vel = frame['rows'][:, [frame['cols'].index(c) for c in ('vx', 'vy', 'vz')]]
            force = frame['rows'][:, [frame['cols'].index(c) for c in ('fx', 'fy', 'fz')]]
            fn = np.linalg.norm(force, axis=1)
            kf = int(np.argmax(fn))
            frame_force_peaks.append((time_ps, float(fn[kf]), int(ids[kf]), int(typ[kf])))
            if tracer_ids is None:
                rho = np.linalg.norm(min_image_xy(xyz[:, :2], frame), axis=1)
                mask = (typ == 4) & (ids <= 4776) & (rho < 8.0) & (xyz[:, 2] >= 12.0) & (xyz[:, 2] < 55.0)
                tracer_ids = ids[mask]
                if len(tracer_ids) < 10:
                    raise RuntimeError(f'only {len(tracer_ids)} starting Ag tracers in diagnostic region')
            loc = np.searchsorted(ids, tracer_ids)
            if np.any(loc >= len(ids)) or not np.array_equal(ids[loc], tracer_ids):
                raise RuntimeError('a tracked Ag atom is missing from a dump frame')
            wrapped = xyz[loc]
            tracer_wrapped.append(wrapped.copy())
            if previous_tracer_wrapped is None:
                cumulative = wrapped.copy()
            else:
                cumulative = cumulative + unwrap_increment(previous_tracer_wrapped, wrapped, frame)
            tracer_unwrapped.append(cumulative.copy())
            previous_tracer_wrapped = wrapped.copy()

            if expected_index % SAMPLE_EVERY_FRAMES == 0:
                order = local_ag_order(frame, time_ps)
                nn_all.extend(order.pop('nn_distances').tolist())
                order_rows.append(order)

            rho = np.linalg.norm(min_image_xy(xyz[:, :2], frame), axis=1)
            zwin = (xyz[:, 2] >= 12.0) & (xyz[:, 2] < 55.0) & (rho < 8.0)
            agt, nag = local_temperature(typ, ids, vel, zwin & (typ == 4) & (ids <= 4305))
            tsct, ntsc = local_temperature(typ, ids, vel, zwin & (typ != 4) & (ids <= 4305))
            local_temps.append((time_ps, agt, nag, tsct, ntsc))
            expected_index += 1
        previous_segment_last = frames[-1]

    times = np.asarray(frame_times)
    if len(times) != 201 or not np.allclose(np.diff(times), 0.005, atol=1e-9):
        raise RuntimeError(f'expected 201 uniform frames at 5 fs, got {len(times)}')
    pos = np.asarray(tracer_unwrapped)
    # Remove collective drift of the tagged Ag cohort before calculating self-MSD.
    pos_rel = pos - np.mean(pos, axis=1, keepdims=True)
    enddisp = np.linalg.norm(pos_rel[-1] - pos_rel[0], axis=1)
    last_frame = previous_segment_last
    last_ids = column(last_frame, 'id').astype(np.int64)
    lastloc = np.searchsorted(last_ids, tracer_ids)
    lastxyz = last_frame['rows'][lastloc][:, [last_frame['cols'].index(c) for c in ('x', 'y', 'z')]]
    last_rho = np.linalg.norm(min_image_xy(lastxyz[:, :2], last_frame), axis=1)
    kept_local = (last_rho < 8.0) & (lastxyz[:, 2] >= 12.0) & (lastxyz[:, 2] < 55.0)
    msd_rows = []
    for k in LAG_STEPS:
        if k >= len(pos_rel):
            continue
        delta = pos_rel[k:] - pos_rel[:-k]
        msd = np.mean(np.sum(delta * delta, axis=2))
        msd_rows.append((k * 0.005, float(msd)))

    with (ROOT / 'analysis_local_order.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['time_ps', 'n_center', 'mean_cn35', 'mean_q6', 'median_q6', 'q6_n'])
        w.writeheader()
        w.writerows(order_rows)
    with (ROOT / 'analysis_msd.csv').open('w', newline='') as f:
        w = csv.writer(f); w.writerow(['lag_ps', 'drift_corrected_MSD_A2']); w.writerows(msd_rows)
    with (ROOT / 'analysis_local_temperature.csv').open('w', newline='') as f:
        w = csv.writer(f); w.writerow(['time_ps', 'Ag_matrix_T_K', 'Ag_matrix_N', 'TSC_T_K', 'TSC_N']); w.writerows(local_temps)
    bins = np.arange(0.0, 4.5001, 0.05)
    hist, edges = np.histogram(np.asarray(nn_all), bins=bins)
    with (ROOT / 'analysis_Ag_neighbor_hist.csv').open('w', newline='') as f:
        w = csv.writer(f); w.writerow(['r_lo_A', 'r_hi_A', 'pair_count'])
        for i, n in enumerate(hist): w.writerow([edges[i], edges[i+1], int(n)])

    order_q = np.array([x['mean_q6'] for x in order_rows])
    order_cn = np.array([x['mean_cn35'] for x in order_rows])
    lag = np.array([x[0] for x in msd_rows]); msdval = np.array([x[1] for x in msd_rows])
    fit = (lag >= 0.2) & (lag <= 1.0)
    slope = float(np.polyfit(lag[fit], msdval[fit], 1)[0]) if fit.sum() >= 2 else float('nan')
    temp_array = np.asarray(local_temps, dtype=float)
    force_peak = max(frame_force_peaks, key=lambda row: row[1])
    summary = {
        'start_time_ps': 105.0, 'end_time_ps': 106.0, 'frames': int(len(times)),
        'frame_spacing_fs': 5.0, 'starting_Ag_tracers': int(len(tracer_ids)),
        'tracer_ids': tracer_ids.astype(int).tolist(),
        'tracer_end_to_end_displacement_A_quantiles': {
            str(q): float(np.percentile(enddisp, q)) for q in (0, 25, 50, 75, 90, 95, 99, 100)},
        'tracer_fraction_moved_over_A': {str(a): float(np.mean(enddisp > a)) for a in (1, 3, 5, 10)},
        'fraction_of_initial_tracers_still_in_r8_z12_55_region': float(np.mean(kept_local)),
        'msd_drift_corrected_A2': {f'{t:.3f}_ps': v for t, v in msd_rows},
        'linear_fit_msd_slope_0p2_to_1ps_A2_per_ps': slope,
        'apparent_D_slope_over_6_A2_per_ps': slope / 6.0,
        'local_matrix_Ag_temperature_K_median_p10_p90': [float(x) for x in np.nanpercentile(temp_array[:, 1], [50, 10, 90])],
        'local_TSC_temperature_K_median_p10_p90': [float(x) for x in np.nanpercentile(temp_array[:, 3], [50, 10, 90])],
        'maximum_force_in_5fs_dump_frames': {'time_ps': force_peak[0], 'force_eV_per_A': force_peak[1],
                                             'atom_id': force_peak[2], 'type': force_peak[3]},
        'local_order_samples': int(len(order_rows)),
        'mean_q6_range': [float(np.nanmin(order_q)), float(np.nanmax(order_q))],
        'mean_q6_first_last': [float(order_q[0]), float(order_q[-1])],
        'mean_coordination_r_le_3p5_range': [float(np.nanmin(order_cn)), float(np.nanmax(order_cn))],
        'mean_coordination_first_last': [float(order_cn[0]), float(order_cn[-1])],
        'boundary_unwrapped_jumps': boundary_jumps,
        'interpretation_limits': [
            'Local Q6 and coordination are structural indicators, not standalone phase classifiers.',
            'The impact-driven trajectory is non-equilibrium; MSD slope is an apparent mobility measure, not an equilibrium self-diffusion coefficient.',
            'The Ag-Ti3SiC2 cross potential remains unvalidated.'
        ]
    }
    (ROOT / 'analysis_summary.json').write_text(json.dumps(summary, indent=2) + '\n')

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), dpi=170, sharex=False, constrained_layout=True)
    ot = np.array([x['time_ps'] for x in order_rows])
    axes[0].plot(ot, order_q, color='#6841a5', lw=1.7)
    axes[0].set_ylabel('mean Ag Q₆')
    label = 'Ag impact-zone dynamics' if DUMP_PREFIX == 'impact' else 'Ag no-impact control'
    axes[0].set_title(f'{label}, cumulative time window 105–106 ps')
    axes[0].grid(alpha=.25)
    axes[1].plot(ot, order_cn, color='#c18c22', lw=1.7)
    axes[1].set_ylabel('mean Ag neighbors\nwithin 3.5 Å')
    axes[1].grid(alpha=.25)
    axes[2].plot(lag, msdval, 'o-', color='#217b8d', lw=1.7, ms=4)
    axes[2].set_xlabel('time lag (ps)')
    axes[2].set_ylabel('drift-corrected MSD (Å²)')
    axes[2].set_xscale('log')
    axes[2].set_yscale('log')
    axes[2].grid(alpha=.25, which='both')
    (ROOT / 'images').mkdir(exist_ok=True)
    fig.savefig(ROOT / 'images/liquid_diagnostic_metrics_5ps.png', bbox_inches='tight')
    plt.close(fig)
    print(json.dumps(summary, indent=2))
    print(f'plot: {ROOT / "images/liquid_diagnostic_metrics_5ps.png"}')


if __name__ == '__main__':
    main()
