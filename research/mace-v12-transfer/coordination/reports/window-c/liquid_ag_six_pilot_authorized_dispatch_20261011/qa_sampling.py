#!/usr/bin/env python3
"""Read-only diagnostics of an owned, genuinely generated liquid-Ag trajectory.

No calculator or dynamics is constructed. A provisional geometry is never an
accepted liquid reference: STRUCTURE_PASS remains false until a separate C
signature. --help needs only Python's standard library.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

# This shared policy imports only standard-library modules, not ASE/calculators.
from run_sampling import RESOURCE_PARAMETER_KEYS


class Hold(RuntimeError):
    pass


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def geometry_identity(atoms, np):
    ids = atoms.arrays.get('atom_ids', atoms.info.get('atom_ids'))
    if ids is None:
        raise Hold('ATOM_IDS_MISSING')
    return {'numbers': [int(x) for x in atoms.numbers],
            'positions': np.round(atoms.positions, 12).tolist(),
            'cell': np.round(atoms.cell.array, 12).tolist(),
            'pbc': [bool(x) for x in atoms.pbc],
            'atom_ids': [int(x) for x in ids]}


def save_json(path, data):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(data, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def save_csv(path, rows, fields):
    with Path(path).open('x', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def pair_vectors(positions, edge, np):
    vectors = positions[None, :, :] - positions[:, None, :]
    vectors -= np.rint(vectors / edge) * edge
    distances = np.linalg.norm(vectors, axis=2)
    return vectors, distances


def liquid_rdf(positions, edge, dr, np):
    """Finite-N normalized all-pair RDF, restricted to half the cubic edge."""
    bins = np.arange(0, edge / 2 + dr * 0.1, dr)
    if len(bins) < 10:
        raise Hold('RDF_GRID_TOO_SHORT')
    histogram = np.zeros(len(bins) - 1)
    upper = np.triu_indices(positions.shape[1], 1)
    minimum = math.inf
    for frame in positions:
        _, distances = pair_vectors(frame, edge, np)
        pair = distances[upper]
        minimum = min(minimum, float(pair.min()))
        histogram += np.histogram(pair, bins=bins)[0]
    shell = 4 * math.pi / 3 * (bins[1:]**3 - bins[:-1]**3)
    pairs = positions.shape[1] * (positions.shape[1] - 1) / 2
    expected = len(positions) * pairs * shell / edge**3
    values = histogram / expected
    radii = (bins[1:] + bins[:-1]) / 2
    smoothed = np.convolve(values, np.ones(5) / 5, mode='same')
    peak_allowed = (radii > 1.8) & (radii < min(3.6, edge / 2 - 0.5))
    peak_indices = np.flatnonzero(peak_allowed)
    if not len(peak_indices):
        raise Hold('RDF_FIRST_PEAK_NOT_RESOLVED')
    peak = int(peak_indices[np.argmax(smoothed[peak_indices])])
    candidates = [i for i in range(peak + max(3, int(0.35 / dr)), len(radii) - 3)
                  if smoothed[i - 1] > smoothed[i] <= smoothed[i + 1]]
    cutoff = float(radii[candidates[0]]) if candidates else None
    return radii, values, cutoff, minimum


def local_order(positions, edge, cutoff, np, sph_harm_y, dot_threshold=0.7,
                min_coherent_bonds=7):
    """Local l=6 qlm; coherent bonds use Re(q_i dot conj(q_j)) > 0.7.

    sph_harm_y(l, m, theta_polar, phi_azimuth) follows SciPy's current API.
    Connected crystal-like clusters require at least seven coherent bonds.
    """
    vectors, distances = pair_vectors(positions, edge, np)
    neighbor = (distances > 1e-12) & (distances < cutoff)
    counts = neighbor.sum(axis=1)
    qlm = np.zeros((len(positions), 13), dtype=complex)
    for atom in range(len(positions)):
        bonded = vectors[atom, neighbor[atom]]
        if not len(bonded):
            continue
        radius = np.linalg.norm(bonded, axis=1)
        theta = np.arccos(np.clip(bonded[:, 2] / radius, -1, 1))
        phi = np.arctan2(bonded[:, 1], bonded[:, 0])
        for column, m in enumerate(range(-6, 7)):
            qlm[atom, column] = sph_harm_y(6, m, theta, phi).mean()
    norm = np.linalg.norm(qlm, axis=1)
    unit = np.zeros_like(qlm)
    nonzero = norm > 1e-14
    unit[nonzero] = qlm[nonzero] / norm[nonzero, None]
    correlations = (unit @ unit.conj().T).real
    coherent = neighbor & (correlations > dot_threshold)
    crystal = coherent.sum(axis=1) >= min_coherent_bonds
    adjacency = coherent & crystal[:, None] & crystal[None, :]
    seen, largest = set(), 0
    for atom in np.flatnonzero(crystal):
        atom = int(atom)
        if atom in seen:
            continue
        queue, size = [atom], 0
        seen.add(atom)
        while queue:
            point = queue.pop()
            size += 1
            for next_point in np.flatnonzero(adjacency[point]):
                next_point = int(next_point)
                if next_point not in seen:
                    seen.add(next_point)
                    queue.append(next_point)
        largest = max(largest, size)
    q6 = np.sqrt(4 * math.pi / 13) * norm
    return counts, q6, int(crystal.sum()), largest


def multiple_origin_msd(positions, times, np):
    """Unwrapped coordinates, instantaneous COM removal, all time origins."""
    coordinates = positions - positions.mean(axis=1, keepdims=True)
    intervals = np.diff(times)
    if len(intervals) < 3 or not np.allclose(intervals, intervals[0], atol=1e-8):
        raise Hold('MSD_REQUIRES_REGULAR_TIME_GRID')
    maximum = len(coordinates) // 2
    rows = []
    for lag in range(1, maximum + 1):
        delta = coordinates[lag:] - coordinates[:-lag]
        values = (delta * delta).sum(axis=2).mean(axis=1)
        rows.append({'lag_fs': float(lag * intervals[0]),
                     'msd_A2': float(values.mean()),
                     'time_origins': int(len(values)),
                     'origin_std_A2': float(values.std())})
    return rows


def fit_msd_windows(rows, np):
    windows = []
    for fraction_start, fraction_end in ((0.1, 0.4), (0.25, 0.65), (0.5, 1.0)):
        start = int(len(rows) * fraction_start)
        end = int(len(rows) * fraction_end)
        selected = rows[start:max(start + 3, end)]
        if len(selected) < 3:
            continue
        x = np.asarray([row['lag_fs'] for row in selected])
        y = np.asarray([row['msd_A2'] for row in selected])
        slope, intercept = np.polyfit(x, y, 1)
        prediction = slope * x + intercept
        denominator = ((y - y.mean())**2).sum()
        r2 = float(1 - ((y - prediction)**2).sum() / denominator) if denominator > 0 else 0.0
        positive = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
        alpha = None
        log_r2 = None
        if np.count_nonzero(positive) >= 3:
            log_x, log_y = np.log(x[positive]), np.log(y[positive])
            exponent, log_intercept = np.polyfit(log_x, log_y, 1)
            log_prediction = exponent * log_x + log_intercept
            log_denominator = ((log_y - log_y.mean())**2).sum()
            alpha = float(exponent)
            log_r2 = float(1 - ((log_y - log_prediction)**2).sum() / log_denominator) if log_denominator > 0 else 0.0
        windows.append({'start_fs': float(x[0]), 'end_fs': float(x[-1]),
                        'loglog_alpha': alpha, 'loglog_R2': log_r2,
                        'slope_A2_fs': float(slope), 'intercept_A2': float(intercept),
                        'R2': r2, 'D_A2_fs_slope_over_6': float(slope / 6)})
    return windows


def diffusion_regime(fits, thresholds=None):
    """Screen ballistic/caged motion using exponent and cross-window slopes.

    These are finite-trajectory screening conventions requiring C-frozen
    bounds. They do not establish a thermodynamic diffusion coefficient.
    No defaults are permitted to turn an unfrozen diagnostic into PASS.
    """
    late_alpha = fits[-1].get('loglog_alpha') if fits else None
    slopes = [fit['slope_A2_fs'] for fit in fits]
    slope_ratio = None
    if len(slopes) == 3 and all(math.isfinite(value) and value > 0 for value in slopes):
        slope_ratio = max(slopes) / min(slopes)
    metrics = {'late_loglog_alpha': late_alpha,
               'window_slope_ratio_max_over_min': slope_ratio,
               'late_loglog_window_start_fs': fits[-1]['start_fs'] if fits else None,
               'late_loglog_window_end_fs': fits[-1]['end_fs'] if fits else None,
               'screen_is_physical_diffusion_validation': False}
    if thresholds is None:
        return metrics, None
    names = ('min_late_msd_loglog_alpha', 'max_late_msd_loglog_alpha',
             'max_msd_window_slope_ratio')
    if any(name not in thresholds for name in names):
        raise Hold('DIFFUSION_REGIME_THRESHOLDS_NOT_EXPLICIT')
    minimum, maximum, ratio_limit = [float(thresholds[name]) for name in names]
    if (not all(math.isfinite(value) for value in (minimum, maximum, ratio_limit)) or
            not 0 < minimum < 1 < maximum < 2 or ratio_limit < 1):
        raise Hold('DIFFUSION_REGIME_THRESHOLDS_INVALID')
    checks = {'late_loglog_diffusive_regime': late_alpha is not None and
                  math.isfinite(late_alpha) and minimum <= late_alpha <= maximum,
              'MSD_window_slope_consistency': slope_ratio is not None and slope_ratio <= ratio_limit}
    return metrics, checks


def thermal_blocks(rows, number, np):
    output = []
    for index, block in enumerate(np.array_split(np.arange(len(rows)), number)):
        if not len(block):
            continue
        selected = [rows[int(i)] for i in block]
        values = {'block': index + 1, 'start_fs': selected[0]['time_fs'],
                  'end_fs': selected[-1]['time_fs'], 'samples': len(selected)}
        for name in ('T_K', 'sampler_energy_eV', 'pressure_eV_A3'):
            data = np.asarray([row[name] for row in selected])
            values[name + '_mean'] = float(data.mean())
            values[name + '_std'] = float(data.std())
        output.append(values)
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--owner', required=True, choices=['D', 'E', 'F'])
    parser.add_argument('--task-id', required=True)
    parser.add_argument('--task-root', type=Path, required=True,
                        help='Only the owned task folder; no root recursion or sealed data')
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True, help='A new QA folder; never overwrite')
    args = parser.parse_args(argv)
    if args.out.exists():
        parser.error('QA output already exists; choose a new folder')
    if args.out.resolve() == args.task_root.resolve() or args.out.resolve() in args.task_root.resolve().parents:
        parser.error('QA output cannot replace or contain task input folder')
    summary = {'task_id': args.task_id, 'owner': args.owner,
               'AUTO_DIAGNOSTIC_PASS': False, 'manualC_REQUIRED': True,
               'STRUCTURE_PASS': False, 'DFT_LABELS': 'ABSENT',
               'status': 'HOLD_FOR_C_REVIEW', 'warnings': [], 'reasons': []}
    args.out.mkdir(parents=True, exist_ok=False)
    try:
        if not args.config.is_file():
            raise Hold('RELEASE_CONFIG_MISSING')
        config = json.loads(args.config.read_text())
        for name in ('density_reference_pass', 'source_release', 'storage_release', 'continuation_release'):
            if config.get(name) is not True:
                raise Hold(name.upper() + '_HOLD')
        host = config.get('hosts', {}).get(args.owner, {})
        if args.task_id not in host.get('allowed_task_ids', []):
            raise Hold('TASK_NOT_OWNED_BY_HOST')
        expected_root = Path(host.get('sampler_work_root', '/UNVERIFIED')).resolve() / args.task_id
        if args.task_root.resolve() != expected_root:
            raise Hold('TASK_ROOT_DIFFERS_FROM_RELEASED_OWNER_ROOT')
        if (args.task_root / 'RUNNING.lock').exists() or (expected_root.parent / 'HOST_SAMPLER_RUNNING.lock').exists():
            raise Hold('SAMPLER_STILL_RUNNING_OR_LOCK_UNREVIEWED')
        names = ('trajectory.traj', 'thermal.csv', 'source_manifest.json', 'sampler_state.json')
        for name in names:
            if not (args.task_root / name).is_file():
                raise Hold('MISSING_REAL_SAMPLER_OUTPUT:' + name)
        summary['input_sha256'] = {name: sha256_file(args.task_root / name) for name in names}
        summary['config_sha256'] = sha256_file(args.config)
        summary['code_sha256'] = sha256_file(Path(__file__))
        manifest = json.loads((args.task_root / 'source_manifest.json').read_text())
        state = json.loads((args.task_root / 'sampler_state.json').read_text())
        matrix_path = args.config.parent / config.get('tasks_file', 'AUTHORIZED_SIX_TASKS.csv')
        if not matrix_path.is_file():
            raise Hold('TASK_MATRIX_MISSING')
        with matrix_path.open(newline='') as stream:
            tasks = list(csv.DictReader(stream))
        expected_matrix_sha = config.get('tasks_sha256')
        if not expected_matrix_sha or sha256_file(matrix_path) != expected_matrix_sha:
            raise Hold('TASK_MATRIX_SHA_MISMATCH')
        selected = [task for task in tasks if task['task_id'] == args.task_id and task['owner'] == args.owner]
        if len(selected) != 1 or len(tasks) != 6:
            raise Hold('TASK_IDENTITY_NOT_IN_SIX_TASK_MATRIX')
        task = selected[0]
        for key in ('task_id', 'owner'):
            if manifest.get(key) != getattr(args, key.replace('-', '_')) or state.get(key) != manifest[key]:
                raise Hold('MANIFEST_STATE_TASK_IDENTITY_MISMATCH:' + key)
        if (manifest.get('RNG_seed') != int(task['RNG_seed']) or
                manifest.get('trajectory_parent') != task['trajectory_parent']):
            raise Hold('MANIFEST_SEED_OR_PARENT_MISMATCH')
        source_root = config.get('lineage', {}).get('source_root_by_task', {}).get(args.task_id, task.get('source_root'))
        if not source_root or 'TO_BE_FROZEN' in source_root or manifest.get('source_root') != source_root:
            raise Hold('SOURCE_ROOT_NOT_FROZEN_OR_MISMATCHED')
        if manifest.get('potential', {}).get('potential_sha256') != config.get('potential', {}).get('sha256'):
            raise Hold('SAMPLER_POTENTIAL_SHA_MISMATCH')
        for key in ('scientific_parameters_sha256', 'task_matrix_sha256'):
            if not manifest.get(key) or state.get(key) != manifest[key]:
                raise Hold('MANIFEST_CHECKPOINT_IDENTITY_MISMATCH:' + key)
        if state.get('source_manifest_sha256') != summary['input_sha256']['source_manifest.json']:
            raise Hold('SOURCE_MANIFEST_DIFFERS_FROM_CHECKPOINT')
        fixed = manifest.get('fixed_identity', {})
        current_scientific_sampling = {key: value for key, value in config.get('sampling', {}).items()
                                       if key not in RESOURCE_PARAMETER_KEYS}
        if fixed.get('sampling') != current_scientific_sampling or fixed.get('potential') != config.get('potential'):
            raise Hold('CURRENT_CONFIG_SCIENTIFIC_PARAMETERS_DIFFER_FROM_SAMPLER')
        if canonical_sha(fixed) != manifest['scientific_parameters_sha256']:
            raise Hold('SCIENTIFIC_PARAMETER_HASH_INVALID')
        if manifest['task_matrix_sha256'] != expected_matrix_sha:
            raise Hold('SAMPLER_TASK_MATRIX_SHA_MISMATCH')
        bundle = (args.task_root / state.get('restart_bundle', '')).resolve()
        if args.task_root.resolve() not in bundle.parents:
            raise Hold('CHECKPOINT_BUNDLE_ESCAPES_TASK_ROOT')
        for name, expected in state.get('bundle_sha256', {}).items():
            target = (bundle / name).resolve()
            if bundle not in target.parents or not target.is_file() or sha256_file(target) != expected:
                raise Hold('CHECKPOINT_BUNDLE_SHA_MISMATCH:' + name)
        for name in ('trajectory.traj', 'thermal.csv'):
            if state.get('bundle_sha256', {}).get(name) != summary['input_sha256'][name]:
                raise Hold('LIVE_OUTPUT_DIFFERS_FROM_FINAL_CHECKPOINT:' + name)
        if state.get('status') != 'TRAJECTORY_COMPLETE_QA_REQUIRED':
            raise Hold('ONLY_WARMUP_OR_INCOMPLETE_TRAJECTORY_EXISTS')
        sampling = config.get('sampling', {})
        expected_steps = sum(int(sampling[key]) for key in ('hot_steps', 'cool_steps', 'equil_steps', 'production_steps'))
        if int(state.get('step', -1)) != expected_steps:
            raise Hold('FINAL_STEP_DIFFERS_FROM_APPROVED_SCHEDULE')
        import numpy as np
        from ase.io import read, write
        from scipy.special import sph_harm_y
        trajectory = read(args.task_root / 'trajectory.traj', index=':')
        if not trajectory:
            raise Hold('EMPTY_TRAJECTORY')
        edge = float(task['box_edge_A'])
        expected_cell = np.eye(3) * edge
        recorded_steps = []
        all_frame_minimum = math.inf
        for atoms in trajectory:
            identity = geometry_identity(atoms, np)
            if (len(atoms) != 64 or not np.all(atoms.numbers == 47) or
                    not np.all(atoms.pbc) or not np.allclose(atoms.cell.array, expected_cell, atol=1e-10, rtol=0)):
                raise Hold('TRAJECTORY_COMPOSITION_CELL_OR_PBC_MISMATCH')
            if identity['atom_ids'] != list(range(1, 65)):
                raise Hold('TRAJECTORY_ATOM_ORDER_MISMATCH')
            if 'momenta' not in atoms.arrays:
                raise Hold('TRAJECTORY_MOMENTA_MISSING')
            if not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.get_momenta()).all():
                raise Hold('NONFINITE_POSITION_OR_MOMENTA')
            if atoms.calc is not None:
                raise Hold('TRAJECTORY_HAS_CALCULATOR_LABELS; geometry-only archive required')
            if (atoms.info.get('task_id') != args.task_id or
                    atoms.info.get('trajectory_parent') != manifest['trajectory_parent'] or
                    atoms.info.get('source_root') != source_root or
                    atoms.info.get('positions_are_unwrapped') is not True):
                raise Hold('FRAME_PROVENANCE_OR_UNWRAPPED_CONTRACT_MISMATCH')
            step = int(atoms.info.get('step', -1))
            if not math.isclose(float(atoms.info.get('time_fs', -1)), step * float(sampling['dt_fs']), abs_tol=1e-8):
                raise Hold('FRAME_STEP_TIME_MISMATCH')
            vectors, distances = pair_vectors(atoms.positions, edge, np)
            all_frame_minimum = min(all_frame_minimum, float(distances[np.triu_indices(64, 1)].min()))
            recorded_steps.append(step)
        if any(b <= a for a, b in zip(recorded_steps, recorded_steps[1:])):
            raise Hold('DUPLICATED_OR_NONMONOTONIC_TRAJECTORY_STEPS')
        if recorded_steps[-1] != expected_steps:
            raise Hold('FINAL_TRAJECTORY_FRAME_MISSING')
        density = 64 * 107.8682 / 6.02214076e23 / (edge**3 * 1e-24)
        if abs(density - float(task['density_g_cm3'])) > 1e-9:
            raise Hold('DENSITY_CELL_MASS_MISMATCH')
        with (args.task_root / 'thermal.csv').open(newline='') as stream:
            thermal = list(csv.DictReader(stream))
        for row in thermal:
            for key in ('step', 'time_fs', 'target_temperature_K', 'T_K', 'sampler_energy_eV', 'pressure_eV_A3', 'COMspeed_A_fs'):
                row[key] = float(row[key])
                if not math.isfinite(row[key]):
                    raise Hold('NONFINITE_THERMAL_RECORD')
            if not math.isclose(row['time_fs'], row['step'] * float(sampling['dt_fs']), abs_tol=1e-8):
                raise Hold('THERMAL_STEP_TIME_MISMATCH')
        if not thermal or any(b['step'] <= a['step'] for a, b in zip(thermal, thermal[1:])):
            raise Hold('EMPTY_OR_DUPLICATED_THERMAL_RECORDS')
        production_start = expected_steps - int(sampling['production_steps'])
        production = [(i, atoms) for i, atoms in enumerate(trajectory)
                      if atoms.info.get('phase') in ('production', 'production_candidate') and int(atoms.info['step']) > production_start]
        thermal_prod = [row for row in thermal if row['phase'] in ('production', 'production_candidate') and row['step'] > production_start]
        if len(production) < 30 or len(thermal_prod) < 30:
            raise Hold('PRODUCTION_WINDOW_TOO_SHORT')
        if any(abs(row['target_temperature_K'] - 1250) > 1e-8 for row in thermal_prod):
            raise Hold('PRODUCTION_TARGET_TEMPERATURE_DIFFERS_FROM_FROZEN_1250K')
        positions = np.asarray([atoms.positions for _, atoms in production])
        times = np.asarray([float(atoms.info['time_fs']) for _, atoms in production])
        qa = config.get('qa', {})
        dr = float(qa.get('rdf_bin_A', 0.05))
        if not math.isfinite(dr) or not 0.01 <= dr <= 0.2:
            raise Hold('RDF_BIN_INVALID')
        radii, rdf, auto_cutoff, minimum = liquid_rdf(positions, edge, dr, np)
        cutoff = qa.get('neighbor_cutoff_A', auto_cutoff)
        if cutoff is None or not 2.5 < float(cutoff) < edge / 2:
            raise Hold('FIRST_RDF_MINIMUM_NOT_RESOLVED; C must freeze neighbor cutoff')
        cutoff = float(cutoff)
        save_csv(args.out / 'rdf.csv', [{'r_A': float(r), 'g_r': float(g)} for r, g in zip(radii, rdf)], ['r_A', 'g_r'])
        frame_rows, local_rows, cn_hist = [], [], {}
        for (index, atoms), position in zip(production, positions):
            counts, q6, crystal_count, largest = local_order(position, edge, cutoff, np, sph_harm_y)
            frame_rows.append({'frame_index': index, 'step': int(atoms.info['step']), 'time_fs': float(atoms.info['time_fs']),
                               'local_q6_mean': float(q6.mean()), 'local_q6_std': float(q6.std()),
                               'crystal_like_atoms': crystal_count, 'largest_crystal_cluster': largest,
                               'coordination_mean': float(counts.mean())})
            for atom, (cn, value) in enumerate(zip(counts, q6), 1):
                local_rows.append({'frame_index': index, 'atom_id': atom, 'coordination': int(cn), 'local_q6': float(value)})
                cn_hist[int(cn)] = cn_hist.get(int(cn), 0) + 1
        save_csv(args.out / 'local_order_frames.csv', frame_rows, list(frame_rows[0]))
        save_csv(args.out / 'local_order_atoms.csv', local_rows, list(local_rows[0]))
        save_csv(args.out / 'coordination_histogram.csv', [{'coordination': key, 'count': value,
                  'fraction': value / sum(cn_hist.values())} for key, value in sorted(cn_hist.items())], ['coordination', 'count', 'fraction'])
        msd = multiple_origin_msd(positions, times, np)
        fits = fit_msd_windows(msd, np)
        regime_metrics, _ = diffusion_regime(fits)
        save_csv(args.out / 'msd.csv', msd, list(msd[0]))
        save_csv(args.out / 'msd_fits.csv', fits, list(fits[0]))
        blocks = thermal_blocks(thermal_prod, 4, np)
        save_csv(args.out / 'thermal_blocks.csv', blocks, list(blocks[0]))
        temperatures = [row['T_K_mean'] for row in blocks]
        energies = [row['sampler_energy_eV_mean'] for row in blocks]
        pressures = [row['pressure_eV_A3_mean'] for row in blocks]
        summary['diagnostics'] = {'production_frames': len(production), 'production_duration_fs': float(times[-1] - times[0]),
            'rdf_radius_max_A': float(radii[-1]), 'neighbor_cutoff_A': cutoff,
            'neighbor_cutoff_source': 'frozen_config' if 'neighbor_cutoff_A' in qa else 'inferred_first_RDF_minimum_for_C_review',
            'minimum_periodic_pair_distance_A': minimum,
            'minimum_pair_distance_all_archived_phases_A': all_frame_minimum,
            'rdf_finite_N_normalization': True,
            'coherent_dot_threshold': 0.7, 'minimum_coherent_bonds': 7,
            'max_largest_crystal_cluster': max(row['largest_crystal_cluster'] for row in frame_rows),
            'mean_crystal_like_fraction': float(np.mean([row['crystal_like_atoms'] / 64 for row in frame_rows])),
            'late_MSD_A2': msd[-1]['msd_A2'], 'MSD_windows': fits,
            'MSD_diffusion_regime': regime_metrics,
            'block_mean_temperatures_K': temperatures, 'block_mean_sampler_energies_eV': energies,
            'block_mean_pressure_eV_A3': pressures,
            'sampler_energy_or_force_is_DFT_truth': False,
            'pressure_not_forced_to_experimental_value': True}
        required_thresholds = ('min_runtime_distance_A', 'temperature_mean_tolerance_K',
            'temperature_block_range_K', 'sampler_energy_block_range_eV_per_atom',
            'pressure_block_range_eV_A3', 'max_crystal_fraction', 'max_cluster_atoms',
            'min_late_msd_A2', 'min_msd_R2', 'min_decorrelation_msd_A2',
            'min_late_msd_loglog_alpha', 'max_late_msd_loglog_alpha',
            'max_msd_window_slope_ratio')
        if qa.get('thresholds_frozen_by_C') is not True or any(key not in qa for key in required_thresholds):
            summary['reasons'].append('QA_THRESHOLDS_NOT_FROZEN_BY_C; diagnostics alone cannot accept liquid')
        else:
            threshold = {key: float(qa[key]) for key in required_thresholds}
            if not all(math.isfinite(value) and value >= 0 for value in threshold.values()):
                raise Hold('QA_THRESHOLD_NOT_FINITE_OR_NEGATIVE')
            checks = {
                'minimum_distance': all_frame_minimum >= threshold['min_runtime_distance_A'],
                'mean_temperature': abs(np.mean(temperatures) - 1250) <= threshold['temperature_mean_tolerance_K'],
                'temperature_blocks': max(temperatures) - min(temperatures) <= threshold['temperature_block_range_K'],
                'energy_blocks': (max(energies) - min(energies)) / 64 <= threshold['sampler_energy_block_range_eV_per_atom'],
                'pressure_blocks': max(pressures) - min(pressures) <= threshold['pressure_block_range_eV_A3'],
                'crystal_fraction': summary['diagnostics']['mean_crystal_like_fraction'] <= threshold['max_crystal_fraction'],
                'crystal_cluster': summary['diagnostics']['max_largest_crystal_cluster'] <= threshold['max_cluster_atoms'],
                'diffusion_amplitude': msd[-1]['msd_A2'] >= threshold['min_late_msd_A2'],
                'diffusion_windows': len(fits) == 3 and all(fit['slope_A2_fs'] > 0 and fit['R2'] >= threshold['min_msd_R2'] for fit in fits)}
            _, regime_checks = diffusion_regime(fits, threshold)
            checks.update(regime_checks)
            centered = positions - positions.mean(axis=1, keepdims=True)
            decorrelation = ((centered[-1] - centered[0])**2).sum(axis=1).mean()
            checks['within_run_decorrelation'] = float(decorrelation) >= threshold['min_decorrelation_msd_A2']
            summary['diagnostics']['endpoint_decorrelation_MSD_A2'] = float(decorrelation)
            summary['checks'] = checks
            summary['reasons'] += [name + '_DIAGNOSTIC_FAILED' for name, passed in checks.items() if not passed]
            summary['AUTO_DIAGNOSTIC_PASS'] = all(checks.values())
        if summary['AUTO_DIAGNOSTIC_PASS']:
            # Select the last production frame only under preregistered diagnostics.
            # Wrap coordinates for DFT; raw unwrapped source trajectory is retained.
            index, candidate = production[-1]
            candidate = candidate.copy()
            candidate.calc = None
            candidate.set_constraint()
            candidate.wrap()
            candidate.set_momenta(np.zeros((64, 3)))
            candidate.info = {'task_id': args.task_id, 'trajectory_parent': manifest['trajectory_parent'],
                              'source_root': source_root, 'atom_ids': list(range(1, 65)),
                              'source_frame_index': index, 'source_time_fs': float(times[-1]),
                              'STRUCTURE_PASS': False, 'DFT_LABELS': 'ABSENT'}
            candidate.set_array('atom_ids', np.arange(1, 65, dtype=int))
            output_path = args.out / 'selected_candidate.extxyz'
            write(output_path, candidate, format='extxyz',
                  columns=['symbols', 'positions', 'atom_ids'], write_results=False)
            # Freeze the serialized geometry actually read by the DFT runner; extxyz
            # uses finite text precision, so the in-memory geometry hash differs.
            serialized = read(output_path)
            identity = geometry_identity(serialized, np)
            freeze = {'task_id': args.task_id, 'owner': args.owner, 'RNG_seed': manifest['RNG_seed'],
                      'trajectory_parent': manifest['trajectory_parent'], 'source_root': source_root,
                      'geometry_identity': identity, 'geometry_sha256': canonical_sha(identity),
                      'candidate_file_sha256': sha256_file(output_path),
                      'candidate_relative_path': output_path.name,
                      'source_trajectory_sha256': summary['input_sha256']['trajectory.traj'],
                      'source_manifest_sha256': summary['input_sha256']['source_manifest.json'],
                      'source_config_sha256': summary['config_sha256'],
                      'source_frame_index': index, 'source_time_fs': float(times[-1]),
                      'selection_rule': 'last production frame after C-frozen diagnostic/decorrelation gates',
                      'AUTO_DIAGNOSTIC_PASS': True, 'manualC_REQUIRED': True,
                      'STRUCTURE_PASS': False, 'data_role': manifest['source_role'],
                      'independent_test_eligible': False,
                      'sampler_labels_are_reference': False, 'DFT_LABELS': 'ABSENT'}
            save_json(args.out / 'candidate_manifest.json', freeze)
            summary['candidate_geometry_sha256'] = freeze['geometry_sha256']
            summary['candidate_file_sha256'] = freeze['candidate_file_sha256']
            summary['status'] = 'AUTO_DIAGNOSTIC_PASS_C_SIGNOFF_REQUIRED'
        summary['warnings'] += ['64-atom local diagnostics do not validate melting or solidification physics',
            'finite-size RDF is limited to half the box edge; no thermodynamic limit claim',
            'Q6/coherent-bond thresholds are screening conventions, not universal liquid criteria',
            'source-root family grouping is preserved; independent seed is not an independent test family']
    except Exception as error:
        summary['AUTO_DIAGNOSTIC_PASS'] = False
        summary['status'] = 'HOLD_FOR_C_REVIEW'
        summary['reasons'].append(type(error).__name__ + ':' + str(error))
    save_json(args.out / 'qa_report.json', summary)
    files = sorted(path for path in args.out.iterdir() if path.is_file())
    with (args.out / 'SHA256SUMS.txt').open('x') as stream:
        for path in files:
            stream.write(f'{sha256_file(path)}  {path.name}\n')
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0 if summary['AUTO_DIAGNOSTIC_PASS'] else 2


if __name__ == '__main__':
    sys.exit(main())
