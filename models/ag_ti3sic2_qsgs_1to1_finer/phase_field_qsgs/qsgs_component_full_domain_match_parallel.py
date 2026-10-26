from __future__ import annotations

import contextlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
MODEL = REPO / 'models/ag_ti3sic2_qsgs_1to1_finer'
PHASE = MODEL / 'phase_field_qsgs'
OUTDIR = PHASE / 'nested_grid_diagnostics'
TMP_FIELDS = Path('/tmp/qsgs50_fields')
TMP_FIELDS.mkdir(parents=True, exist_ok=True)
OUTDIR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(MODEL))
sys.path.insert(0, str(PHASE))
from qsgs_component_full_domain_match import (  # noqa: E402
    COARSE_GROUP,
    NATIVE_GROUP,
    GEOMETRY,
    ELECTRON_TABLE,
    _run_one,
    _spatial_comparison,
    _temporal_comparison,
)


def run_case(dt_fs: float, model: str, group: tuple[int, int, int]):
    label = f'{dt_fs:g}fs_{model}'
    log_path = OUTDIR / f'qsgs_component_full_domain_50fs_{label}.log'
    field_path = TMP_FIELDS / f'{label}.npz'
    with log_path.open('w', encoding='utf-8') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        metrics, fields = _run_one(
            GEOMETRY,
            ELECTRON_TABLE,
            tile_xy=3,
            group=group,
            preconditioner='amg',
            dt_ps=dt_fs * 1.0e-3,
            duration_ps=0.05,
            label=label,
        )
        np.savez(field_path, **fields)
    return label, metrics, str(field_path), str(log_path)


def main():
    cases = [
        (0.3125, 'component_2x1x1', COARSE_GROUP),
        (0.3125, 'native_voxel', NATIVE_GROUP),
        (0.15625, 'component_2x1x1', COARSE_GROUP),
        (0.15625, 'native_voxel', NATIVE_GROUP),
    ]
    print('Launching four independent full-domain runs on four CPU workers.', flush=True)
    completed = {}
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(run_case, *case) for case in cases]
        for future in as_completed(futures):
            label, metrics, field_path, log_path = future.result()
            completed[label] = {'metrics': metrics, 'fields': field_path, 'log': log_path}
            print(
                f'Completed {label}: wall={metrics["run_seconds"]:.1f}s, '
                f'Tl,max={metrics["max_lattice_temperature_K"]:.5f} K, '
                f'liquid={metrics["ag_liquid_volume_fraction"]:.8g}',
                flush=True,
            )

    by_dt = {}
    source_powers = []
    for dt_fs in (0.3125, 0.15625):
        tag = f'{dt_fs:g}fs'
        comp_key = f'{tag}_component_2x1x1'
        native_key = f'{tag}_native_voxel'
        comp = completed[comp_key]
        native = completed[native_key]
        with np.load(comp['fields']) as a, np.load(native['fields']) as b:
            spatial = _spatial_comparison(
                {k: a[k] for k in a.files}, {k: b[k] for k in b.files}
            )
        by_dt[tag] = {
            'component_2x1x1': comp['metrics'],
            'native_voxel': native['metrics'],
            'component_vs_native_spatial_error': spatial,
        }
        source_powers += [comp['metrics']['source_power_W'], native['metrics']['source_power_W']]

    refinement = {}
    for model in ('component_2x1x1', 'native_voxel'):
        old = completed[f'0.3125fs_{model}']
        new = completed[f'0.15625fs_{model}']
        with np.load(old['fields']) as a, np.load(new['fields']) as b:
            refinement[model] = _temporal_comparison(
                {k: b[k] for k in b.files}, {k: a[k] for k in a.files}
            )
    by_dt['0.15625fs']['time_step_refinement_vs_previous'] = {
        'dt_fs_previous': 0.3125,
        'dt_fs_current': 0.15625,
        **refinement,
    }

    source_relative_difference = (max(source_powers) - min(source_powers)) / max(abs(source_powers[0]), 1.0e-30)
    spatial_gate = all(
        row['component_vs_native_spatial_error'][key] <= 0.05
        for row in by_dt.values()
        for key in (
            'electron_excess_energy_eV_A3_relative_L1',
            'lattice_enthalpy_eV_A3_relative_L1',
            'electron_to_lattice_energy_eV_A3_relative_L1',
        )
    )
    energy_gate = all(
        result['metrics']['energy_residual_relative_to_source'] <= 1.0e-10
        for result in completed.values()
    )
    timestep_gate = all(
        refinement[model][key] <= 0.05
        for model in ('component_2x1x1', 'native_voxel')
        for key in (
            'electron_excess_energy_eV_A3_relative_L1',
            'lattice_enthalpy_eV_A3_relative_L1',
            'electron_to_lattice_energy_eV_A3_relative_L1',
            'electron_temperature_K_peak_relative_change',
            'lattice_temperature_K_peak_relative_change',
            'ag_liquid_fraction_mean_absolute_change',
        )
    )
    source_gate = source_relative_difference <= 1.0e-12
    output = {
        'status': 'FULL_DOMAIN_QSGS_COMPONENT_TTM_MATCHED_0.05PS',
        'baseline_commit': '3feb5bc',
        'acceptance': {
            'source_power_match_passed': source_gate,
            'energy_ledger_passed': energy_gate,
            'spatial_L1_error_le_5_percent_passed': spatial_gate,
            'time_step_refinement_le_5_percent_passed': timestep_gate,
            'numerical_50fs_gate_passed': source_gate and energy_gate and spatial_gate and timestep_gate,
            'stable_melt_pool_observed': False,
            'note': 'Numerical transport acceptance does not establish a stable melt pool or validate scaffold/interface heat transfer.',
        },
        'geometry_file': GEOMETRY.name,
        'electron_table': ELECTRON_TABLE.name,
        'native_voxel_resolution_A': 'approximately 3 Å source voxels',
        'tile_xy': 3,
        'duration_ps': 0.05,
        'component_coarse_voxels_xyz': list(COARSE_GROUP),
        'reference_coarse_voxels_xyz': list(NATIVE_GROUP),
        'heat_source': {
            'peak_flux_eV_ps_A2': 5.0,
            'gaussian_radius_A': 30.0,
            'gaussian_definition': 'exp(-r^2/R^2)',
            'absorption_depth_A': 20.0,
            'source_mask': 'native binary Ag voxels; exact integrals conservatively aggregated to components',
            'source_center': 'center of 3x3 periodic domain in x/y; depth from top z surface',
            'source_power_W_by_case': source_powers,
            'relative_power_difference': source_relative_difference,
        },
        'boundary_conditions': 'periodic x/y including outer supercell seams; insulated z; boundary A; no scaffold heat storage or cross-phase transport',
        'time_steps_fs': [0.3125, 0.15625],
        'linear_preconditioner': 'amg',
        'results_by_time_step': by_dt,
        'run_logs': {k: v['log'] for k, v in completed.items()},
        'interpretation': 'Full 3x3 QSGS continuum transport comparison only; candidate Ag properties; Ti3SiC2 is a geometric mask.',
    }
    output_path = OUTDIR / 'qsgs_component_full_domain_50fs_match.json'
    output_path.write_text(json.dumps(output, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(output['acceptance'], indent=2))
    print('source relative difference:', source_relative_difference)
    for tag, row in by_dt.items():
        print(tag, row.get('component_vs_native_spatial_error', {}))
    print('wrote', output_path)


if __name__ == '__main__':
    main()
