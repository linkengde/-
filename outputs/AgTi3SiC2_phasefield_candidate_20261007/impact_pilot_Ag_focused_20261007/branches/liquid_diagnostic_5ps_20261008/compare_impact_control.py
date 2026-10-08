#!/usr/bin/env python3
from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
CONTROL = ROOT.parent / 'liquid_control_noimpact_5ps_20261008'


def read_csv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def load(root, filename, fields):
    rows = read_csv(root / filename)
    return {name: np.array([float(row[col]) for row in rows]) for name, col in fields.items()}


def smooth(values, width=41):
    return np.convolve(values, np.ones(width) / width, mode='same')


def main():
    impact_temp = load(ROOT, 'analysis_local_temperature.csv', {'t':'time_ps', 'Ag':'Ag_matrix_T_K', 'TSC':'TSC_T_K'})
    control_temp = load(CONTROL, 'analysis_local_temperature.csv', {'t':'time_ps', 'Ag':'Ag_matrix_T_K', 'TSC':'TSC_T_K'})
    impact_order = load(ROOT, 'analysis_local_order.csv', {'t':'time_ps', 'q6':'mean_q6', 'cn':'mean_cn35'})
    control_order = load(CONTROL, 'analysis_local_order.csv', {'t':'time_ps', 'q6':'mean_q6', 'cn':'mean_cn35'})
    impact_msd = load(ROOT, 'analysis_msd.csv', {'lag':'lag_ps', 'msd':'drift_corrected_MSD_A2'})
    control_msd = load(CONTROL, 'analysis_msd.csv', {'lag':'lag_ps', 'msd':'drift_corrected_MSD_A2'})
    colors = {'impact':'#c45035', 'control':'#337f91'}
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), dpi=180, constrained_layout=True)
    ax = axes[0, 0]
    for dat, label, c in ((impact_temp, 'continued Ag impacts', colors['impact']),
                          (control_temp, 'no-impact control', colors['control'])):
        ax.plot(dat['t'], dat['Ag'], color=c, alpha=.22, lw=.7)
        ax.plot(dat['t'], smooth(dat['Ag']), color=c, lw=1.7, label=label)
    ax.set_ylabel('local Ag kinetic temperature (K)')
    ax.set_title('Ag zone temperature (0.2 ps rolling mean)')
    ax.grid(alpha=.25); ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    ax.plot(impact_order['t'], impact_order['q6'], color=colors['impact'], lw=1.5, label='continued Ag impacts')
    ax.plot(control_order['t'], control_order['q6'], color=colors['control'], lw=1.5, label='no-impact control')
    ax.set_ylabel('mean local Ag Q₆')
    ax.set_title('Local orientational order (sampled every 50 fs)')
    ax.grid(alpha=.25); ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    ax.plot(impact_order['t'], impact_order['cn'], color=colors['impact'], lw=1.5, label='continued Ag impacts')
    ax.plot(control_order['t'], control_order['cn'], color=colors['control'], lw=1.5, label='no-impact control')
    ax.set_ylabel('mean Ag neighbors within 3.5 Å')
    ax.set_xlabel('cumulative time (ps)')
    ax.set_title('Local Ag coordination')
    ax.grid(alpha=.25); ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    ax.plot(impact_msd['lag'], impact_msd['msd'], 'o-', ms=4, color=colors['impact'], label='continued Ag impacts')
    ax.plot(control_msd['lag'], control_msd['msd'], 'o-', ms=4, color=colors['control'], label='no-impact control')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('time lag (ps)')
    ax.set_ylabel('drift-corrected MSD (Å²)')
    ax.set_title('MSD of 227 Ag atoms tagged at t=100 ps')
    ax.grid(alpha=.25, which='both'); ax.legend(frameon=False, fontsize=8)
    fig.suptitle('Impact and no-impact trajectories from the same 100 ps state', fontsize=14, weight='bold')
    out = ROOT / 'images/impact_vs_control_5ps.png'
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)

    si = json.loads((ROOT / 'analysis_summary.json').read_text())
    sc = json.loads((CONTROL / 'analysis_summary.json').read_text())
    comparison = {
        'starting_state': 'same hit-446 (100 ps) coordinates and velocities',
        'window_ps': [100.0, 105.0],
        'frames_per_branch': 1001,
        'sampling_interval_fs': 5.0,
        'tagged_initial_Ag_atoms': 227,
        'drift_corrected_MSD_A2': {
            '1ps': {'impact': si['msd_drift_corrected_A2']['1.000_ps'], 'control': sc['msd_drift_corrected_A2']['1.000_ps']},
            '5ps': {'impact': si['msd_drift_corrected_A2']['5.000_ps'], 'control': sc['msd_drift_corrected_A2']['5.000_ps']},
        },
        'apparent_MSD_slope_0p2_to_1ps_A2_per_ps': {
            'impact': si['linear_fit_msd_slope_0p2_to_1ps_A2_per_ps'],
            'control': sc['linear_fit_msd_slope_0p2_to_1ps_A2_per_ps'],
        },
        'tracer_fraction_still_in_r8_z12_55': {
            'impact': si['fraction_of_initial_tracers_still_in_r8_z12_55_region'],
            'control': sc['fraction_of_initial_tracers_still_in_r8_z12_55_region'],
        },
        'mean_Q6_first_last': {
            'impact': si['mean_q6_first_last'],
            'control': sc['mean_q6_first_last'],
        },
        'local_Ag_temperature_median_p10_p90_K': {
            'impact': si['local_matrix_Ag_temperature_K_median_p10_p90'],
            'control': sc['local_matrix_Ag_temperature_K_median_p10_p90'],
        },
        'impact_branch_largest_sampled_force': si['maximum_force_in_5fs_dump_frames'],
        'control_branch_largest_sampled_force': sc['maximum_force_in_5fs_dump_frames'],
        'limits': [
            'The high-temperature impact trajectory is non-equilibrium; fitted MSD slope is apparent mobility, not equilibrium diffusion.',
            'Q6 and coordination are local structural indicators; they do not alone establish a liquid phase.',
            'Ag-Ti3SiC2 cross interactions are not independently validated.'
        ]
    }
    (ROOT / 'comparison_impact_control.json').write_text(json.dumps(comparison, indent=2) + '\n')
    print(json.dumps(comparison, indent=2))
    print(out)


if __name__ == '__main__':
    main()
