#!/usr/bin/env python3
from __future__ import annotations
import csv
import gzip
import json
import math
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/agti-dt-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/agti-dt-cache')
os.environ.setdefault('FONTCONFIG_PATH', '/etc/fonts')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
BASE = HERE.parent / 'liquid_diagnostic_5ps_20261008'
PAIR_IDS = (4508, 4771)
TIME_ORIGIN_PS = 103.8


def frames(path: Path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        while True:
            line = f.readline()
            if not line:
                return
            if line.strip() != 'ITEM: TIMESTEP':
                continue
            step = int(f.readline())
            f.readline()
            natoms = int(f.readline())
            f.readline()
            bounds = [f.readline().split() for _ in range(3)]
            cols = f.readline().split()[2:]
            ix = {name: cols.index(name) for name in ('id', 'type', 'xu', 'yu', 'zu', 'fx', 'fy', 'fz')}
            selected = {}
            max_force = (-math.inf, -1, -1)
            finite = True
            for _ in range(natoms):
                row = f.readline().split()
                aid, typ = int(row[ix['id']]), int(row[ix['type']])
                fx, fy, fz = (float(row[ix[k]]) for k in ('fx', 'fy', 'fz'))
                force = math.sqrt(fx*fx + fy*fy + fz*fz)
                finite &= math.isfinite(force)
                if force > max_force[0]:
                    max_force = (force, aid, typ)
                if aid in PAIR_IDS:
                    selected[aid] = {
                        'type': typ,
                        'xyz': tuple(float(row[ix[k]]) for k in ('xu', 'yu', 'zu')),
                        'force': force,
                    }
            a, b = (selected[aid] for aid in PAIR_IDS)
            delta = [a['xyz'][k] - b['xyz'][k] for k in range(3)]
            yield {
                'step': step,
                'time_local_ps': step * (0.000025 if 'dt025' in path.name else 0.00005),
                'time_cumulative_ps': TIME_ORIGIN_PS + step * (0.000025 if 'dt025' in path.name else 0.00005),
                'natoms': natoms,
                'finite_forces': finite,
                'max_force_eV_A': max_force[0],
                'max_force_atom_id': max_force[1],
                'pair_distance_A': math.sqrt(sum(v*v for v in delta)),
                'pair_4508_force_eV_A': a['force'],
                'pair_4771_force_eV_A': b['force'],
                'bounds': bounds,
            }


def thermo(path: Path):
    rows = []
    with path.open(errors='replace') as f:
        active = False
        columns = []
        for line in f:
            parts = line.split()
            if parts and parts[0] == 'Step' and 'Time' in parts and 'PotEng' in parts:
                columns = parts
                active = True
                continue
            if not active:
                continue
            if len(parts) == len(columns):
                try:
                    values = [float(v) for v in parts]
                except ValueError:
                    active = False
                    continue
                rows.append(dict(zip(columns, values)))
            elif rows:
                active = False
    return rows


def summarize(rows, dt_fs):
    for row in rows:
        row['time_cumulative_ps'] = TIME_ORIGIN_PS + row['step'] * dt_fs / 1000.0
    pairmin = min(rows, key=lambda r: r['pair_distance_A'])
    forcepeak = max(rows, key=lambda r: r['max_force_eV_A'])
    pairforcepeak = max(rows, key=lambda r: max(r['pair_4508_force_eV_A'], r['pair_4771_force_eV_A']))
    around = [r for r in rows if 103.88 <= r['time_cumulative_ps'] <= 103.93]
    return {
        'frames': len(rows),
        'atom_counts': sorted(set(r['natoms'] for r in rows)),
        'all_sampled_forces_finite': all(r['finite_forces'] for r in rows),
        'minimum_target_pair_distance': {k: pairmin[k] for k in ('time_cumulative_ps', 'pair_distance_A', 'pair_4508_force_eV_A', 'pair_4771_force_eV_A')},
        'maximum_sampled_force': {k: forcepeak[k] for k in ('time_cumulative_ps', 'max_force_eV_A', 'max_force_atom_id')},
        'maximum_target_pair_force': {
            'time_cumulative_ps': pairforcepeak['time_cumulative_ps'],
            'pair_4508_force_eV_A': pairforcepeak['pair_4508_force_eV_A'],
            'pair_4771_force_eV_A': pairforcepeak['pair_4771_force_eV_A'],
            'pair_distance_A': pairforcepeak['pair_distance_A'],
        },
        'minimum_pair_distance_103p88_to_103p93_ps': min((r['pair_distance_A'] for r in around), default=None),
    }


def energy_balance(rows):
    corrected = [r['TotEng'] + r['f_bathfix'] for r in rows]
    return {
        'quantity': 'Etot + cumulative Langevin energy exchange',
        'start_eV': corrected[0],
        'end_eV': corrected[-1],
        'range_eV': max(corrected)-min(corrected),
        'start_to_end_delta_eV': corrected[-1]-corrected[0],
    }


def paired_plot(baseline, refined):
    fig, (axr, axf) = plt.subplots(2, 1, figsize=(9.2, 7.0), sharex=True, constrained_layout=True)
    for name, rows, color, marker in [
        ('0.05 fs timestep (5 fs samples)', baseline, '#d1495b', 'o'),
        ('0.025 fs timestep (2.5 fs samples)', refined, '#00798c', 's'),
    ]:
        t = [r['time_cumulative_ps'] for r in rows]
        d = [r['pair_distance_A'] for r in rows]
        fp = [max(r['pair_4508_force_eV_A'], r['pair_4771_force_eV_A']) for r in rows]
        axr.plot(t,d,marker=marker,ms=3.3,lw=1.4,color=color,label=name)
        axf.plot(t,fp,marker=marker,ms=3.3,lw=1.4,color=color,label=name)
    axr.set_ylabel('Ag–Ag separation (Å)')
    axr.set_title('Incoming Ag projectile (4771) and target Ag atom (4508)')
    axr.legend(frameon=False)
    axr.grid(alpha=.25)
    axf.set_ylabel('Pair force magnitude (eV/Å)')
    axf.set_xlabel('Cumulative time (ps)')
    axf.grid(alpha=.25)
    axf.legend(frameon=False)
    fig.suptitle('Hit 466 timestep sensitivity',fontsize=15)
    out=HERE/'images/dt_sensitivity_hit466.png'
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out,dpi=170)
    plt.close(fig)


def main():
    baseline = list(frames(BASE/'dumps/impact_466_62eV.dump.gz'))
    refined = list(frames(HERE/'dumps/impact_466_dt025.dump'))
    for r in baseline:
        r['time_cumulative_ps'] = TIME_ORIGIN_PS + r['time_local_ps']
    for r in refined:
        r['time_cumulative_ps'] = TIME_ORIGIN_PS + r['time_local_ps']
    outcsv = HERE/'analysis/pair_force_timeseries.csv'
    with outcsv.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['branch','step','time_cumulative_ps','natoms','pair_distance_A','pair_4508_force_eV_A','pair_4771_force_eV_A','max_force_eV_A','max_force_atom_id','finite_forces'])
        writer.writeheader()
        for name, rows in [('0.05_fs', baseline), ('0.025_fs', refined)]:
            for r in rows:
                writer.writerow({'branch':name, **{k:r[k] for k in writer.fieldnames if k != 'branch'}})
    result = {
        'same_input_state_sha256': '0de755108b45f7c53371a96e5774fbf546b0ec2beeffcfc0b357095dbe7f3e0e',
        'injection': {'hit':466,'Ag_projectile_id':4771,'target_Ag_id':4508,'energy_eV':62.0,'velocity_z_A_per_ps':-105.3164,'position_A':[6.6810405777,20.9416452647,60.0],'langevin_seed':624166},
        'comparison': {
            'baseline': {'timestep_fs':0.05,'steps':4000,'duration_ps':0.2,**summarize(baseline,0.05),'thermo_energy_balance':energy_balance(thermo(BASE/'logs/impact_466_62eV.lammps'))},
            'half_timestep': {'timestep_fs':0.025,'steps':8000,'duration_ps':0.2,**summarize(refined,0.025),'thermo_energy_balance':energy_balance(thermo(HERE/'logs/impact_466_dt025.lammps'))},
        },
        'thermo_baseline': thermo(BASE/'logs/impact_466_62eV.lammps'),
        'thermo_half_timestep': thermo(HERE/'logs/impact_466_dt025.lammps'),
        'limits':['The Langevin thermostat is stochastic; same seed and initial state do not make trajectories identical at different dt.','The largest sampled force can depend on dump-frame spacing; it is not a continuous-time maximum.','This is a one-event timestep sensitivity check, not a validation of the Ag-Ti3SiC2 cross potential.'],
    }
    at_event=[]
    for label, rows in [('baseline',baseline),('half_timestep',refined)]:
        row=min(rows,key=lambda r:abs(r['time_cumulative_ps']-103.905))
        at_event.append({'branch':label,**{k:row[k] for k in ('time_cumulative_ps','pair_distance_A','pair_4508_force_eV_A','pair_4771_force_eV_A')}})
    result['shared_time_103p905_ps']=at_event
    outjson = HERE/'analysis/dt_sensitivity_summary.json'
    outjson.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    for name in ('baseline','half_timestep'):
        s=result['comparison'][name]
        print(name, json.dumps(s,ensure_ascii=False))
    paired_plot(baseline,refined)
    print('thermo rows',len(result['thermo_baseline']),len(result['thermo_half_timestep']))


if __name__=='__main__':
    main()
