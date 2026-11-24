#!/workspace/venvs/gpaw/bin/python
"""Summarize the archived 52-atom EAM/TSCBOP/Morse sampler cycle."""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path('/workspace/potential-fit/round7_pod_snap/acquisition/md')
LOG = ROOT / 'interface_heat_cool_morse.log'
DUMP = ROOT / 'interface_heat_cool_morse.dump'
OUT = Path(__file__).resolve().parent / 'v14_stage2_p0_ag_dft_diagnostics' / 'legacy_hybrid_cycle_audit.json'


def main():
    text = LOG.read_text(errors='replace')
    thermo = []
    for line in text.splitlines():
        m = re.match(r'^\s*(\d+)\s+([\d.eE+-]+)\s+(\d+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s*$', line)
        if m:
            s, t, n, temp, pe, ke, etot, fmax_s = m.groups()
            fmax = float(fmax_s)
            thermo.append({'step': int(s), 'time_ps': float(t), 'natoms': int(n),
                           'temperature_K': float(temp), 'pe_eV': float(pe),
                           'ke_eV': float(ke), 'etotal_eV': float(etot),
                           'fmax_eV_A': fmax})
    frames = read(DUMP, index=':', format='lammps-dump-text')
    frame_rows = []
    all_counts = None
    for a in frames:
        typ = np.asarray(a.arrays['type'], dtype=int)
        symbols = [{1: 'Ti', 2: 'Si', 3: 'C', 4: 'Ag'}[int(x)] for x in typ]
        a.set_chemical_symbols(symbols)
        if all_counts is None:
            all_counts = {s: symbols.count(s) for s in sorted(set(symbols))}
        d = a.get_all_distances(mic=True)
        s = np.asarray(symbols)
        ag = np.flatnonzero(s == 'Ag')
        sk = np.flatnonzero(s != 'Ag')
        ti = np.flatnonzero(s == 'Ti')
        c = np.flatnonzero(s == 'C')
        agag = float(d[np.ix_(ag, ag)][np.triu_indices(len(ag), 1)].min()) if len(ag) > 1 else None
        cross = float(d[np.ix_(ag, sk)].min()) if len(ag) and len(sk) else None
        tic = float(d[np.ix_(ti, c)].min()) if len(ti) and len(c) else None
        frame_rows.append({'step': len(frame_rows) * 500, 'time_ps': len(frame_rows) * 0.05,
                           'natoms': len(a), 'min_Ag_Ag_A': agag,
                           'min_Ag_skeleton_A': cross, 'min_Ti_C_A': tic})
    if not thermo or not frame_rows:
        raise SystemExit('Missing thermo or trajectory rows')
    peak = max(thermo, key=lambda x: x['temperature_K'])
    last_window = [r for r in thermo if r['time_ps'] >= 3.0]
    out = {
        'status': 'completed exploratory code-path test; not evidence of physical Ag melting/freezing',
        'potential': 'Ag_u3.eam + TSCBOP.tersoff + calibration-v0 Ag/skeleton Morse cross',
        'potential_use': 'exploratory sampler only; Morse cross is not a validated reference',
        'input': str(ROOT / 'in.interface_cycle_morse_exploratory'),
        'duration_ps': 3.5, 'timestep_ps': 0.0001,
        'atom_counts': all_counts,
        'thermo_rows': len(thermo), 'peak_temperature': peak,
        'final_thermo': thermo[-1],
        'last_0p5ps_temperature_min_K': min(r['temperature_K'] for r in last_window),
        'last_0p5ps_temperature_max_K': max(r['temperature_K'] for r in last_window),
        'all_atoms_retained': all(r['natoms'] == thermo[0]['natoms'] for r in thermo),
        'dump_frames': len(frame_rows),
        'minima_over_dump': {k: min(r[k] for r in frame_rows if r[k] is not None)
                             for k in ['min_Ag_Ag_A', 'min_Ag_skeleton_A', 'min_Ti_C_A']},
        'initial_distances': frame_rows[0], 'final_dump_distances': frame_rows[-1],
        'limitations': [
            'Only four Ag atoms are present; no meaningful liquid-pool size or connectivity can be established.',
            'The whole 52-atom system is thermostatted, and this is not an electron-lattice TTM test.',
            'Global temperature oscillates during the nominal heat/hold schedule.',
            'The run does not include Q6/PTM/MSD phase classification or independent DFT checks of the cross interaction.',
            'This archived test does not establish a validated cooling/solidification result.'
        ]
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
