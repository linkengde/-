#!/usr/bin/env python3
"""Check LAMMPS energy-force consistency on one four-element interface frame."""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path('/workspace/potential-fit/round7_pod_snap')
MODEL = ROOT / 'models/final_v10'
SOURCE = ROOT / 'dataset_v10_v25_supplement/validation/frame_0000.xyz'
WORK = MODEL / 'finite_difference_check'
LMP = Path('/workspace/software/lammps-pod-snap-clean/build/lmp')
ELEMENTS = ['Ag', 'C', 'Si', 'Ti']
H = 1e-4


def pair_config(name: str) -> tuple[str, str]:
    if name == 'POD':
        return ('pod', f"pair_coeff * * {MODEL/'POD/AgTiSiC_POD.pod'} "
                f"{MODEL/'POD/pod_coefficients.pod'} Ag C Si Ti")
    return ('snap', f"pair_coeff * * {MODEL/'SNAP/AgTiSiC_SNAP_final_v10.snapcoeff'} "
            f"{MODEL/'SNAP/AgTiSiC_SNAP_final_v10.snapparam'} Ag C Si Ti")


def energy_force(atoms, name: str, case: Path, tag: str) -> tuple[float, np.ndarray]:
    case.mkdir(parents=True, exist_ok=True)
    data = case / f'{tag}.data'
    out = case / f'{tag}.energy'
    dump = case / f'{tag}.dump'
    inp = case / f'{tag}.in'
    write(data, atoms, format='lammps-data', atom_style='atomic',
          specorder=ELEMENTS, masses=True, reduce_cell=False, force_skew=True)
    boundary = ' '.join('p' if x else 'f' for x in atoms.pbc)
    style, coeff = pair_config(name)
    inp.write_text(f'''clear
units metal
atom_style atomic
boundary {boundary}
read_data {data}
pair_style {style}
{coeff}
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
thermo 1
thermo_style custom step atoms pe
thermo_modify lost error flush yes format float %.16g
run 0
print "ENERGY $(pe:%.17g)" file {out} screen no
write_dump all custom {dump} id fx fy fz modify sort id
''')
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = ('/workspace/software/openmpi-dev-5.0.7/usr/lib/x86_64-linux-gnu/openmpi/lib:'
                              '/workspace/software/lammps-pod-snap-clean/build:' + env.get('LD_LIBRARY_PATH', ''))
    env['OMP_NUM_THREADS'] = '1'
    result = subprocess.run([str(LMP), '-in', str(inp), '-log', 'none'], cwd=case,
                            env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True)
    if result.returncode:
        raise RuntimeError(f'{name} {tag} failed: {result.stdout[-3000:]}')
    energy = float(out.read_text().strip().split()[1])
    rows = np.loadtxt(dump, skiprows=9)
    forces = rows[np.argsort(rows[:, 0].astype(int)), 1:4]
    return energy, forces


def main():
    atoms = read(SOURCE, format='extxyz')
    rows = {}
    for name in ('POD', 'SNAP'):
        folder = WORK / name
        e0, f0 = energy_force(atoms.copy(), name, folder, 'base')
        tests = []
        # One atom of every element; perturb all three Cartesian components.
        selected = [next(i for i, sym in enumerate(atoms.get_chemical_symbols()) if sym == el)
                    for el in ('Ag', 'C', 'Si', 'Ti')]
        for i in selected:
            for axis in range(3):
                plus, minus = atoms.copy(), atoms.copy()
                plus.positions[i, axis] += H
                minus.positions[i, axis] -= H
                ep, _ = energy_force(plus, name, folder, f'i{i}_a{axis}_plus')
                em, _ = energy_force(minus, name, folder, f'i{i}_a{axis}_minus')
                fd_force = -(ep - em) / (2 * H)
                model_force = float(f0[i, axis])
                tests.append({'atom_index_1based': i + 1,
                              'element': atoms[i].symbol,
                              'axis': 'xyz'[axis],
                              'force_lammps_eV_A': model_force,
                              'force_finite_difference_eV_A': fd_force,
                              'abs_error_eV_A': abs(model_force - fd_force)})
        err = np.asarray([x['abs_error_eV_A'] for x in tests])
        rows[name] = {'source': str(SOURCE), 'natoms': len(atoms), 'step_A': H,
                      'selected_elements': ['Ag', 'C', 'Si', 'Ti'],
                      'components_checked': len(tests),
                      'max_abs_error_eV_A': float(err.max()),
                      'rms_abs_error_eV_A': float(np.sqrt(np.mean(err ** 2))),
                      'energy_eV': e0, 'checks': tests}
    output = WORK / 'summary.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(rows, indent=2) + '\n')
    print(json.dumps({k: {key: v for key, v in r.items() if key != 'checks'}
                      for k, r in rows.items()}, indent=2))


if __name__ == '__main__':
    main()
