#!/workspace/venvs/gpaw/bin/python
"""Exploratory fixed-geometry Ag FCC DFT equation-of-state diagnostic.

These points are diagnostics only. They must not be silently added to the
formal V14 training set until the DFT convergence protocol is approved.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write
from gpaw import GPAW, FermiDirac, Mixer, PW

WORK = Path(__file__).resolve().parent
OUT = WORK / 'v14_stage2_p0_ag_dft_diagnostics' / 'ag_fcc_eos_retry'
LATTICE_PARAMETERS_A = [3.90, 3.60, 3.30]
ECUT_EV = 350
KPTS = (3, 3, 3)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    progress = OUT / 'progress.json'
    rows = json.loads(progress.read_text()).get('rows', []) if progress.exists() else []
    completed_a = {float(r['lattice_parameter_A']) for r in rows if r.get('status') == 'complete'}
    for a in LATTICE_PARAMETERS_A:
        if a in completed_a:
            continue
        tag = f'a_{a:.2f}A'
        atoms = bulk('Ag', 'fcc', a=a, cubic=True)
        log = OUT / f'{tag}.gpaw.txt'
        calc = GPAW(
            mode=PW(ECUT_EV), xc='PBE', kpts=KPTS, symmetry='off',
            occupations=FermiDirac(0.10),
            mixer=Mixer(beta=0.05, nmaxold=8, weight=100),
            convergence={'density': 1e-5, 'eigenstates': 1e-5, 'energy': 1e-5},
            maxiter=180, txt=str(log),
        )
        atoms.calc = calc
        started = time.time()
        row = {'lattice_parameter_A': a, 'volume_per_atom_A3': atoms.get_volume()/len(atoms),
               'nearest_neighbor_A': a/np.sqrt(2), 'natoms': len(atoms),
               'ecut_eV': ECUT_EV, 'kpts': list(KPTS), 'sigma_eV': 0.10,
               'symmetry': 'off', 'status': 'running'}
        try:
            e_extrap = float(atoms.get_potential_energy(force_consistent=False))
            e_free = float(atoms.get_potential_energy(force_consistent=True))
            forces = np.asarray(atoms.get_forces(), dtype=float)
            stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
            log_text = log.read_text(errors='replace') if log.exists() else ''
            converged = (bool(getattr(getattr(calc, 'scf', None), 'converged', False))
                         or ('Converged in ' in log_text and 'Did not converge!' not in log_text))
            row.update({
                'energy_free_eV': e_free, 'energy_extrapolated_eV': e_extrap,
                'force_rmse_vector_eV_A': float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
                'max_force_eV_A': float(np.linalg.norm(forces, axis=1).max()),
                'stress_eV_A3': stress.tolist(), 'scf_converged': converged,
                'scf_convergence_source': 'gpaw_flag_or_log',
                'scf_iterations': int(calc.get_number_of_iterations()),
                'status': 'complete' if converged and np.isfinite([e_free, e_extrap]).all() else 'not_converged',
            })
            if converged:
                atoms.info['diagnostic_only'] = True
                atoms.info['energy_free_eV'] = e_free
                atoms.info['energy_extrapolated_eV'] = e_extrap
                atoms.info['stress_eV_A3'] = stress.tolist()
                atoms.calc = SinglePointCalculator(atoms, energy=e_free, forces=forces, stress=stress.reshape(-1))
                xyz = OUT / f'{tag}.extxyz'
                write(xyz, atoms, format='extxyz')
                row['extxyz'] = str(xyz)
                row['extxyz_sha256'] = sha(xyz)
        except Exception as exc:  # keep failed points visible in the audit
            row.update({'status': 'failed', 'error': f'{type(exc).__name__}: {exc}'})
        row['elapsed_seconds'] = round(time.time()-started, 3)
        row['log'] = str(log)
        row['log_sha256'] = sha(log) if log.exists() else None
        rows.append(row)
        (OUT / 'progress.json').write_text(json.dumps({'rows': rows}, indent=2)+'\n')
        print(json.dumps({k: row.get(k) for k in ('lattice_parameter_A', 'status', 'energy_free_eV', 'scf_iterations', 'elapsed_seconds')}), flush=True)
    good = [r for r in rows if r.get('status') == 'complete']
    if good:
        ref = next((r for r in good if r['lattice_parameter_A'] == 4.09), None)
        if ref:
            for r in good:
                r['relative_energy_to_a4p09_eV_per_atom'] = (r['energy_free_eV']-ref['energy_free_eV'])/r['natoms']
    summary = {'purpose': 'diagnostic_only_not_formal_V14_training', 'code': 'GPAW/PBE',
               'settings': {'ecut_eV': ECUT_EV, 'kpts': list(KPTS), 'sigma_eV': 0.10,
                            'max_scf_iterations': 180,
                            'symmetry': 'off', 'scf_tolerances': {'density':1e-5,'eigenstates':1e-5,'energy':1e-5}},
               'rows': rows}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({'summary':str(OUT/'summary.json'), 'n_complete':len(good), 'n_total':len(rows)}, indent=2))


if __name__ == '__main__':
    main()
