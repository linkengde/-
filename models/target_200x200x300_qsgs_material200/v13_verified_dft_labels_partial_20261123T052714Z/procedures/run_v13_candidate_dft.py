#!/workspace/venvs/gpaw/bin/python
"""Generate a provenance-checked fixed-geometry GPAW energy/force label."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
import traceback
from pathlib import Path

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from gpaw import GPAW, FermiDirac, Mixer, PW
from gpaw.mpi import world

BASE = Path('/workspace/potential-fit/round7_pod_snap/models/target_200x200x300_qsgs_material200/v13_candidate_sources_20261123')
MANIFEST = BASE / 'training_candidate_manifest_v13.json'
LABEL_ROOT = BASE / 'dft_labels'
TOL = {'density': 1e-5, 'eigenstates': 1e-5, 'energy': 1e-5}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def canonical_hash(atoms):
    frac = atoms.get_scaled_positions(wrap=True)
    rows = sorted((s, *(round(float(x), 8) for x in q)) for s, q in zip(atoms.get_chemical_symbols(), frac))
    obj = {'cell': np.round(np.asarray(atoms.cell), 8).tolist(),
           'pbc': np.asarray(atoms.pbc, bool).tolist(), 'atoms': rows}
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def scf_converged(calc, log):
    value = bool(getattr(getattr(calc, 'scf', None), 'converged', False))
    if not value and log.exists():
        t = log.read_text(errors='replace')
        value = 'Converged in ' in t and 'Did not converge!' not in t
    return value


def run_one(row):
    cid = row['candidate_id']
    out = LABEL_ROOT / cid
    out.mkdir(parents=True, exist_ok=True)
    geometry = Path(row['geometry_file'])
    geometry_sha = sha(geometry)
    atoms = read(str(geometry), format='extxyz')
    chash = canonical_hash(atoms)
    if geometry_sha != row['geometry_file_sha256'] or chash != row['canonical_geometry_sha256']:
        raise RuntimeError(f'{cid}: geometry hash does not match frozen manifest')
    kpts = row.get('dft_kpoints') or {
        'Ag_solid': (2, 2, 2), 'Ag_liquid': (2, 2, 2),
        'Ti3SiC2_solid_thermal': (3, 3, 1),
        'Ag_solid_liquid_coexistence': (3, 3, 1),
        'Ag_Ti3SiC2_interface': (4, 4, 1),
        'Ag_rich_Si_local_environment': (4, 4, 1),
    }[row['category']]
    kpts = tuple(int(x) for x in kpts)
    ecut = int(row.get('dft_ecut_eV', 350))
    log = out / 'gpaw_scf.txt'
    label_json = out / 'label.json'
    label_extxyz = out / 'label.extxyz'
    record = {
        'candidate_id': cid, 'category': row['category'], 'geometry_file': str(geometry),
        'geometry_sha256': geometry_sha, 'canonical_geometry_sha256': chash,
        'natoms': len(atoms), 'formula': atoms.get_chemical_formula(),
        'cell_A': np.asarray(atoms.cell).tolist(), 'pbc': np.asarray(atoms.pbc, bool).tolist(),
        'source_parent_id': row['parent_source_id'], 'random_seed': row['random_seed'],
        'target_temperature_K': row['target_temperature_K'], 'observed_temperature_K': row['observed_temperature_K'],
        'xc': 'PBE', 'code': 'GPAW', 'gpaw_version': '26.7.0', 'ase_version': '3.29.0',
        'paw_data_version': 'gpaw PAW dataset 1.2.1', 'mode': 'plane waves',
        'ecut_eV': ecut, 'kpts': list(kpts), 'symmetry': 'off', 'sigma_eV': 0.10,
        'scf_convergence': TOL, 'max_scf_iterations': 250,
        'training_energy_quantity': 'force-consistent finite-smearing free energy at sigma=0.10 eV',
        'extrapolated_energy_also_saved': True, 'stress_is_fit_target': False,
        'log': str(log), 'label_extxyz': str(label_extxyz),
        'run_status': 'running', 'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    }
    started = time.time()
    try:
        config = {
            'mode': PW(ecut), 'xc': 'PBE', 'kpts': kpts, 'symmetry': 'off',
            'occupations': FermiDirac(0.10), 'mixer': Mixer(beta=0.05, nmaxold=8, weight=100),
            'convergence': TOL, 'maxiter': 250,
            'parallel': {'domain': max(1, world.size)}, 'txt': str(log),
        }
        if not np.asarray(atoms.pbc, bool).all():
            config['poissonsolver'] = {'dipolelayer': 'xy'}
        atoms.calc = GPAW(**config)
        e_extrap = float(atoms.get_potential_energy(force_consistent=False))
        e_free = float(atoms.get_potential_energy(force_consistent=True))
        forces = np.asarray(atoms.get_forces(), dtype=float)
        scf_ok = scf_converged(atoms.calc, log)
        stress = None
        try:
            stress = np.asarray(atoms.get_stress(voigt=False), dtype=float).tolist()
        except Exception:
            pass
        finite = bool(np.isfinite(e_extrap) and np.isfinite(e_free) and np.isfinite(forces).all())
        record.update({
            'energy_extrapolated_eV': e_extrap,
            'energy_force_consistent_free_eV': e_free,
            'energy_fit_target_eV': e_free,
            'free_energy_eV_per_atom': e_free / len(atoms),
            'extrapolated_energy_eV_per_atom': e_extrap / len(atoms),
            'forces_eV_A': forces.tolist(),
            'force_rmse_component_eV_A': float(np.sqrt(np.mean(forces**2))),
            'force_rmse_vector_norm_eV_A': float(np.sqrt(np.mean(np.sum(forces**2, axis=1)))),
            'max_force_norm_eV_A': float(np.linalg.norm(forces, axis=1).max()),
            'stress_eV_A3': stress,
            'scf_converged': scf_ok, 'finite_energy_force': finite,
            'scf_iterations': int(atoms.calc.get_number_of_iterations()),
            'run_status': 'complete' if scf_ok and finite else 'scf_or_finite_check_failed',
        })
        if scf_ok and finite:
            atoms.info['energy_force_consistent_free_eV'] = e_free
            atoms.info['energy_extrapolated_eV'] = e_extrap
            atoms.info['energy_fit_target_kind'] = record['training_energy_quantity']
            atoms.info['dft_candidate_id'] = cid
            atoms.calc = SinglePointCalculator(atoms, energy=e_free, forces=forces,
                                               **({'stress': np.asarray(stress).reshape(-1)} if stress is not None else {}))
            write(str(label_extxyz), atoms, format='extxyz')
    except Exception as exc:
        record.update({'run_status': 'failed', 'error': f'{type(exc).__name__}: {exc}',
                       'traceback': traceback.format_exc()})
    record['elapsed_seconds'] = round(time.time() - started, 3)
    record['finished_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    label_json.write_text(json.dumps(record, indent=2) + '\n')
    if world.rank == 0:
        print(json.dumps({k: record.get(k) for k in ('candidate_id', 'run_status', 'scf_iterations',
              'elapsed_seconds', 'max_force_norm_eV_A', 'energy_fit_target_eV')}), flush=True)
    return record


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidate-id', required=True)
    args = ap.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    rows = [x for x in manifest['candidates'] if x['candidate_id'] == args.candidate_id]
    if len(rows) != 1:
        raise SystemExit(f"candidate id {args.candidate_id} not found exactly once")
    run_one(rows[0])


if __name__ == '__main__':
    main()
