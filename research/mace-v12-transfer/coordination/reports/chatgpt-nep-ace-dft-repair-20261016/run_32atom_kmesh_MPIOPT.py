#!/usr/bin/env python3
"""Fixed-geometry Ag10Ti12Si2C8 GPAW k-mesh pilot; dry-run by default.

This code is NOT a DFT result. Real SCF requires --execute and validated PAW.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'Ti_gap_2p3_frozen_32atoms_UNLABELED.extxyz'
EXPECTED_SHA = '575686f5c9749c23c2258e7491633cd265367bad00432d0b9078feb89fe0f6b8'
KPOINTS = {'1x4x2': (1, 4, 2), '1x6x4': (1, 6, 4), '1x8x6': (1, 8, 6)}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kmesh', choices=KPOINTS, required=True)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--paw-manifest', type=Path, help='JSON keyed Ag/Ti/Si/C, each with path and sha256')
    args = parser.parse_args()

    source_sha = sha256(INPUT)
    if source_sha != EXPECTED_SHA:
        raise SystemExit('FROZEN_INPUT_SHA_MISMATCH: stop, no calculations')
    lines = INPUT.read_text().splitlines()
    if int(lines[0]) != 32 or 'pbc="T T T"' not in lines[1] or len(lines) != 34:
        raise SystemExit('FROZEN_STRUCTURE_INVALID')
    if not args.execute:
        print(json.dumps({'status': 'INPUT_PREFLIGHT_PASS', 'SCF_STARTED': False,
                          'mesh': KPOINTS[args.kmesh], 'atoms': 32,
                          'input_sha256': source_sha}, indent=2))
        return
    if args.output is None or args.paw_manifest is None:
        raise SystemExit('EXECUTION_BLOCKED: require --output and --paw-manifest')
    output = args.output.resolve()
    if output.exists():
        raise SystemExit('EXECUTION_BLOCKED: refusing to overwrite existing output')

    import gpaw_data  # only on the verified compute node
    actual_setup_root = Path(gpaw_data.datapath()).resolve()
    manifest = json.loads(args.paw_manifest.read_text())
    checked = {}
    for symbol in ('Ag', 'Ti', 'Si', 'C'):
        item = manifest.get(symbol, {})
        p = Path(item.get('path', '/nonexistent')).resolve()
        expected = str(item.get('sha256', '')).lower()
        if not re.fullmatch(r'[0-9a-f]{64}', expected):
            raise SystemExit(f'PAW_MANIFEST_INVALID: {symbol}')
        if p != actual_setup_root / (symbol + '.PBE.gz') or not p.is_file():
            raise SystemExit(f'PAW_ACTUAL_RUNTIME_PATH_MISMATCH: {symbol}')
        if sha256(p) != expected:
            raise SystemExit(f'PAW_SHA_MISMATCH: {symbol}')
        checked[symbol] = {'path': str(p), 'sha256': expected}
    os.environ['GPAW_SETUP_PATH'] = str(actual_setup_root)

    import ase
    import gpaw
    import numpy as np
    from ase.io import read
    from gpaw import GPAW, PW, FermiDirac, Mixer
    from gpaw.mpi import world
    if not ase.__version__.startswith('3.29.') or not gpaw.__version__.startswith('26.7.'):
        raise SystemExit(f'VERSION_MISMATCH: ASE {ase.__version__}, GPAW {gpaw.__version__}')
    if world.size != 4:
        raise SystemExit(f'MPI_RANK_MISMATCH: expected 4, got {world.size}')
    atoms = read(str(INPUT), format='extxyz')
    if len(atoms) != 32 or not np.all(atoms.pbc):
        raise SystemExit('ASE_ATOMIC_MAPPING_OR_PBC_FAILED')
    if {s: atoms.symbols.count(s) for s in ['Ag','Ti','Si','C']} != {'Ag':10,'Ti':12,'Si':2,'C':8}:
        raise SystemExit('ASE_COMPOSITION_FAILED')
    if world.rank == 0:
        output.mkdir(parents=True, exist_ok=False)
    world.barrier()
    calc = GPAW(mode=PW(500), xc='PBE', kpts={'size': KPOINTS[args.kmesh], 'gamma': True},
                occupations=FermiDirac(0.10), mixer=Mixer(beta=0.05, nmaxold=8, weight=100),
                convergence={'energy': 1e-5, 'density': 1e-5, 'eigenstates': 1e-5},
                maxiter=160, spinpol=False, txt=str(output / 'gpaw.log'))
    atoms.calc = calc
    energy = float(atoms.get_potential_energy(force_consistent=False))
    free_energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = np.asarray(atoms.get_forces(), dtype=float)
    if not np.isfinite([energy, free_energy]).all() or not np.isfinite(forces).all():
        raise SystemExit('NONFINITE_DFT_LABEL')
    if world.rank == 0:
        result = {'status': 'SCF_COMPLETED_NOT_KPOINT_VALIDATED', 'atom_count': 32,
                  'kmesh': KPOINTS[args.kmesh], 'input_sha256': source_sha,
                  'native_energy_eV': energy, 'force_consistent_free_energy_eV': free_energy,
                  'forces_eV_A': forces.tolist(), 'max_force_eV_A': float(np.linalg.norm(forces, axis=1).max()),
                  'GPAW_version': gpaw.__version__, 'ASE_version': ase.__version__,
                  'mpi_ranks': world.size, 'paw_setup_files': checked,
                  'method': 'PBE PW500 Fermi0.10 gamma-centered spin-unpolarized fixed geometry'}
        (output / 'RESULT.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
        print(json.dumps({'status': result['status'], 'directory': str(output)}), flush=True)
    world.barrier()


if __name__ == '__main__':
    main()
