#!/workspace/venvs/dft-labels/bin/python
"""DFT single point for new training candidate similar to reserved SI01."""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
from ase.io import read, write
from gpaw import FermiDirac, GPAW, Mixer, PW
from gpaw.mpi import world

HERE = Path(__file__).resolve().parent
source = HERE / 'si_hotspot_04_candidate.extxyz'
outdir = HERE / 'labels'
outdir.mkdir(exist_ok=True)
atoms = read(source)
center_id = int(atoms.info['center_atom_id'])
ids = np.asarray(atoms.arrays['source_atom_id'], dtype=int)
matched = np.flatnonzero(ids == center_id)
if len(matched) != 1:
    raise RuntimeError(f'Expected one center atom ID {center_id}; found {len(matched)}')
ci = int(matched[0])
start = time.time()
stem = 'si_hotspot_04_pw350_gamma'
log = outdir / f'{stem}.txt'
atoms.calc = GPAW(mode=PW(350), xc='PBE', kpts=(1, 1, 1), symmetry='off',
                  occupations=FermiDirac(0.10),
                  mixer=Mixer(beta=0.05, nmaxold=8, weight=100),
                  convergence={'density': 1e-5, 'eigenstates': 1e-5, 'energy': 1e-5},
                  maxiter=160, parallel={'domain': max(1, world.size)}, txt=str(log))
energy = float(atoms.get_potential_energy())
free_energy = float(atoms.get_potential_energy(force_consistent=True))
forces = np.asarray(atoms.get_forces(), dtype=float)
iterations = int(atoms.calc.get_number_of_iterations())
if not np.isfinite(forces).all() or not np.isfinite(energy):
    raise RuntimeError('DFT returned a non-finite energy or force')
ref_center = forces[ci]
pod_center = np.asarray(atoms.info['center_POD_v11_force_eV_A'], dtype=float)
atoms.info.update({
    'REF_energy': energy,
    'REF_energy_unshifted_PBE_PW350_Gamma_eV': energy,
    'free_energy_eV': free_energy,
    'label_status': 'GPAW PBE/PW350 Gamma single-point label for new Si training candidate',
    'DFT_method': 'GPAW PBE/PW350 Gamma, Fermi-Dirac sigma 0.10 eV',
    'DFT_cutoff_eV': 350,
    'DFT_kpts': [1, 1, 1],
    'DFT_center_force_eV_A': ref_center.tolist(),
    'DFT_center_force_norm_eV_A': float(np.linalg.norm(ref_center)),
    'DFT_minus_POD_v11_center_force_eV_A': (ref_center - pod_center).tolist(),
    'DFT_minus_POD_v11_center_force_error_eV_A': float(np.linalg.norm(ref_center - pod_center)),
    'DFT_edge_force_note': 'Only the center Si force is used in fitting; local crop edge forces are excluded.',
    'DFT_energy_note': 'Finite-crop energy has zero fit weight.',
    'model_data_role': 'training_candidate; SI01 center atom 104359 remains excluded from training',
})
atoms.arrays['REF_forces'] = forces
atoms.calc = None
xyz = outdir / f'{stem}.extxyz'
write(xyz, atoms, format='extxyz')
record = {
    'status': 'COMPUTED_AUDIT_PENDING', 'index': 4, 'center_atom_id': center_id,
    'center_species': 'Si', 'natoms': len(atoms), 'formula': atoms.get_chemical_formula(),
    'energy_eV': energy, 'free_energy_eV': free_energy,
    'center_DFT_force_eV_A': ref_center.tolist(),
    'center_DFT_force_norm_eV_A': float(np.linalg.norm(ref_center)),
    'center_POD_v11_force_eV_A': pod_center.tolist(),
    'center_force_error_eV_A': float(np.linalg.norm(ref_center - pod_center)),
    'max_DFT_force_norm_eV_A': float(np.linalg.norm(forces, axis=1).max()),
    'iterations': iterations, 'elapsed_s': time.time() - start,
    'role': 'new training candidate; SI01 remains development holdout',
    'output': str(xyz), 'log': str(log),
}
(outdir / f'{stem}.json').write_text(json.dumps(record, indent=2) + '\n')
if world.rank == 0:
    print(json.dumps(record), flush=True)
