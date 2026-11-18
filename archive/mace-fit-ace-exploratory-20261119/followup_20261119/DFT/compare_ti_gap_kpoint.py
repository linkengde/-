#!/usr/bin/env python3
"""Compare the converged 1x6x3 GPAW label against 1x4x2 and Gamma."""
import json,hashlib
from pathlib import Path
import numpy as np
from ase.io import read
root=Path(__file__).resolve().parent
label=read(root/'ti_gap_1x6x3.extxyz')
coarse=read(root.parent/'ti_gap_1x4x2.extxyz')
ref=json.loads((root/'gamma_reference.json').read_text())
forces=np.asarray(label.arrays['REF_forces'],float)
gamma=np.asarray(ref['gamma_force_vectors_eV_A'],float)
delta=forces-gamma
coarse_forces=np.asarray(coarse.arrays['REF_forces'],float)
if (label.get_chemical_symbols()!=coarse.get_chemical_symbols() or
    not np.allclose(label.positions,coarse.positions,atol=1e-10,rtol=0) or
    not np.allclose(label.cell,coarse.cell,atol=1e-10,rtol=0)):
 raise ValueError('The 1x4x2 and 1x6x3 labels do not use identical geometry')
mesh_delta=forces-coarse_forces
energy=float(label.info['REF_energy'])
coarse_energy=float(coarse.info['REF_energy'])
result={
 'status':'SCF_CONVERGED',
 'natoms':len(label),
 'kpoints':[1,6,3],
 'energy_eV':energy,
 'energy_delta_vs_gamma_eV':energy-ref['gamma_energy_eV'],
 'energy_delta_vs_gamma_meV_atom':1000*(energy-ref['gamma_energy_eV'])/len(label),
 'energy_1x4x2_eV':coarse_energy,
 'energy_delta_vs_1x4x2_eV':energy-coarse_energy,
 'energy_delta_vs_1x4x2_meV_atom':1000*(energy-coarse_energy)/len(label),
 'force_component_rmse_vs_1x4x2_eV_A':float(np.sqrt(np.mean(mesh_delta**2))),
 'force_vector_rmse_vs_1x4x2_eV_A':float(np.sqrt(np.mean(np.sum(mesh_delta**2,axis=1)))),
 'force_vector_max_error_vs_1x4x2_eV_A':float(np.linalg.norm(mesh_delta,axis=1).max()),
 'force_component_rmse_vs_gamma_eV_A':float(np.sqrt(np.mean(delta**2))),
 'force_vector_rmse_vs_gamma_eV_A':float(np.sqrt(np.mean(np.sum(delta**2,axis=1)))),
 'force_vector_max_error_vs_gamma_eV_A':float(np.linalg.norm(delta,axis=1).max()),
 'absolute_force_vector_max_eV_A':float(np.linalg.norm(forces,axis=1).max()),
 'source_method':label.info['source_method'],
 'input_sha256':hashlib.sha256((root/'input.extxyz').read_bytes()).hexdigest(),
 'label_sha256':hashlib.sha256((root/'ti_gap_1x6x3.extxyz').read_bytes()).hexdigest(),
}
(root/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
