"""Read-only independent reproduction of the three published mesh comparisons."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from ase.io import read
from ase.geometry import find_mic

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[3] / 'periodic_interface_v4/pbe_interface_v19_cif_convergence_pilot'
loaded = {}
for k in (4, 5, 6):
    label = f'AgSi_COD9009647_pilot_k{k}x{k}'
    folder = ROOT / 'calculations' / label
    subprocess.run([sys.executable, str(ROOT/'verify_result.py'), label, str(folder)], check=True, stdout=subprocess.DEVNULL)
    for line in (folder/'SHA256SUMS.txt').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        assert hashlib.sha256((folder/name).read_bytes()).hexdigest() == digest
    loaded[k] = (json.loads((folder/'summary.json').read_text()), read(folder/(label+'_PW_PBE.extxyz')))
results = []
for k0, k1 in ((4, 5), (4, 6), (5, 6)):
    s0, a0 = loaded[k0]; s1, a1 = loaded[k1]
    for x, y in ((a0.positions,a1.positions), (a0.cell.array,a1.cell.array), (a0.numbers,a1.numbers), (a0.pbc,a1.pbc), (a0.arrays['lammps_id'],a1.arrays['lammps_id']), (a0.arrays['central_pair'],a1.arrays['central_pair'])):
        assert np.array_equal(x, y)
    err = a1.arrays['PW_PBE_forces'] - a0.arrays['PW_PBE_forces']
    symbols = np.array(a0.get_chemical_symbols()); marked = a0.arrays['central_pair'].astype(bool)
    ag = np.flatnonzero(marked & (symbols=='Ag')); si = np.flatnonzero(marked & (symbols=='Si'))
    assert len(ag)==len(si)==1
    ag, si = int(ag[0]), int(si[0])
    v = find_mic(a0.positions[si]-a0.positions[ag], a0.cell, a0.pbc)[0]
    metrics = {
        'native_energy_delta_meV_atom': (s1['energy_eV_cell']-s0['energy_eV_cell'])*1000/len(a0),
        'free_energy_delta_meV_atom': (s1['free_energy_eV_cell']-s0['free_energy_eV_cell'])*1000/len(a0),
        'force_vector_difference_RMSE_eV_A': float(np.sqrt(np.mean(np.sum(err**2,axis=1)))),
        'maximum_atom_force_difference_eV_A': float(np.linalg.norm(err,axis=1).max()),
        'pair_signed_projection_delta_eV_A': float(np.dot(err[si]-err[ag],v/np.linalg.norm(v))),
    }
    published = ROOT/f'mesh_comparison_AgSi_COD9009647_pilot_k{k0}x{k0}__AgSi_COD9009647_pilot_k{k1}x{k1}.json'
    reference = json.loads(published.read_text())
    for key, value in metrics.items():
        assert np.isclose(value, reference['metrics'][key], rtol=1e-12, atol=1e-12), key
    results.append({'meshes':[k0,k1], 'metrics':metrics, 'reproduces_published':True, 'published_sha256':hashlib.sha256(published.read_bytes()).hexdigest()})
result = {'status':'REPRODUCED_FORCE_GATE_FAIL', 'archive_inventories_verified':[4,5,6], 'comparisons':results, 'limits':{'energy_meV_atom':2, 'vector_eV_A':0.01, 'pair_eV_A':0.02}, 'limitations':['Numerical pilots only; no model accuracy or dataset integration claim.', 'Fixed width 0.1 eV does not isolate smearing causality.', 'No new DFT or sealed-label access.']}
(OUT/'review.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
