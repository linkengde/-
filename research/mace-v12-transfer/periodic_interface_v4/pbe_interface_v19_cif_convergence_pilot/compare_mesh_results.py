"""Compare verified identical-geometry pilot archives; no model inference."""
import hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
from ase.io import read
from ase.geometry import find_mic
ROOT=Path(__file__).resolve().parent
records=json.loads((ROOT/'input_manifest.json').read_text())['records'];loaded=[]
for r in records:
 d=ROOT/'calculations'/r['label']
 if not (d/'summary.json').exists():raise SystemExit('BLOCKED: both verified mesh results are required')
 subprocess.run([sys.executable,str(ROOT/'verify_result.py'),r['label'],str(d)],check=True,stdout=subprocess.DEVNULL)
 s=json.loads((d/'summary.json').read_text());a=read(d/(r['label']+'_PW_PBE.extxyz'));loaded.append((s,a));assert r['role']=='numerical_convergence_pilot_not_training_or_validation_label'
(s0,a0),(s1,a1)=loaded
assert np.array_equal(a0.positions,a1.positions) and np.array_equal(a0.cell,a1.cell) and np.array_equal(a0.numbers,a1.numbers) and np.array_equal(a0.arrays['lammps_id'],a1.arrays['lammps_id'])
f0=a0.arrays['PW_PBE_forces'];f1=a1.arrays['PW_PBE_forces'];err=f1-f0;symbols=np.array(a0.get_chemical_symbols());marked=a0.arrays['central_pair'].astype(bool);ia=np.flatnonzero(marked&(symbols=='Ag'));ix=np.flatnonzero(marked&(symbols=='Si'));assert len(ia)==len(ix)==1;ia,ix=int(ia[0]),int(ix[0]);v=find_mic(a0.positions[ix]-a0.positions[ia],a0.cell,a0.pbc)[0];unit=v/np.linalg.norm(v)
metrics={'native_energy_delta_meV_atom':(s1['energy_eV_cell']-s0['energy_eV_cell'])*1000/len(a0),'free_energy_delta_meV_atom':(s1['free_energy_eV_cell']-s0['free_energy_eV_cell'])*1000/len(a0),'force_vector_difference_RMSE_eV_A':float(np.sqrt(np.mean(np.sum(err**2,axis=1)))),'maximum_atom_force_difference_eV_A':float(np.linalg.norm(err,axis=1).max()),'pair_signed_projection_delta_eV_A':float(np.dot(err[ix]-err[ia],unit))}
passes={'energy':abs(metrics['native_energy_delta_meV_atom'])<=2,'vector':metrics['force_vector_difference_RMSE_eV_A']<=.01,'pair':abs(metrics['pair_signed_projection_delta_eV_A'])<=.02}
result={'status':'PILOT_WITHIN_PROVISIONAL_BUDGET' if all(passes.values()) else 'PILOT_MESH_SENSITIVITY_EXCEEDS_BUDGET','metrics':metrics,'checks':passes,'input_sha256':records[0]['input_sha256'],'result_files_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for r in records for p in (ROOT/'calculations'/r['label']).iterdir() if p.is_file() and p.suffix!='.gpw'},'limitations':['Two mesh points do not establish mesh convergence.','Vacuum/dipole/slab/coherent-strain checks remain separate.','No MACE error conclusion or dataset label integration from numerical pilots.']}
(ROOT/'mesh_comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
