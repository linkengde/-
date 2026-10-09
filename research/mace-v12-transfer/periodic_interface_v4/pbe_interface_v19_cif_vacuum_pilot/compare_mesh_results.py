"""Compare verified identical-geometry pilot archives; no model inference."""
import hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
from ase.io import read
from ase.geometry import find_mic
ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'input_manifest.json').read_text())
BASELINE=ROOT.parent/'pbe_interface_v19_cif_convergence_pilot'
assert hashlib.sha256((BASELINE/'input_manifest.json').read_bytes()).hexdigest()==manifest['baseline_manifest_sha256']
new_records=manifest['records']
baseline_records=json.loads((BASELINE/'input_manifest.json').read_text())['records']
WIDTH005=ROOT.parent/'pbe_interface_v19_cif_smearing_pilot'
assert hashlib.sha256((WIDTH005/'input_manifest.json').read_bytes()).hexdigest()==manifest['width005_manifest_sha256']
width005_records=json.loads((WIDTH005/'input_manifest.json').read_text())['records']
all_records=new_records+baseline_records+width005_records
def record_root(r):return ROOT if r in new_records else (WIDTH005 if r in width005_records else BASELINE)
selected=sys.argv[1:] or ['AgSi_COD9009647_pilot_gamma','AgSi_COD9009647_pilot_k2x2']
assert len(selected)==2 and selected[0]!=selected[1]
records=[next(r for r in all_records if r['label']==name) for name in selected];loaded=[]
for r in records:
 d=record_root(r)/'calculations'/r['label']
 if not (d/'summary.json').exists():raise SystemExit('BLOCKED: both verified mesh results are required')
 subprocess.run([sys.executable,str(record_root(r)/'verify_result.py'),r['label'],str(d)],check=True,stdout=subprocess.DEVNULL)
 s=json.loads((d/'summary.json').read_text());a=read(d/(r['label']+'_PW_PBE.extxyz'));loaded.append((s,a));assert r['role']=='numerical_convergence_pilot_not_training_or_validation_label'
(s0,a0),(s1,a1)=loaded
assert records[0]['label']=='AgSi_COD9009647_pilot_k6x6' and records[1]['label']=='AgSi_COD9009647_vacuum26_k6x6_sigma0p10'
assert np.array_equal(a0.numbers,a1.numbers) and np.array_equal(a0.arrays['lammps_id'],a1.arrays['lammps_id']) and np.array_equal(a0.arrays['central_pair'],a1.arrays['central_pair']) and np.array_equal(a0.pbc,a1.pbc)
assert np.allclose(a1.positions-a0.positions,[0,0,5],rtol=0,atol=1e-8) and np.allclose(a1.cell.array-a0.cell.array,[[0,0,0],[0,0,0],[0,0,10]],rtol=0,atol=1e-8)
assert s0['method']==s1['method'], 'Only vacuum geometry may change'
f0=a0.arrays['PW_PBE_forces'];f1=a1.arrays['PW_PBE_forces'];err=f1-f0;symbols=np.array(a0.get_chemical_symbols());marked=a0.arrays['central_pair'].astype(bool);ia=np.flatnonzero(marked&(symbols=='Ag'));ix=np.flatnonzero(marked&(symbols=='Si'));assert len(ia)==len(ix)==1;ia,ix=int(ia[0]),int(ix[0]);v=find_mic(a0.positions[ix]-a0.positions[ia],a0.cell,a0.pbc)[0];unit=v/np.linalg.norm(v)
metrics={'native_energy_delta_meV_atom':(s1['energy_eV_cell']-s0['energy_eV_cell'])*1000/len(a0),'free_energy_delta_meV_atom':(s1['free_energy_eV_cell']-s0['free_energy_eV_cell'])*1000/len(a0),'force_vector_difference_RMSE_eV_A':float(np.sqrt(np.mean(np.sum(err**2,axis=1)))),'maximum_atom_force_difference_eV_A':float(np.linalg.norm(err,axis=1).max()),'pair_signed_projection_delta_eV_A':float(np.dot(err[ix]-err[ia],unit))}
passes={'energy':abs(metrics['native_energy_delta_meV_atom'])<=2,'vector':metrics['force_vector_difference_RMSE_eV_A']<=.01,'pair':abs(metrics['pair_signed_projection_delta_eV_A'])<=.02}
result={'control':'expanded vacuum; fixed relative ionic geometry/method','smearing_eV':[s0['method']['smearing_eV'],s1['method']['smearing_eV']],'meshes':[r['kpts'] for r in records],'labels':selected,'status':'PILOT_WITHIN_PROVISIONAL_BUDGET' if all(passes.values()) else 'PILOT_MESH_SENSITIVITY_EXCEEDS_BUDGET','metrics':metrics,'checks':passes,'input_sha256':records[0]['input_sha256'],'result_files_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for r in records for p in (record_root(r)/'calculations'/r['label']).iterdir() if p.is_file() and p.suffix!='.gpw'},'limitations':['Two mesh points do not establish mesh convergence.','Vacuum/dipole/slab/coherent-strain checks remain separate.','No MACE error conclusion or dataset label integration from numerical pilots.']}
(ROOT/('mesh_comparison_'+selected[0]+'__'+selected[1]+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
