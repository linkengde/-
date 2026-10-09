import hashlib,json,tempfile
from pathlib import Path
import numpy as np
from ase.io import read,write
from ase.geometry import find_mic
from lammps import lammps
ROOT=Path(__file__).resolve().parent;REPO=ROOT.parents[5];base=REPO/'research/mace-v12-transfer';listing=base/'coordination/reports/window-b/v19_Ag_symmetric_displacement_volume_controls/eam_static_input_list.json';pot=base/'coordination/reports/window-b/current_heating_potential_handoff_20261009/potentials/Ag_u3.eam';rows=[]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for item in json.loads(listing.read_text())['inputs']:
 src=REPO/item['path'];assert sha(src)==item['sha256'];a=read(src);assert set(a.get_chemical_symbols())=={'Ag'} and a.pbc.all();ref=a.arrays['PW_PBE_forces'];summary=json.loads((src.parent/'summary.json').read_text());assert summary['scf_converged']
 with tempfile.TemporaryDirectory() as tmp:
  data=Path(tmp)/'in.data';write(data,a,format='lammps-data',specorder=['Ag'],atom_style='atomic',masses=True)
  l=lammps(cmdargs=['-log','none','-screen','none'])
  try:
   for cmd in ['units metal','atom_style atomic','boundary p p p','read_data '+str(data),'pair_style eam','pair_coeff 1 1 '+str(pot),'run 0']:l.command(cmd)
   ids=l.numpy.extract_atom('id').copy()[:len(a)];order=np.argsort(ids);assert np.array_equal(ids[order],np.arange(1,len(a)+1));f=l.numpy.extract_atom('f').copy()[:len(a)][order];x=l.numpy.extract_atom('x').copy()[:len(a)][order];energy=float(l.get_thermo('pe'));assert find_mic(x-a.positions,a.cell,a.pbc)[1].max()<1e-10 and np.isfinite(f).all() and np.isfinite(energy)
  finally:l.close()
 rows.append({'label':item['label'],'input_sha256':item['sha256'],'volume_A3':float(a.get_volume()),'EAM_energy_eV':energy,'DFT_native_eV':summary['energy_eV_cell'],'DFT_free_eV':summary['free_energy_eV_cell'],'force_vector_rmse_eV_A':float(np.sqrt(np.mean(np.sum((f-ref)**2,axis=1)))),'max_atom_force_error_eV_A':float(np.linalg.norm(f-ref,axis=1).max()),'atom1_EAM_force_eV_A':f[0].tolist(),'atom1_DFT_force_eV_A':ref[0].tolist(),'EAM_forces_eV_A':f.tolist()})
b=next(r for r in rows if r['label']=='Ag_baseline_k6x6x6_v19')
for r in rows:
 for field in ['EAM_energy_eV','DFT_native_eV','DFT_free_eV']:r[field+'_delta_meV_atom']=1000*(r[field]-b[field])/32
 r['relative_native_energy_error_meV_atom']=r['EAM_energy_eV_delta_meV_atom']-r['DFT_native_eV_delta_meV_atom'];r['relative_free_energy_error_meV_atom']=r['EAM_energy_eV_delta_meV_atom']-r['DFT_free_eV_delta_meV_atom']
report={'status':'UNCHANGED_EAM_LOCAL_STATIC_DIAGNOSTIC','EAM_sha256':sha(pot),'cases':rows,'limitations':['Only five related near-equilibrium Ag controls; no independent validation or MD','Absolute energy zeros not compared; baseline-subtracted energies only','Finite-width native/free conventions separate','Volume k6 numerical convergence provisional','Single sampled volume trend is not optimized lattice or300K stability','No interface or separation-force gate implied']}
(ROOT/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps([{k:v for k,v in r.items() if k!='EAM_forces_eV_A'} for r in rows],indent=2))
