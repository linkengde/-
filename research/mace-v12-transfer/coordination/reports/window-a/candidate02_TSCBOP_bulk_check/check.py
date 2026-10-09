import hashlib,json,tempfile
from pathlib import Path
import numpy as np
from ase.io import read,write
from ase.geometry import find_mic
from lammps import lammps
ROOT=Path(__file__).resolve().parent;REPO=ROOT.parents[5];BASE=REPO/'research/mace-v12-transfer';root=BASE/'periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations';potdir=BASE/'coordination/reports/window-b/current_heating_potential_handoff_20261009/potentials';pot=potdir/'TSCBOP.tersoff';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();rows=[]
for label in ['Ti3SiC2_baseline_k4x4x2_v19','Ti3SiC2_baseline_k6x6x2_v19']:
 p=root/label/(label+'_PW_PBE.extxyz');a=read(p);summary=json.loads((root/label/'summary.json').read_text());assert summary['scf_converged'] and set(a.get_chemical_symbols())=={'Ti','Si','C'}
 with tempfile.TemporaryDirectory() as td:
  data=Path(td)/'in.data';write(data,a,format='lammps-data',specorder=['Ti','Si','C'],atom_style='atomic',masses=True);l=lammps(cmdargs=['-log','none','-screen','none'])
  try:
   for c in ['units metal','atom_style atomic','boundary p p p','read_data '+str(data),'pair_style tersoff','pair_coeff * * '+str(pot)+' Ti Si C','run 0']:l.command(c)
   ids=l.numpy.extract_atom('id').copy()[:len(a)];order=np.argsort(ids);assert np.array_equal(ids[order],np.arange(1,len(a)+1));f=l.numpy.extract_atom('f').copy()[:len(a)][order];x=l.numpy.extract_atom('x').copy()[:len(a)][order];energy=float(l.get_thermo('pe'));assert find_mic(x-a.positions,a.cell,a.pbc)[1].max()<1e-10 and np.isfinite(f).all()
  finally:l.close()
 ref=a.arrays['PW_PBE_forces'];delta=f-ref;rows.append({'label':label,'input_sha256':sha(p),'TSCBOP_energy_eV':energy,'DFT_native_eV':summary['energy_eV_cell'],'DFT_free_eV':summary['free_energy_eV_cell'],'force_vector_rmse_eV_A':float(np.sqrt(np.mean(np.sum(delta*delta,axis=1)))),'max_atom_force_error_eV_A':float(np.linalg.norm(delta,axis=1).max()),'per_species_force_rmse_eV_A':{s:float(np.sqrt(np.mean(np.sum(delta[np.array(a.get_chemical_symbols())==s]**2,axis=1)))) for s in ['C','Si','Ti']}})
report={'status':'STATIC_TSCBOP_FIXED_BULK_DIAGNOSTIC','potential_sha256':sha(pot),'cases':rows,'limitations':['Equilibrium bulk only; no strained/thermal validation yet','Absolute energy zeros unmatched; no raw energy score','Tersoff result is one term in the mixed potential; interface data include cross interactions','Stage69 configuration unconfirmed; no fit or dynamics']};(ROOT/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(rows,indent=2))
