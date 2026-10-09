import hashlib,json,tempfile
from pathlib import Path
import numpy as np
from ase.io import read,write
from ase.geometry import find_mic
from lammps import lammps
ROOT=Path(__file__).resolve().parent;REPO=ROOT.parents[5];BASE=REPO/'research/mace-v12-transfer';source=BASE/'coordination/reports/window-b/v19_Ag_shear_k6_k8_controls/report.json';report=json.loads(source.read_text());pot=BASE/'coordination/reports/window-b/current_heating_potential_handoff_20261009/potentials/Ag_u3.eam';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();rows=[]
items=[{'label':'Ag_baseline_k8x8x8_v19','path':'research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations/Ag_baseline_k8x8x8_v19/Ag_baseline_k8x8x8_v19_PW_PBE.extxyz'}]
items += [{'label':x['label'],'path':x['path'],'sha256':x['sha256']} for x in report['EAM_static_input_list'] if x['label']!='Ag_baseline_k8x8x8_v19']
for item in items:
 p=REPO/item['path'];assert not item.get('sha256') or sha(p)==item['sha256'];a=read(p);assert set(a.get_chemical_symbols())=={'Ag'} and a.pbc.all();ref=a.arrays['PW_PBE_forces']
 with tempfile.TemporaryDirectory() as t:
  data=Path(t)/'in.data';write(data,a,format='lammps-data',specorder=['Ag'],atom_style='atomic',masses=True);l=lammps(cmdargs=['-log','none','-screen','none'])
  try:
   for cmd in ['units metal','atom_style atomic','boundary p p p','read_data '+str(data),'pair_style eam','pair_coeff 1 1 '+str(pot),'run 0']:l.command(cmd)
   ids=l.numpy.extract_atom('id').copy()[:len(a)];order=np.argsort(ids);assert np.array_equal(ids[order],np.arange(1,len(a)+1));f=l.numpy.extract_atom('f').copy()[:len(a)][order];x=l.numpy.extract_atom('x').copy()[:len(a)][order];assert find_mic(x-a.positions,a.cell,a.pbc)[1].max()<1e-10 and np.isfinite(f).all();energy=float(l.get_thermo('pe'))
  finally:l.close()
 eam_energy=float(l.get_thermo('pe')) if False else None
 delta=f-ref;rows.append({'label':item['label'],'input_sha256':sha(p),'EAM_energy_eV':energy,'force_vector_rmse_eV_A':float(np.sqrt(np.mean(np.sum(delta*delta,axis=1)))),'max_atom_force_error_eV_A':float(np.linalg.norm(delta,axis=1).max()),'DFT_energy_eV':a.info['PW_PBE_energy_eV'],'EAM_forces_eV_A':f.tolist()})
base=rows[0]
for row in rows:
 row['EAM_delta_from_baseline_meV_atom']=1000*(row['EAM_energy_eV']-base['EAM_energy_eV'])/32
 row['DFT_native_delta_from_baseline_meV_atom']=1000*(row['DFT_energy_eV']-base['DFT_energy_eV'])/32
 row['EAM_minus_DFT_native_response_meV_atom']=row['EAM_delta_from_baseline_meV_atom']-row['DFT_native_delta_from_baseline_meV_atom']
report_out={'status':'STATIC_AG_EAM_DIAGNOSTIC','EAM_sha256':sha(pot),'cases':rows,'limits':['Raw potential and DFT energies use different zeros and are not scored','Static pure Ag shear only; not interface validation or thermal/elastic-stability proof','No fitting, dynamics or parameter changes']};(ROOT/'report.json').write_text(json.dumps(report_out,indent=2)+'\n');print(json.dumps([{k:v for k,v in x.items() if k!='EAM_forces_eV_A'} for x in rows],indent=2))
