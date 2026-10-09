import json,hashlib
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parents[5]/'research/mace-v12-transfer/periodic_interface_v4'
entries=[('periodic','pbe_interface_v19_cif_convergence_pilot','AgSi_COD9009647_pilot_k6x6'),('slab_off','pbe_interface_v19_cif_slab_controls','AgSi_originalvac_slab_dipole_off_v19'),('slab_on','pbe_interface_v19_cif_slab_controls','AgSi_originalvac_slab_dipole_on_v19')];data={}
for tag,folder,label in entries:
 p=BASE/folder/'calculations'/label;s=json.loads((p/'summary.json').read_text());assert s['scf_converged'];a=read(p/(label+'_PW_PBE.extxyz'));data[tag]=(a,s)
 for line in (p/'SHA256SUMS.txt').read_text().splitlines():
  digest,name=line.split('  ',1);assert hashlib.sha256((p/name).read_bytes()).hexdigest()==digest
results=[]
for left,right in [('periodic','slab_off'),('slab_off','slab_on')]:
 a,s=data[left];b,t=data[right];assert a.get_chemical_symbols()==b.get_chemical_symbols() and np.allclose(a.positions,b.positions,rtol=0,atol=1e-10) and np.allclose(a.cell,b.cell,rtol=0,atol=1e-10)
 assert np.array_equal(a.arrays['lammps_id'],b.arrays['lammps_id']) and np.array_equal(a.arrays['central_pair'],b.arrays['central_pair'])
 for k in ['xc','cutoff_eV','smearing_eV','kpts']:assert s['method'][k]==t['method'][k],k
 if left=='periodic':assert a.pbc.tolist()==[True]*3 and b.pbc.tolist()==[True,True,False]
 else:assert np.array_equal(a.pbc,b.pbc) and s['method']['poissonsolver']=={} and t['method']['poissonsolver']=={'dipolelayer':'xy'}
 f=a.arrays['PW_PBE_forces'];g=b.arrays['PW_PBE_forces'];df=g-f;i,j=np.flatnonzero(a.arrays['central_pair'])
 if a[i].symbol!='Ag':i,j=j,i
 # Contact vector uses shared slab XY periodic boundary; separation lies inside slab.
 temp=a.copy();temp.set_pbc([True,True,False]);v=temp.get_distance(int(i),int(j),mic=True,vector=True);u=v/np.linalg.norm(v)
 energy=1000*(t['energy_eV_cell']-s['energy_eV_cell'])/len(a);free=1000*(t['free_energy_eV_cell']-s['free_energy_eV_cell'])/len(a);rms=float(np.sqrt(np.mean(np.sum(df*df,axis=1))));pair=float(np.dot(df[j]-df[i],u))
 results.append({'pair':[left,right],'native_energy_delta_meV_atom':energy,'free_energy_delta_meV_atom':free,'force_vector_rms_eV_A':rms,'maximum_atom_difference_eV_A':float(np.linalg.norm(df,axis=1).max()),'marked_pair_projection_delta_eV_A':pair,'checks':{'native_energy':abs(energy)<=2,'free_energy':abs(free)<=2,'vector':rms<=.01,'marked_pair':abs(pair)<=.02},'left_pbc':a.pbc.tolist(),'right_pbc':b.pbc.tolist(),'left_poisson':s['method'].get('poissonsolver',{}),'right_poisson':t['method'].get('poissonsolver',{})})
report={'status':'SAME_GEOMETRY_BOUNDARY_AND_DIPOLE_DIAGNOSTIC','comparisons':results,'limits':['Boundary control and dipole control separated; not kmesh repeats','Single original-vacuum slab only, not expanded slab/dipole convergence','No independent interface/model force pass','Native/free energies separate; finite-smearing convention unchanged']}
(ROOT/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
