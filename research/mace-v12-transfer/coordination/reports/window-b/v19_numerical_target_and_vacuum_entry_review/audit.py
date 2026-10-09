"""Verify numerical archives and prepare disabled vacuum controls; no GPAW calculation."""
import hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
from ase.io import read,write
OUT=Path(__file__).resolve().parent
BASE=OUT.parents[3]/'periodic_interface_v4'
GPAW=Path('/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
records=[]
for dirname,width in [('pbe_interface_v19_cif_convergence_pilot',0.10),('pbe_interface_v19_cif_smearing_pilot',0.05),('pbe_interface_v19_cif_smearing020_pilot',0.20)]:
 root=BASE/dirname
 for mesh in [5,6]:
  label=f'AgSi_COD9009647_pilot_k{mesh}x{mesh}'+('' if width==.10 else f'_sigma0p{int(width*100):02}')
  folder=root/'calculations'/label
  v=json.loads(subprocess.check_output([sys.executable,str(root/'verify_result.py'),label,str(folder)]))
  assert v['status']=='PASS'
  for line in (folder/'SHA256SUMS.txt').read_text().splitlines():
   h,n=line.split(maxsplit=1);assert sha(folder/n)==h
  s=json.loads((folder/'summary.json').read_text())
  assert s['method']['smearing_eV']==width
  records.append({'label':label,'mesh':mesh,'width_eV':width,'summary':s,'inventory_sha256':sha(folder/'SHA256SUMS.txt')})
# Reproduce available comparison results using read-only evaluation of production comparator.
root=BASE/'pbe_interface_v19_cif_smearing020_pilot'
comparisons=[]
for m in [5,6]:
 b=f'AgSi_COD9009647_pilot_k{m}x{m}';pairs=[(b,b+'_sigma0p05'),(b,b+'_sigma0p20'),(b+'_sigma0p05',b+'_sigma0p20')]
 for a,b in pairs:
  source=(root/'compare_mesh_results.py').read_text()
  source=source[:source.index("(ROOT/('mesh_comparison_")]
  old=sys.argv;sys.argv=['compare',a,b]
  scope={'__file__':str(root/'compare_mesh_results.py')};exec(compile(source,str(root/'compare_mesh_results.py'),'exec'),scope);sys.argv=old
  comparisons.append(scope['result'])
for suffix in ['_sigma0p05','_sigma0p20']:
 source=(root/'compare_mesh_results.py').read_text();source=source[:source.index("(ROOT/('mesh_comparison_")]
 old=sys.argv;sys.argv=['compare','AgSi_COD9009647_pilot_k5x5'+suffix,'AgSi_COD9009647_pilot_k6x6'+suffix]
 scope={'__file__':str(root/'compare_mesh_results.py')};exec(compile(source,str(root/'compare_mesh_results.py'),'exec'),scope);sys.argv=old;comparisons.append(scope['result'])
input_path=BASE/'pbe_interface_v19_cif_convergence_pilot/inputs/AgSi_cod9009647_slab_v19_dev_proposal.extxyz'
a=read(input_path);assert np.all(a.pbc) and np.allclose(a.cell.array[2,:2],0) and np.allclose(a.cell.array[:2,2],0)
zmin,zmax=a.positions[:,2].min(),a.positions[:,2].max();height=float(a.cell[2,2]);span=float(zmax-zmin)
inputs=OUT/'inputs';inputs.mkdir(exist_ok=True)
proposals=[]
for increment in [0.,10.]:
 for boundary,dipole in [(True,False),(False,False),(False,True)]:
  b=a.copy();b.cell[2,2]=height+increment;b.positions[:,2]+=increment/2
  if not boundary:b.pbc=[True,True,False]
  assert np.max(np.abs((b.positions-b.positions[0])-(a.positions-a.positions[0])))<1e-12
  for name in ['lammps_id','central_pair']:assert np.array_equal(b.arrays[name],a.arrays[name])
  label=f'AgSi_vacuum_plus{int(increment)}_periodicZ{int(boundary)}_dipole{int(dipole)}_PROPOSAL'
  file=inputs/(label+'.extxyz');write(file,b,format='extxyz');c=read(file)
  assert np.max(np.abs(c.positions-b.positions))<1e-8
  assert np.max(np.abs(c.cell.array-b.cell.array))<1e-8
  proposals.append({'label':label,'input':str(file.relative_to(OUT)),'sha256':sha(file),'launch_enabled':False,'owner_instance':None,'method':{'xc':'PBE','cutoff_eV':500,'kpts':None,'smearing_eV':None,'poissonsolver':{'dipolelayer':'xy'} if dipole else None},'pbc':b.pbc.tolist(),'vacuum_gap_A':height+increment-span,'lower_margin_A':float(b.positions[:,2].min()),'upper_margin_A':float(b.cell[2,2]-b.positions[:,2].max()),'cell_height_A':height+increment,'translation_z_A':increment/2,'blocked_reason':'Await fixed-width/mesh numerical decision after A6/8, exact owner registration and reviewed per-record runner/verifier. Dipole toggle requires explicit pbc_z=False control; original all-periodic input is incompatible with this API.'})
manifest={'launch_enabled':False,'source_sha256':sha(input_path),'original_vacuum_gap_A':height-span,'ionic_span_A':span,'serialization_position_tolerance_A':1e-8,'relative_geometry_tolerance_A':1e-12,'records':proposals,'note':'No relaxation. Width and mesh deliberately unresolved, not silently inherited or selected by pass. Paired controls separate vacuum change from dipole/boundary change.'}
(OUT/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
sources=['new/ase_interface.py','new/calculation.py','new/pw/poisson.py','new/energies.py','occupations.py','old/calculator.py','old/pw/hamiltonian.py','dipole_correction.py']
result={'status':'REVIEW_COMPLETE_VACUUM_LAUNCH_BLOCKED','archives_verified':records,'comparisons':comparisons,'vacuum_manifest':manifest,'installed_source_evidence':{name:sha(GPAW/name) for name in sources},'conclusions':['SCF pass does not imply numerical mesh convergence.','0.20 eV 5/6 vector0.0099259526 passes locally; finite-width smoothing is not proof of zero-temperature accuracy.','Native extrapolated energy and force-consistent free energy are distinct; changing energy key does not change returned forces.','Preserve legacy targets as historical; choose future common width/mesh/energy/boundary policy explicitly and relabel incompatible rows before mixing.']}
(OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n')
print('Verified six archives, eight comparisons and six disabled vacuum/dipole proposals; original gap',height-span)
