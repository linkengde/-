"""NONPHYSICAL fixtures exercising fail-closed metadata/identity checks; no DFT."""
import copy,json,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
from ase.io import write
from entry import ROOT,load,validate,sha
results=[]
m=json.loads((ROOT/'input_manifest.json').read_text());r=m['records'][0]
# Construct no calculator; all static inputs validate.
for rec in m['records']:validate(m,rec)
for command in [[sys.executable,str(ROOT/'run_queue.py'),r['label']],[sys.executable,str(ROOT/'run_pw_reference.py'),r['label'],str(ROOT/r['input']),'/tmp/forbidden-dft-draft-output']]:
 p=subprocess.run(command,capture_output=True,text=True);assert p.returncode!=0 and 'Launch disabled' in p.stderr;results.append('disabled '+Path(command[1]).name+' rejected')
for field,value in [('pbc',[True,True,True]),('poissonsolver',{'dipolelayer':'xz'}),('input_sha256','0'*64),('ordered_ids',[0]*len(r['ordered_ids']))]:
 wrong=copy.deepcopy(r);wrong[field]=value
 try:validate(m,wrong)
 except (AssertionError,ValueError):results.append('input '+field+' rejected')
 else:raise AssertionError('accepted '+field)
a=validate(m,r)
with tempfile.TemporaryDirectory(prefix='NONPHYSICAL-slab-verifier-') as tmp:
 folder=Path(tmp);b=a.copy();b.arrays['PW_PBE_forces']=np.zeros((len(a),3));b.info['PW_PBE_energy_eV']=0.
 output=folder/(r['label']+'_PW_PBE.extxyz');write(output,b,format='extxyz')
 summary={'label':r['label'],'source_sha256':r['input_sha256'],'atoms':len(a),'formula':a.get_chemical_formula(),'scf_converged':True,'scf_iterations':1,'mpi_ranks':4,'energy_eV_cell':0.,'free_energy_eV_cell':0.,'method':r['method'],'energy_convention':m['energy_convention'],'dataset_role':r['role'],'fmax_eV_A':0.,'rms_force_eV_A':0.}
 progress={'label':r['label'],'status':'complete','scf_converged':True,'iteration':1,'source_sha256':r['input_sha256'],'mpi_ranks':4}
 (folder/'gpaw.log').write_text('NONPHYSICAL TEST FIXTURE: NOT A CALCULATION\nConverged in 1 steps\n')
 def fixture(s,p,b):
  (folder/'summary.json').write_text(json.dumps(s));(folder/'progress.json').write_text(json.dumps(p));write(output,b,format='extxyz')
  (folder/'SHA256SUMS.txt').write_text(''.join(f'{sha(f)}  {f.name}\n' for f in sorted(folder.iterdir()) if f.name!='SHA256SUMS.txt'))
 def verify():return subprocess.run([sys.executable,str(ROOT/'verify_result.py'),r['label'],str(folder)],capture_output=True,text=True)
 fixture(summary,progress,b);assert verify().returncode==0;results.append('nonphysical valid-shape fixture accepted: not physical evidence')
 for name,mutate in [('method',lambda s,p,b:s['method'].update(smearing_eV=.2)),('Poisson',lambda s,p,b:s['method'].update(poissonsolver={'dipolelayer':'xy'})),('PBC',lambda s,p,b:b.set_pbc([True,True,True])),('source_hash',lambda s,p,b:s.update(source_sha256='0'*64)),('ids',lambda s,p,b:b.arrays['lammps_id'].__setitem__(0,-1)),('unconverged',lambda s,p,b:s.update(scf_converged=False)),('energy_convention',lambda s,p,b:s.update(energy_convention='free only')),('nonfinite',lambda s,p,b:b.arrays['PW_PBE_forces'].__setitem__((0,0),float('nan')))]:
  s=copy.deepcopy(summary);p=copy.deepcopy(progress);bb=b.copy();mutate(s,p,bb);fixture(s,p,bb);assert verify().returncode!=0,name;results.append('result '+name+' rejected despite recomputed fixture inventory')
 fixture(summary,progress,b);(folder/'unexpected.txt').write_text('extra');assert verify().returncode!=0;results.append('extraneous archive file rejected')
( ROOT/'validation.json').write_text(json.dumps({'status':'STATIC_REJECTION_TESTS_PASS','physical_calculations':0,'fixture_warning':'All synthetic energies/forces/logs are nonphysical, temporary and not archived as DFT.','checks':results},indent=2)+'\n')
print('PASS',len(results),'static/fixture cases; no DFT')
