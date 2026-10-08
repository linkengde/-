"""Build hash-pinned training36 only; never reuse acquisition parents as validation."""
import argparse,hashlib,importlib.util,json,subprocess,sys
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
REPO=ROOT.parents[3]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');args=ap.parse_args()
old=BASE/'mace_periodic_v18_clean_core/data';om=json.loads((old/'dataset_manifest.json').read_text());assert sha(old/'train.extxyz')==om['train_input_sha256']
acq=BASE/'pbe_interface_v19_residual_acquisition';im=json.loads((acq/'input_manifest.json').read_text());assert len(im['records'])==6
frames=read(old/'train.extxyz',':');assert len(frames)==30
payload=(old/'train.extxyz').read_bytes();assert payload.endswith(b'\n');newframes=[];rows=[];pins={str((old/'train.extxyz').relative_to(REPO)):sha(old/'train.extxyz'),str((old/'dataset_manifest.json').relative_to(REPO)):sha(old/'dataset_manifest.json'),str((acq/'input_manifest.json').relative_to(REPO)):sha(acq/'input_manifest.json')}
for r in im['records']:
 d=acq/'calculations'/r['label'];subprocess.run([sys.executable,str(acq/'verify_result.py'),r['label'],str(d)],check=True,stdout=subprocess.DEVNULL)
 for line in (d/'SHA256SUMS.txt').read_text().splitlines():
  h,n=line.split('  ',1);assert sha(d/n)==h
 s=json.loads((d/'summary.json').read_text());source=d/(r['label']+'_PW_PBE.extxyz');a=read(source);energy=s['energy_eV_cell'];forces=a.arrays['PW_PBE_forces'];assert a.info['PW_PBE_energy_eV']==energy and np.isfinite(forces).all()
 ids=a.arrays['lammps_id'];pair=a.arrays['central_pair'];header='Lattice="'+' '.join(format(v,'.17g') for v in a.cell.array.ravel())+'" Properties=species:S:1:pos:R:3:REF_forces:R:3:lammps_id:I:1:central_pair:I:1 pbc="T T T"'
 header+=' REF_energy='+format(energy,'.17g')+' free_energy_eV='+format(s['free_energy_eV_cell'],'.17g')+' config_type='+r['label']+' dataset_role=training_acquisition_reused_development_family source_method="GPAW 26.7.0 PW-PBE 500 eV Gamma single point" source_input_sha256='+r['input_sha256']
 text=str(len(a))+'\n'+header+'\n'+''.join(sym+' '+' '.join(format(v,'.17g') for v in list(x)+list(f))+f' {int(i)} {int(p)}\n' for sym,x,f,i,p in zip(a.get_chemical_symbols(),a.positions,forces,ids,pair));payload+=text.encode();newframes.append(a)
 for p in [source,d/'summary.json',d/'verification.json',d/'SHA256SUMS.txt']:pins[str(p.relative_to(REPO))]=sha(p)
 rows.append({'frame_index':30+len(rows),'config_type':r['label'],'input_sha256':r['input_sha256'],'label_sha256':sha(source),'lineage_role':r['role'],'candidate_id':r['candidate_id'],'source_method':s['method'],'energy_convention':s['energy_convention']})
elements=['Ag','C','Si','Ti'];matrix=[[a.get_chemical_symbols().count(e) for e in elements] for a in frames+newframes]
spec=importlib.util.spec_from_file_location('rank',BASE/'mace_periodic_v17_interface_energy/preflight_v17.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);assert mod.exact_rank(np.asarray(matrix,dtype=int))==4
manifest={'status':'TRAINING_CANDIDATE_ONLY_BLOCKED_PENDING_INDEPENDENT_DEVELOPMENT','roles':{'train_candidate':36,'development_validation':0,'test':0},'train_sha256':hashlib.sha256(payload).hexdigest(),'source_files_sha256':pins,'base30_membership_and_bytes_preserved':True,'base_rows':om['rows'],'new_rows':rows,'composition_matrix':matrix,'composition_elements':elements,'composition_rank_exact':4,'withheld_labels_opened':False,'training_authorized':False,'energy_convention':om['energy_convention'],'limitations':['Reused V17 development acquisition families are not fresh validation.','Inherited frame0 Gamma convergence and method-uniformity limitations remain.','Training recipe will remain fixed versus V18; no accuracy improvement asserted.']}
data=ROOT/'training_candidate'
if args.check:assert (data/'train.extxyz').read_bytes()==payload and json.loads((data/'manifest.json').read_text())==manifest
else:
 assert not data.exists(),'Existing candidate refused'
 data.mkdir();(data/'train.extxyz').write_bytes(payload);(data/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(data/'SHA256SUMS.txt').write_text(''.join(f'{sha(data/n)}  {n}\n' for n in ['train.extxyz','manifest.json']))
result=read(data/'train.extxyz',':');assert len(result)==36
for b,a in zip(result[30:],newframes):assert np.array_equal(b.positions,a.positions) and np.array_equal(b.arrays['REF_forces'],a.arrays['PW_PBE_forces']) and b.info['REF_energy']==a.info['PW_PBE_energy_eV']
print(manifest['status'],'36 frames; exact rank4; six labels reproduced exactly')
