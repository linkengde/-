#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, json, os, subprocess, shutil
import numpy as np
from ase.io import read

PKG=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
POT=Path('/workspace/-/research/mace-v12-transfer/coordination/reports/window-e/ti3sic2_backbone_e5_force_jacobian_20261014/potentials/TSCBOP.source.tersoff')
POT_SHA='729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7'
LMP=Path('/workspace/.venvs/lammps-run0/bin/lmp')
TASKS=[
 ('TSC96_DEV_PRISM_10m10_COD9009647_8L_REF','REF'),
 ('TSC96_DEV_PRISM_10m10_COD9009647_8L_RESP_TIZ020','RESP_TIZ020'),
]
TYPE={'Ti':1,'Si':2,'C':3}; MASS={'Ti':47.867,'Si':28.085,'C':12.011}

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def geom(a):
 h=hashlib.sha256(); h.update(' '.join(a.get_chemical_symbols()).encode())
 h.update(np.asarray(a.cell,dtype='<f8').tobytes()); h.update(np.asarray(np.round(a.positions,8),dtype='<f8').tobytes())
 h.update(np.asarray(a.pbc,dtype=np.uint8).tobytes()); h.update(np.asarray(a.arrays['source_atom_id'],dtype='<i8').tobytes())
 return h.hexdigest()
def make_run(atoms, run_dir):
 run_dir.mkdir(parents=True,exist_ok=True)
 cell=np.asarray(atoms.cell,dtype=float)
 if atoms.pbc.tolist()!=[True,True,True]: raise ValueError('expected fully periodic slab box')
 if np.max(np.abs([cell[0,1],cell[0,2],cell[1,2]]))>1e-12: raise ValueError('unexpected nonrestricted triclinic cell')
 lx,ly,lz=cell[0,0],cell[1,1],cell[2,2]; xy,xz,yz=cell[1,0],cell[2,0],cell[2,1]
 frac=np.asarray(atoms.positions)@np.linalg.inv(cell); img=np.floor(frac+1e-12).astype(int); xyz=(frac-img)@cell
 ids=np.asarray(atoms.arrays['source_atom_id'],int); syms=atoms.get_chemical_symbols()
 if len(set(ids.tolist()))!=len(atoms) or ids.min()<1: raise ValueError('source IDs must be unique positive integers')
 lines=['Static TSCBOP diagnostic; source IDs retained','','%d atoms'%len(atoms),'3 atom types','',f'0 {lx:.16g} xlo xhi',f'0 {ly:.16g} ylo yhi',f'0 {lz:.16g} zlo zhi',f'{xy:.16g} {xz:.16g} {yz:.16g} xy xz yz','','Masses','']
 for el,t in TYPE.items(): lines.append(f'{t} {MASS[el]}')
 lines+=['','Atoms # atomic','']
 with open(run_dir/'atoms.data','w') as f:
  f.write('\n'.join(lines)+'\n')
  for i,e,r in zip(ids,syms,xyz): f.write(f'{i} {TYPE[e]} {r[0]:.16g} {r[1]:.16g} {r[2]:.16g}\n')
 with open(run_dir/'coordinate_map.csv','w',newline='') as f:
  w=csv.writer(f,lineterminator='\n'); w.writerow(['source_atom_id','element','x_source_A','y_source_A','z_source_A','x_wrapped_A','y_wrapped_A','z_wrapped_A'])
  for i,e,p,q in zip(ids,syms,atoms.positions,xyz): w.writerow([int(i),e,*map(float,p),*map(float,q)])
 shutil.copy2(POT,run_dir/'TSCBOP.tersoff')
 (run_dir/'run.in').write_text(f'''clear\nunits metal\natom_style atomic\nboundary p p p\nread_data {(run_dir/'atoms.data').resolve()}\npair_style tersoff\npair_coeff * * {(run_dir/'TSCBOP.tersoff').resolve()} Ti Si C\nneighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes\nthermo 1\nthermo_style custom step atoms pe\ndump forces all custom 1 forces.dump id type x y z fx fy fz\ndump_modify forces sort id format line "%d %d %.16g %.16g %.16g %.16g %.16g %.16g"\nrun 0 post no\nundump forces\n''')
 return ids,syms,xyz

def parse_dump(p):
 lines=p.read_text().splitlines(); sec=next(i for i,x in enumerate(lines) if x.startswith('ITEM: ATOMS'))
 cols=lines[sec].split()[2:]; ix={k:i for i,k in enumerate(cols)}
 out={}
 inv={v:k for k,v in TYPE.items()}
 for line in lines[sec+1:]:
  if not line.strip(): continue
  r=line.split(); i=int(r[ix['id']]); vals=np.array([float(r[ix[k]]) for k in ('x','y','z','fx','fy','fz')])
  if not np.isfinite(vals).all(): raise ValueError(f'nonfinite LAMMPS record: {i}')
  out[i]=(inv[int(r[ix['type']])],vals[:3],vals[3:])
 return out

def metrics(err, ids, mask):
 e=err[mask]; norms=np.linalg.norm(e,axis=1); local=np.flatnonzero(mask)
 return {'atom_count':int(mask.sum()),'component_rmse_eV_A':float(np.sqrt(np.mean(e**2))),
         'vector_rmse_eV_A':float(np.sqrt(np.mean(np.sum(e*e,axis=1)))),
         'p95_vector_error_eV_A':float(np.percentile(norms,95)),
         'max_vector_error_eV_A':float(norms.max()),'max_error_source_atom_id':int(ids[local[np.argmax(norms)]])}

def main():
 if sha(POT)!=POT_SHA: raise ValueError('frozen TSCBOP SHA mismatch')
 if not LMP.is_file(): raise FileNotFoundError('approved LAMMPS run0 executable is missing')
 rows={}
 for sid,suffix in TASKS:
  inp=PKG/'inputs'/f'{sid}.extxyz'; result_dir=PKG/'results'/sid/'attempt-0001'
  r=json.loads((result_dir/'RESULT.json').read_text()); status=json.loads((result_dir/'RUN_STATUS.json').read_text())
  if status.get('status')!='SUCCESS' or not status.get('scf_converged'): raise ValueError(f'{sid}: DFT labels are not converged')
  a=read(inp); lab=read(result_dir/'labelled.extxyz')
  if len(a)!=96 or not np.array_equal(a.numbers,lab.numbers) or not np.array_equal(a.positions,lab.positions) or not np.array_equal(a.cell.array,lab.cell.array) or not np.array_equal(a.pbc,lab.pbc): raise ValueError(f'{sid}: DFT input/output geometry mismatch')
  if sha(inp)!=r['input_sha256'] or geom(a)!=r['geometry_sha256']: raise ValueError(f'{sid}: geometry/input SHA mismatch')
  ids,syms,xyz=make_run(a,OUT/'runs'/suffix)
  dft=np.asarray(r['forces_eV_per_angstrom'],float)
  if dft.shape!=(len(a),3) or not np.isfinite(dft).all(): raise ValueError(f'{sid}: invalid DFT forces')
  run=OUT/'runs'/suffix; env=os.environ.copy(); env['OMP_NUM_THREADS']='1'; env['OMP_DYNAMIC']='FALSE'
  if not (run/'forces.dump').is_file():
   p=subprocess.run([str(LMP),'-in','run.in','-log','lammps.log','-screen','screen.log','-nocite'],cwd=run,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
   (run/'launcher.log').write_text(p.stdout)
   if p.returncode!=0 or not (run/'forces.dump').is_file(): raise RuntimeError(f'{sid}: LAMMPS run0 failed ({p.returncode})')
  elif not (run/'launcher.log').is_file() or 'run 0' not in (run/'run.in').read_text():
   raise ValueError(f'{sid}: existing run0 output lacks a matching launcher/input record')
  lmp=parse_dump(run/'forces.dump'); byid={int(i):n for n,i in enumerate(ids)}
  if set(lmp)!=set(ids.tolist()): raise ValueError(f'{sid}: LAMMPS IDs differ')
  ff=np.zeros_like(dft); maxcoord=0.
  for i in ids:
   j=byid[int(i)]; el,coord,force=lmp[int(i)]
   if el!=syms[j]: raise ValueError(f'{sid}: element mismatch at source ID {i}')
   maxcoord=max(maxcoord,float(np.linalg.norm(coord-xyz[j])))
   if np.linalg.norm(coord-xyz[j])>1e-8: raise ValueError(f'{sid}: coordinate mismatch at source ID {i}')
   ff[j]=force
  err=ff-dft; rows[suffix]={'id':ids,'syms':np.array(syms),'dft':dft,'lmp':ff,'err':err,'maxcoord':maxcoord,'sha':sha(inp),'geom':r['geometry_sha256'],'scf':r['scf_steps'],'run_sha':sha(run/'run.in'),'data_sha':sha(run/'atoms.data')}
  with open(OUT/f'{suffix}_FORCE_COMPARISON.csv','w',newline='') as f:
   cols=['source_atom_id','element','x_A','y_A','z_A','dft_fx_eV_A','dft_fy_eV_A','dft_fz_eV_A','tscbop_fx_eV_A','tscbop_fy_eV_A','tscbop_fz_eV_A','error_fx_eV_A','error_fy_eV_A','error_fz_eV_A','force_error_norm_eV_A']
   w=csv.DictWriter(f,fieldnames=cols,lineterminator='\n'); w.writeheader()
   for j,i in enumerate(ids):
    w.writerow({'source_atom_id':int(i),'element':syms[j],'x_A':a.positions[j,0],'y_A':a.positions[j,1],'z_A':a.positions[j,2],**{f'dft_f{ax}_eV_A':dft[j,k] for k,ax in enumerate('xyz')},**{f'tscbop_f{ax}_eV_A':ff[j,k] for k,ax in enumerate('xyz')},**{f'error_f{ax}_eV_A':err[j,k] for k,ax in enumerate('xyz')},'force_error_norm_eV_A':np.linalg.norm(err[j])})
 with open(OUT/'RESPONSE_DIFFERENCE.csv','w',newline='') as f:
  cols=['source_atom_id','element','dft_delta_fx','dft_delta_fy','dft_delta_fz','tscbop_delta_fx','tscbop_delta_fy','tscbop_delta_fz','response_error_norm_eV_A','within_3p3A_displaced']
  w=csv.DictWriter(f,fieldnames=cols,lineterminator='\n');w.writeheader()
  ref,resp=rows['REF'],rows['RESP_TIZ020']; a=read(PKG/'inputs/TSC96_DEV_PRISM_10m10_COD9009647_8L_REF.extxyz'); r=read(PKG/'inputs/TSC96_DEV_PRISM_10m10_COD9009647_8L_RESP_TIZ020.extxyz')
  displaced=np.asarray([3,90],dtype=int)
  pos={int(i):a.positions[j] for j,i in enumerate(ref['id'])}; disp_pos=[r.positions[np.where(r.arrays['source_atom_id']==i)[0][0]] for i in displaced]
  n=0; errs=[]
  for j,i in enumerate(ref['id']):
   dfd=resp['dft'][j]-ref['dft'][j]; dlp=resp['lmp'][j]-ref['lmp'][j]; e=dlp-dfd
   near=any(np.linalg.norm(pos[int(i)]-p)<3.3 for p in disp_pos); n+=int(near); errs.append(e)
   w.writerow({'source_atom_id':int(i),'element':ref['syms'][j],'dft_delta_fx':dfd[0],'dft_delta_fy':dfd[1],'dft_delta_fz':dfd[2],'tscbop_delta_fx':dlp[0],'tscbop_delta_fy':dlp[1],'tscbop_delta_fz':dlp[2],'response_error_norm_eV_A':np.linalg.norm(e),'within_3p3A_displaced':near})
  errs=np.asarray(errs); near_mask=np.array([any(np.linalg.norm(pos[int(i)]-p)<3.3 for p in disp_pos) for i in ref['id']]); allmask=np.ones(len(ref['id']),bool)
 summary={'scope':'Frozen-coordinate LAMMPS run 0 only; no minimization, MD, or parameter changes. Diagnostic candidate-family analysis; roles unchanged.',
  'potential':{'path':str(POT),'sha256':POT_SHA,'pair_style':'tersoff','pair_coeff_mapping':['Ti','Si','C'],'units':'metal; forces eV/Angstrom'},
  'lammps':{'executable':str(LMP),'version':'22 Jul 2025 Update 4','run_steps':0,'OMP_NUM_THREADS':1},
  'structures':{},'response_error':{},'interpretation':'These geometries share the COD 9009647 parent lineage. This evaluates one geometry family, not an independent crystallographic-source blind test.'}
 for suffix,row in rows.items():
  mask=np.ones(len(row['id']),bool)
  summary['structures'][suffix]={'input_sha256':row['sha'],'geometry_sha256':row['geom'],'scf_steps':row['scf'],'atom_count':96,'wrapped_coordinate_max_error_A':row['maxcoord'],'run_input_sha256':row['run_sha'],'lammps_data_sha256':row['data_sha'],'all':metrics(row['err'],row['id'],mask),'by_element':{e:metrics(row['err'],row['id'],row['syms']==e) for e in ('Ti','Si','C')}}
 e=np.asarray(errs)
 summary['response_error']={'definition':'(F_TSCBOP_RESP-F_TSCBOP_REF)-(F_DFT_RESP-F_DFT_REF)','all':metrics(e,ref['id'],allmask),'within_3p3A_displaced':metrics(e,ref['id'],near_mask),'displaced_source_atom_ids':[int(x) for x in displaced],'near_count':int(near_mask.sum())}
 (OUT/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
