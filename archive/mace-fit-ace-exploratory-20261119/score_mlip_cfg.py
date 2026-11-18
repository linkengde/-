#!/usr/bin/env python3
import json,re,sys
from pathlib import Path
import numpy as np

def parse(path):
 text=Path(path).read_text()
 out=[]
 for b in re.split(r'(?=^BEGIN_CFG\s*$)',text,flags=re.M):
  if not b.strip(): continue
  m=re.search(r'^\s*Size\s*\n\s*(\d+)',b,re.M); assert m
  n=int(m.group(1)); a=re.search(r'^\s*AtomData:.*\n((?:.*\n){%d})'%n,b,re.M); assert a
  rows=[]
  for line in a.group(1).splitlines():
   x=line.split(); assert len(x)>=8,(path,line)
   rows.append([int(x[0]),int(x[1]),*map(float,x[2:8])])
  e=re.search(r'^\s*Energy\s*\n\s*([-+0-9.eE]+)',b,re.M); assert e
  out.append({'n':n,'types':np.array([x[1] for x in rows]),'forces':np.array([x[5:8] for x in rows]),'energy':float(e.group(1))})
 return out

def main(refp,predp,outp):
 ref,pred=parse(refp),parse(predp); assert len(ref)==len(pred)
 de=[];df=[];si_df=[];si_max_site=[];maxforce=[]
 per=[]
 for r,p in zip(ref,pred):
  assert r['n']==p['n'] and np.array_equal(r['types'],p['types'])
  fdelta=p['forces']-r['forces'];de.append((p['energy']-r['energy'])/r['n']);df.extend(fdelta.reshape(-1))
  si=(r['types']==2)
  if si.any():
   se=np.linalg.norm(fdelta[si],axis=1);si_df.extend(se)
   si_max_site.append(float(se.max()))
   id_hot=np.argmax(np.linalg.norm(r['forces'][si],axis=1))
   maxforce.append(float(se[id_hot]))
  per.append({'natoms':r['n'],'energy_error_eV_atom':(p['energy']-r['energy'])/r['n'],'force_component_rmse_eV_A':float(np.sqrt(np.mean(fdelta**2))),'force_vector_rmse_eV_A':float(np.sqrt(np.mean(np.sum(fdelta**2,axis=1)))),'max_force_vector_error_eV_A':float(np.linalg.norm(fdelta,axis=1).max()),'max_reference_force_Si_site_error_eV_A':maxforce[-1] if si.any() else None})
 out={'energy_rmse_meV_atom':float(1000*np.sqrt(np.mean(np.square(de)))),'energy_mae_meV_atom':float(1000*np.mean(np.abs(de))),'force_component_rmse_eV_A':float(np.sqrt(np.mean(np.square(df)))),'force_vector_rmse_eV_A':float(np.sqrt(np.mean([v['force_vector_rmse_eV_A']**2 for v in per]))),'max_force_vector_error_eV_A':float(max(v['max_force_vector_error_eV_A'] for v in per)),'si_site_vector_error_rmse_eV_A':float(np.sqrt(np.mean(np.square(si_df)))) if si_df else None,'si_site_vector_error_max_eV_A':float(max(si_max_site)) if si_max_site else None,'max_reference_force_Si_site_error_eV_A':float(max(maxforce)) if maxforce else None,'si_center_metric_status':'proxy only: maximum-force Si site in each development frame; no registered hotspot center IDs were found','per_configuration':per}
 Path(outp).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__': main(*sys.argv[1:])
