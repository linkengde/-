#!/usr/bin/env python3
from pathlib import Path
import json,time
import numpy as np,pandas as pd
from ase.io import read
from pyace import create_multispecies_basis_config,PyACECalculator
from pyace.linearacefit import LinearACEFit,LinearACEDataset
from sklearn.linear_model import Ridge
ROOT=Path('/workspace/potential-comparison/mace_fit_data_20261118');OUT=ROOT/'ace_mace_cutoff6';OUT.mkdir(exist_ok=True)
RCUT=6.0

def ds(path):
 rows=[]
 for a in read(path,':'):
  e=float(a.info['REF_energy']);f=np.asarray(a.arrays['REF_forces'],float).copy();a=a.copy()
  if a.pbc.any() and not a.pbc.all():
   old=a.get_all_distances(mic=True);cell=np.asarray(a.cell).copy()
   for axis in range(3):
    if not a.pbc[axis]:
     h=float(np.linalg.norm(cell[axis]));span=float(np.ptp(a.get_scaled_positions(wrap=False)[:,axis])*h)
     if h-span<RCUT+2:cell[axis]*=(span+RCUT+2)/h
   a.set_cell(cell,scale_atoms=False);a.pbc=True;new=a.get_all_distances(mic=True)
   mask=(old<RCUT+1e-8)|(new<RCUT+1e-8);assert np.allclose(old[mask],new[mask],atol=1e-8)
  rows.append({'ase_atoms':a,'energy':e,'energy_corrected':e,'forces':f,'NUMBER_OF_ATOMS':len(a)})
 return pd.DataFrame(rows)
config={'deltaSplineBins':0.001,'elements':['Ag','C','Si','Ti'],'embeddings':{'ALL':{'npot':'FinnisSinclairShiftedScaled','fs_parameters':[1,1],'ndensity':1}},'bonds':{'ALL':{'radbase':'SBessel','radparameters':[6.25],'rcut':RCUT,'dcut':0.01}},'functions':{'ALL':{'nradmax_by_orders':[8,3,1],'lmax_by_orders':[0,3,1]}}}
(OUT/'basis_config.json').write_text(json.dumps(config,indent=2)+'\n')
basis=create_multispecies_basis_config(potential_config=config)
train=LinearACEDataset(basis,ds(ROOT/'train.extxyz'));valid=LinearACEDataset(basis,ds(ROOT/'dev.extxyz'))
start=time.monotonic();train.get_design_matrix(max_workers=2,verbose=True);valid.get_design_matrix(max_workers=2,verbose=True)
y=train.get_target_vector();nE=len(train.df)
def evalfit(fit,df):
 calc=PyACECalculator(fit.get_bbasis());de=[];dc=[];si=[];per=[]
 for _,row in df.iterrows():
  a=row['ase_atoms'].copy();a.calc=calc;ep=float(a.get_potential_energy());fp=a.get_forces();f=np.asarray(row['forces']);d=fp-f;v=np.linalg.norm(d,axis=1);de.append((ep-float(row['energy']))/len(a));dc.extend(d.reshape(-1));mask=np.array(a.get_chemical_symbols())=='Si';si.extend(v[mask]);per.append({'config_type':a.info.get('config_type'),'energy_error_eV_atom':de[-1],'force_component_rmse_eV_A':float(np.sqrt(np.mean(d**2)))})
 return {'energy_rmse_meV_atom':float(1000*np.sqrt(np.mean(np.square(de)))),'force_component_rmse_eV_A':float(np.sqrt(np.mean(np.square(dc)))),'si_site_vector_error_max_eV_A':float(max(si)),'per_configuration':per}
results=[]
for alpha in [1e-4,1e-3,1e-2,1e-1,1.0]:
 w=np.ones(len(y));w[:nE]=10000.;w[nE:]=1.0
 fit=LinearACEFit(train_dataset=train,model=Ridge(alpha=alpha));fit.fit(sample_weight=w)
 name=f'AgCSiTi_B1616_E10000_F1_a{alpha:g}.ace.yaml';fit.get_bbasis().save(str(OUT/name))
 item={'basis_functions':int(basis.total_number_of_functions),'energy_weight':10000,'force_weight':1,'alpha':alpha,'train':evalfit(fit,train.df),'development':evalfit(fit,valid.df),'model':name}
 (OUT/f'{name}.metrics.json').write_text(json.dumps(item,indent=2)+'\n');results.append(item)
 print(json.dumps({'alpha':alpha,'train':{k:v for k,v in item['train'].items() if k!='per_configuration'},'development':{k:v for k,v in item['development'].items() if k!='per_configuration'}},separators=(',',':')),flush=True)
summary={'basis_functions':int(basis.total_number_of_functions),'cutoff_A':RCUT,'wall_s':time.monotonic()-start,'results':results,'weighting':{'energy':10000,'force':1},'alpha_values':[1e-4,1e-3,1e-2,1e-1,1.0]}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
