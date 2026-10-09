import csv,json,hashlib
from pathlib import Path
import numpy as np
from ase.io import read
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[5]
BUNDLE=REPO/'research/mace-v12-transfer/exports/cross_potential_dft_reference_bundle_20261009'
m=json.loads((BUNDLE/'bundle_manifest.json').read_text())
rows=[]
for e in m['entries']:
 a=read(BUNDLE/e['input']);b=read(BUNDLE/e['result_extxyz'])
 assert a.get_chemical_symbols()==b.get_chemical_symbols() and np.allclose(a.positions,b.positions,rtol=0,atol=1e-12)
 assert np.allclose(a.cell,b.cell,rtol=0,atol=1e-12) and np.array_equal(a.pbc,b.pbc)
 assert np.array_equal(a.arrays['lammps_id'],b.arrays['lammps_id'])
 pair=np.flatnonzero(a.arrays['central_pair']);assert len(pair)==2
 i,j=pair
 if a[i].symbol!='Ag':i,j=j,i
 assert a[i].symbol=='Ag' and a[j].symbol in ['Ti','Si','C']
 v=a.get_distance(int(i),int(j),mic=True,vector=True);d=np.linalg.norm(v);f=b.arrays['PW_PBE_forces'];assert np.isfinite(f).all()
 rows.append(dict(label=e['label'],role=e['role'],contact='Ag-'+a[j].symbol,distance_A=float(d),pair_projection_eV_A=float(np.dot(f[j]-f[i],v/d)),vector_rms_eV_A=float(np.sqrt(np.mean(np.sum(f*f,axis=1)))),native_eV=e['native_energy_eV'],free_eV=e['free_energy_eV'],sigma_eV=e['method']['smearing_eV'],kmesh='x'.join(map(str,e['method']['kpts'])),input_sha256=e['input_sha256']))
with (ROOT/'reference_ledger.csv').open('w') as out:
 w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
report={'status':'REFERENCE_INVENTORY_ONLY_NOT_POTENTIAL_VALIDATION','cases':len(rows),'mixed_potential_predictions':0,'mixed_potential_parameters_available':False,'units':m['units'],'local_diagnostic_rows':[r for r in rows if r['role']=='reused_training_family_local_reference'],'blockers':['Actual Stage69 pair_style/pair_coeff, EAM/Tersoff files, Morse parameters/cutoff/shift and element-type mapping absent in current checkout','No matched separated-phase/slab reference establishing absolute interface energy','Six signed residual perturbations are not systematic normal-separation scans or independent validation','Finite-width native energy is not exactly force-consistent with archived forces; preserve free energy separately'],'next':'Obtain exact unchanged mixed-potential files and calculate static predictions on identical inputs; fit energy differences and forces jointly only within matched method/reference groups. Supplement distance coverage after numerical settings review.'}
(ROOT/'report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'cases':len(rows),'diagnostics':report['local_diagnostic_rows']},ensure_ascii=False,indent=2))
