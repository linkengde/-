from pathlib import Path
import csv, math, json, hashlib
import numpy as np
from ase.io import read

root=Path('/workspace/-/research/mace-v12-transfer/coordination/reports/window-e/tsc_surface_dev_prismatic_20261020')
out=root/'analysis/ENVIRONMENT_ERROR_ANALYSIS_20261021'
ref_path=root/'inputs/TSC96_DEV_PRISM_10m10_COD9009647_8L_REF.extxyz'
ref_csv=root/'analysis/DEV_TSCBOP_RUN0_PAIR_20261021/REF_FORCE_COMPARISON.csv'
resp_csv=root/'analysis/DEV_TSCBOP_RUN0_PAIR_20261021/RESP_TIZ020_FORCE_COMPARISON.csv'
resp_path=root/'inputs/TSC96_DEV_PRISM_10m10_COD9009647_8L_RESP_TIZ020.extxyz'
delta_csv=root/'analysis/DEV_TSCBOP_RUN0_PAIR_20261021/RESPONSE_DIFFERENCE.csv'
atoms=read(ref_path)
resp_atoms=read(resp_path)
ids=np.asarray(atoms.arrays['source_atom_id'],int)
symbols=np.asarray(atoms.get_chemical_symbols())
pos=np.asarray(atoms.positions,float)
assert len(atoms)==96 and np.array_equal(ids,np.arange(1,97))

# Read and bind force rows to source IDs and coordinates.
def load_force(path, geometry):
    rows={}
    with open(path,newline='') as f:
        for row in csv.DictReader(f):
            aid=int(row['source_atom_id'])
            xyz=np.array([float(row[k]) for k in ('x_A','y_A','z_A')])
            rows[aid]=(row,xyz)
    assert sorted(rows)==ids.tolist()
    for i,aid in enumerate(ids):
        assert rows[int(aid)][0]['element']==geometry.get_chemical_symbols()[i]
        assert np.max(np.abs(rows[int(aid)][1]-geometry.positions[i]))<1e-7
    return rows
ref=load_force(ref_csv,atoms); resp=load_force(resp_csv,resp_atoms)
# Force error vector, surface-normal depth, and first-shell counts.
err={}; dft={}; tsc={}; depth=[]; side=[]; c_ti=[]; si_ti=[]; ti_ti=[]
for i,aid in enumerate(ids):
    row=ref[int(aid)][0]
    err[int(aid)]=np.array([float(row[f'error_f{axis}_eV_A']) for axis in 'xyz'])
    dft[int(aid)]=np.array([float(row[f'dft_f{axis}_eV_A']) for axis in 'xyz'])
    tsc[int(aid)]=np.array([float(row[f'tscbop_f{axis}_eV_A']) for axis in 'xyz'])
z=pos[:,2]; zlo=float(z.min()); zhi=float(z.max())
for i,sym in enumerate(symbols):
    dep=min(float(z[i]-zlo),float(zhi-z[i]))
    depth.append(dep); side.append('bottom' if z[i]-zlo <= zhi-z[i] else 'top')
    cnt_ti_c=cnt_ti_si=cnt_ti_ti=cnt_si_ti=0
    for j,other in enumerate(symbols):
        if i==j: continue
        dist=float(atoms.get_distance(i,j,mic=True))
        pair={sym,other}
        if pair=={'Ti','C'} and dist<2.50:
            cnt_ti_c+=1
        if pair=={'Ti','Si'} and dist<3.10:
            if sym=='Ti': cnt_ti_si+=1
            else: cnt_si_ti+=1
        if sym=='Ti' and other=='Ti' and dist<3.25:
            cnt_ti_ti+=1
    c_ti.append(cnt_ti_c); si_ti.append(cnt_ti_si if sym=='Ti' else cnt_si_ti); ti_ti.append(cnt_ti_ti)
depth=np.array(depth); side=np.array(side)
# Set structural classes from distance to the nearest slab edge; cutoffs descriptive.
def depth_bin(d):
    if d<2.5: return 'outer_<2.5A'
    if d<5.0: return 'subsurface_2.5-5A'
    return 'inner_>=5A'
classes=np.array([depth_bin(d) for d in depth])

# Make response rows indexed by atom ID.
response={}
with open(delta_csv,newline='') as f:
    for row in csv.DictReader(f):
        aid=int(row['source_atom_id'])
        response[aid]=row
assert sorted(response)==ids.tolist()

def stats(vectors, indices):
    v=np.asarray([vectors[int(ids[i])] for i in indices],float)
    norms=np.linalg.norm(v,axis=1)
    return {'n_atoms':int(len(indices)), 'vector_rmse_eV_A':float(np.sqrt(np.mean(np.sum(v*v,axis=1)))) if len(indices) else None,
            'p95_error_norm_eV_A':float(np.quantile(norms,.95)) if len(indices) else None,
            'max_error_norm_eV_A':float(np.max(norms)) if len(indices) else None,
            'mean_error_norm_eV_A':float(np.mean(norms)) if len(indices) else None}

rows=[]
for i,aid in enumerate(ids):
    q=response[int(aid)]
    de=np.array([float(q[f'dft_delta_f{axis}']) for axis in 'xyz'])
    te=np.array([float(q[f'tscbop_delta_f{axis}']) for axis in 'xyz'])
    re=te-de
    rows.append({'source_atom_id':int(aid),'element':str(symbols[i]),'x_A':pos[i,0],'y_A':pos[i,1],'z_A':pos[i,2],
                 'nearest_surface':side[i],'depth_A':depth[i],'depth_bin':classes[i],
                 'n_TiC_partners_first_shell':int(c_ti[i]),'n_TiSi_partners_first_shell':int(si_ti[i]),
                 'ref_error_fx_eV_A':err[int(aid)][0],'ref_error_fy_eV_A':err[int(aid)][1],'ref_error_fz_eV_A':err[int(aid)][2],
                 'ref_error_norm_eV_A':float(np.linalg.norm(err[int(aid)])),
                 'dft_force_norm_eV_A':float(np.linalg.norm(dft[int(aid)])),
                 'tscbop_force_norm_eV_A':float(np.linalg.norm(tsc[int(aid)])),
                 'response_error_fx_eV_A':re[0],'response_error_fy_eV_A':re[1],'response_error_fz_eV_A':re[2],
                 'response_error_norm_eV_A':float(np.linalg.norm(re))})

with open(out/'FORCE_ERROR_BY_ENVIRONMENT.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

summary={'source_ref_input_sha256':hashlib.sha256(ref_path.read_bytes()).hexdigest(),
         'source_ref_force_csv_sha256':hashlib.sha256(ref_csv.read_bytes()).hexdigest(),
         'response_force_csv_sha256':hashlib.sha256(resp_csv.read_bytes()).hexdigest(),
         'response_difference_csv_sha256':hashlib.sha256(delta_csv.read_bytes()).hexdigest(),
         'geometry':{'atom_count':len(atoms),'composition':{s:int(np.sum(symbols==s)) for s in sorted(set(symbols))},'surface_span_A':zhi-zlo,
                     'z_cell_A':float(atoms.cell[2,2]),'vacuum_A':float(atoms.cell[2,2]-(zhi-zlo)),'depth_class_cutoffs_A':[2.5,5.0]},
         'coordination_cutoffs_A':{'Ti-C':2.50,'Ti-Si':3.10,'Ti-C':2.50},
         'reference_by_element':{},'reference_by_depth':{},'reference_by_element_and_depth':{},
         'response_by_element':{},'response_by_depth':{},'response_by_element_and_depth':{},
         'coordination_error_summary':{}}
for sym in sorted(set(symbols)):
    ix=np.flatnonzero(symbols==sym)
    summary['reference_by_element'][sym]=stats(err,ix)
    vec={int(ids[i]):np.array([float(response[int(ids[i])][f'tscbop_delta_f{a}'])-float(response[int(ids[i])][f'dft_delta_f{a}']) for a in 'xyz']) for i in ix}
    summary['response_by_element'][sym]=stats(vec,ix)
for b in sorted(set(classes),key=lambda x:['outer_<2.5A','subsurface_2.5-5A','inner_>=5A'].index(x)):
    ix=np.flatnonzero(classes==b)
    summary['reference_by_depth'][b]=stats(err,ix)
    vec={int(ids[i]):np.array([float(response[int(ids[i])][f'tscbop_delta_f{a}'])-float(response[int(ids[i])][f'dft_delta_f{a}']) for a in 'xyz']) for i in ix}
    summary['response_by_depth'][b]=stats(vec,ix)
    for sym in sorted(set(symbols)):
        sub=np.array([i for i in ix if symbols[i]==sym],int)
        if len(sub):
            summary['reference_by_element_and_depth'][f'{sym}/{b}']=stats(err,sub)
            summary['response_by_element_and_depth'][f'{sym}/{b}']=stats(vec,sub)
# Group exact first-shell counts, showing force residuals rather than asserting ideal coordination.
for sym in sorted(set(symbols)):
    for n in sorted(set(c_ti[i] if sym in ('Ti','C') else si_ti[i] for i in range(len(ids)) if symbols[i]==sym)):
        ix=np.array([i for i in range(len(ids)) if symbols[i]==sym and (c_ti[i] if sym in ('Ti','C') else si_ti[i])==n],int)
        summary['coordination_error_summary'][f'{sym}_n_relevant_neighbors_{n}']=stats(err,ix)
with open(out/'SUMMARY.json','w') as f: json.dump(summary,f,indent=2,sort_keys=True); f.write('\n')
print(json.dumps(summary,indent=2))
