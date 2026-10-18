from pathlib import Path
from ase.io import read,write
import csv,json,sys,hashlib
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parent
TYPE_Z={'Ti':22,'Si':14,'C':6,'Ag':47}
PAIR={'Ti_C':(22,6),'Ag_C':(47,6),'Ag_Ti':(47,22)}

def finalize(model,path,model_path,complete,stop_reason):
    frames=read(path,index=':',format='extxyz')
    rows=[]; highs=[]; ids=np.asarray(frames[0].arrays['source_id'],dtype=int); mobile=np.asarray(frames[0].arrays['mobile'],dtype=bool)
    for atoms in frames:
        pos=atoms.positions.copy(); forces=np.asarray(atoms.arrays['mace_force_eV_A']); fm=np.linalg.norm(forces,axis=1)
        active_f=fm[mobile]; L=atoms.cell.lengths(); p=pos.copy();p[:,1]%=L[1];p[:,2]%=L[2]
        tree=cKDTree(p,boxsize=[1e6,L[1],L[2]]); pp=tree.query_pairs(3.4,output_type='ndarray')
        if len(pp):
            dv=p[pp[:,0]]-p[pp[:,1]]
            dv[:,1]-=np.rint(dv[:,1]/L[1])*L[1];dv[:,2]-=np.rint(dv[:,2]/L[2])*L[2]
            dd=np.linalg.norm(dv,axis=1)
        else: dd=np.array([])
        row={'model':model,'step':int(atoms.info['step']),'time_ps':float(atoms.info['time_ps']),'T_mobile_K':float(atoms.info['T_mobile_K']),'PE_eV':float(atoms.info['PE_eV']),
             'N_mobile':int(mobile.sum()),'Fmax_mobile_eV_A':float(active_f.max()),'P95_mobile_eV_A':float(np.percentile(active_f,95)),'P99_mobile_eV_A':float(np.percentile(active_f,99)),
             'Fmax_mobile_source_id':int(ids[np.flatnonzero(mobile)[np.argmax(active_f)]])}
        z=np.array([TYPE_Z[s] for s in atoms.get_chemical_symbols()])
        for label,(z1,z2) in PAIR.items():
            m=((z[pp[:,0]]==z1)&(z[pp[:,1]]==z2))|((z[pp[:,0]]==z2)&(z[pp[:,1]]==z1)) if len(pp) else np.zeros(0,bool)
            m &= (mobile[pp[:,0]]|mobile[pp[:,1]]) if len(pp) else m
            row['min_'+label+'_A']=float(dd[m].min()) if m.any() else float('nan')
        imap={int(k):i for i,k in enumerate(ids)}
        for label,a,b in [('C19219_Ti25627',19219,25627),('C19219_Ti25628',19219,25628),('C19219_Ag38116',19219,38116)]:
            i,j=imap[a],imap[b]; d=pos[i]-pos[j]
            d[1]-=np.rint(d[1]/L[1])*L[1];d[2]-=np.rint(d[2]/L[2])*L[2]
            row[label+'_A']=float(np.linalg.norm(d))
        rows.append(row)
        for i in np.flatnonzero(mobile & (fm>2.0)):
            highs.append({'model':model,'step':row['step'],'time_ps':row['time_ps'],'source_id':int(ids[i]),'element':atoms[i].symbol,'fmag_eV_A':float(fm[i]),'x_A':pos[i,0],'y_A':pos[i,1],'z_A':pos[i,2]})
    out=path.parent
    with (out/'METRICS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with (out/'HIGH_FORCE_ACTIVE_gt2.csv').open('w',newline='') as f:
        cols=list(highs[0]) if highs else ['model','step','time_ps','source_id','element','fmag_eV_A','x_A','y_A','z_A']
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(highs)
    write(out/'FINAL_LOCAL_PATCH.extxyz',frames[-1],format='extxyz')
    if model_path.exists(): modelsha=hashlib.sha256(model_path.read_bytes()).hexdigest()
    else: modelsha=None
    report={'model':model,'model_sha256':modelsha,'device':'CPU','dtype':'float64','MACE_version':'0.3.16','steps_completed':int(rows[-1]['step']),
            'simulated_time_ps':rows[-1]['time_ps'],'requested_steps':int(sys.argv[6]) if len(sys.argv)>6 else (500 if model=='MACE_MP' else 50),'dt_fs':0.2,'thermostat':'ASE NoseHooverChainNVT, 300 K, tdamp 100 fs',
            'local_patch_atoms':len(frames[-1]),'mobile_atoms':int(mobile.sum()),'fixed_buffer_atoms':int((~mobile).sum()),
            'Fmax_mobile_eV_A':max(r['Fmax_mobile_eV_A'] for r in rows),'P99_mobile_eV_A_max':max(r['P99_mobile_eV_A'] for r in rows),
            'minimum_Ti_C_A':min(r['min_Ti_C_A'] for r in rows),'T_mobile_range_K':[min(r['T_mobile_K'] for r in rows),max(r['T_mobile_K'] for r in rows)],
            'status':'completed' if complete else 'partial_smoke_test','stop_reason':stop_reason,
            'interpretation':'Local fixed-buffer diagnostic only; not a full-model MD or whole-backbone stability validation.'}
    (out/'RUN_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report)

if __name__=='__main__':
    # model key, extxyz path, weight path, complete flag, stop reason
    spec=sys.argv[1:]
    model,rel,weights,complete,reason=spec
    finalize(model,ROOT/rel,ROOT.parent.parent.parent.parent.parent.parent/weights if False else Path(weights),complete=='yes',reason)
