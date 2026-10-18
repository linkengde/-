from pathlib import Path
import csv, json
import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list
from scipy import sparse
from scipy.sparse.csgraph import connected_components

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'input/NEW_MODEL_ATOMS.data'
TYPE_SYMBOL={1:'Ti',2:'Si',3:'C',4:'Ag'}
SYMBOL_Z={'Ti':22,'Si':14,'C':6,'Ag':47}

def read_data(path):
    lines=path.read_text().splitlines()
    n=None; bounds={}; atom_start=None
    for k,line in enumerate(lines):
        s=line.strip()
        if s.endswith('atoms') and s.split()[0].isdigit(): n=int(s.split()[0])
        t=s.split()
        if len(t)>=4 and t[-2:] == ['xlo','xhi']: bounds['x']=(float(t[0]),float(t[1]))
        if len(t)>=4 and t[-2:] == ['ylo','yhi']: bounds['y']=(float(t[0]),float(t[1]))
        if len(t)>=4 and t[-2:] == ['zlo','zhi']: bounds['z']=(float(t[0]),float(t[1]))
        if s.startswith('Atoms'):
            atom_start=k+1
            break
    while atom_start<len(lines) and not lines[atom_start].strip(): atom_start+=1
    rows=[]
    for line in lines[atom_start:]:
        s=line.strip()
        if not s: break
        t=s.split()
        if len(t)<5: break
        aid,typ=int(t[0]),int(t[1]); rows.append((aid,typ,*map(float,t[2:5])))
        if n and len(rows)==n: break
    a=np.array(rows,dtype=float)
    ids=a[:,0].astype(int); types=a[:,1].astype(int); pos=a[:,2:5]
    cell=np.array([bounds['x'][1]-bounds['x'][0],bounds['y'][1]-bounds['y'][0],bounds['z'][1]-bounds['z'][0]])
    return ids,types,pos,cell

def read_dump(path):
    text=path.read_text().splitlines(); frames=[]; i=0
    while i<len(text):
        if text[i].strip()!='ITEM: TIMESTEP': i+=1; continue
        step=int(text[i+1]); i+=2
        assert text[i].startswith('ITEM: NUMBER OF ATOMS'); n=int(text[i+1]); i+=2
        assert text[i].startswith('ITEM: BOX BOUNDS')
        bounds=[]
        for row in text[i+1:i+4]:
            v=row.split(); bounds.append((float(v[0]),float(v[1])))
        i+=4
        assert text[i].startswith('ITEM: ATOMS')
        cols=text[i].split()[2:]; i+=1
        data=np.loadtxt(text[i:i+n],dtype=float,ndmin=2); i+=n
        ix={c:k for k,c in enumerate(cols)}
        order=np.argsort(data[:,ix['id']])
        data=data[order]
        frames.append((step,data,cols,np.array([b[1]-b[0] for b in bounds])))
    return frames

def make_atoms(types,pos,cell):
    syms=[TYPE_SYMBOL[int(t)] for t in types]
    return Atoms(numbers=[SYMBOL_Z[s] for s in syms],positions=pos,cell=cell,pbc=[False,True,True])

def force_stats(data,cols):
    ix={c:k for k,c in enumerate(cols)}
    f=data[:,[ix['fx'],ix['fy'],ix['fz']]]
    norm=np.linalg.norm(f,axis=1)
    return norm,f,{
      'n':int(len(norm)),'finite':bool(np.isfinite(norm).all()),
      'Fmax_eV_A':float(norm.max()),'P95_eV_A':float(np.percentile(norm,95)),
      'P99_eV_A':float(np.percentile(norm,99)),'mean_eV_A':float(norm.mean()),
      'rms_eV_A':float(np.sqrt(np.mean(norm**2))),
      'n_gt_2':int((norm>2).sum()),'n_gt_5':int((norm>5).sum()),'n_gt_10':int((norm>10).sum())}

def graph_connectivity(atoms,types,phase_types,cutoff):
    mask=np.isin(types,phase_types); glob=np.flatnonzero(mask); loc=np.full(len(types),-1,int); loc[glob]=np.arange(len(glob))
    ii,jj=neighbor_list('ij',atoms,cutoff,self_interaction=False)
    keep=(loc[ii]>=0)&(loc[jj]>=0)&(ii<jj)
    edges=np.vstack([loc[ii[keep]],loc[jj[keep]]])
    mat=sparse.coo_matrix((np.ones(edges.shape[1],dtype=np.int8),(edges[0],edges[1])),shape=(len(glob),len(glob)))
    mat=mat+mat.T
    nc,labels=connected_components(mat,directed=False,return_labels=True)
    counts=np.bincount(labels) if len(labels) else np.array([],int)
    return {'n_atoms':int(len(glob)),'cutoff_A':cutoff,'components':int(nc),'largest_component':int(counts.max()) if len(counts) else 0,'largest_fraction':float(counts.max()/len(glob)) if len(counts) else 0.0}

def distance_audit(atoms,types,limit=5.5):
    ii,jj,dd=neighbor_list('ijd',atoms,limit,self_interaction=False)
    keep=ii<jj; ii,jj,dd=ii[keep],jj[keep],dd[keep]
    pairs={}
    for typ1,typ2,name in [(1,3,'Ti-C'),(1,2,'Ti-Si'),(4,1,'Ag-Ti'),(4,2,'Ag-Si'),(4,3,'Ag-C'),(4,4,'Ag-Ag')]:
        mask=((types[ii]==typ1)&(types[jj]==typ2))|((types[ii]==typ2)&(types[jj]==typ1))
        vals=dd[mask]
        pairs[name]={'count_within_5p5A':int(len(vals)),'minimum_A':float(vals.min()) if len(vals) else None,
                     'n_lt_1p8A':int((vals<1.8).sum()),'n_lt_2p0A':int((vals<2.0).sum()),'n_lt_2p2A':int((vals<2.2).sum())}
    return pairs,ii,jj,dd

orig_ids,orig_types,orig_pos,cell=read_data(DATA)
pre=read_dump(ROOT/'pre_relax.dump')[-1]
post=read_dump(ROOT/'post_relax.dump')[-1]
frames=read_dump(ROOT/'minimization_frames.dump')
assert len(orig_ids)==39486 and len(post[1])==39486
assert np.array_equal(orig_ids,pre[1][:,pre[2].index('id')].astype(int))
assert np.array_equal(orig_ids,post[1][:,post[2].index('id')].astype(int))
assert np.array_equal(orig_types,post[1][:,post[2].index('type')].astype(int))
predata,precols=pre[1],pre[2]; postdata,postcols=post[1],post[2]
pre_pos=predata[:,[precols.index('x'),precols.index('y'),precols.index('z')]]
post_pos=postdata[:,[postcols.index('x'),postcols.index('y'),postcols.index('z')]]
pre_force,pre_vec,pre_stats=force_stats(predata,precols)
post_force,post_vec,post_stats=force_stats(postdata,postcols)
# Minimum-image displacement in periodic Y and Z only.
disp=post_pos-pre_pos
for ax in [1,2]: disp[:,ax]-=np.rint(disp[:,ax]/cell[ax])*cell[ax]
dnorm=np.linalg.norm(disp,axis=1)
pre_atoms=make_atoms(orig_types,pre_pos,cell); post_atoms=make_atoms(orig_types,post_pos,cell)
pre_dist,_,_,_=distance_audit(pre_atoms,orig_types)
post_dist,ii,jj,dd=distance_audit(post_atoms,orig_types)
pre_conn={'Ag':graph_connectivity(pre_atoms,orig_types,[4],3.05),'TSC':graph_connectivity(pre_atoms,orig_types,[1,2,3],3.30)}
post_conn={'Ag':graph_connectivity(post_atoms,orig_types,[4],3.05),'TSC':graph_connectivity(post_atoms,orig_types,[1,2,3],3.30)}
# nearest unlike-phase distances for high-force atoms from within the audited 5.5-A pair set
idxmap={int(a):k for k,a in enumerate(orig_ids)}
high=np.argsort(post_force)[::-1][:50]
rows=[]
for k in high:
    aid=int(orig_ids[k]); typ=int(orig_types[k]); sym=TYPE_SYMBOL[typ]
    other=(orig_types[jj] if typ==4 else orig_types[ii])
    edge_atoms=(ii if typ==4 else jj)
    m=edge_atoms==k
    opposite=np.isin(other[m],[1,2,3]) if typ==4 else (other[m]==4)
    dmin=float(dd[m][opposite].min()) if opposite.any() else ''
    rows.append([aid,sym,*post_pos[k].tolist(),float(post_force[k]),dmin,float(dnorm[k])])
with (ROOT/'HIGH_FORCE_POST_RELAX.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.writer(f,lineterminator='\n'); w.writerow(['atom_id','element','x_A','y_A','z_A','force_eV_A','nearest_opposite_phase_within_5p5A_A','displacement_from_run0_A']); w.writerows(rows)
steps=[]
for step,data,cols,_ in frames:
  _,_,s=force_stats(data,cols)
  idx={c:i for i,c in enumerate(cols)}; norms=np.linalg.norm(data[:,[idx['fx'],idx['fy'],idx['fz']]],axis=1)
  top=int(np.argmax(norms)); steps.append([step,float(norms[top]),int(data[top,idx['id']]),float(np.percentile(norms,95)),float(np.percentile(norms,99)),float(norms.mean()),int((norms>2).sum()),int((norms>5).sum())])
with (ROOT/'MINIMIZATION_STEP_SUMMARY.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.writer(f,lineterminator='\n'); w.writerow(['step','Fmax_eV_A','Fmax_atom_id','P95_eV_A','P99_eV_A','mean_force_eV_A','n_F_gt_2','n_F_gt_5']); w.writerows(steps)
by_element={}
for typ,sym in TYPE_SYMBOL.items():
 m=orig_types==typ; v=post_force[m]
 by_element[sym]={'n':int(m.sum()),'Fmax':float(v.max()),'P95':float(np.percentile(v,95)),'P99':float(np.percentile(v,99)),'n_gt_2':int((v>2).sum()),'n_gt_5':int((v>5).sum())}
result={'atoms':len(orig_ids),'cell_A':cell.tolist(),'pbc':'f p p','ids_types_preserved':True,'finite_forces':pre_stats['finite'] and post_stats['finite'],'run0_force':pre_stats,'post_relax_force':post_stats,'post_relax_force_by_element':by_element,'displacement_A':{'P50':float(np.percentile(dnorm,50)),'P95':float(np.percentile(dnorm,95)),'P99':float(np.percentile(dnorm,99)),'max':float(dnorm.max()),'n_gt_1A':int((dnorm>1).sum()),'n_gt_3A':int((dnorm>3).sum()),'atom_id_max':int(orig_ids[np.argmax(dnorm)])},'pair_distances_pre_relax':pre_dist,'pair_distances_post_relax':post_dist,'connectivity_pre_relax':pre_conn,'connectivity_post_relax':post_conn,'minimization_frame_count':len(frames),'minimization_steps_saved':[int(x[0]) for x in steps]}
(ROOT/'RELAXATION_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
with (ROOT/'RELAXATION_AUDIT.md').open('w') as f:
 f.write('# 39,486-atom exploratory short-relaxation audit\n\n')
 f.write('## Identity and numerical checks\n\n')
 f.write(f"- Atom count: {len(orig_ids):,}; IDs and types preserved in final data: **{result['ids_types_preserved']}**. PBC: `f p p`; cell: {cell[0]:.6f} × {cell[1]:.6f} × {cell[2]:.6f} Å.\n")
 f.write(f"- Finite force vectors: **{result['finite_forces']}**. Initial Fmax/mean: {pre_stats['Fmax_eV_A']:.4f}/{pre_stats['mean_eV_A']:.4f} eV/Å. Final Fmax/P95/P99/mean: {post_stats['Fmax_eV_A']:.4f}/{post_stats['P95_eV_A']:.4f}/{post_stats['P99_eV_A']:.4f}/{post_stats['mean_eV_A']:.4f} eV/Å.\n")
 f.write(f"- Maximum minimum-image displacement: {dnorm.max():.4f} Å (atom {int(orig_ids[np.argmax(dnorm)])}); P95/P99: {np.percentile(dnorm,95):.4f}/{np.percentile(dnorm,99):.4f} Å; atoms displaced >1 Å: {int((dnorm>1).sum())}.\n\n")
 f.write('## Connectivity at archived geometry cutoffs\n\n')
 for phase,cut in [('Ag',3.05),('TSC',3.30)]:
  a=pre_conn[phase]; b=post_conn[phase]
  f.write(f"- {phase}, cutoff {cut:.2f} Å: before {a['components']} component(s), largest {a['largest_fraction']:.6f}; after {b['components']} component(s), largest {b['largest_fraction']:.6f}.\n")
 f.write('\n## Pair-distance minima\n\n| pair | before (Å) | after (Å) | count after <2.0 Å |\n|---|---:|---:|---:|\n')
 for name,v in post_dist.items(): f.write(f"| {name} | {pre_dist[name]['minimum_A']:.4f} | {v['minimum_A']:.4f} | {v['n_lt_2p0A']} |\n")
 f.write('\n## Decision\n\n')
 f.write('The run reached its 500-iteration cap; it did not meet the requested 0.5 eV/Å Fmax tolerance. The minimization trace includes a large transient force spike. Review `MINIMIZATION_STEP_SUMMARY.csv` and `HIGH_FORCE_POST_RELAX.csv`; this result alone does not authorize heating. The next decision is based on whether the final high-force atoms and short contacts are isolated surface/interface rearrangements and whether both networks remain connected.\n')
print(json.dumps(result,indent=2))
print('minimization frames',len(frames),'steps', [s[0] for s in steps])
