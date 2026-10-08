#!/usr/bin/env python3
from __future__ import annotations
import importlib.util, json, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
AUDIT_LIB=ROOT/'audit'
sys.path.insert(0,str(AUDIT_LIB))
import audit_crater as ac
conn_spec=importlib.util.spec_from_file_location('audit_connectivity', AUDIT_LIB/'audit_connectivity.py')
conn=importlib.util.module_from_spec(conn_spec); conn_spec.loader.exec_module(conn)
START=ROOT/'structures/start_after40.data'
FINAL=ROOT/'structures/cooled_after45_0p5ps.data'
DUMP=ROOT/'dumps/cool_after_45_0p5ps.dump'
CENTER=np.array([9.9573898076,24.2578814220])
MASS={1:47.867,2:28.0855,3:12.0107,4:107.87}
M_VV2E=1.0364269e-4
KB=8.617333262e-5

def frames(path):
    lines=path.read_text().splitlines(); i=0
    while i<len(lines):
        if lines[i]!='ITEM: TIMESTEP': i+=1; continue
        step=int(lines[i+1]); n=int(lines[i+3]); j=i+4
        while not lines[j].startswith('ITEM: ATOMS'): j+=1
        cols=lines[j].split()[2:]
        a=np.array([[float(x) for x in row.split()] for row in lines[j+1:j+1+n]])
        yield step,{c:a[:,k] for k,c in enumerate(cols)}
        i=j+1+n

def frame_positions(d):
    return np.column_stack([d[k] for k in ('x','y','z')])

def local_r(atoms, xyz=None):
    if xyz is None: xyz=atoms['xyz']
    # Exact minimum image in the 2-D triclinic plane.
    origin=np.array([atoms['bounds']['x'][0],atoms['bounds']['y'][0]])
    H=atoms['H']
    frac=(np.linalg.inv(H)@(xyz[:,:2]-origin).T).T
    fc=np.linalg.inv(H)@(CENTER-origin)
    df=frac-fc
    cand=np.stack([(df+[i,j])@H.T for i in (-1,0,1) for j in (-1,0,1)])
    return np.sqrt(np.min(np.sum(cand*cand,axis=2),axis=0))

def temperature(types, vel, mask):
    ix=np.flatnonzero(mask); n=len(ix)
    if n<2: return {'N':n,'T_K':None}
    v=vel[ix].copy()
    m=np.array([MASS[int(t)] for t in types[ix]])
    v-=np.average(v,axis=0,weights=m)
    ke=.5*M_VV2E*np.sum(m[:,None]*v*v)
    return {'N':int(n),'T_K':float(2*ke/((3*n-3)*KB))}

def displacements(a,b,sel):
    ia={int(i):j for j,i in enumerate(a['ids'])}; ib={int(i):j for j,i in enumerate(b['ids'])}
    ids=np.array(sorted(set(ia)&set(ib)),dtype=int)
    ai=np.array([ia[int(i)] for i in ids]); bi=np.array([ib[int(i)] for i in ids])
    use=sel(a,ai,ids)
    ai=ai[use];bi=bi[use]
    df=b['frac'][bi]-a['frac'][ai]
    cand=np.stack([(df+[i,j])@a['H'].T for i in (-1,0,1) for j in (-1,0,1)])
    image=np.argmin(np.sum(cand*cand,axis=2),axis=0)
    dxy=cand[image,np.arange(len(ai))]
    dz=b['xyz'][bi,2]-a['xyz'][ai,2]
    ds=np.sqrt(np.sum(dxy*dxy,axis=1)+dz*dz)
    return {'N':int(len(ds)),'rms_A':float(np.sqrt(np.mean(ds**2))),'median_A':float(np.median(ds)),
      'p90_A':float(np.percentile(ds,90)),'fraction_gt_1A':float(np.mean(ds>1)),
      'fraction_gt_2A':float(np.mean(ds>2)),'fraction_gt_5A':float(np.mean(ds>5))}

def main():
    s=ac.read_data(START); f=ac.read_data(FINAL)
    # final phase-specific temperature statistics in impact region and radial shells
    summaries=[]; last=None; peak={}
    for step,d in frames(DUMP):
        xyz=frame_positions(d); typ=d['type'].astype(int); ids=d['id'].astype(int)
        vel=np.column_stack([d[k] for k in ('vx','vy','vz')]); r=local_r(f,xyz)
        allhot=(r<8)&(xyz[:,2]>=12)&(xyz[:,2]<55)
        row={'step':int(step),'time_ps':float(step*0.00005)}
        for rad,label in [(4,'r_lt_4'),(8,'r_lt_8')]:
            region=(r<rad)&(xyz[:,2]>=12)&(xyz[:,2]<55)
            row[label]={
              'TSC':temperature(typ,vel,region&(ids<=4305)&(typ!=4)),
              'Ag_matrix':temperature(typ,vel,region&(ids<=4305)&(typ==4)),
              'Ag_all':temperature(typ,vel,region&(typ==4)),
            }
        summaries.append(row); last=(step,d,xyz,typ,ids,vel,r,allhot)
    assert last is not None
    step,d,xyz,typ,ids,vel,r,hot=last
    # max instantaneous local phase temperature over all recorded samples
    for row in summaries:
        for shell in ('r_lt_4','r_lt_8'):
            for phase in ('TSC','Ag_matrix','Ag_all'):
                value=row[shell][phase]['T_K']
                if value is not None:
                    key=f'{shell}_{phase}'; old=peak.get(key,{'T_K':-1})
                    if value>old['T_K']: peak[key]={'T_K':value,'step':row['step'],'N':row[shell][phase]['N']}
    # Displacements from quenched 30-shot structure to the final relaxed state, restricted to initially local material.
    rs=ac.min_image_xy(s); local0=(rs<8)&(s['xyz'][:,2]>=12)&(s['xyz'][:,2]<55)
    mov={
      'TSC':displacements(s,f,lambda a,ix,ids: (ids<=4305)&local0[ix]&(a['typ'][ix]!=4)),
      'Ag_matrix':displacements(s,f,lambda a,ix,ids: (ids<=4305)&local0[ix]&(a['typ'][ix]==4)),
    }
    ff=np.column_stack([d[k] for k in ('fx','fy','fz')]); norms=np.linalg.norm(ff,axis=1); k=int(np.argmax(norms))
    out={
      'total_impacts':45,'relaxation_after_hit45_ps':0.5,'timestep':int(step),'time_ps':float(step*0.00005),
      'atom_count_start':len(s['ids']),'atom_count_final':len(f['ids']),
      'local_phase_temperature_final':summaries[-1], 'maximum_sampled_local_phase_temperature':peak,
      'displacement_from_start_after40_state':mov,
      'max_force_eV_per_A':float(norms[k]),'max_force_atom_id':int(ids[k]),'max_force_type':int(typ[k]),
      'max_force_xyz_A':xyz[k].tolist(),
      'note':'Local kinetic temperatures and displacement are indicators, not proof of liquid. Ag-Ti3SiC2 cross-potential is exploratory and unvalidated.'}
    (ROOT/'audit/audit_cool_after45.json').write_text(json.dumps({'summary':out,'samples':summaries},indent=2)+'\n')
    print(json.dumps(out,indent=2))
    print('\nSURFACE PROFILE')
    ac.summary(f,'40 impacts + five 200 fs pulses + 0.5 ps heat-sink relaxation')
    print('\nNEAREST PAIRS')
    ac.pair_audit(f,2.5)
    print('\nCONNECTIVITY')
    a=np.array([f['lx'],0,0]); b=np.array([f['tilt'][0],f['ly'],0])
    for name,mask,cut in [('Ti3SiC2',f['typ']!=4,3.1),('Ag',f['typ']==4,3.35)]:
        largest,ncomp,frac=conn.largest_component(f['xyz'][mask],f['ids'][mask],cut,a,b)
        print(f'{name}: cutoff={cut:.2f} A largest={largest}/{mask.sum()} ({frac:.6f}); clusters={ncomp}')
if __name__=='__main__': main()
