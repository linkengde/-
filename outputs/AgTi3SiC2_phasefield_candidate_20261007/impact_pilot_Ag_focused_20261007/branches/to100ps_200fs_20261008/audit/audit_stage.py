#!/usr/bin/env python3
from __future__ import annotations
import json, re, sys
from pathlib import Path
import importlib.util
import numpy as np
import audit_crater as ac
ROOT=Path(__file__).resolve().parents[1]
conn_spec=importlib.util.spec_from_file_location('audit_connectivity',ROOT/'audit/audit_connectivity.py')
conn=importlib.util.module_from_spec(conn_spec); conn_spec.loader.exec_module(conn)
MASS={1:47.867,2:28.0855,3:12.0107,4:107.87}; KB=8.617333262e-5; M_VV2E=1.0364269e-4

def read_last_dump(path):
    lines=path.read_text(errors='replace').splitlines(); i=0; result=None
    while i<len(lines):
        if lines[i]!='ITEM: TIMESTEP': i+=1; continue
        step=int(lines[i+1]); n=int(lines[i+3]); j=i+4
        while not lines[j].startswith('ITEM: ATOMS'): j+=1
        cols=lines[j].split()[2:]
        a=np.array([[float(v) for v in line.split()] for line in lines[j+1:j+1+n]])
        result=(step,cols,a); i=j+1+n
    if result is None: raise RuntimeError(f'no dump frame: {path}')
    return result

def local_temp(types,vel,mask):
    ix=np.flatnonzero(mask)
    if len(ix)<2:return {'N':len(ix),'T_K':None}
    m=np.array([MASS[int(t)] for t in types[ix]])
    v=vel[ix]-np.average(vel[ix],axis=0,weights=m)
    ke=.5*M_VV2E*np.sum(m[:,None]*v*v)
    return {'N':int(len(ix)),'T_K':float(2*ke/((3*len(ix)-3)*KB))}

def main(hit):
    data_path=ROOT/'structures'/f'impact_{hit:02d}_62eV.data'
    dump_path=ROOT/'dumps'/f'impact_{hit:02d}_62eV.dump'
    log_path=ROOT/'logs'/f'impact_{hit:02d}_62eV.lammps'
    a=ac.read_data(data_path)
    expected=4341+(hit-36)
    if len(a['ids'])!=expected: raise RuntimeError(f'atom count {len(a["ids"])} != expected {expected}')
    log=log_path.read_text(errors='replace')
    if 'ERROR:' in log or 'Lost atoms' in log or re.search(r'\b(?:nan|inf)\b',log): raise RuntimeError('LAMMPS error, lost atoms, or lowercase nan/inf in thermo/log')
    if 'Total wall time:' not in log: raise RuntimeError('LAMMPS log did not reach normal completion')
    step,cols,rows=read_last_dump(dump_path)
    typ=rows[:,cols.index('type')].astype(int); ids=rows[:,cols.index('id')].astype(int)
    xyz=rows[:,[cols.index('x'),cols.index('y'),cols.index('z')]]
    vel=rows[:,[cols.index('vx'),cols.index('vy'),cols.index('vz')]]
    frc=rows[:,[cols.index('fx'),cols.index('fy'),cols.index('fz')]]
    fn=np.linalg.norm(frc,axis=1); k=int(np.argmax(fn))
    if not np.isfinite(rows).all(): raise RuntimeError('non-finite dump values')
    dx=xyz[:,0]-9.9573898076; dy=xyz[:,1]-24.2578814220; rr=np.hypot(dx,dy)
    w=(xyz[:,2]>=12)&(xyz[:,2]<55)
    temp={}
    for rad,key in ((4,'r_lt_4'),(8,'r_lt_8')):
        reg=w&(rr<rad)
        temp[key]={'TSC':local_temp(typ,vel,reg&(typ!=4)&(ids<=4305)),
                   'Ag_matrix':local_temp(typ,vel,reg&(typ==4)&(ids<=4305))}
    H1=np.array([a['lx'],0.,0.]); H2=np.array([a['tilt'][0],a['ly'],0.])
    con={}
    for name,mask,cut in (('TSC',a['typ']!=4,3.1),('Ag',a['typ']==4,3.35)):
        nlargest,ncomp,frac=conn.largest_component(a['xyz'][mask],a['ids'][mask],cut,H1,H2)
        con[name]={'largest':int(nlargest),'total':int(mask.sum()),'clusters':int(ncomp),'fraction':float(frac),'cutoff_A':cut}
    out={'hit':hit,'step':step,'time_ps':step*0.00005,'atoms':len(a['ids']),
         'type_counts':{str(int(t)):int((a['typ']==t).sum()) for t in np.unique(a['typ'])},
         'local_kinetic_temperature_K':temp,
         'max_force_eV_per_A':float(fn[k]),'max_force_atom_id':int(ids[k]),'max_force_type':int(typ[k]),
         'atoms_above_z55':int(np.sum(xyz[:,2]>55)), 'connectivity':con}
    (ROOT/'audit'/f'audit_hit_{hit:02d}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    ac.summary(a,f'hit {hit}')
    ac.pair_audit(a,2.5)
if __name__=='__main__':
    if len(sys.argv)!=2: raise SystemExit('usage: audit_stage.py HIT')
    main(int(sys.argv[1]))
