#!/usr/bin/env python3
"""Fit a small, explicitly idealized Ag/MAX cross Morse model to DFT interface forces.

The five V13 training interfaces calibrate the three Ag-C/Si/Ti pair terms;
V13_VAL_INTERFACE_01 is held out.  This is an exploratory calibration, not a
transferable interface potential or an MD-use certification.
"""
from __future__ import annotations
import json, os, subprocess
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from ase.io import read, write
from ase.neighborlist import neighbor_list

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'v14_stage2_p0_ag_dft_diagnostics' / 'ideal_cross_morse_fit'
OUT.mkdir(parents=True, exist_ok=True)
POT = Path(os.environ.get('V14_POT_DIR', '/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/potentials'))
LMP = os.environ.get('LAMMPS_BIN', '/workspace/venvs/podfit/bin/lmp')
SOURCES = [ROOT/'fit_v13'/'training'/f'V13_IFACE_{i:02d}.xyz' for i in range(1,6)]
SOURCES += [ROOT/'fit_v13'/'validation'/'V13_VAL_INTERFACE_01.xyz']
OLD = np.array([1.35077230, 2.0, 2.15255924,
                0.93944948, 2.0, 2.40407837,
                0.94814194, 2.0, 2.33708864])
PAIR_ORDER = [('Ag','C'), ('Ag','Si'), ('Ag','Ti')]
SPECORDER = ['Ag','C','Si','Ti']

def lmp_base_forces(atoms, tag):
    data, inp, log, dump = [OUT / f'{tag}.{x}' for x in ('data','in','log','dump')]
    write(data, atoms, format='lammps-data', atom_style='atomic', units='metal', specorder=SPECORDER)
    inp.write_text(f'''clear
units metal
atom_style atomic
boundary p p p
atom_modify sort 0 0
read_data {data}
mass 1 107.8682
mass 2 12.011
mass 3 28.0855
mass 4 47.867
pair_style hybrid eam tersoff morse 8.0
pair_coeff 1 1 eam {POT/'Ag_u3.eam'}
pair_coeff * * tersoff {POT/'TSCBOP.tersoff'} NULL C Si Ti
pair_coeff 1 2 morse 0.0 2.0 2.15 8.0
pair_coeff 1 3 morse 0.0 2.0 2.40 8.0
pair_coeff 1 4 morse 0.0 2.0 2.34 8.0
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
thermo_style custom step atoms pe press
thermo 1
run 0
write_dump all custom {dump} id type x y z fx fy fz modify sort id
''')
    env=os.environ.copy(); env.update({'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    p=subprocess.run([LMP,'-in',str(inp),'-log',str(log),'-screen','none'],cwd=OUT,env=env,capture_output=True,text=True)
    if p.returncode: raise RuntimeError(f'{tag}: {p.stdout[-1000:]} {p.stderr[-1000:]}')
    lines=dump.read_text().splitlines(); k=lines.index('ITEM: ATOMS id type x y z fx fy fz'); n=int(lines[3])
    a=np.array([[float(x) for x in row.split()] for row in lines[k+1:k+1+n]])
    return a[np.argsort(a[:,0].astype(int))][:,5:8]

def pair_geometry(atoms):
    syms=np.asarray(atoms.get_chemical_symbols())
    ii,jj,ss,rr=neighbor_list('ijSd',atoms,8.0,self_interaction=False)
    geom={}
    for sa,sb in PAIR_ORDER:
        mask=(syms[ii]==sa)&(syms[jj]==sb)
        i,j,s,r=ii[mask],jj[mask],ss[mask],rr[mask]
        dv=atoms.positions[j]+s@np.asarray(atoms.cell)-atoms.positions[i]
        geom[(sa,sb)]=(i,j,dv,r)
    return geom

def cross_force(natoms, geom, params):
    f=np.zeros((natoms,3),dtype=float)
    for pidx,(sa,sb) in enumerate(PAIR_ORDER):
        D,alpha,r0=params[pidx*3:pidx*3+3]
        ii,jj,dv,rr=geom[(sa,sb)]
        if not len(ii): continue
        x=np.exp(np.clip(-alpha*(rr-r0),-40,40))
        dedr=2*D*alpha*x*(1-x)
        ff=dedr[:,None]*dv/rr[:,None]
        np.add.at(f,ii,ff); np.add.at(f,jj,-ff)
    return f

frames=[]
for src in SOURCES:
    atoms=read(src,format='extxyz'); fref=np.asarray(atoms.get_forces())
    tag=src.stem
    fbase=lmp_base_forces(atoms,tag+'_base')
    frames.append({'atoms':atoms,'fref':fref,'fbase':fbase,'geom':pair_geometry(atoms),'tag':tag})

def residual(p, frame_indices):
    out=[]
    for k in frame_indices:
        fr=frames[k]
        out.extend((fr['fbase']+cross_force(len(fr['atoms']),fr['geom'],p)-fr['fref']).ravel())
    # Weak ridge keeps the underdetermined interface fit near the stable
    # exploratory starting values without dominating force matching.
    out.extend(0.015*(p-OLD))
    return np.asarray(out)

train=np.arange(5); test=[5]
lower=np.tile([0.02,0.5,1.4],3); upper=np.tile([5.0,8.0,4.0],3)
fit=least_squares(lambda p:residual(p,train),OLD,bounds=(lower,upper),
                  max_nfev=500,ftol=1e-10,xtol=1e-10,gtol=1e-10,verbose=0)
params=fit.x

def metrics(fr,p):
    d=fr['fbase']+cross_force(len(fr['atoms']),fr['geom'],p)-fr['fref']
    ev=np.linalg.norm(d,axis=1); by={}
    syms=np.asarray(fr['atoms'].get_chemical_symbols())
    for s in ('Ag','C','Si','Ti'):
        ix=syms==s
        by[s]={'n':int(ix.sum()),'force_rmse_vector_eV_A':float(np.sqrt(np.mean(np.sum(d[ix]**2,axis=1)))),
               'error_p95_eV_A':float(np.percentile(ev[ix],95)),'error_max_eV_A':float(ev[ix].max())}
    return {'name':fr['tag'],'force_rmse_component_eV_A':float(np.sqrt(np.mean(d**2))),
            'force_rmse_vector_eV_A':float(np.sqrt(np.mean(np.sum(d**2,axis=1)))),
            'force_error_p95_eV_A':float(np.percentile(ev,95)),
            'force_error_max_eV_A':float(ev.max()),'by_element':by}

result={'purpose':'ideal-contact exploratory fit; five interface training geometries, one parent-disjoint held-out validation geometry',
        'model':'Foiles Ag_u3 EAM + Ti3SiC2 TSCBOP + pairwise Ag-C/Ag-Si/Ag-Ti Morse cross terms',
        'independent_transferability_established':False,
        'fit_success':bool(fit.success),'fit_message':fit.message,'nfev':int(fit.nfev),
        'parameters':{f'{a}-{b}':{'D_eV':float(params[i*3]),'alpha_Ainv':float(params[i*3+1]),'r0_A':float(params[i*3+2])} for i,(a,b) in enumerate(PAIR_ORDER)},
        'initial_parameters':{f'{a}-{b}':{'D_eV':float(OLD[i*3]),'alpha_Ainv':float(OLD[i*3+1]),'r0_A':float(OLD[i*3+2])} for i,(a,b) in enumerate(PAIR_ORDER)},
        'train': [metrics(frames[i],params) for i in train],
        'heldout': [metrics(frames[i],params) for i in test],
        'limitations':['force-only calibration; does not fit cohesive/interface energies',
                       'three pairwise Morse terms cannot represent angular or charge-transfer bonding',
                       'ideal interface labels are small and sparse',
                       'held-out result is a small-set check, not broad validation']}
(OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
(OUT/'cross_morse_parameters.json').write_text(json.dumps(result['parameters'],indent=2)+'\n')
print(json.dumps(result,indent=2))
