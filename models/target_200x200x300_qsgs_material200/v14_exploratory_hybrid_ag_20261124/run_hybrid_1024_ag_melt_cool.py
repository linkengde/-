#!/usr/bin/env python3
"""Small, imposed-temperature Ag melt/cool screen on a 1024-atom ideal interface.

Ag uses Foiles EAM, Ti3SiC2 uses TSCBOP, cross interactions use a DFT-force-
calibrated Morse candidate plus a cross-only ZBL short-range guard. This is
not a TTM simulation or validation of the full porous target model.
"""
import json, os, re, subprocess
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.neighborlist import neighbor_list
from scipy.special import sph_harm_y

ROOT=Path(__file__).resolve().parent
TARGET=int(os.environ.get('V14_AG_TARGET_K','1800'))
HEAT_PS=float(os.environ.get('V14_AG_HEAT_RAMP_PS','0.5'))
HOLD_PS=float(os.environ.get('V14_AG_HOLD_PS','1.0'))
COOL_PS=float(os.environ.get('V14_AG_COOL_PS','1.5'))
FINAL_PS=float(os.environ.get('V14_AG_FINAL_HOLD_PS','0.5'))
RUN_TAG=os.environ.get('V14_RUN_TAG',f'{TARGET}K')
HEAT_STEPS=round(HEAT_PS*1000); HOLD_STEPS=round(HOLD_PS*1000)
COOL_STEPS=round(COOL_PS*1000); FINAL_STEPS=round(FINAL_PS*1000)
OUT=ROOT/'v14_stage2_p0_ag_dft_diagnostics'/f'hybrid_1024_ag_melt_cool_{RUN_TAG}'
OUT.mkdir(parents=True,exist_ok=True)
SRC=ROOT/'fit_v13'/'validation'/'V13_VAL_INTERFACE_01.xyz'
POT=Path(os.environ.get('V14_POT_DIR', '/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer/potentials'))
LMP=os.environ.get('LAMMPS_BIN', '/workspace/venvs/podfit/bin/lmp')
CELL=read(SRC,format='extxyz').repeat((4,4,1))
DATA=OUT/'ideal_interface_4x4x1_1024.data'
write(DATA,CELL,format='lammps-data',atom_style='atomic',units='metal',specorder=['Ag','C','Si','Ti'])
INIT=OUT/'minimized_initial.dump'; TRAJ=OUT/'cycle.dump'; LOG=OUT/'log.lammps'; IN=OUT/'in.cycle'
IN.write_text(f'''clear
units metal
atom_style atomic
boundary p p p
atom_modify sort 0 0
read_data {DATA}
mass 1 107.8682
mass 2 12.011
mass 3 28.0855
mass 4 47.867
group silver type 1
group skeleton type 2 3 4
pair_style hybrid/overlay eam tersoff morse 8.0 zbl 1.5 2.4
pair_coeff 1 1 eam {POT/'Ag_u3.eam'}
pair_coeff * * tersoff {POT/'TSCBOP.tersoff'} NULL C Si Ti
pair_coeff 1 2 morse 0.020000 7.804451 1.400000 8.0
pair_coeff 1 3 morse 0.020000 3.107766 1.400605 8.0
pair_coeff 1 4 morse 0.678644 1.155035 2.989514 8.0
pair_coeff 1 2 zbl 47 6
pair_coeff 1 3 zbl 47 14
pair_coeff 1 4 zbl 47 22
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
variable fmag atom sqrt(fx*fx+fy*fy+fz*fz)
compute fmax all reduce max v_fmag
thermo 50
thermo_style custom step time temp pe press c_fmax
min_style fire
min_modify dmax 0.05
minimize 1.0e-8 1.0e-3 2000 20000
write_dump all custom {INIT} id type x y z fx fy fz modify sort id
write_data {OUT/'after_fire.data'}
reset_timestep 0
compute TAg silver temp
compute TSkel skeleton temp
compute MSDAg silver msd com yes
thermo_style custom step time temp c_TAg c_TSkel c_MSDAg[4] pe ke etotal press c_fmax
velocity silver create 300.0 24601 mom yes rot yes dist gaussian
velocity skeleton create 300.0 24602 mom yes rot yes dist gaussian
fix sink skeleton nvt temp 300.0 300.0 0.1
fix equil silver nvt temp 300.0 300.0 0.1
timestep 0.001
dump tr all custom 50 {TRAJ} id type x y z vx vy vz fx fy fz
dump_modify tr sort id
run 200
unfix equil
fix heat silver nvt temp 300.0 {TARGET}.0 0.1
run {HEAT_STEPS}
unfix heat
fix hold silver nvt temp {TARGET}.0 {TARGET}.0 0.1
run {HOLD_STEPS}
unfix hold
fix cool silver nvt temp {TARGET}.0 300.0 0.1
run {COOL_STEPS}
unfix cool
fix final silver nvt temp 300.0 300.0 0.1
run {FINAL_STEPS}
write_data {OUT/'after_cycle.data'}
''')
env=os.environ.copy();env.update({'OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
p=subprocess.run([LMP,'-sf','omp','-pk','omp','4','-in',str(IN),'-log',str(LOG),'-screen','none'],cwd=OUT,env=env,capture_output=True,text=True)
if p.returncode: raise RuntimeError((p.stdout+'\n'+p.stderr)[-5000:])

def q6_ag(a,cutoff=3.8):
    sy=np.asarray(a.get_chemical_symbols());ii,jj,S=neighbor_list('ijS',a,cutoff)
    m=(sy[ii]=='Ag')&(sy[jj]=='Ag');ii=ii[m];jj=jj[m];S=S[m]
    v=a.positions[jj]+S@np.asarray(a.cell)-a.positions[ii];r=np.linalg.norm(v,axis=1)
    good=r>1e-8;ii=ii[good];v=v[good];r=r[good]
    theta=np.arccos(np.clip(v[:,2]/r,-1,1));phi=np.arctan2(v[:,1],v[:,0])
    q=np.zeros((len(a),13),complex)
    for k,mv in enumerate(range(-6,7)):np.add.at(q[:,k],ii,sph_harm_y(6,mv,theta,phi))
    cnt=np.bincount(ii,minlength=len(a));nz=cnt>0;q[nz]/=cnt[nz,None]
    qi=np.sqrt(4*np.pi/13*np.sum(np.abs(q)**2,axis=1));return qi[sy=='Ag']

def analyze(frames,offset=0,spacing=50):
    rows=[];prev=None;cum=None
    for k,a in enumerate(frames):
        sy=np.asarray(a.get_chemical_symbols());ai=np.where(sy=='Ag')[0]
        ci=np.where(sy=='C')[0];si=np.where(sy=='Si')[0];ti=np.where(sy=='Ti')[0]
        ii,jj,dd=neighbor_list('ijd',a,5.0)
        def md(x,y,same=False):
            m=(sy[ii]==x)&(sy[jj]==y)
            if same:m &= ii<jj
            val=dd[m]
            return float(val.min()) if len(val) else None
        if prev is None:cum=np.zeros((len(ai),3))
        else:
            dv=a.positions[ai]-prev[ai];frac=dv@np.linalg.inv(np.asarray(a.cell));frac-=np.rint(frac);cum+=frac@np.asarray(a.cell)
        prev=a.positions.copy();q=q6_ag(a)
        step=offset+k*spacing
        rows.append({'step':step,'time_ps':step*0.001,'natoms':len(a),'Q6_Ag_mean':float(q.mean()),
            'Q6_Ag_p50':float(np.median(q)),'Q6_Ag_fraction_gt_0p4':float(np.mean(q>0.4)),
            'MSD_Ag_A2':float(np.mean(np.sum(cum*cum,axis=1))),
            'min_Ag_Ag_A':md('Ag','Ag',True),'min_Ag_C_A':md('Ag','C'),'min_Ag_Si_A':md('Ag','Si'),
            'min_Ag_Ti_A':md('Ag','Ti'),'min_Ti_C_A':md('Ti','C')})
    return rows

text=LOG.read_text();thermo=[];active=False
for line in text.splitlines():
    if all(x in line for x in ['Step','Time','Temp','c_TAg','c_TSkel','c_MSDAg[4]','PotEng','KinEng','TotEng','Press','c_fmax']):active=True;continue
    if active:
        v=line.split()
        if len(v)==11 and v[0].isdigit():
            try:thermo.append([float(x) for x in v])
            except ValueError:pass
        elif thermo:active=False
relax=re.search(r'Energy initial, next-to-last, final =\s*([+\-0-9.eE]+)\s+([+\-0-9.eE]+)\s+([+\-0-9.eE]+)',text)
frelax=re.search(r'Force max component initial, final =\s*([+\-0-9.eE]+)\s+([+\-0-9.eE]+)',text)
sym=np.asarray(CELL.get_chemical_symbols());counts={s:int(np.sum(sym==s)) for s in sorted(set(sym))}
m={'Ag':107.8682,'C':12.011,'Si':28.0855,'Ti':47.867};mag=counts['Ag']*m['Ag'];msk=sum(counts[s]*m[s] for s in ('C','Si','Ti'))
frames=[read(INIT,index=0,format='lammps-dump-text',specorder=['Ag','C','Si','Ti'])]+read(TRAJ,index=':',format='lammps-dump-text',specorder=['Ag','C','Si','Ti'])
result={'status':'completed','purpose':'exploratory imposed-temperature Ag melt/cool screening; not spontaneous laser/TTM',
 'source_geometry':str(SRC),'replication':[4,4,1],'atoms_initial':len(CELL),'element_counts':counts,'Ag_to_skeleton_mass_ratio':mag/msk,
 'dt_fs':1.0,'potential':'Foiles Ag_u3 EAM + TSCBOP Ti3SiC2 + force-calibrated Ag/skeleton Morse + cross-only ZBL inner1.5 outer2.4 A',
 'schedule_ps':{'300K_equilibration':0.2,'Ag_ramp_300_to_target':HEAT_PS,'target_Ag_K':TARGET,'Ag_target_hold':HOLD_PS,
                'Ag_cool_to_300K':COOL_PS,'Ag_final_300K_hold':FINAL_PS,
                'total':0.2+HEAT_PS+HOLD_PS+COOL_PS+FINAL_PS},
 'skeleton_thermostat':'NVT at 300 K throughout as an ideal cold sink; Ag receives a separate imposed NVT schedule',
 'FIRE_energy_initial_nexttolast_final_eV':list(map(float,relax.groups())) if relax else None,
 'FIRE_fmax_component_initial_final_eV_A':list(map(float,frelax.groups())) if frelax else None,
 'thermo_columns':['step','time_ps','T_all_K','T_Ag_K','T_skeleton_K','MSD_Ag_A2','PE_eV','KE_eV','Etot_eV','pressure_bar','Fmax_eV_A'],
 'thermo_rows':thermo,'structure_rows':analyze(frames,0,50),
 'limitations':['fixed volume','1x ideal periodic interface repeated 4x4, not the target porous model','direct Ag NVT ramp is only a potential check and replaces TTM','Q6 uses Ag-Ag neighbors within 3.8 A','long cooling and independent validation remain outstanding']}
(OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'atoms':len(CELL),'counts':counts,'mass_ratio':mag/msk,
 'FIRE_energy':result['FIRE_energy_initial_nexttolast_final_eV'],'FIRE_Fmax':result['FIRE_fmax_component_initial_final_eV_A'],
 'thermo_rows':len(thermo),'structure_rows':len(result['structure_rows']),'summary':str(OUT/'summary.json')},indent=2))
