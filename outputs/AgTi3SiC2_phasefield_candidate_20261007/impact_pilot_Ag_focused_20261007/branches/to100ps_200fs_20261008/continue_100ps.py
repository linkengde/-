#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math, os, re, shutil, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
LMP=Path('/workspace/.local/lammps-static-30Sep2026/lammps-static/bin/lmp')
SOURCE_HIGHFLUX=ROOT.parent/'highflux_200fs_20261008'
TEMPLATE=ROOT/'inputs/impact_template.lmp'
CENTER=np.array([9.9573898076,24.2578814220])
KB=8.617333262e-5
MIN_LAUNCH_CLEARANCE_A=2.5
sys.path.insert(0,str(ROOT/'audit'))
import audit_crater

def positions(first,last,seed=20261012):
    rng=np.random.default_rng(seed)
    out=[]
    R=7.5; lmax=math.log1p((R/6.0)**2)
    # Use one global seeded stream across batch invocations; advancing through
    # skipped shots prevents every 25-shot batch from repeating the same map.
    for _ in range(first-1):
        rng.random(); rng.uniform(0,2*math.pi)
    for hit in range(first,last+1):
        u=float(rng.random()); r=6.0*math.sqrt(math.exp(u*lmax)-1.0)
        theta=float(rng.uniform(0,2*math.pi))
        x,y=CENTER+r*np.array([math.cos(theta),math.sin(theta)])
        out.append((hit,r,theta,x,y,60.0,62,-105.3164,200))
    return out

def launch_clearance(row):
    d=audit_crater.read_data(ROOT/'structures/state_current.data')
    _,_,_,x,y,z,_,_,_=row
    target=np.linalg.solve(d['H'],np.array([x-d['bounds']['x'][0],y-d['bounds']['y'][0]]))
    df=d['frac']-target
    df-=np.round(df)
    dxy=df@d['H'].T
    dist=np.sqrt(np.sum(dxy*dxy,axis=1)+(d['xyz'][:,2]-z)**2)
    k=int(np.argmin(dist))
    return float(dist[k]),int(d['ids'][k]),int(d['typ'][k])

def choose_safe_launch(row):
    hit,r,theta,x,y,z,energy,vz,spacing=row
    clearance,near_id,near_type=launch_clearance(row)
    if clearance>=MIN_LAUNCH_CLEARANCE_A:
        return row,clearance,0,near_id,near_type
    rng=np.random.default_rng(np.random.SeedSequence([20261012,hit,3232026]))
    R=7.5; lmax=math.log1p((R/6.0)**2)
    for attempt in range(1,10001):
        u=float(rng.random()); r=6.0*math.sqrt(math.exp(u*lmax)-1.0)
        theta=float(rng.uniform(0,2*math.pi))
        x,y=CENTER+r*np.array([math.cos(theta),math.sin(theta)])
        candidate=(hit,r,theta,x,y,z,energy,vz,spacing)
        clearance,near_id,near_type=launch_clearance(candidate)
        if clearance>=MIN_LAUNCH_CLEARANCE_A:
            return candidate,clearance,attempt,near_id,near_type
    raise RuntimeError(f'could not find a non-overlapping launch point for hit {hit}')

def make_input(hit,row):
    row,clearance,resamples,near_id,near_type=choose_safe_launch(row)
    _,r,theta,x,y,z,energy,vz,spacing=row
    s=TEMPLATE.read_text()
    repl={
      'read_data structures/impact_54_62eV.data':'read_data structures/state_current.data',
      'group projectiles id 4306:4359':f'group projectiles id 4306:{4304+hit}',
      'group new_projectile id 4360':f'group new_projectile id {4305+hit}',
      '623755':str(623700+hit),
      'dump tr all custom 1000 ': 'dump tr all custom 4000 ',
      'dumps/impact_55_62eV.dump':f'dumps/impact_{hit}_62eV.dump',
      'write_data structures/impact_55_62eV.data':'write_data structures/state_current.data',
      'write_restart restart/impact_55_62eV.restart':'write_restart restart/state_current.restart',
    }
    for a,b in repl.items():
        if a not in s: raise RuntimeError(f'template token missing: {a}')
        s=s.replace(a,b)
    s=re.sub(r'^create_atoms 4 single .* units box$',f'create_atoms 4 single {x:.10f} {y:.10f} {z:.1f} units box',s,flags=re.M)
    if 'impact_55_62eV' in s: raise RuntimeError('old stage number remains in generated input')
    path=ROOT/'inputs'/f'impact_{hit}_62eV.lmp'
    path.write_text(s)
    return path,row,clearance,resamples,near_id,near_type

def parse_pairs(text):
    d={}
    for phase,val in re.findall(r'^\s*(Ag-Ag|Ag-TSC|C-C): min=([0-9.]+) A',text,re.M): d[phase]=float(val)
    return d

def checkpoint(hit,dump_path):
    tag=f'after_hit_{hit:03d}'
    shutil.copy2(ROOT/'structures/state_current.data',ROOT/'checkpoints'/f'{tag}.data')
    shutil.copy2(ROOT/'restart/state_current.restart',ROOT/'checkpoints'/f'{tag}.restart')
    shutil.copy2(dump_path,ROOT/'checkpoints'/f'{tag}.dump')

def one_hit(hit,row,threads):
    inp,row,clearance,resamples,near_id,near_type=make_input(hit,row)
    name=f'impact_{hit}_62eV'
    _,r,theta,x,y,z,energy,vz,spacing=row
    used_path=ROOT/'impact_positions_used_323_446.csv'
    if not used_path.exists():
        with used_path.open('w',newline='') as f:
            csv.writer(f).writerow(['hit','radius_A','theta_rad','x_A','y_A','launch_z_A','energy_eV','Ag_velocity_z_A_per_ps','spacing_fs','min_initial_distance_A','resamples','nearest_initial_atom_id','nearest_initial_atom_type'])
    with used_path.open('a',newline='') as f:
        csv.writer(f).writerow([hit,f'{r:.10f}',f'{theta:.12f}',f'{x:.10f}',f'{y:.10f}',f'{z:.3f}',energy,f'{vz:.4f}',spacing,f'{clearance:.6f}',resamples,near_id,near_type]); f.flush(); os.fsync(f.fileno())
    log=ROOT/'logs'/f'{name}.lammps'
    console=ROOT/'logs'/f'{name}.console.txt'
    dump=ROOT/'dumps'/f'{name}.dump'
    cmd=[str(LMP),'-sf','omp','-pk','omp',str(threads),'-in',f'inputs/{name}.lmp','-log',f'logs/{name}.lammps']
    env=os.environ.copy(); env['OMP_NUM_THREADS']=str(threads)
    with console.open('w') as f:
        p=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,env=env)
    if p.returncode!=0: raise RuntimeError(f'LAMMPS failed at hit {hit}: {console}')
    if not log.is_file() or not dump.is_file() or not (ROOT/'structures/state_current.data').is_file() or not (ROOT/'restart/state_current.restart').is_file():
        raise RuntimeError(f'missing stage output at hit {hit}')
    logtext=log.read_text(errors='replace')
    if 'ERROR:' in logtext or 'Lost atoms' in logtext or re.search(r'\b(?:nan|inf)\b',logtext) or 'Total wall time:' not in logtext:
        raise RuntimeError(f'LAMMPS log check failed at hit {hit}')
    stage_data=ROOT/'structures'/f'{name}.data'
    shutil.copy2(ROOT/'structures/state_current.data',stage_data)
    env['MPLCONFIGDIR']='/tmp/agti-mpl/config'; env['XDG_CACHE_HOME']='/tmp/agti-mpl/cache'
    audit_txt=ROOT/'audit'/f'audit_hit_{hit}.txt'
    with audit_txt.open('w') as f:
        q=subprocess.run([sys.executable,'audit/audit_stage.py',str(hit)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,env=env)
    if q.returncode!=0: raise RuntimeError(f'audit failed at hit {hit}: {audit_txt}')
    js=json.loads((ROOT/'audit'/f'audit_hit_{hit:02d}.json').read_text())
    pair=parse_pairs(audit_txt.read_text(errors='replace'))
    if js['atoms'] != 4305+hit: raise RuntimeError(f'atom count mismatch at {hit}: {js["atoms"]}')
    fmax=js['max_force_eV_per_A']
    tc=js['connectivity']['TSC']['fraction']; ac=js['connectivity']['Ag']['fraction']
    if fmax>100.0: raise RuntimeError(f'max force exceeded 100 eV/A at hit {hit}: {fmax:.3f}')
    if tc<0.99 or ac<0.95: raise RuntimeError(f'connectivity threshold crossed at hit {hit}: TSC={tc:.5f}, Ag={ac:.5f}')
    for key,limit in (('Ag-Ag',1.0),('Ag-TSC',1.0),('C-C',1.10)):
        if pair.get(key,999.0)<limit: raise RuntimeError(f'{key} distance below {limit} A at hit {hit}: {pair.get(key)}')
    ta=js['local_kinetic_temperature_K']['r_lt_8']['Ag_matrix']['T_K']
    tt=js['local_kinetic_temperature_K']['r_lt_8']['TSC']['T_K']
    exposure=21.8+0.2*(hit-55)
    rowout=[hit,f'{exposure:.3f}',js['atoms'],f'{fmax:.6f}',f'{ta or 0:.3f}',f'{tt or 0:.3f}',f'{tc:.6f}',f'{ac:.6f}',f'{pair.get("Ag-Ag",float("nan")):.6f}',f'{pair.get("Ag-TSC",float("nan")):.6f}',f'{pair.get("C-C",float("nan")):.6f}']
    with (ROOT/'run_progress.csv').open('a',newline='') as f:
        csv.writer(f).writerow(rowout); f.flush(); os.fsync(f.fileno())
    do_checkpoint=((hit-55)%25==0 or hit==446)
    if do_checkpoint:
        checkpoint(hit,dump)
    print(f'PASS hit={hit} exposure={exposure:.1f}ps atoms={js["atoms"]} maxF={fmax:.1f} TAg={ta:.0f}K TTSC={tt:.0f}K conn={tc:.4f}/{ac:.4f} minCC={pair.get("C-C",float("nan")):.3f}A launch_clearance={clearance:.2f}A resamples={resamples} checkpoint={do_checkpoint}',flush=True)
    if not do_checkpoint:
        stage_data.unlink(missing_ok=True); dump.unlink(missing_ok=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--first',type=int,required=True)
    ap.add_argument('--last',type=int,required=True)
    ap.add_argument('--threads',type=int,default=5)
    args=ap.parse_args()
    if not LMP.is_file(): raise FileNotFoundError(LMP)
    if args.first<56 or args.last>446 or args.first>args.last: raise ValueError('requested range must be within 56..446')
    pp=positions(args.first,args.last)
    csvpath=ROOT/f'impact_positions_candidates_{args.first}_{args.last}.csv'
    with csvpath.open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['hit','radius_A','theta_rad','x_A','y_A','launch_z_A','energy_eV','Ag_velocity_z_A_per_ps','spacing_fs']); w.writerows(pp)
    for row in pp:
        one_hit(row[0],row,args.threads)
    print(f'BATCH_COMPLETE {args.first}-{args.last}; cumulative bombardment time {21.8+0.2*(args.last-55):.1f} ps',flush=True)
if __name__=='__main__': main()
