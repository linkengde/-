#!/usr/bin/env python3
"""Summarize the second seed41 localized-pulse experiment."""
from pathlib import Path
import csv, json, math
import numpy as np

ROOT = Path(__file__).resolve().parent
LXY = 200.0

def dump_frames(path, type_filter=4):
    out=[]
    with open(path) as f:
        while True:
            line=f.readline()
            if not line: break
            if not line.startswith('ITEM: TIMESTEP'): continue
            step=int(f.readline().strip())
            f.readline(); n=int(f.readline().strip())
            f.readline(); bounds=[tuple(map(float,f.readline().split()[:2])) for _ in range(3)]
            cols=f.readline().strip().split()[2:]
            idx={k:i for i,k in enumerate(cols)}
            data=np.empty((n,len(cols)),dtype=np.float64)
            for i in range(n): data[i]=np.fromstring(f.readline(),sep=' ')
            if type_filter is not None and 'type' in idx:
                data=data[data[:,idx['type']]==type_filter]
            out.append((step,idx,data,bounds))
    return out

def surface_metrics(frame, bw=4.0):
    step,idx,a,_=frame
    x=a[:,idx['x']]-100.; y=a[:,idx['y']]-100.; z=a[:,idx['z']]
    r=np.hypot(x,y)
    keep=(np.abs(x)<50)&(np.abs(y)<50)&(z>=150)&(z<=230)
    n=int(round(100/bw)); h=np.full((n,n),np.nan)
    ix=np.floor((x[keep]+50)/bw).astype(int); iy=np.floor((y[keep]+50)/bw).astype(int)
    for bx,by,bz in zip(ix,iy,z[keep]):
        if 0<=bx<n and 0<=by<n and (np.isnan(h[by,bx]) or bz>h[by,bx]): h[by,bx]=bz
    center=[]; ring=[]
    for by,bx in zip(*np.where(np.isfinite(h))):
        rr=math.hypot((bx+.5)*bw-50,(by+.5)*bw-50)
        if rr<10: center.append(h[by,bx])
        elif 20<=rr<40: ring.append(h[by,bx])
    return {
        'step':step,'time_ps':sim_time(step),'n_Ag':len(a),
        'center_top_median_A':float(np.median(center)) if center else None,
        'ring_top_median_A':float(np.median(ring)) if ring else None,
        'center_minus_ring_A':float(np.median(center)-np.median(ring)) if center and ring else None,
        'center_r20_z150_220_Ag':int(np.count_nonzero((r<20)&(z>=150)&(z<=220))),
        'n_z_gt_220':int(np.count_nonzero(z>220)),
        'n_z_gt_250':int(np.count_nonzero(z>250)),
        'n_z_gt_300':int(np.count_nonzero(z>300)),
        'n_z_gt_400':int(np.count_nonzero(z>400)),
        'zmax_A':float(z.max()),
    }

def sim_time(step):
    # LAMMPS reset the step origin when the second NVE continuation began.
    if step <= 28400: return 5.8+(step-20400)*0.00025
    return 7.8+(step-28400)*0.00025

def main():
    frames={}
    for name in ('seed41_pulse2_ag.dump','seed41_afterpulse2_NVE_ag_nomsd.dump'):
        for fr in dump_frames(ROOT/name): frames[fr[0]]=fr
    traj=[frames[k] for k in sorted(frames)]
    # Integrate wrapped x/y displacements frame-by-frame to estimate MSD.
    first=traj[0]; c0=first[1]; a0=first[2]
    start_ids=a0[:,c0['id']].astype(np.int64)
    prev={int(row[c0['id']]):row[[c0['x'],c0['y'],c0['z']]].copy() for row in a0}
    disp={i:np.zeros(3,dtype=np.float64) for i in start_ids}
    x0=a0[:,c0['x']]-100.; y0=a0[:,c0['y']]-100.
    local_ids=set(start_ids[np.hypot(x0,y0)<25.0].tolist())
    qrows=[]; trajectory=[]
    for j,fr in enumerate(traj):
        step,idx,a,_=fr
        ids=a[:,idx['id']].astype(np.int64)
        pos=a[:,[idx['x'],idx['y'],idx['z']]]
        if j:
            newprev={}
            for atomid,p in zip(ids,pos):
                atomid=int(atomid)
                if atomid not in prev: continue
                delta=p-prev[atomid]
                delta[:2]=(delta[:2]+LXY/2)%LXY-LXY/2
                disp[atomid]+=delta
                newprev[atomid]=p.copy()
            prev=newprev
        present=set(ids.tolist())
        valid=[i for i in start_ids if int(i) in present]
        local=[i for i in valid if int(i) in local_ids]
        msd_all=float(np.mean([disp[int(i)]@disp[int(i)] for i in valid])) if valid else None
        msd_local=float(np.mean([disp[int(i)]@disp[int(i)] for i in local])) if local else None
        metric=surface_metrics(fr)
        trajectory.append(metric)
        qkey=next((k for k in idx if k.startswith('c_q6Ag')),None)
        if qkey:
            q=a[:,idx[qkey]]
            cx=a[:,idx['x']]-100.; cy=a[:,idx['y']]-100.
            local_now=np.hypot(cx,cy)<25.0
            qrows.append({'step':step,'time_ps':sim_time(step),
                'Q6mean':float(np.nanmean(q)),
                'Q6_gt_0p30_fraction':float(np.nanmean(q>0.30)),
                'Q6mean_center_r25':float(np.nanmean(q[local_now])) if np.any(local_now) else None,
                'MSD_Ag_A2':msd_all,'MSD_initial_center_r25_A2':msd_local,
                'surviving_initial_Ag':len(valid),'surviving_initial_center_Ag':len(local)})
    # Stable endpoint states for surface change comparisons.
    p580=ROOT.parent/'seed41_local_pulse_crater_20261126'/'seed41_local_pulse_final_5p8ps_all.dump'
    p780=ROOT/'seed41_pulse2_final_7p8ps_all.dump'
    p1080=ROOT/'seed41_afterpulse2_final_10p8ps_nomsd_all.dump'
    endpoints=[]
    for p in (p580,p780,p1080): endpoints.append(surface_metrics(dump_frames(p)[0]))
    # Collect thermo rows from the pulse2 and NVE logs.
    thermo=[]
    for fname in ('log.seed41_pulse2_1479eV_2ps.lammps','log.seed41_afterpulse2_nve_7p8_to10p8ps_nomsd.lammps'):
        has_msd='pulse2_1479eV' in fname
        for line in (ROOT/fname).read_text(errors='replace').splitlines():
            w=line.split()
            if not w or not w[0].isdigit() or len(w)<18: continue
            try:
                zidx=11 if has_msd else 10
                nidx=zidx+1
                row={'step':int(w[0]),'time_ps':float(w[1]),'TAg_K':float(w[2]),'Tspot_K':float(w[3]),
                     'PE_eV':float(w[4]),'KE_eV':float(w[5]),'Etot_eV':float(w[6]),'FmaxAg_eV_A':float(w[7]),
                     'Q6mean':float(w[8]),'Q6_gt_0p30_fraction':float(w[9]),'zmax_A':float(w[zidx]),
                     'N_gt220':int(float(w[nidx])),'N_gt250':int(float(w[nidx+1])),'N_gt300':int(float(w[nidx+2])),
                     'N_gt400':int(float(w[nidx+3])) if not has_msd else None,
                     'N_center_r25_gt220':int(float(w[nidx+4])) if not has_msd else int(float(w[nidx+3])),
                     'N_center_r25_gt300':int(float(w[nidx+5])) if not has_msd else int(float(w[nidx+5])),
                     'atoms':int(float(w[-1]))}
                thermo.append(row)
            except (ValueError,IndexError): pass
    thermo=list({r['step']:r for r in thermo}.values()); thermo.sort(key=lambda x:x['step'])
    report={'run':'seed41 second local velocity impulse + NVE, exploratory; no TTM',
        'pulse_energy_eV':1479.2492958829243,'timestep_fs':0.25,
        'box_A':[200,200,500],'initial_atoms':469040,'final_atoms':469038,
        'lost_Ag_at_zhi':2,'pulse2_seed':20261142,
        'pulse2_profile':{'center_A':[100,100,180.09],'sigma_xy_A':12,'attenuation_depth_A':6},
        'endpoints':endpoints,'trajectory_frames':len(traj),'thermo_records':len(thermo),
        'thermo_first':thermo[0] if thermo else None,'thermo_last':thermo[-1] if thermo else None,
        'etotal_change_including_atom_losses_eV':(thermo[-1]['Etot_eV']-thermo[0]['Etot_eV']) if thermo else None,
        'max_Q6_ag_mean':max((x['Q6mean'] for x in qrows),default=None),
        'max_Q6_gt_0p30_fraction':max((x['Q6_gt_0p30_fraction'] for x in qrows),default=None),
        'max_MSD_Ag_A2':max((x['MSD_Ag_A2'] for x in qrows if x['MSD_Ag_A2'] is not None),default=None),
        'max_MSD_initial_center_r25_A2':max((x['MSD_initial_center_r25_A2'] for x in qrows if x['MSD_initial_center_r25_A2'] is not None),default=None),
        'caveat':'Two atoms crossed the finite zhi=500 A boundary and were removed; these are ejection/escape events in this finite-cell model, not quantitative evaporation. Skeleton is frozen; Ag-skeleton interaction is repulsive LJ, and heating is a localized velocity impulse, not TTM. Q6 and endpoint MSD are exploratory disorder/mobility indicators; not standalone proof of a liquid phase.'}
    (ROOT/'seed41_double_pulse_summary.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    with (ROOT/'seed41_double_pulse_trajectory.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=trajectory[0].keys()); w.writeheader(); w.writerows(trajectory)
    with (ROOT/'seed41_double_pulse_q6_msd.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=qrows[0].keys()); w.writeheader(); w.writerows(qrows)
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=='__main__': main()
