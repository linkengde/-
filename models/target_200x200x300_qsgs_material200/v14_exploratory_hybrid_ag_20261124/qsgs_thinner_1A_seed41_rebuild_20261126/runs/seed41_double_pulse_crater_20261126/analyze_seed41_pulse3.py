#!/usr/bin/env python3
"""Analyze third same-dose impulse and its 2 ps NVE follow-up."""
from pathlib import Path
import csv,json
import numpy as np
from analyze_seed41_double_pulse import dump_frames,surface_metrics

ROOT=Path(__file__).resolve().parent
def sim_time(step):
    if step<=48400: return 10.8+(step-40400)*0.00025
    return 12.8+(step-48400)*0.00025

def main():
    frames={}
    for name in ('seed41_pulse3_ag.dump','seed41_afterpulse3_NVE_ag_nomsd.dump'):
        for fr in dump_frames(ROOT/name): frames[fr[0]]=fr
    traj=[]; qrows=[]
    first=frames[min(frames)]; c0=first[1]; a0=first[2]
    ids0=a0[:,c0['id']].astype(np.int64)
    x0=a0[:,c0['x']]-100.; y0=a0[:,c0['y']]-100.
    local_ids=set(ids0[np.hypot(x0,y0)<25].tolist())
    prev={int(row[c0['id']]):row[[c0['x'],c0['y'],c0['z']]].copy() for row in a0}
    disp={int(i):np.zeros(3) for i in ids0}
    for j,step in enumerate(sorted(frames)):
        fr=frames[step]; idx=fr[1]; a=fr[2]; ids=a[:,idx['id']].astype(np.int64)
        pos=a[:,[idx['x'],idx['y'],idx['z']]]
        if j:
            newprev={}
            for atomid,p in zip(ids,pos):
                atomid=int(atomid)
                if atomid not in prev: continue
                d=p-prev[atomid]
                d[:2]=(d[:2]+100.)%200.-100.
                disp[atomid]+=d
                newprev[atomid]=p.copy()
            prev=newprev
        present=set(ids.tolist())
        valid=[int(i) for i in ids0 if int(i) in present]
        local=[i for i in valid if i in local_ids]
        row=surface_metrics(fr); row['time_ps']=sim_time(step); traj.append(row)
        qkey=next((k for k in idx if k.startswith('c_q6Ag')),None)
        if qkey:
            q=a[:,idx[qkey]]; x=a[:,idx['x']]-100.; y=a[:,idx['y']]-100.
            now_local=np.hypot(x,y)<25
            qrows.append({'step':step,'time_ps':sim_time(step),'Q6mean':float(np.nanmean(q)),
                'Q6_gt_0p30_fraction':float(np.nanmean(q>0.30)),
                'Q6mean_center_r25':float(np.nanmean(q[now_local])) if np.any(now_local) else None,
                'MSD_Ag_A2':float(np.mean([disp[i]@disp[i] for i in valid])) if valid else None,
                'MSD_initial_center_r25_A2':float(np.mean([disp[i]@disp[i] for i in local])) if local else None,
                'surviving_initial_Ag':len(valid),'surviving_initial_center_Ag':len(local)})
    endpoints=[]
    for name in ('seed41_afterpulse2_final_10p8ps_nomsd_all.dump','seed41_pulse3_final_12p8ps_all.dump','seed41_afterpulse3_final_14p8ps_nomsd_all.dump'):
        fr=dump_frames(ROOT/name)[0]; m=surface_metrics(fr); m['time_ps']=sim_time(m['step']); endpoints.append(m)
    thermo=[]
    for fname in ('log.seed41_pulse3_1479eV_2ps.lammps','log.seed41_afterpulse3_nve_12p8_to14p8ps_nomsd.lammps'):
        for line in (ROOT/fname).read_text(errors='replace').splitlines():
            w=line.split()
            if not w or not w[0].isdigit() or len(w)<18: continue
            try:
                thermo.append({'step':int(w[0]),'time_ps':float(w[1]),'TAg_K':float(w[2]),'Tspot_K':float(w[3]),
                    'PE_eV':float(w[4]),'KE_eV':float(w[5]),'Etot_eV':float(w[6]),'FmaxAg_eV_A':float(w[7]),
                    'Q6mean':float(w[8]),'Q6_gt_0p30_fraction':float(w[9]),'zmax_A':float(w[10]),
                    'N_gt220':int(float(w[11])),'N_gt250':int(float(w[12])),'N_gt300':int(float(w[13])),
                    'N_gt400':int(float(w[14])),'N_center_r25_gt220':int(float(w[15])),
                    'N_center_r25_gt300':int(float(w[16])),'atoms':int(float(w[-1]))})
            except (ValueError,IndexError): pass
    thermo=list({r['step']:r for r in thermo}.values()); thermo.sort(key=lambda r:r['time_ps'])
    report={'run':'seed41 third localized velocity impulse + 2 ps NVE; exploratory, no TTM',
        'pulse_energy_eV':1479.2492958829243,'seed':20261143,'time_step_fs':0.25,
        'pulse_profile':{'center_A':[100,100,180.09],'sigma_xy_A':12,'attenuation_depth_A':6},
        'start_atoms':469038,'after_pulse3_atoms':469035,
        'after_NVE_atoms':endpoints[-1]['n_Ag']+358428,
        'atom_losses_during_pulse3':endpoints[0]['n_Ag']-endpoints[1]['n_Ag'],
        'atom_losses_during_NVE3':endpoints[1]['n_Ag']-endpoints[2]['n_Ag'],
        'endpoints':endpoints,'trajectory_frames':len(traj),'thermo_records':len(thermo),
        'thermo_first':thermo[0] if thermo else None,'thermo_last':thermo[-1] if thermo else None,
        'max_FmaxAg_eV_A':max((r['FmaxAg_eV_A'] for r in thermo),default=None),
        'max_MSD_Ag_A2':max((r['MSD_Ag_A2'] for r in qrows if r['MSD_Ag_A2'] is not None),default=None),
        'max_MSD_initial_center_r25_A2':max((r['MSD_initial_center_r25_A2'] for r in qrows if r['MSD_initial_center_r25_A2'] is not None),default=None),
        'caveat':'Ag atoms crossing zhi=500 A are removed; count indicates finite-cell ejection/escape, not evaporation yield. Frozen skeleton, repulsive Ag-skeleton LJ cross interactions, localized velocity impulse rather than TTM.'}
    (ROOT/'seed41_pulse3_summary.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    with (ROOT/'seed41_pulse3_trajectory.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=traj[0].keys()); w.writeheader(); w.writerows(traj)
    with (ROOT/'seed41_pulse3_q6_msd.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=qrows[0].keys()); w.writeheader(); w.writerows(qrows)
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=='__main__': main()
