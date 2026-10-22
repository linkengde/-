#!/usr/bin/env python3
"""MACE V17 CPU-only single-local-source exploratory Ag/Ti3SiC2 run."""
from __future__ import annotations
import csv, hashlib, json, os, resource, sys, time
from pathlib import Path
import numpy as np
import torch
from ase import units
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import FIRE
from ase.neighborlist import neighbor_list
from mace.calculators import MACECalculator
from scipy.special import sph_harm_y
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT=Path(__file__).resolve().parent
REPO=next(p for p in (ROOT,*ROOT.parents) if (p/".git").exists())
MODEL=REPO/"research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_03/models/MACE_periodic_v17_provisional_screening.model"
EXPECTED="559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3"
DT_FS=0.5
FRICTION_FS=0.02
HOT_K=1900.0
COLD_K=300.0
SEED=202610212
MAX_HEAT_PS=4.0
MAX_HEAT_WALL_S=12*3600

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def attach(a):
    got=sha(MODEL)
    if got!=EXPECTED:raise RuntimeError(f"V17 SHA mismatch: {got}")
    torch.set_num_threads(4)
    a.calc=MACECalculator(model_paths=str(MODEL),device="cpu",default_dtype="float64",enable_cueq=False)
    return got
def masks(a):
    ph=np.asarray(a.arrays["phase"])
    return ph=="Ag",ph=="Ti3SiC2",np.asarray(a.arrays["hot_source"]).astype(bool),np.asarray(a.arrays["tsc_heat_sink"]).astype(bool)
def q6_and_coord(a,cutoff=3.7):
    ag,_,_,_=masks(a); idx=np.flatnonzero(ag); slot=np.full(len(a),-1,int);slot[idx]=np.arange(len(idx))
    i,j,S,D=neighbor_list("ijSD",a,cutoff)
    m=ag[i]&ag[j]
    ii,jj,dd=i[m],j[m],D[m]
    coord=np.bincount(ii,minlength=len(a))[idx]
    q=np.full(len(a),np.nan,float)
    by=[[] for _ in range(len(a))]
    for u,v,d in zip(ii,jj,dd):by[int(u)].append(d)
    l=6; fac=4*np.pi/(2*l+1)
    for u in idx:
        vec=np.asarray(by[int(u)])
        if len(vec)==0:continue
        norm=np.linalg.norm(vec,axis=1); z=np.clip(vec[:,2]/norm,-1,1)
        theta=np.arccos(z);phi=np.arctan2(vec[:,1],vec[:,0])
        qlm=[np.mean(sph_harm_y(l,m,theta,phi)) for m in range(-l,l+1)]
        q[u]=np.sqrt(fac*np.sum(np.abs(qlm)**2))
    co=np.zeros(len(a),int);co[idx]=coord
    return q,co
def cluster_stats(a,cand,cut=3.6):
    idx=np.flatnonzero(cand);n=len(a)
    if not len(idx):return 0,0
    i,j=neighbor_list("ij",a,cut)
    m=cand[i]&cand[j]&(i!=j)
    mp=np.full(n,-1,int);mp[idx]=np.arange(len(idx));u,v=mp[i[m]],mp[j[m]]
    mat=coo_matrix((np.ones(len(u),int),(u,v)),shape=(len(idx),len(idx))).tocsr()
    nc,lab=connected_components(mat,directed=False)
    sizes=np.bincount(lab)
    return int(sizes.max()),int(nc)
def min_pairs(a):
    sym=np.asarray(a.get_chemical_symbols());i,j,d=neighbor_list("ijd",a,4.6);m=i<j
    pairs={"Ti-C":((sym[i]=="Ti")&(sym[j]=="C"))|((sym[i]=="C")&(sym[j]=="Ti")),
           "Ti-Si":((sym[i]=="Ti")&(sym[j]=="Si"))|((sym[i]=="Si")&(sym[j]=="Ti")),
           "Ag-Ti":((sym[i]=="Ag")&(sym[j]=="Ti"))|((sym[i]=="Ti")&(sym[j]=="Ag")),
           "Ag-Si":((sym[i]=="Ag")&(sym[j]=="Si"))|((sym[i]=="Si")&(sym[j]=="Ag")),
           "Ag-C":((sym[i]=="Ag")&(sym[j]=="C"))|((sym[i]=="C")&(sym[j]=="Ag")),
           "Ag-Ag":(sym[i]=="Ag")&(sym[j]=="Ag")}
    out={}
    for name,p in pairs.items():
        mm=m&p
        out["min_"+name.replace("-","_")+"_A"]=float(d[mm].min()) if np.any(mm) else None
    return out
def temps(a,mask):
    ids=np.flatnonzero(mask)
    if len(ids)<2:return float("nan")
    v=a.get_velocities()[ids];m=a.get_masses()[ids,None];vcm=(m*v).sum(0)/m.sum();dof=max(1,3*len(ids)-3)
    return float(np.sum(m*(v-vcm)**2)/(dof*units.kB))
def metrics(a,refpos,step,stage):
    ag,tsc,hot,sink=masks(a);forces=np.asarray(a.get_forces());fn=np.linalg.norm(forces,axis=1)
    if not np.isfinite(fn).all() or not np.isfinite(a.positions).all():raise FloatingPointError("non-finite coordinate/force")
    q,coord=q6_and_coord(a)
    disp=np.linalg.norm(a.positions-refpos,axis=1)
    # Q6 threshold is paired with substantial per-atom migration and a usable
    # Ag coordination, preventing surface under-coordination alone counting.
    cand=ag&(coord>=8)&(q<0.38)&(disp>=1.4)
    maxc,nc=cluster_stats(a,cand)
    cand_hot=cand&hot; maxhot,nchot=cluster_stats(a,cand_hot)
    outside=ag&~hot
    cand_out=cand&outside;maxout,nout=cluster_stats(a,cand_out)
    row={"stage":stage,"step":int(step),"time_ps":float(step*DT_FS/1000),"T_system_K":float(a.get_temperature()),
      "T_Ag_K":temps(a,ag),"T_hot_Ag_K":temps(a,hot),"T_TSC_K":temps(a,tsc),"T_TSC_sink_K":temps(a,sink),
      "E_potential_eV":float(a.get_potential_energy()),"Fmax_eV_A":float(fn.max()),"P99_eV_A":float(np.percentile(fn,99)),
      "Fmax_Ag_eV_A":float(fn[ag].max()),"Fmax_TSC_eV_A":float(fn[tsc].max()),
      "hot_Ag_count":int(hot.sum()),"liquid_candidate_Ag":int(cand.sum()),"max_liquid_cluster_Ag":maxc,"liquid_cluster_count":nc,
      "hot_liquid_candidate_Ag":int(cand_hot.sum()),"hot_max_liquid_cluster_Ag":maxhot,"hot_liquid_cluster_count":nchot,
      "outside_hot_liquid_candidate_Ag":int(cand_out.sum()),"outside_hot_max_cluster_Ag":maxout,"outside_hot_cluster_count":nout,
      "q6_hot_median":float(np.nanmedian(q[hot])),"q6_hot_p25":float(np.nanpercentile(q[hot],25)),
      "coord_Ag_hot_median":float(np.median(coord[hot])),"hot_Ag_displacement_rms_A":float(np.sqrt(np.mean(disp[hot]**2))),
      "n_hot_displacement_ge_1p4_A":int(np.sum(hot&(disp>=1.4))),**min_pairs(a)}
    return row,forces,fn,q,coord,disp,cand
def write_frame(path,a,row,forces,fn,q,coord,disp,cand):
    a.arrays["mace_force_eV_A"]=forces.copy();a.arrays["mace_force_norm_eV_A"]=fn.copy()
    a.arrays["Ag_q6"]=np.nan_to_num(q,nan=-1.0);a.arrays["Ag_coordination"]=coord.astype(int)
    a.arrays["displacement_from_reference_A"]=disp
    a.arrays["liquid_candidate"]=cand.astype(np.int8)
    a.info.update({k:(v.item() if isinstance(v,np.generic) else v) for k,v in row.items() if np.isscalar(v)})
    write(path,a,format="extxyz",append=Path(path).exists())
def rss_gib():return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024/1024
def severe(row):
    return row["Fmax_eV_A"]>100 or row["min_Ti_C_A"] is None or row["min_Ti_C_A"]<1.45 or row["min_Ag_Ag_A"]<1.70 or min([row[k] for k in ("min_Ag_Ti_A","min_Ag_Si_A","min_Ag_C_A") if row[k] is not None],default=10)<1.70

def persistent_melt(rows, hot_count):
    if len(rows)<3:return False
    need=max(12,int(np.ceil(0.32*hot_count)))
    tail=rows[-3:]
    return all(r["hot_max_liquid_cluster_Ag"]>=need and r["q6_hot_median"]<0.42 and
               r["hot_Ag_displacement_rms_A"]>=1.5 and r["n_hot_displacement_ge_1p4_A"]>=need
               for r in tail)

def persistent_solidification(rows):
    if len(rows)<3:return False
    tail=rows[-3:]
    return all(r["T_hot_Ag_K"]<1235.0 and r["q6_hot_median"]>0.48 and
               r["hot_max_liquid_cluster_Ag"]<10 for r in tail)
def main():
    stage=sys.argv[1] if len(sys.argv)>1 else "all"
    if stage not in {"all","static","baseline","heat","cool"}:raise SystemExit("stage: all|static|baseline|heat|cool")
    atomfile=ROOT/"COMPACT_BICONTINUOUS.extxyz"
    a=read(atomfile);modelsha=attach(a)
    refpos=a.positions.copy()
    t0=time.perf_counter();row,f,fn,q,coord,disp,cand=metrics(a,refpos,0,"static")
    row["static_wall_seconds"]=time.perf_counter()-t0
    row["model_sha256"]=modelsha;row["torch_version"]=torch.__version__;row["mace_version"]="0.3.16";row["dtype"]="float64";row["device"]="CPU"
    write_frame(ROOT/"STATIC_RUN0.extxyz",a,row,f,fn,q,coord,disp,cand)
    (ROOT/"STATIC_RUN0_REPORT.json").write_text(json.dumps(row,indent=2,allow_nan=True)+"\n")
    with (ROOT/"STATIC_FORCES.csv").open("w",newline="") as fp:
        w=csv.writer(fp);w.writerow(["atom_id","element","phase","x_A","y_A","z_A","force_norm_eV_A","Fx_eV_A","Fy_eV_A","Fz_eV_A","Ag_q6","Ag_coordination","interface_atom","hot_source","TSC_sink"])
        sy=a.get_chemical_symbols();ph=a.arrays["phase"]
        for i in range(len(a)):w.writerow([int(a.arrays["atom_id"][i]),sy[i],ph[i],*a.positions[i],fn[i],*f[i],q[i],coord[i],int(a.arrays["interface_atom"][i]),int(a.arrays["hot_source"][i]),int(a.arrays["tsc_heat_sink"][i])])
    print("STATIC",json.dumps(row,allow_nan=True),flush=True)
    if stage in {"static"}:return
    if severe(row):
        print("STOP: severe force/contact at static run0",flush=True);return
    if stage=="all":
        # Bounded 25-step whole-cell FIRE, used only to remove initial contact
        # stress. No atoms are deleted or constrained.
        opt=FIRE(a,dt=0.05,maxstep=0.035,downhill_check=False,logfile=None)
        relax0=time.perf_counter();rrows=[]
        for k in range(26):
            rr,ff,nn,qq,cc,dd,ca=metrics(a,refpos,k,"relax")
            rrows.append(rr)
            if k%5==0 or k==25:write_frame(ROOT/"RELAXATION.extxyz",a,rr,ff,nn,qq,cc,dd,ca);print("FIRE",k,"Fmax",rr["Fmax_eV_A"],flush=True)
            if severe(rr):print("STOP: severe during FIRE",json.dumps(rr),flush=True);return
            if k==25:break
            opt.step(ff)
        with (ROOT/"RELAXATION_METRICS.csv").open("w",newline="") as fp:
            w=csv.DictWriter(fp,fieldnames=list(rrows[0]));w.writeheader();w.writerows(rrows)
        (ROOT/"RELAXATION_REPORT.json").write_text(json.dumps({"steps":25,"wall_seconds":time.perf_counter()-relax0,"initial_Fmax":row["Fmax_eV_A"],"final":rrows[-1],"note":"short, bounded all-atom FIRE; no atoms deleted/fixed"},indent=2,allow_nan=True)+"\n")
        row=rrows[-1]
        # Establish realistic thermal velocities and a 300 K baseline.
        rng=np.random.default_rng(SEED);MaxwellBoltzmannDistribution(a,temperature_K=300.0,rng=rng)
        for arr in ("STATIC_RUN0.extxyz",):pass
        baseline=ROOT/"HEAT_COOL_OVITO.extxyz"
        if baseline.exists():baseline.unlink()
        write(baseline,a,format="extxyz")
        refpos=a.positions.copy()
        t_start=time.perf_counter();dyn=Langevin(a,timestep=DT_FS*units.fs,temperature_K=300.0,friction=FRICTION_FS/units.fs,fixcm=False,rng=rng)
        brow=[]
        for step in range(101):
            r,ff,nn,qq,cc,dd,ca=metrics(a,refpos,step,"baseline_300K")
            if step%10==0 or step==100:write_frame(baseline,a,r,ff,nn,qq,cc,dd,ca);print("BASE",step,"T",round(r["T_system_K"]),"F",round(r["Fmax_eV_A"],1),flush=True)
            brow.append(r)
            if severe(r):print("STOP: severe 300K baseline",json.dumps(r),flush=True);return
            if step==100:break
            dyn.run(1)
        with (ROOT/"BASELINE_METRICS.csv").open("w",newline="") as fp:
            w=csv.DictWriter(fp,fieldnames=list(brow[0]));w.writeheader();w.writerows(brow)
        (ROOT/"BASELINE_REPORT.json").write_text(json.dumps({"steps":100,"dt_fs":DT_FS,"physical_time_ps":0.05,"wall_seconds":time.perf_counter()-t_start,"final":brow[-1],"status":"short_baseline_completed"},indent=2,allow_nan=True)+"\n")
        write(ROOT/"BASELINE_CHECKPOINT.extxyz",a,format="extxyz")
        row=brow[-1]
        # State transitions continuously from baseline. One and only one Ag
        # Langevin source; mobile TSC heat sink at the opposite free surface.
        ag,tsc,hot,sink=masks(a)
        temp=np.full((len(a),1),COLD_K);fric=np.zeros((len(a),1))
        temp[hot]=HOT_K;fric[hot]=FRICTION_FS/units.fs;fric[sink]=FRICTION_FS/units.fs
        dyn=Langevin(a,timestep=DT_FS*units.fs,temperature_K=temp,friction=fric,fixcm=False,rng=rng)
        # Continue in 0.20 ps chunks, with a four-ps / twelve-hour safety cap.
        # The Langevin state and RNG remain live across chunk boundaries.
        nheat_max=int(round(MAX_HEAT_PS*1000/DT_FS));hotrows=[];ht0=time.perf_counter();heat_path=baseline
        melt=False; melt_start_ps=None; stop_reason="safety_limit"
        for step in range(1,nheat_max+1):
            if time.perf_counter()-ht0>=MAX_HEAT_WALL_S:
                stop_reason="12h_wall_limit";break
            dyn.run(1)
            if step%10==0 or step==nheat_max:
                rr,ff,nn,qq,cc,dd,ca=metrics(a,refpos,100+step,"heating")
                write_frame(heat_path,a,rr,ff,nn,qq,cc,dd,ca);hotrows.append(rr)
                if step%50==0 or step==nheat_max:
                    print("HEAT",step,"t_ps",round(step*DT_FS/1000,4),"T_hot",round(rr["T_hot_Ag_K"]),"hot_liq",rr["hot_liquid_candidate_Ag"],"hot_cluster",rr["hot_max_liquid_cluster_Ag"],"outside_cluster",rr["outside_hot_max_cluster_Ag"],"F",round(rr["Fmax_eV_A"],1),flush=True)
                if severe(rr):stop_reason="severe_force_or_short_contact";print("STOP: severe heating condition",json.dumps(rr),flush=True);break
                if persistent_melt(hotrows,int(hot.sum())):
                    melt=True;stop_reason="persistent_localized_Ag_liquid_cluster"
                    melt_start_ps=float((hotrows[-3]["step"]-100)*DT_FS/1000);print("MELTING_ESTABLISHED at",melt_start_ps,"ps",flush=True);break
                if rr["outside_hot_max_cluster_Ag"]>=10:
                    stop_reason="remote_liquid_cluster_competes_with_single_pool";print("STOP: remote candidate Ag disorder",flush=True);break
                if step%400==0:
                    write(ROOT/f"HEATING_CHECKPOINT_{step:05d}.extxyz",a,format="extxyz")
                    (ROOT/"HEATING_CONTINUATION_STATE.json").write_text(json.dumps({"completed_heating_steps":step,"heating_time_ps":step*DT_FS/1000,"random_generator_state":rng.bit_generator.state},indent=2,default=lambda x:x.tolist() if hasattr(x,"tolist") else int(x))+"\n")
            if step%400==0:continue
        heat_done=int(hotrows[-1]["step"]-100) if hotrows else 0
        with (ROOT/"HEATING_METRICS.csv").open("w",newline="") as fp:
            w=csv.DictWriter(fp,fieldnames=list(hotrows[0]));w.writeheader();w.writerows(hotrows)
        (ROOT/"HEATING_REPORT.json").write_text(json.dumps({"dt_fs":DT_FS,"requested_steps":nheat_max,"completed_steps":int(hotrows[-1]["step"]-100) if hotrows else 0,"physical_time_ps":float((hotrows[-1]["step"]-100)*DT_FS/1000) if hotrows else 0,"wall_seconds":time.perf_counter()-ht0,"melting_established":melt,"melting_start_ps":melt_start_ps,"stop_reason":stop_reason,"final":hotrows[-1] if hotrows else None,"status":"melted" if melt else "not_melted_before_stop","hot_source_count":int(hot.sum()),"hot_source_K":HOT_K,"sink_count":int(sink.sum()),"sink_K":COLD_K,"max_heat_ps":MAX_HEAT_PS,"max_heat_wall_seconds":MAX_HEAT_WALL_S},indent=2,allow_nan=True)+"\n")
        write(ROOT/"HEATING_CHECKPOINT.extxyz",a,format="extxyz")
        # Only cool after sustained local melting evidence. Heat is switched off;
        # the cold boundary is the mobile Ti3SiC2 sink, with all Ag unthermostatted.
        coolrows=[];cool_ok=False;cool_reason="no_melt_established"
        if melt:
            t_cool=time.perf_counter();temps_arr=np.full((len(a),1),COLD_K);fric_arr=np.zeros((len(a),1));fric_arr[sink]=FRICTION_FS/units.fs
            cdyn=Langevin(a,timestep=DT_FS*units.fs,temperature_K=temps_arr,friction=fric_arr,fixcm=False,rng=rng)
            coolpath=heat_path;write(ROOT/"COOLING_START_CHECKPOINT.extxyz",a,format="extxyz")
            n_cool=4000
            for cs in range(1,n_cool+1):
                cdyn.run(1)
                if cs%10==0 or cs==n_cool:
                    rr,ff,nn,qq,cc,dd,ca=metrics(a,refpos,100+heat_done+cs,"cooling")
                    write_frame(coolpath,a,rr,ff,nn,qq,cc,dd,ca);coolrows.append(rr)
                    if cs%200==0:print("COOL",cs,"t_ps",round(cs*DT_FS/1000,4),"T_hot",round(rr["T_hot_Ag_K"]),"q6",round(rr["q6_hot_median"],3),"cluster",rr["hot_max_liquid_cluster_Ag"],flush=True)
                    if severe(rr):cool_reason="severe_force_or_short_contact";break
                    if persistent_solidification(coolrows):cool_ok=True;cool_reason="cooling_and_reordering_criteria_persisted";break
            if coolrows:
                with (ROOT/"COOLING_METRICS.csv").open("w",newline="") as fp:
                    w=csv.DictWriter(fp,fieldnames=list(coolrows[0]));w.writeheader();w.writerows(coolrows)
                (ROOT/"COOLING_REPORT.json").write_text(json.dumps({"dt_fs":DT_FS,"physical_time_ps":float((coolrows[-1]["step"]-(100+heat_done))*DT_FS/1000),"wall_seconds":time.perf_counter()-t_cool,"resolidification_established":cool_ok,"stop_reason":cool_reason,"final":coolrows[-1],"sink_count":int(sink.sum()),"sink_target_K":COLD_K},indent=2,allow_nan=True)+"\n")
                write(ROOT/"COOLING_CHECKPOINT.extxyz",a,format="extxyz")
        (ROOT/"PERFORMANCE_REPORT.json").write_text(json.dumps({"static_force_seconds":row.get("static_wall_seconds"),"total_wall_seconds":time.perf_counter()-t0,"peak_RSS_GiB":rss_gib(),"cpu_threads":torch.get_num_threads(),"device":"CPU","dt_fs":DT_FS},indent=2)+"\n")
    elif stage=="baseline":
        print("Use stage=all for a continuous baseline and local heat run.")

if __name__=="__main__":main()
