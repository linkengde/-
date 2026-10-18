from pathlib import Path
import csv, json
import numpy as np
from ase import Atoms
from ase.io import write
from ase.neighborlist import neighbor_list
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent
SYMBOL={1:"Ti",2:"Si",3:"C",4:"Ag"}
Z={"Ti":22,"Si":14,"C":6,"Ag":47}

def read_data(path):
    lines=path.read_text().splitlines(); n=None; b={}; start=None
    for i,line in enumerate(lines):
        f=line.split()
        if len(f)==2 and f[1]=="atoms": n=int(f[0])
        if len(f)>=4 and f[-2:] in (["xlo","xhi"],["ylo","yhi"],["zlo","zhi"]): b[f[-2][0]]=(float(f[0]),float(f[1]))
        if line.strip().startswith("Atoms"): start=i+1; break
    while not lines[start].strip(): start+=1
    rows=[x.split() for x in lines[start:start+n]]
    ids=np.array([int(x[0]) for x in rows]); typ=np.array([int(x[1]) for x in rows]); pos=np.array([[float(a) for a in x[2:5]] for x in rows])
    order=np.argsort(ids); ids=ids[order]; typ=typ[order]; pos=pos[order]
    cell=np.array([b[a][1]-b[a][0] for a in "xyz"])
    return ids,typ,pos,cell

def read_dump(path):
    lines=path.read_text().splitlines(); out=[]; i=0
    while i<len(lines):
        if lines[i]!="ITEM: TIMESTEP": i+=1; continue
        step=int(lines[i+1]); n=int(lines[i+3]); cols=lines[i+8].split()[2:]
        a=np.loadtxt(lines[i+9:i+9+n],ndmin=2); a=a[np.argsort(a[:,cols.index("id")])]
        out.append((step,a,cols)); i+=9+n
    return out

def stats(x):
    return {"Fmax":float(x.max()),"P95":float(np.percentile(x,95)),"P99":float(np.percentile(x,99)),"mean":float(x.mean()),"n_gt_2":int((x>2).sum()),"n_gt_5":int((x>5).sum()),"n_gt_10":int((x>10).sum())}

ids,typ,relaxed_pos,cell=read_data(ROOT/"post_relax.data")
frames=read_dump(ROOT/"trajectory_300K.dump")
assert len(frames)>=2 and [x[0] for x in frames]==[0,100,200]
assert all(len(x[1])==len(ids) for x in frames)
assert all(np.array_equal(x[1][:,x[2].index("id")].astype(int),ids) for x in frames)
assert all(np.array_equal(x[1][:,x[2].index("type")].astype(int),typ) for x in frames)

# Current nearest opposite-phase distance using a 3x3 Y/Z image expansion.
phase_groups={}
with (ROOT/"input/PHASE_ATOM_GROUPS.csv").open() as f:
    phase_groups={int(r["atom_id"]):r["geometry_region"] for r in csv.DictReader(f)}

results=[]; high_rows=[]; pair_rows=[]; ext_frames=[]
for step,a,cols in frames:
    ix={k:i for i,k in enumerate(cols)}
    pos=a[:,[ix["x"],ix["y"],ix["z"]]]
    forces=a[:,[ix["fx"],ix["fy"],ix["fz"]]]
    fn=np.linalg.norm(forces,axis=1)
    atoms=Atoms(numbers=[Z[SYMBOL[int(t)]] for t in typ],positions=pos,cell=np.diag(cell),pbc=[False,True,True])
    ii,jj,dd=neighbor_list("ijd",atoms,3.4,self_interaction=False)
    keep=ii<jj; ii,jj,dd=ii[keep],jj[keep],dd[keep]
    pair_values={}
    pairs=[(1,3,"Ti-C"),(1,2,"Ti-Si"),(4,1,"Ag-Ti"),(4,2,"Ag-Si"),(4,3,"Ag-C"),(4,4,"Ag-Ag")]
    for ta,tb,name in pairs:
        m=((typ[ii]==ta)&(typ[jj]==tb))|((typ[ii]==tb)&(typ[jj]==ta)); d=dd[m]
        pair_values[name]={"minimum_A":float(d.min()) if len(d) else None,"n_lt_1p8":int((d<1.8).sum()),"n_lt_1p9":int((d<1.9).sum()),"n_lt_2p0":int((d<2.0).sum()),"n_lt_2p1":int((d<2.1).sum()),"n_lt_2p2":int((d<2.2).sum())}
        pair_rows.append([step,name,pair_values[name]["minimum_A"],*(pair_values[name][k] for k in ("n_lt_1p8","n_lt_1p9","n_lt_2p0","n_lt_2p1","n_lt_2p2"))])
    conn={}
    for name,phase,cut in [("TSC",{1,2,3},2.70),("Ag",{4},2.90),("Ag",{4},3.05)]:
        glob=np.flatnonzero(np.isin(typ,list(phase))); loc=np.full(len(typ),-1,int); loc[glob]=np.arange(len(glob))
        m=(loc[ii]>=0)&(loc[jj]>=0)&(dd<=cut)
        mat=sparse.coo_matrix((np.ones(int(m.sum()),dtype=np.int8),(loc[ii[m]],loc[jj[m]])),shape=(len(glob),len(glob))); mat=mat+mat.T
        nc,labs=connected_components(mat,directed=False,return_labels=True); c=np.bincount(labs)
        conn[f"{name}_{cut:.2f}"]={"components":int(nc),"largest_atoms":int(c.max()),"largest_fraction":float(c.max()/len(glob))}
    # Phase-resolved forces and local IDs over 5 eV/A.
    by_element={}
    for t,s in SYMBOL.items():
        m=typ==t; by_element[s]={**stats(fn[m]),"n_atoms":int(m.sum())}
    ag=pos[typ==4]
    ag_img=np.concatenate([ag+np.array([0,sy*cell[1],sz*cell[2]]) for sy in (-1,0,1) for sz in (-1,0,1)])
    nearest=np.full(len(ids),np.inf); tsc=np.flatnonzero(typ!=4)
    nearest[tsc]=cKDTree(ag_img).query(pos[tsc],k=1)[0]
    top=np.flatnonzero(fn>5.0)
    for k in top:
        d=pos-pos[k]; d[:,1:]-=np.rint(d[:,1:]/cell[1:])*cell[1:]; rr=np.linalg.norm(d,axis=1); rr[k]=np.inf; js=np.argsort(rr)[:8]
        for j in js:
            high_rows.append([step,int(ids[k]),SYMBOL[int(typ[k])],float(fn[k]),*pos[k],float(nearest[k]),phase_groups.get(int(ids[k]),""),int(ids[j]),SYMBOL[int(typ[j])],float(rr[j])])
    results.append({"step":step,"time_ps":step*0.001,"finite_forces":bool(np.isfinite(fn).all()),"force":stats(fn),"force_by_element":by_element,"pair_distances":pair_values,"connectivity":conn,"n_atoms":len(ids)})
    vel=a[:,[ix["vx"],ix["vy"],ix["vz"]]]
    ext=Atoms(numbers=[Z[SYMBOL[int(t)]] for t in typ],positions=pos,cell=np.diag(cell),pbc=[False,True,True])
    ext.arrays["id"]=ids.copy(); ext.arrays["forces"]=forces.copy(); ext.arrays["velocities"]=vel.copy(); ext.info["step"]=step; ext.info["time_ps"]=step*0.001
    ext_frames.append(ext)

# Displacement at 200 fs relative to the post-minimization state.
last=frames[-1]; ix={k:i for i,k in enumerate(last[2])}; p200=last[1][:,[ix["x"],ix["y"],ix["z"]]]
delta=p200-relaxed_pos; delta[:,1:]-=np.rint(delta[:,1:]/cell[1:])*cell[1:]; dn=np.linalg.norm(delta,axis=1)
displacement={}
for t,s in SYMBOL.items():
    v=dn[typ==t]; displacement[s]={"P50_A":float(np.percentile(v,50)),"P95_A":float(np.percentile(v,95)),"P99_A":float(np.percentile(v,99)),"max_A":float(v.max()),"n_gt_1A":int((v>1).sum())}

with (ROOT/"300K_FRAME_AUDIT.csv").open("w",newline="") as f:
    w=csv.writer(f,lineterminator="\n"); w.writerow(["step","time_ps","pair","minimum_A","n_lt_1p8","n_lt_1p9","n_lt_2p0","n_lt_2p1","n_lt_2p2"]); w.writerows(pair_rows)
with (ROOT/"300K_HIGH_FORCE_ATOMS.csv").open("w",newline="") as f:
    w=csv.writer(f,lineterminator="\n"); w.writerow(["step","atom_id","element","force_eV_A","x_A","y_A","z_A","nearest_opposite_phase_A","initial_geometry_region","neighbor_id","neighbor_element","distance_A"]); w.writerows(high_rows)
write(ROOT/"TRAJECTORY_300K_OVITO.extxyz",ext_frames,format="extxyz")
min_frames=[]
for step,a,cols in read_dump(ROOT/"minimization_frames.dump"):
    ix={k:i for i,k in enumerate(cols)}
    pos=a[:,[ix["x"],ix["y"],ix["z"]]]; forces=a[:,[ix["fx"],ix["fy"],ix["fz"]]]
    ext=Atoms(numbers=[Z[SYMBOL[int(t)]] for t in typ],positions=pos,cell=np.diag(cell),pbc=[False,True,True])
    ext.arrays["id"]=ids.copy(); ext.arrays["forces"]=forces.copy()
    ext.info["stage"]="FIRE_minimization"; ext.info["step"]=step
    min_frames.append(ext)
for ext in ext_frames: ext.info["stage"]="NVT_300K_interrupted"
write(ROOT/"MINIMIZATION_OVITO.extxyz",min_frames,format="extxyz")

result={"scope":"Interrupted 300 K NVT screening only; no heating/cooling/TTM/DFT/training","source_data_sha256":"e2b62a68c46e23037370a6ca20951aa0774edea2d728d3c3aadbcba4b3755674","potential":"calibration-v0, unchanged copies","frames":results,"displacement_from_relaxed_at_0p2ps":displacement,"stop_reason":"Repeated localized C/Ti interface force outliers >60 eV/A appeared during the first 0.2 ps.","source_model_unchanged":True}
(ROOT/"SCREEN_300K_AUDIT.json").write_text(json.dumps(result,indent=2)+"\n")
with (ROOT/"SCREEN_300K_AUDIT.md").open("w") as f:
    f.write("# 300 K baseline screen: stopped after 0.2 ps\n\n")
    f.write("The test was run on `post_relax.data`, an independent copy of the user-supplied 39,486-atom model. The original model and potential files remain unchanged. Calibration-v0 mapping was preserved. The run used NVT at 300 K, 1 fs timestep, and was interrupted after the 200 fs frame; no heat pulse or TTM was applied.\n\n")
    f.write("## Observed force spike\n\n| time | temperature from LAMMPS | Fmax (eV/Å) | mean force | atoms F>5 | atoms F>10 |\n|---:|---:|---:|---:|---:|---:|\n")
    lmp={0:(300.0,-227196.0481),100:(198.9974283,-226598.9318),200:(263.0009004,-226778.0677)}
    for x in results:
        step=x["step"]; force=x["force"]
        f.write(f"| {step/1000:.3f} ps | {lmp[step][0]:.1f} K | {force['Fmax']:.3f} | {force['mean']:.3f} | {force['n_gt_5']} | {force['n_gt_10']} |\n")
    f.write("\nAt 0.1 ps, the maximum force was 66.49 eV/Å on C ID 19219; at 0.2 ps the LAMMPS log reports Fmax=67.01 eV/Å. The 0.1 ps neighborhood contains Ti ID 25627 at 1.809 Å, Ti ID 25628 at 1.914 Å, and Ag ID 38116 at 1.981 Å from that carbon. The same C ID had already reached 50.83 eV/Å near minimization iteration 475. Its source-region label is a geometric `interior` candidate but it is also an initial Ag-interface candidate (nearest Ag 2.851 Å); it is not evidence of a bulk-only force-field test. This recurring local response is the immediate stop condition.\n\n")
    f.write("## Skeleton and contact summary\n\n")
    for x in results:
        f.write(f"- step {x['step']}: TSC graph at 2.70 Å has {x['connectivity']['TSC_2.70']['components']} component(s); Ag graph at 3.05 Å has {x['connectivity']['Ag_3.05']['components']} component(s), largest fraction {x['connectivity']['Ag_3.05']['largest_fraction']:.5f}.\n")
    f.write("\nThe test was intentionally stopped before equilibration because the repeated localized force spikes are inconsistent with proceeding to a driven melt-pool calculation under this calibration-v0 mapping. The partial trajectory frames (0, 100, 200 fs) are provided for OVITO inspection; they are not a completed 300 K stability test. No heating, cooling, TTM, DFT, or model/potential modification occurred.\n")
local=ROOT/"RELAXATION_LOCAL_GEOMETRY_AUDIT.md"
text=local.read_text()
text=text.replace("No 300 K, heating/cooling, TTM, DFT, or training run has been started.","A separate 300 K NVT screen was started on the relaxed copy and interrupted after 0.2 ps because a repeated localized interface force spike appeared; see `SCREEN_300K_AUDIT.md`. No heating/cooling, TTM, DFT, or training run was started.")
local.write_text(text)
print(json.dumps({"steps":[x["step"] for x in results],"forces":[x["force"] for x in results],"disp":displacement,"high_force_neighbor_rows":len(high_rows)},indent=2))
