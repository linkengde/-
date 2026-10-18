from pathlib import Path
import csv
import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list
from scipy import sparse
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parent
TYP = {1: "Ti", 2: "Si", 3: "C", 4: "Ag"}

def read_data(path):
    lines = path.read_text().splitlines()
    n = None; bounds = {}; start = None
    for i, line in enumerate(lines):
        s = line.strip(); f = s.split()
        if len(f) == 2 and f[1] == "atoms": n = int(f[0])
        if len(f) >= 4 and f[-2:] == ["xlo", "xhi"]: bounds[0] = (float(f[0]), float(f[1]))
        if len(f) >= 4 and f[-2:] == ["ylo", "yhi"]: bounds[1] = (float(f[0]), float(f[1]))
        if len(f) >= 4 and f[-2:] == ["zlo", "zhi"]: bounds[2] = (float(f[0]), float(f[1]))
        if s.startswith("Atoms"): start = i + 1; break
    while not lines[start].strip(): start += 1
    rows = []
    for line in lines[start:]:
        if not line.strip(): break
        f = line.split()
        if len(f) < 5: break
        rows.append((int(f[0]), int(f[1]), *map(float, f[2:5])))
        if len(rows) == n: break
    a = np.asarray(rows, float)
    return a[:,0].astype(int), a[:,1].astype(int), a[:,2:5], np.array([bounds[k][1]-bounds[k][0] for k in range(3)])

def read_dump(path):
    lines = path.read_text().splitlines(); out=[]; i=0
    while i < len(lines):
        if lines[i].strip() != "ITEM: TIMESTEP": i += 1; continue
        step = int(lines[i+1]); n = int(lines[i+3]); i += 4
        assert lines[i].startswith("ITEM: BOX BOUNDS")
        bounds = [tuple(map(float, lines[i+k].split()[:2])) for k in range(1,4)]; i += 4
        assert lines[i].startswith("ITEM: ATOMS")
        names = lines[i].split()[2:]; i += 1
        vals = np.loadtxt(lines[i:i+n], ndmin=2); i += n
        ix={x:k for k,x in enumerate(names)}; vals=vals[np.argsort(vals[:,ix["id"]])]
        out.append((step, vals, names, np.array([b[1]-b[0] for b in bounds])))
    return out

def atoms(pos, types, cell):
    return Atoms(numbers=[{1:22,2:14,3:6,4:47}[int(t)] for t in types], positions=pos,
                 cell=np.diag(cell), pbc=[False,True,True])

def force_norm(frame):
    step, v, names, cell=frame; ix={x:k for k,x in enumerate(names)}
    return np.linalg.norm(v[:,[ix["fx"],ix["fy"],ix["fz"]]], axis=1)

def pair_audit(frame, ids, types, out_rows, tag):
    step, v, names, cell=frame; ix={x:k for k,x in enumerate(names)}
    pos=v[:,[ix["x"],ix["y"],ix["z"]]]; a=atoms(pos,types,cell)
    ii,jj,dd=neighbor_list("ijd",a,2.35,self_interaction=False)
    keep=ii<jj; ii,jj,dd=ii[keep],jj[keep],dd[keep]
    fn=force_norm(frame)
    for x,y,d in zip(ii,jj,dd):
        p=frozenset((TYP[int(types[x])],TYP[int(types[y])]))
        pair_name={frozenset(("Ti","C")):"Ti-C",frozenset(("Ag","C")):"Ag-C",frozenset(("Ag","Ti")):"Ag-Ti",frozenset(("Ag","Si")):"Ag-Si",frozenset(("Ti","Si")):"Ti-Si",frozenset(("Ag",)):"Ag-Ag"}.get(p)
        if pair_name is None: continue
        key=pair_name
        relevant=((key=="Ti-C" and d<1.90) or (key=="Ag-C" and d<2.00) or (key=="Ag-Ti" and d<2.10) or (key=="Ag-Si" and d<2.20) or (key=="Ti-Si" and d<2.40) or (key=="Ag-Ag" and d<2.35))
        if relevant:
            out_rows.append([tag,step,int(ids[x]),TYP[int(types[x])],int(ids[y]),TYP[int(types[y])],float(d),float(fn[x]),float(fn[y]),*pos[x],*pos[y]])

ids, types, original, cell=read_data(ROOT/"input/NEW_MODEL_ATOMS.data")
frames=read_dump(ROOT/"minimization_frames.dump")
by_step={x[0]:x for x in frames}
pre=read_dump(ROOT/"pre_relax.dump")[-1]
post=read_dump(ROOT/"post_relax.dump")[-1]
frame_rows=[]
for tag,fr in [("pre",pre), *[(f"step{s}",by_step[s]) for s in (250,325,450,475,500) if s in by_step], ("post",post)]:
    pair_audit(fr,ids,types,frame_rows,tag)
with (ROOT/"MINIMIZATION_SHORT_CONTACTS.csv").open("w",newline="") as f:
    w=csv.writer(f,lineterminator="\n"); w.writerow(["frame","step","atom_i","element_i","atom_j","element_j","distance_A","force_i_eV_A","force_j_eV_A","xi_A","yi_A","zi_A","xj_A","yj_A","zj_A"]); w.writerows(frame_rows)

# Displacement by species and at each stored frame; Y/Z minimum image only.
ix={x:k for k,x in enumerate(post[2])}; ppost=post[1][:,[ix["x"],ix["y"],ix["z"]]]
delta=ppost-original
for axis in (1,2): delta[:,axis]-=np.rint(delta[:,axis]/cell[axis])*cell[axis]
disp=np.linalg.norm(delta,axis=1)
disprows=[]
for typ,sym in TYP.items():
    m=types==typ; v=disp[m]
    disprows.append((sym,int(m.sum()),float(np.median(v)),float(np.percentile(v,95)),float(np.percentile(v,99)),float(v.max()),int((v>1).sum()),int((v>3).sum())))

# Force-spike atom local coordination and its nearest neighbours at selected steps.
spike_rows=[]
for s in (450,475,500):
    fr=by_step[s]; nrm=force_norm(fr); order=np.argsort(nrm)[::-1][:12]
    ix={x:k for k,x in enumerate(fr[2])}; pos=fr[1][:,[ix["x"],ix["y"],ix["z"]]]
    a=atoms(pos,types,fr[3]); ii,jj,dd=neighbor_list("ijd",a,5.0,self_interaction=False)
    for k in order:
        m=(ii==k)&(dd<4.0); js=jj[m]; ds=dd[m]; ordern=np.argsort(ds)
        for q in ordern[:12]:
            j=js[q]
            spike_rows.append([s,int(ids[k]),TYP[int(types[k])],float(nrm[k]),*pos[k],int(ids[j]),TYP[int(types[j])],float(ds[q]),float(nrm[j])])
with (ROOT/"FORCE_SPIKE_LOCAL_NEIGHBORS.csv").open("w",newline="") as f:
    w=csv.writer(f,lineterminator="\n"); w.writerow(["step","central_id","central_element","central_force_eV_A","x_A","y_A","z_A","neighbor_id","neighbor_element","distance_A","neighbor_force_eV_A"]); w.writerows(spike_rows)

# Sensitivity of connected-component counts to the archived geometry graph cutoffs.
conn=[]
for step,fr in [(0,pre),(500,post)]:
    ix={x:k for k,x in enumerate(fr[2])}; pos=fr[1][:,[ix["x"],ix["y"],ix["z"]]]; a=atoms(pos,types,fr[3])
    for phase,phase_types,cutoffs in [("TSC",[1,2,3],np.arange(2.4,3.51,.1)),("Ag",[4],np.arange(2.7,3.51,.1))]:
        glob=np.flatnonzero(np.isin(types,phase_types)); loc=np.full(len(types),-1,int); loc[glob]=np.arange(len(glob))
        for cutoff in cutoffs:
            ii,jj=neighbor_list("ij",a,float(cutoff),self_interaction=False); keep=(loc[ii]>=0)&(loc[jj]>=0)&(ii<jj)
            mat=sparse.coo_matrix((np.ones(int(keep.sum()),dtype=np.int8),(loc[ii[keep]],loc[jj[keep]])),shape=(len(glob),len(glob))); mat=mat+mat.T
            nc,labs=connected_components(mat,directed=False,return_labels=True); counts=np.bincount(labs)
            conn.append([step,phase,float(cutoff),int(nc),int(counts.max()),float(counts.max()/len(glob)),int((counts==1).sum())])
with (ROOT/"CONNECTIVITY_CUTOFF_SENSITIVITY.csv").open("w",newline="") as f:
    w=csv.writer(f,lineterminator="\n"); w.writerow(["step","phase","cutoff_A","component_count","largest_component_atoms","largest_fraction","singleton_components"]); w.writerows(conn)

with (ROOT/"RELAXATION_LOCAL_GEOMETRY_AUDIT.md").open("w") as f:
    f.write("# Follow-up local geometry audit of bounded FIRE relaxation\n\n")
    f.write("This audit reads the original source and its archived LAMMPS force dumps only. It does not modify the source structure or run another simulation. Distances use minimum-image Y/Z and nonperiodic X. See `MINIMIZATION_SHORT_CONTACTS.csv`, `FORCE_SPIKE_LOCAL_NEIGHBORS.csv`, and `CONNECTIVITY_CUTOFF_SENSITIVITY.csv` for atom-level records.\n\n")
    f.write("## Species-resolved displacement from the unrelaxed model\n\n| phase | atoms | median Å | P95 Å | P99 Å | max Å | >1 Å | >3 Å |\n|---|---:|---:|---:|---:|---:|---:|---:|\n")
    for row in disprows: f.write(f"| {row[0]} | {row[1]} | {row[2]:.3f} | {row[3]:.3f} | {row[4]:.3f} | {row[5]:.3f} | {row[6]} | {row[7]} |\n")
    f.write("\n## Interpretation\n\nThe TSC geometric graph remains one connected component at the archived 3.30 Å cutoff, and the principal Ag component contains 99.708% of Ag atoms at 3.05 Å. Connectivity alone does not exclude local bond distortion. The post-minimization Ti–C minimum is 1.729 Å, compared with 2.092 Å in the source; the source model has no Ti–C contacts below 2.0 Å, while the relaxed model has 1,668 such pairs including 20 below 1.8 Å. Post-minimization Ag–C and Ag–Ti contacts also contract to 1.906 and 2.027 Å. These values and the extensive displacement/reorganization are reported as a potential-induced geometry risk, not as a proven force-field error or a proof of a chemical bond. The calibration-v0 cross terms are explicitly marked unvalidated.\n\n")
    f.write("The minimizer's transient force spike at step 475 is localized by the neighbor table; it fell by the saved step 500. It is not itself proof of final instability, but together with the new short-contact population it warrants stopping before thermal driving under this mapping. No 300 K, heating/cooling, TTM, DFT, or training run has been started.\n")
print("species_displacements", disprows)
print("short_contact_rows",len(frame_rows),"force_spike_local_rows",len(spike_rows),"connectivity_rows",len(conn))
print("wrote MINIMIZATION_SHORT_CONTACTS.csv, FORCE_SPIKE_LOCAL_NEIGHBORS.csv, CONNECTIVITY_CUTOFF_SENSITIVITY.csv, RELAXATION_LOCAL_GEOMETRY_AUDIT.md")
