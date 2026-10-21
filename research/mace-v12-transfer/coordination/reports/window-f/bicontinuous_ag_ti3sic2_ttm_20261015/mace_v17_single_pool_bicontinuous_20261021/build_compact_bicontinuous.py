#!/usr/bin/env python3
"""Build a compact, ideal-lattice-derived 3-D bicontinuous Ag/Ti3SiC2 slab.

The model is an exploratory CPU-sized representative, not a coordinate crop or
trajectory restart. X is free with vacuum at both ends; Y/Z are periodic.
"""
from __future__ import annotations

import csv, hashlib, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase import Atoms
from ase.build import make_supercell
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write
from ase.neighborlist import neighbor_list
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT=Path(__file__).resolve().parent
REPO=next(p for p in (ROOT,*ROOT.parents) if (p/".git").exists())
CIF=REPO/"research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif"
SEED=202610211
AG_A0=4.09
MIN_CROSS=2.20
TARGET_WT=0.5995
PRE_EXCLUSION_TARGET_WT=0.6170
VAC=8.0
SIGMA_A=(3.0,3.8,3.8)
GRID=(40,64,64)

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def parent_lattices():
    cif=read(CIF)
    if cif.get_chemical_formula()!="C4Si2Ti6": raise RuntimeError("Unexpected Ti3SiC2 CIF")
    pmat=np.array([[0,0,1],[1,0,0],[1,2,0]],int)
    t=make_supercell(cif,pmat)
    rot=np.array([[0,1,0],[0,0,1],[1,0,0]],float)
    t.positions=t.positions@rot; t.cell.array=t.cell.array@rot; t.wrap()
    ax,ay,az=t.cell.lengths()
    ny_t,nz_t=8,5
    ly_t,lz_t=ny_t*ay,nz_t*az
    agp=AG_A0/np.sqrt(2)
    ny_a,nz_a=9,9
    ly_a,lz_a=ny_a*agp,nz_a*agp
    ly=np.sqrt(ly_t*ly_a); lz=np.sqrt(lz_t*lz_a)
    sty,stz=ly/ly_t,lz/lz_t; say,saz=ly/ly_a,lz/lz_a
    # One complete 17.6235-A c repeat; periodic translations in Y/Z preserve
    # exact parent crystal topology, with only the measured matching strains.
    tpos=[]; tsym=[]
    for iy in range(ny_t):
      for iz in range(nz_t):
        q=t.positions.copy(); q[:,1]=(q[:,1]+iy*ay)*sty; q[:,2]=(q[:,2]+iz*az)*stz
        tpos.append(q); tsym.extend(t.get_chemical_symbols())
    tpos=np.vstack(tpos); tsym=np.array(tsym)
    # Orthorhombic FCC Ag with [001] along X and [110]/[1-10] along Y/Z.
    ay_a=az_a=agp; basis=np.array([[0,0,0],[AG_A0/2,ay_a/2,az_a/2]])
    nx_a=int(np.ceil(ax/AG_A0));
    ix,iy,iz,ib=np.meshgrid(np.arange(nx_a),np.arange(ny_a),np.arange(nz_a),np.arange(2),indexing="ij")
    idx=np.column_stack([ix.ravel(),iy.ravel(),iz.ravel(),ib.ravel()])
    apos=basis[idx[:,3]].copy(); apos[:,0]+=idx[:,0]*AG_A0
    apos[:,1]=(apos[:,1]+idx[:,1]*ay_a)*say
    apos[:,2]=(apos[:,2]+idx[:,2]*az_a)*saz
    keep=apos[:,0]<ax-1e-9; apos=apos[keep]
    return ax,ly,lz,tpos,tsym,apos,{
      "TSC_Y_strain":sty-1,"TSC_Z_strain":stz-1,"Ag_Y_strain":say-1,"Ag_Z_strain":saz-1,
      "TSC_parent_repeat_A":[float(ax),float(ay),float(az)],"TSC_repeats_YZ":[ny_t,nz_t],
      "Ag_repeat_A":[float(AG_A0),float(agp),float(agp)],"Ag_repeats_YZ":[ny_a,nz_a]}

def field_grid(ax,ly,lz):
    rng=np.random.default_rng(SEED)
    f=gaussian_filter(rng.normal(size=GRID),
       sigma=(SIGMA_A[0]/(ax/(GRID[0]-1)),SIGMA_A[1]/(ly/GRID[1]),SIGMA_A[2]/(lz/GRID[2])),
       mode=("reflect","wrap","wrap"))
    f-=f.mean(); f/=f.std()
    gx=np.linspace(0,ax,GRID[0])[:,None,None]
    gy=np.arange(GRID[1])[None,:,None]*ly/GRID[1]
    gz=np.arange(GRID[2])[None,None,:]*lz/GRID[2]
    dy=np.minimum(np.abs(gy-ly/2),ly-np.abs(gy-ly/2))
    dz=np.minimum(np.abs(gz-lz/2),lz-np.abs(gz-lz/2))
    f += 1.15*np.exp(-0.5*(((ax-gx)/4.2)**2+(dy/6.5)**2+(dz/6.5)**2))
    f-=f.mean(); f/=f.std()
    return f

def sample(f,pos,ax,ly,lz):
    c=np.vstack((pos[:,0]/ax*(GRID[0]-1),pos[:,1]%ly/ly*GRID[1],pos[:,2]%lz/lz*GRID[2]))
    return map_coordinates(f,c,order=1,mode="wrap",prefilter=False)

def choose_threshold(a,t,ts):
    ma=107.8682; masses={"Ti":47.867,"Si":28.085,"C":12.011}
    lo,hi=min(float(a.min()),float(t.min()))-1e-8,max(float(a.max()),float(t.max()))+1e-8
    for _ in range(48):
        mid=(lo+hi)/2; nag=np.count_nonzero(a>mid)
        keep=t<=mid
        mt=sum(np.count_nonzero(keep & (ts==e)) * masses[e] for e in masses)
        wt=nag*ma/(nag*ma+mt)
        if wt>PRE_EXCLUSION_TARGET_WT:lo=mid
        else:hi=mid
    return (lo+hi)/2

def remove_cross_clashes(atoms,phase,scores,threshold):
    alive=np.ones(len(atoms),bool); depth=np.abs(scores-threshold); removed=[]
    for _ in range(20):
        ids=np.flatnonzero(alive); a=atoms[ids]; p=phase[ids]
        i,j,d=neighbor_list("ijd",a,MIN_CROSS)
        cross=(p[i]!=p[j])&(i<j)
        i,j,d=i[cross],j[cross],d[cross]
        if not len(d):break
        # Remove the more interfacial site; ties reproducibly remove Ag.
        deg=np.bincount(np.r_[i,j],minlength=len(a));
        for edge in np.argsort(d,kind="stable"):
            u,v=ids[int(i[edge])],ids[int(j[edge])]
            if not alive[u] or not alive[v]:continue
            if depth[u]<depth[v]-1e-10: victim=u
            elif depth[v]<depth[u]-1e-10: victim=v
            elif deg[i[edge]]>=deg[j[edge]]:victim=u if phase[u]=="Ag" else v
            else:victim=v if phase[v]=="Ag" else u
            alive[victim]=False; removed.append(int(victim))
    else: raise RuntimeError("interface conflict resolution did not converge")
    return atoms[alive],phase[alive],scores[alive],removed

def graph(atoms,phase_name,cutoff):
    p=np.asarray(atoms.arrays["phase"]); ids=np.flatnonzero(p==phase_name)
    i,j,s=neighbor_list("ijS",atoms,cutoff)
    m=(p[i]==phase_name)&(p[j]==phase_name)
    u,v,ss=i[m],j[m],s[m].astype(int)
    n=len(atoms); adj=[[] for _ in range(n)]
    for a,b,shift in zip(u,v,ss):adj[a].append((b,shift));adj[b].append((a,-shift))
    seen=set(); wind=np.zeros(3,bool); comp=[]
    for root in ids:
      if int(root) in seen:continue
      seen.add(int(root)); stack=[int(root)]; image={int(root):np.zeros(3,int)}; nodes=[]; w=np.zeros(3,bool)
      while stack:
        a=stack.pop();nodes.append(a)
        for b,sh in adj[a]:
          new=image[a]+sh
          if b not in seen:seen.add(b); image[b]=new;stack.append(b)
          else:w|=new!=image[b]
      comp.append((len(nodes),nodes,w));
    comp.sort(key=lambda x:x[0],reverse=True)
    largest=comp[0]
    wind=largest[2]
    xx=atoms.positions[np.asarray(largest[1]),0]
    return {"phase":phase_name,"cutoff_A":cutoff,"components":len(comp),"largest_atoms":largest[0],
      "largest_fraction":largest[0]/len(ids),"wind_Y":bool(wind[1]),"wind_Z":bool(wind[2]),
      "X_span_A":[float(xx.min()),float(xx.max())]}

def main():
    ax,ly,lz,tp,ts,ap,strain=parent_lattices(); f=field_grid(ax,ly,lz)
    at=sample(f,ap,ax,ly,lz); tt=sample(f,tp,ax,ly,lz)
    # Adjust pre-removal threshold to compensate for deleted unlike-phase clash
    # sites; a small discrete composition offset is reported, not hidden.
    t=choose_threshold(at,tt,ts)
    tkeep=tt<=t; akeep=at>t
    pos=np.vstack((tp[tkeep],ap[akeep])); sy=np.r_[ts[tkeep],np.full(akeep.sum(),"Ag")]
    ph=np.r_[np.full(tkeep.sum(),"Ti3SiC2"),np.full(akeep.sum(),"Ag")]
    score=np.r_[tt[tkeep],at[akeep]]
    pos[:,0]+=VAC
    cell=np.diag([ax+2*VAC,ly,lz])
    atoms=Atoms(sy,positions=pos,cell=cell,pbc=[False,True,True])
    atoms,ph,score,removed=remove_cross_clashes(atoms,ph,score,t)
    sy=np.asarray(atoms.get_chemical_symbols())
    ids=np.arange(1,len(atoms)+1,dtype=int)
    atoms.arrays["atom_id"]=ids; atoms.arrays["phase"]=ph.astype("U8"); atoms.arrays["field_value"]=score
    # one localized surface heat target: Ag sites in one radial patch at X+.
    xhot=VAC+ax-1.5; y0,z0=ly/2,lz/2
    dy=np.minimum(np.abs(atoms.positions[:,1]-y0),ly-np.abs(atoms.positions[:,1]-y0))
    dz=np.minimum(np.abs(atoms.positions[:,2]-z0),lz-np.abs(atoms.positions[:,2]-z0))
    r=np.sqrt(dy**2+dz**2+(atoms.positions[:,0]-xhot)**2)
    hot=(ph=="Ag")&(r<=7.0)
    sink=(ph=="Ti3SiC2")&(atoms.positions[:,0] <= VAC+5.5)
    interface=np.zeros(len(atoms),bool)
    ii,jj,dd=neighbor_list("ijd",atoms,4.5)
    cross=(ph[ii]!=ph[jj])
    interface[np.unique(np.r_[ii[cross],jj[cross]])]=True
    atoms.arrays["hot_source"]=hot.astype(np.int8);atoms.arrays["tsc_heat_sink"]=sink.astype(np.int8)
    atoms.arrays["interface_atom"]=interface.astype(np.int8)
    atoms.arrays["region"]=np.where(hot,"single_hot_Ag_patch",np.where(sink,"TSC_heat_sink",np.where(interface,"interface",np.where(ph=="Ag","Ag_network","TSC_backbone")))).astype("U24")
    atoms.info.update({"model_status":"COMPACT_3D_BICONTINUOUS_EXPLORATORY_GEOMETRY",
      "construction":"independently generated ideal Ti3SiC2 and FCC Ag lattices classified by a seeded 3D periodic smooth random field; no historical trajectory coordinates",
      "seed":SEED,"phase_field_threshold":float(t),"pre_exclusion_target_Ag_wt_fraction":PRE_EXCLUSION_TARGET_WT,"pbc":"F T T","single_heat_source":"Ag atoms within 7 A spherical patch centered near X+ free surface",
      "heat_sink":"mobile Ti3SiC2 in bottom 5.5 A slab; no frozen backbone"})
    write(ROOT/"COMPACT_BICONTINUOUS.extxyz",atoms,format="extxyz")
    write(ROOT/"COMPACT_BICONTINUOUS.data",atoms,format="lammps-data",atom_style="atomic",specorder=["Ti","Si","C","Ag"])
    with (ROOT/"ATOM_GROUPS.csv").open("w",newline="") as fcsv:
      w=csv.writer(fcsv);w.writerow(["atom_id","element","phase","region","hot_source","tsc_heat_sink","interface_atom"])
      for k in range(len(atoms)):w.writerow([int(ids[k]),sy[k],ph[k],atoms.arrays["region"][k],int(hot[k]),int(sink[k]),int(interface[k])])
    np.savetxt(ROOT/"INTERFACE_SITE_REMOVALS.csv",np.array(removed,dtype=int),fmt="%d",header="original_combined_0based_index_removed_to_resolve_unlike_contacts_below_2.20A")
    masses={e:float(atomic_masses[atomic_numbers[e]]) for e in ["Ag","Ti","Si","C"]}
    counts={e:int(np.count_nonzero(sy==e)) for e in ["Ag","Ti","Si","C"]}
    mag=counts["Ag"]*masses["Ag"];mtsc=sum(counts[e]*masses[e] for e in ["Ti","Si","C"])
    crossd=dd[cross]
    # Per-phase nearest-neighbor graph checks with modest cutoff sensitivity.
    graphs=[graph(atoms,"Ag",3.15),graph(atoms,"Ti3SiC2",2.75)]
    audits={"atom_count":len(atoms),"counts":counts,"Ag_mass_fraction":mag/(mag+mtsc),"TSC_mass_fraction":mtsc/(mag+mtsc),
      "box_A":cell.diagonal().tolist(),"material_X_A":[VAC,VAC+ax],"vacuum_each_X_surface_A":VAC,
      "pbc":[False,True,True],"lattice_matching":strain,"removed_interface_sites":len(removed),
      "minimum_unlike_phase_contact_A":float(crossd.min()) if len(crossd) else None,
      "unlike_pair_percentiles_A":np.percentile(crossd,[1,5,50,95]).tolist() if len(crossd) else [],
      "interface_atom_count":int(interface.sum()),"hot_source_Ag_count":int(hot.sum()),"hot_source_Ag_mass_fraction":float(hot.sum()/np.sum(ph=="Ag")),
      "hot_source_position_A":[float(xhot),float(y0),float(z0)],"hot_source_radius_A":7.0,"TSC_sink_count":int(sink.sum()),
      "phase_connectivity":graphs,"ideal_parent_CIF_sha256":sha(CIF),"field_sigma_A":SIGMA_A,
      "scientific_scope":"thin compact 3-D representative slab; not the full 39,486 atom double-continuous model; geometric connectivity is not force stability or melting evidence"}
    (ROOT/"GEOMETRY_AUDIT.json").write_text(json.dumps(audits,indent=2)+"\n")
    fig,axp=plt.subplots(figsize=(6.8,5.8),dpi=170)
    palette={"Ag":"#e5a129","Ti3SiC2":"#d83dc3"}
    view=np.array([1.7,-1.0,0.8]);view/=np.linalg.norm(view);up0=np.array([0,0,1.]);right=np.cross(up0,view);right/=np.linalg.norm(right);up=np.cross(view,right)
    xy=np.c_[atoms.positions@right,atoms.positions@up];dep=atoms.positions@view;order=np.argsort(dep)
    axp.scatter(xy[order,0],xy[order,1],c=[palette[x] for x in ph[order]],s=4,linewidth=0,alpha=.9)
    axp.set_aspect('equal');axp.axis('off');axp.set_title('Compact 3-D bicontinuous Ag (gold) / Ti₃SiC₂ (magenta)\nred outline = single local heating patch')
    hh=np.flatnonzero(hot)
    if len(hh):axp.scatter(xy[hh,0],xy[hh,1],facecolors='none',edgecolors='red',s=20,linewidth=.4)
    fig.tight_layout();fig.savefig(ROOT/"COMPACT_BICONTINUOUS_VIEW.png",bbox_inches='tight');plt.close(fig)
    print(json.dumps(audits,indent=2))

if __name__=="__main__":main()
