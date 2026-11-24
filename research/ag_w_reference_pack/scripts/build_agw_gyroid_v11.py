#!/usr/bin/env python3
"""Build the V11 ideal bicontinuous Ag-W gyroid candidate and TTM geometry grids."""
from __future__ import annotations
import argparse, gzip, hashlib, json
from pathlib import Path
import numpy as np
from scipy.ndimage import label
from scipy.spatial import cKDTree

L, BOX_Z = 200.0, 300.0
N_AG, N_W = 49, 63
A_AG, A_W = 4.09, 3.16
M_AG, M_W = 107.8682, 183.84
V_AG, V_W = A_AG**3/4.0, A_W**3/2.0
I0, R, SKIN, SURFACE = 5.0, 30.0, 20.0, 200.0
CENTER = (100.0, 100.0)

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()

def lattice(n,a,basis):
    ijk=np.indices((n,n,n),dtype=np.int16).reshape(3,-1).T.astype(float)
    return np.concatenate([(ijk+np.asarray(b))*a for b in basis])

def gyroid(xyz):
    p=2*np.pi*np.asarray(xyz)/L
    x,y,z=p.T
    return np.sin(x)*np.cos(y)+np.sin(y)*np.cos(z)+np.sin(z)*np.cos(x)

def thresholds(ga,gw,gap):
    lo,hi=-0.5,1.0
    for _ in range(60):
        mid=(lo+hi)/2; ta=mid-gap/2; tw=mid+gap/2
        na=np.count_nonzero(ga<ta); nw=np.count_nonzero(gw>tw)
        ratio=na*M_AG/(nw*M_W)
        if ratio>1: hi=mid
        else: lo=mid
    cand=[]
    for mid in (lo,hi,(lo+hi)/2):
        ta,tw=mid-gap/2,mid+gap/2
        na=int(np.count_nonzero(ga<ta)); nw=int(np.count_nonzero(gw>tw))
        ratio=na*M_AG/(nw*M_W)
        cand.append((abs(ratio-1),mid,ta,tw,na,nw,ratio))
    return min(cand)[1:]

def cross_min(ag,w):
    shifts=[[i*L,j*L,0.] for i in (-1,0,1) for j in (-1,0,1)]
    tw=cKDTree(np.concatenate([w+s for s in shifts]))
    ta=cKDTree(np.concatenate([ag+s for s in shifts]))
    da,_=tw.query(ag,k=1,workers=-1); dw,_=ta.query(w,k=1,workers=-1)
    return {"Ag_to_W_min_A":float(da.min()),"Ag_to_W_p01_A":float(np.quantile(da,.01)),
            "W_to_Ag_min_A":float(dw.min()),"W_to_Ag_p01_A":float(np.quantile(dw,.01))}

def connectivity(ta,tw,n=100):
    x=(np.arange(n)+.5)*L/n
    X,Y,Z=np.meshgrid(x,x,x,indexing="ij")
    p=2*np.pi/L
    g=np.sin(p*X)*np.cos(p*Y)+np.sin(p*Y)*np.cos(p*Z)+np.sin(p*Z)*np.cos(p*X)
    out={}
    for name,mask in (("Ag",g<ta),("W",g>tw)):
        lab,_=label(np.tile(mask,(3,3,1)),structure=np.ones((3,3,3),bool))
        c=lab[n:2*n,n:2*n,:]
        ids,cnt=np.unique(c[c>0],return_counts=True)
        best=int(ids[np.argmax(cnt)]); xyz=np.argwhere(lab==best)
        out[name]={"largest_component_share":float(cnt.max()/mask.sum()),
                   "connected":bool(cnt.max()==mask.sum()),
                   "spans_z":bool(xyz[:,2].min()==0 and xyz[:,2].max()==n-1),
                   "voxel_fraction":float(mask.mean())}
    return out

def fractions(ta,tw,h,z0=40.,z1=200.,ns=5,source=False):
    nx=round(L/h); nz=round((z1-z0)/h)
    cx=(np.arange(nx)+.5)*h; cy=cx.copy(); cz=z0+(np.arange(nz)+.5)*h
    sub=((np.arange(ns)+.5)/ns-.5)*h
    a=np.zeros((nx,nx,nz)); w=np.zeros_like(a); v=np.zeros_like(a)
    q=np.zeros_like(a) if source else None
    dv=h**3/ns**3; p=2*np.pi/L
    for dx in sub:
      for dy in sub:
       for dz in sub:
        X,Y,Z=np.meshgrid(cx+dx,cy+dy,cz+dz,indexing="ij")
        g=np.sin(p*X)*np.cos(p*Y)+np.sin(p*Y)*np.cos(p*Z)+np.sin(p*Z)*np.cos(p*X)
        ma=g<ta; mw=g>tw
        a+=ma; w+=mw; v+=~(ma|mw)
        if source:
            r2=(X-CENTER[0])**2+(Y-CENTER[1])**2
            q+=(I0/SKIN)*np.exp(-r2/R**2-np.maximum(0,SURFACE-Z)/SKIN)*ma*dv
    a/=ns**3; w/=ns**3; v/=ns**3
    return a,w,v,q

def write_data(path,ag,w):
    with gzip.open(path,"wt",encoding="ascii",compresslevel=6) as f:
        n=len(ag)+len(w)
        f.write("Ag-W V11 ideal gyroid; type 1=W BCC, type 2=Ag FCC\n\n")
        f.write(f"{n} atoms\n2 atom types\n\n0 {L} xlo xhi\n0 {L} ylo yhi\n0 {BOX_Z} zlo zhi\n\n")
        f.write("Masses\n\n1 183.84\n2 107.8682\n\nAtoms # atomic\n\n")
        i=1
        for xyz in w:
            f.write(f"{i} 1 {xyz[0]:.8f} {xyz[1]:.8f} {xyz[2]:.8f}\n"); i+=1
        for xyz in ag:
            f.write(f"{i} 2 {xyz[0]:.8f} {xyz[1]:.8f} {xyz[2]:.8f}\n"); i+=1

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--gap-level",type=float,default=.12)
    p.add_argument("--fraction-samples",type=int,default=5)
    args=p.parse_args()
    if args.gap_level<=0 or args.fraction_samples<2: raise SystemExit("invalid gap or sampling")
    out=args.output.resolve(); out.mkdir(parents=True,exist_ok=True)
    aa=L/N_AG; aw=L/N_W
    ag0=lattice(N_AG,aa,((0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)))
    w0=lattice(N_W,aw,((0,0,0),(.5,.5,.5)))
    ga,gw=gyroid(ag0),gyroid(w0)
    mid,ta,tw,_,_,ratio=thresholds(ga,gw,args.gap_level)
    ag=ag0[ga<ta]; w=w0[gw>tw]
    ratio=len(ag)*M_AG/(len(w)*M_W)
    if abs(ratio-1)>.001: raise RuntimeError(f"mass ratio failed: {ratio}")
    cross=cross_min(ag,w); dmin=min(cross["Ag_to_W_min_A"],cross["W_to_Ag_min_A"])
    if dmin<2.6: raise RuntimeError(f"Ag-W contact too short: {dmin}")
    conn=connectivity(ta,tw)
    if not all(conn[s]["connected"] and conn[s]["spans_z"] for s in ("Ag","W")):
        raise RuntimeError("one phase is not bicontinuous")
    data=out/"agw_v11_gyroid.data.gz"; write_data(data,ag,w)
    a10,w10,v10,src=fractions(ta,tw,10.,ns=args.fraction_samples,source=True)
    a5,w5,v5,_=fractions(ta,tw,5.,ns=max(3,args.fraction_samples-2),source=False)
    ag_newton=ag[(ag[:,2]>=40)&(ag[:,2]<=200)]
    qa=(I0/SKIN)*np.exp(-((ag_newton[:,0]-100)**2+(ag_newton[:,1]-100)**2)/R**2
                         -np.maximum(0,SURFACE-ag_newton[:,2])/SKIN)*V_AG
    p_atom=float(qa.sum()); p_geom=float(src.sum())
    scale=p_atom/p_geom; src*=scale
    npz=out/"agw_v11_geometry_fractions.npz"
    np.savez_compressed(npz,fraction_ag_10A=a10.astype("f4"),fraction_w_10A=w10.astype("f4"),
        fraction_void_10A=v10.astype("f4"),source_shape_10A=src.astype("f8"),
        fraction_ag_5A=a5.astype("f4"),fraction_w_5A=w5.astype("f4"),fraction_void_5A=v5.astype("f4"),
        grid10_spacing_A=[10.,10.,10.],grid10_origin_A=[0.,0.,40.],
        grid5_spacing_A=[5.,5.,5.],grid5_origin_A=[0.,0.,40.],
        threshold_mid=mid,threshold_ag=ta,threshold_w=tw)
    va=len(ag)*V_AG; vw=len(w)*V_W
    report={"title":"Ag-W V11 ideal bicontinuous gyroid geometry",
      "status":"GEOMETRY_BUILD_PASS_NOT_YET_LAMMPS_VALIDATED",
      "geometry":"single-period gyroid; W positive labyrinth, Ag complementary labyrinth; explicit interface void band",
      "box_A":[L,L,BOX_Z],"material_region_A":[L,L,L],"boundary":"p p f",
      "atoms":{"Ag":len(ag),"W":len(w),"total":len(ag)+len(w)},
      "mass_ratio_Ag_over_W":ratio,
      "lattices":{"Ag":{"structure":"FCC","cells":N_AG,"a_A":aa,"strain_percent":100*(aa/A_AG-1)},
                  "W":{"structure":"BCC","cells":N_W,"a_A":aw,"strain_percent":100*(aw/A_W-1)}},
      "gyroid":{"gap_level_width":args.gap_level,"threshold_mid":mid,"threshold_ag":ta,"threshold_w":tw},
      "cross_phase_distances":cross,"continuous_geometry_connectivity":conn,
      "atomic_volume_check":{"Ag_A3":va,"W_A3":vw,"solid_fraction_using_TTM_volumes":(va+vw)/L**3},
      "ttm_geometry":{"grid_A":10,"shape":list(a10.shape),"max_combined_fraction":float((a10+w10).max()),
        "mean_fractions":{"Ag":float(a10.mean()),"W":float(w10.mean()),"void":float(v10.mean())},
        "fraction_sum_max_error":float(np.max(np.abs(a10+w10+v10-1))),
        "source_atom_quadrature_eV_ps":p_atom,"source_geometry_integral_before_scale_eV_ps":p_geom,
        "source_scale":scale,"source_power_eV_ps":float(src.sum())},
      "limitations":["Ideal gyroid, not an exact QSGS voxel copy.",
        "Initial interface void band must be relaxed with Dai EFS.",
        "No LAMMPS minimization, cold MD, or TTM step is claimed by this geometry builder.",
        "Geometry fractions are frozen for the short TTM software-link pilot; dynamic remapping is needed after substantial melting."]}
    rp=out/"v11_geometry_report.json"; rp.write_text(json.dumps(report,indent=2)+"\n")
    manifest={"title":"Ag-W V11 exploratory geometry manifest","md_start_gate":"EXPLORATORY_ONLY",
      "static_audit_status":"NOT_RUN","data_file":data.name,"data_sha256":sha(data),
      "geometry_fractions_file":npz.name,"geometry_fractions_sha256":sha(npz),
      "atom_count":len(ag)+len(w),"Ag_atoms":len(ag),"W_atoms":len(w),
      "mass_ratio_Ag_over_W":ratio,"box_A":[L,L,BOX_Z],"grid_A":10,"grid_shape":list(a10.shape),
      "fixed_bottom_region_A":[0,20],"thermostat_buffer_A":[20,40],"ttm_region_A":[40,200],
      "geometry_method":"gyroid level-set volume fractions; phase-restricted atom CIC for lattice kinetic-energy handoff",
      "source":{"I0_eV_ps_A2":I0,"R_A":R,"skin_depth_A":SKIN,"center_xy_A":[100,100],"surface_z_A":SURFACE}}
    mp=out/"v11_exploratory_manifest.json"; mp.write_text(json.dumps(manifest,indent=2)+"\n")
    files=[data,npz,rp,mp]
    (out/"SHA256SUMS.txt").write_text("".join(f"{sha(x)}  {x.name}\n" for x in files))
    print(json.dumps({"status":report["status"],"atoms":report["atoms"],"mass_ratio":ratio,
      "cross_min_A":dmin,"connectivity":conn,"grid_max_fill":report["ttm_geometry"]["max_combined_fraction"],
      "source_eV_ps":p_atom,"output":str(out)},indent=2))
if __name__=="__main__": main()
