#!/usr/bin/env python3
"""Compare frozen and final Ag-W snapshots: density, neighbors, Voro, T and forces."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.spatial import cKDTree
from scipy.stats import kurtosis, skew

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
REPO = Path(__file__).resolve().parents[5]
TTM = MODEL / "ttm"
V9 = MODEL / "interface_reconstruction_v9_hybrid"
V6 = MODEL / "interface_reconstruction_v6_exploratory"
DATA = MODEL / "interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz"
PRE = HERE / "v10_voronoi_force_pre.dump.gz"
FINAL = HERE / "v10_voronoi_force_final.dump.gz"
PAIR_CSV = HERE / "v10_pair_distance_hist.csv"
OVERFILL_CSV = HERE / "v10_5A_overfilled_cells.csv"
OUT = HERE / "v10_local_density_coord_voronoi_temperature.json"
CROSSCHECK_LOG = HERE / "v10_temperature_crosscheck.log"
sys.path.insert(0, str(TTM))
sys.path.insert(0, str(V9))
from agw_ttm_reference import (  # noqa: E402
    AG_ATOMIC_VOLUME_A3, W_ATOMIC_VOLUME_A3,
    _raw_cic_volume_power, gaussian_ag_source_power,
)
from diagnose_agw_ttm_mapping_v7 import read_lammps_data_atoms  # noqa: E402
from analyze_v9_cold_baseline import read_thermo  # noqa: E402

MASS = {1: 183.84, 2: 107.8682}
TYPE_NAME = {1: "W", 2: "Ag"}
M_VV2E = 1.0364269e-4
KB = 8.617333262145e-5
PAIR_NAMES = ["W-W", "Ag-W", "Ag-Ag"]
PAIR_CODES = {"W-W": 0, "Ag-W": 1, "Ag-Ag": 2}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_custom_dump(path: Path) -> tuple[list[str], np.ndarray]:
    op = gzip.open if path.suffix == ".gz" else open
    with op(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("ITEM: ATOMS"):
                cols = line.split()[2:]
                break
        else:
            raise ValueError(f"No ATOMS section in {path}")
        table = np.loadtxt(stream)
    order = np.argsort(table[:, cols.index("id")].astype(np.int64))
    return cols, table[order]


def columns(cols: list[str], table: np.ndarray) -> dict[str, np.ndarray]:
    return {name: table[:, cols.index(name)] for name in cols}


def grid10_counts(xyz: np.ndarray, types: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    nx, ny, nz = 20, 20, 30
    x = np.mod(xyz[:, 0], 200.0)
    y = np.mod(xyz[:, 1], 200.0)
    ix = np.floor(x / 10.0).astype(int).clip(0, nx-1)
    iy = np.floor(y / 10.0).astype(int).clip(0, ny-1)
    iz = np.floor(xyz[:, 2] / 10.0).astype(int).clip(0, nz-1)
    cid = (ix * ny + iy) * nz + iz
    count_ag = np.bincount(cid[types == 2], minlength=nx*ny*nz).reshape((nx,ny,nz))
    count_w = np.bincount(cid[types == 1], minlength=nx*ny*nz).reshape((nx,ny,nz))
    return count_ag, count_w


def summarize_density(count_ag: np.ndarray, count_w: np.ndarray) -> dict:
    both = count_ag + count_w
    material = np.zeros_like(both, dtype=bool)
    material[:, :, :20] = True  # z=0..200 A; exclude top vacuum.
    nonempty = material & (both > 0)
    ag_cells=material&(count_ag>0);w_cells=material&(count_w>0)
    densities = {
        "Ag_atoms_per_A3_where_Ag_present": count_ag[ag_cells] / 1000.0,
        "W_atoms_per_A3_where_W_present": count_w[w_cells] / 1000.0,
        "total_atoms_per_A3_where_any_species_present": both[nonempty] / 1000.0,
    }
    def stats(arr: np.ndarray) -> dict:
        return {"mean":float(arr.mean()),"p05":float(np.quantile(arr,.05)),"p50":float(np.quantile(arr,.5)),
                "p95":float(np.quantile(arr,.95)),"max":float(arr.max())}
    phase_cat = np.full(both.shape, -1, dtype=np.int8)
    phase_cat[material & (count_ag>0) & (count_w==0)] = 0
    phase_cat[material & (count_w>0) & (count_ag==0)] = 1
    phase_cat[material & (count_ag>0) & (count_w>0)] = 2
    names = {0:"Ag_only_10A_cells",1:"W_only_10A_cells",2:"mixed_10A_cells"}
    category_counts = {name:int(np.count_nonzero(phase_cat==code)) for code,name in names.items()}
    return {
        "material_z_scope_A":[0,200],
        "nonempty_10A_cells":int(nonempty.sum()),
        "nonempty_cell_density_atoms_per_A3":{name:stats(values) for name,values in densities.items()},
        "phase_cell_counts":category_counts,
        "phase_cell_code":phase_cat,
    }


def pair_code(t1: np.ndarray, t2: np.ndarray) -> np.ndarray:
    return np.where(t1 == t2, np.where(t1 == 1, 0, 2), 1).astype(np.int8)


def build_neighbor_tree(xyz: np.ndarray, box: float=200.0, margin: float=6.0) -> tuple[cKDTree,np.ndarray]:
    base = xyz.copy()
    base[:,0:2] %= box
    clouds=[base]
    owners=[np.arange(len(xyz),dtype=np.int64)]
    choices=(-box,0.0,box)
    for sx in choices:
        if sx < 0: mx = base[:,0] > box-margin
        elif sx > 0: mx = base[:,0] < margin
        else: mx = np.ones(len(base),dtype=bool)
        for sy in choices:
            if sx==0 and sy==0: continue
            if sy < 0: my = base[:,1] > box-margin
            elif sy > 0: my = base[:,1] < margin
            else: my = np.ones(len(base),dtype=bool)
            take=np.flatnonzero(mx & my)
            if not len(take): continue
            cloud=base[take].copy();cloud[:,0]+=sx;cloud[:,1]+=sy
            clouds.append(cloud);owners.append(take.astype(np.int64))
    tree_xyz=np.concatenate(clouds,axis=0)
    owner=np.concatenate(owners)
    return cKDTree(tree_xyz,leafsize=32,compact_nodes=True,balanced_tree=True),owner


def local_neighbor_analysis(xyz: np.ndarray, types: np.ndarray, region_code: np.ndarray) -> tuple[dict,np.ndarray]:
    tree,owner=build_neighbor_tree(xyz)
    query=xyz.copy();query[:,0:2]%=200.0
    k=24
    dist,tree_idx=tree.query(query,k=k,workers=4)
    neighbor=owner[tree_idx]
    center=np.broadcast_to(np.arange(len(xyz))[:,None],neighbor.shape)
    valid=(neighbor!=center)&np.isfinite(dist)&(dist<5.0)
    center_t=types[center]
    neigh_t=types[neighbor]
    code=pair_code(center_t,neigh_t)
    bin_width=.025
    bins=np.arange(0.0,5.0+bin_width,bin_width)
    bin_id=np.minimum((dist/bin_width).astype(np.int32),len(bins)-2)
    joint=code[valid].astype(np.int32)*(len(bins)-1)+bin_id[valid]
    hist=np.bincount(joint,minlength=3*(len(bins)-1)).reshape(3,-1)
    nearest=np.min(np.where(valid,dist,np.inf),axis=1)

    centers=.5*(bins[:-1]+bins[1:])
    cutoffs={}; peaks={}
    for c,name in enumerate(PAIR_NAMES):
        shell=4.0*np.pi/3.0*(bins[1:]**3-bins[:-1]**3)
        profile=hist[c]/np.maximum(shell,1e-12)
        smoothed=gaussian_filter1d(profile.astype(float),sigma=1.5)
        peak_window=(centers>=1.6)&(centers<=3.7)
        peak_i=int(np.argmax(np.where(peak_window,smoothed,-np.inf)))
        minima,_=find_peaks(-smoothed)
        after=minima[(minima>peak_i)&(centers[minima]<4.8)]
        cutoff_i=int(after[0]) if len(after) else int(np.argmin(smoothed[peak_i+1:]/1.0)+peak_i+1)
        cutoff=float(centers[cutoff_i])
        cutoffs[name]=cutoff
        peaks[name]={"first_shell_peak_A":float(centers[peak_i]),"first_minimum_cutoff_A":cutoff}

    # Per-center species-resolved first-shell coordination using the first
    # minima found from the typed neighbor-distance distributions.
    coord=np.zeros((3,len(xyz)),dtype=np.int16)
    for c,name in enumerate(PAIR_NAMES):
        take=valid&(code==c)&(dist<cutoffs[name])
        coord[c]=np.count_nonzero(take,axis=1)
    face_summary={}
    category_names={0:"Ag_only_cells",1:"W_only_cells",2:"mixed_interface_cells"}
    for t in (1,2):
        for cat,catname in category_names.items():
            select=(types==t)&(region_code==cat)&(nearest<np.inf)
            if not np.any(select): continue
            if t==1: pair_codes=(0,1)
            else: pair_codes=(1,2)
            for c in pair_codes:
                cname=PAIR_NAMES[c]
                arr=coord[c,select]
                face_summary[f"{TYPE_NAME[t]}_centers_{catname}_neighbors_{cname}"]={
                    "atoms":int(select.sum()),"mean":float(arr.mean()),"p50":float(np.quantile(arr,.5)),
                    "p95":float(np.quantile(arr,.95)),"zero_fraction":float(np.mean(arr==0)),
                }
    dist_summary={}
    for c,name in enumerate(PAIR_NAMES):
        values=dist[valid&(code==c)]
        dist_summary[name]={"candidate_directed_neighbors_under_5A":int(len(values)),
                            "min_A":float(values.min()),"p01_A":float(np.quantile(values,.01)),
                            "p50_A":float(np.quantile(values,.5)),"p99_A":float(np.quantile(values,.99))}
    result={
        "xy_periodic_z_open":True,"search_radius_A":5.0,"queried_neighbors_per_atom":k-1,
        "first_shell_cutoffs_A":cutoffs,"first_shell_peaks_and_minima":peaks,
        "typed_neighbor_distance_summary":dist_summary,
        "nearest_neighbor_summary_A":{"min":float(nearest.min()),"p01":float(np.quantile(nearest,.01)),
            "p50":float(np.quantile(nearest,.5)),"p95":float(np.quantile(nearest,.95)),"max":float(nearest.max())},
        "pair_coordination_by_first_minimum":face_summary,
        "coordination_definition":"typed neighbor counts below the first minimum in the smooth typed pair-distance distribution; 10 A cell phase categories",
    }
    return result,hist


def temperature_summary(xyz: np.ndarray, vel: np.ndarray, types: np.ndarray, ids: np.ndarray,
                        fixed_ids: np.ndarray, phase_cat: np.ndarray) -> dict:
    nxy,nz=20,30
    x=np.mod(xyz[:,0],200);y=np.mod(xyz[:,1],200)
    ix=np.floor(x/10).astype(int).clip(0,nxy-1);iy=np.floor(y/10).astype(int).clip(0,nxy-1)
    iz=np.floor(xyz[:,2]/10).astype(int).clip(0,nz-1)
    cats=phase_cat.ravel()[(ix*nxy+iy)*nz+iz]
    fixed=np.isin(ids,fixed_ids)
    mobile=~fixed
    summaries={}
    cat_lookup={0:"Ag_only_10A_cells",1:"W_only_10A_cells",2:"mixed_interface_10A_cells"}
    selections={"all_mobile":mobile}
    for t in (2,1):
        selections[f"{TYPE_NAME[t]}_mobile"] = mobile&(types==t)
        for cat,nm in cat_lookup.items():
            selections[f"{TYPE_NAME[t]}_{nm}"]=mobile&(types==t)&(cats==cat)
    for name,sel in selections.items():
        idx=np.flatnonzero(sel)
        if len(idx)<3: continue
        masses=np.asarray([MASS[int(t)] for t in types[idx]])
        vv=vel[idx]
        com=np.average(vv,axis=0,weights=masses)
        rel=vv-com
        ke_lammps=np.sum(.5*masses[:,None]*M_VV2E*vv**2)
        ke_peculiar=np.sum(.5*masses[:,None]*M_VV2E*rel**2)
        dof=max(3*len(idx)-3,1)
        temp_lammps=2*ke_lammps/(dof*KB)
        temp_peculiar=2*ke_peculiar/(dof*KB)
        normalized=rel*np.sqrt(masses[:,None]*M_VV2E/(KB*temp_peculiar))
        flat=normalized.ravel()
        summaries[name]={"atoms":int(len(idx)),
            "temperature_K_LAMMPS_compute_temp_convention_no_velocity_bias":float(temp_lammps),
            "temperature_K_after_removing_group_COM_drift_diagnostic":float(temp_peculiar),
            "kinetic_energy_eV_including_group_drift":float(ke_lammps),
            "group_COM_velocity_A_ps":com.tolist(),
            "normalized_velocity_components":{"mean":float(flat.mean()),"variance":float(flat.var()),
              "skewness":float(skew(flat)),"kurtosis_pearson":float(kurtosis(flat,fisher=False))}}
    # Cellwise species temperatures; only cells with at least six moving atoms.
    cell=(ix*nxy+iy)*nz+iz
    local=[]
    for t in (1,2):
        idx=np.flatnonzero(mobile&(types==t))
        cid=cell[idx];vv=vel[idx];mass=MASS[t]
        count=np.bincount(cid,minlength=nxy*nxy*nz)
        sums=np.stack([np.bincount(cid,weights=vv[:,a],minlength=nxy*nxy*nz) for a in range(3)],axis=1)
        sums2=np.bincount(cid,weights=np.sum(vv*vv,axis=1),minlength=nxy*nxy*nz)
        good=np.flatnonzero(count>=6)
        if len(good):
            centered=sums2[good]-np.sum(sums[good]**2,axis=1)/count[good]
            dof=3*(count[good]-1)
            tl=mass*M_VV2E*np.maximum(centered,0.0)/(dof*KB)
            local.append({"element":TYPE_NAME[t],"cells":int(len(good)),"atoms_per_cell_min":6,
                "temperature_K":{"p05":float(np.quantile(tl,.05)),"p50":float(np.quantile(tl,.5)),
                                 "p95":float(np.quantile(tl,.95)),"max":float(tl.max())}})
    return {"fixed_bottom_ids_excluded":int(fixed.sum()),
        "group_temperature_note":"LAMMPS compute temp removes translational degrees of freedom from the DOF count but does not subtract group COM velocity from kinetic energy; both that value and an explicitly COM-removed diagnostic are reported. Cellwise temperatures remove local mean velocity.",
        "groups":summaries,"cellwise_temperature":local}


def force_summary(force: np.ndarray, types: np.ndarray, phase_cat: np.ndarray, ids: np.ndarray,
                  fixed_ids: np.ndarray, xyz: np.ndarray) -> dict:
    f=np.linalg.norm(force,axis=1)
    mobile=~np.isin(ids,fixed_ids)
    rows=[]
    category_names={0:"Ag_only_cells",1:"W_only_cells",2:"mixed_interface_cells"}
    for t in (1,2):
        s=mobile&(types==t);a=f[s]
        rows.append({"element":TYPE_NAME[t],"atoms":int(s.sum()),"p50_eV_A":float(np.quantile(a,.5)),
                     "p95_eV_A":float(np.quantile(a,.95)),"p99_eV_A":float(np.quantile(a,.99)),
                     "above_2eV_A":int(np.count_nonzero(a>2.0)),"above_3eV_A":int(np.count_nonzero(a>3.0))})
    by_phase=[]
    for t in (1,2):
        for code,name in category_names.items():
            sel=mobile&(types==t)&(phase_cat==code)
            if not np.any(sel): continue
            a=f[sel]
            by_phase.append({"element":TYPE_NAME[t],"10A_cell_region":name,"atoms":int(sel.sum()),
                "p50_eV_A":float(np.quantile(a,.5)),"p95_eV_A":float(np.quantile(a,.95)),
                "p99_eV_A":float(np.quantile(a,.99)),"max_eV_A":float(a.max()),
                "above_2eV_A":int(np.count_nonzero(a>2.0)),"above_3eV_A":int(np.count_nonzero(a>3.0))})
    top=np.argsort(np.where(mobile,f,-1.0))[-25:][::-1]
    return {"mobile_max_eV_A":float(f[mobile].max()),"mobile_rms_eV_A":float(np.sqrt(np.mean(f[mobile]**2))),
            "mobile_above_2eV_A":int(np.count_nonzero(mobile&(f>2))),"mobile_above_3eV_A":int(np.count_nonzero(mobile&(f>3))),
            "species":rows,"species_by_10A_phase_region":by_phase,
            "top_25_mobile_atoms":[{"id":int(ids[i]),"type":TYPE_NAME[int(types[i])],"force_eV_A":float(f[i]),
                "xyz_A":xyz[i].tolist(),"grid_region_code":int(phase_cat[i])} for i in top]}


def state_report(path: Path, label: str, fixed_ids: np.ndarray) -> tuple[dict,np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    cols,table=load_custom_dump(path);d=columns(cols,table)
    ids=d["id"].astype(np.int64);types=d["type"].astype(np.int32)
    xyz=np.column_stack([d[a] for a in ("x","y","z")]);vel=np.column_stack([d[a] for a in ("vx","vy","vz")])
    force=np.column_stack([d[a] for a in ("fx","fy","fz")]);voro=d["c_vcell[1]"];faces=d["c_vcell[2]"]
    if len(ids)!=482601 or len(np.unique(ids))!=len(ids): raise RuntimeError(f"Atom count/ID issue in {path}")
    if int(np.count_nonzero(types==2))!=304046 or int(np.count_nonzero(types==1))!=178555: raise RuntimeError("Species count changed")
    ca,cw=grid10_counts(xyz,types);density=summarize_density(ca,cw);phase=density.pop("phase_cell_code")
    ix=np.floor(np.mod(xyz[:,0],200)/10).astype(int).clip(0,19)
    iy=np.floor(np.mod(xyz[:,1],200)/10).astype(int).clip(0,19)
    iz=np.floor(xyz[:,2]/10).astype(int).clip(0,29)
    cell=(ix*20+iy)*30+iz
    phase_atom=phase.ravel()[cell]
    neighbor,hist=local_neighbor_analysis(xyz,types,phase_atom)
    temps=temperature_summary(xyz,vel,types,ids,fixed_ids,phase)
    forces=force_summary(force,types,phase_atom,ids,fixed_ids,xyz)
    interior=(xyz[:,2]>=20)&(xyz[:,2]<190)
    surface=(xyz[:,2]>=190)&(xyz[:,2]<220)
    bands={"interior_20_to_190A":interior,"surface_sensitive_190_to_220A":surface,"z_above_200A":xyz[:,2]>200}
    voro_summary={}
    for nm,mask in bands.items():
        for t in (1,2):
            m=mask&(types==t)
            if np.any(m):
                voro_summary[f"{nm}_{TYPE_NAME[t]}"]={"atoms":int(m.sum()),
                    "volume_A3":{"p05":float(np.quantile(voro[m],.05)),"p50":float(np.quantile(voro[m],.5)),
                                 "p95":float(np.quantile(voro[m],.95)),"max":float(voro[m].max())},
                    "faces":{"p50":float(np.quantile(faces[m],.5)),"p95":float(np.quantile(faces[m],.95))}}
    for t,ref in ((1,W_ATOMIC_VOLUME_A3),(2,AG_ATOMIC_VOLUME_A3)):
        key=f"interior_20_to_190A_{TYPE_NAME[t]}"
        if key in voro_summary:
            med=voro_summary[key]["volume_A3"]["p50"]
            voro_summary[key]["reference_atomic_volume_A3_for_diagnostic_only"] = float(ref)
            voro_summary[key]["median_minus_reference_A3"] = float(med-ref)
            voro_summary[key]["median_relative_difference_percent"] = float(100.*(med/ref-1.))
    report={"snapshot":label,"path":str(path.relative_to(MODEL)),"sha256":sha(path),
        "atoms":int(len(ids)),"Ag":int(np.count_nonzero(types==2)),"W":int(np.count_nonzero(types==1)),
        "box_A":[200,200,300],"boundary":"p p f","density_10A":density,
        "typed_pair_distances_and_coordination":neighbor,"voronoi":{"method":"LAMMPS compute voronoi/atom",
            "interpretation":"Voronoi cell volumes partition the geometric domain; they are not treated as atomic solid occupied volume. Surface-sensitive atoms are reported separately and excluded from the interior comparison.",
            "region_statistics":voro_summary},
        "temperature":temps,"forces":forces}
    del table,d
    return report,ca,cw,hist,xyz,types


def overfill_by_phase(xyz: np.ndarray,types: np.ndarray,static_newton: np.ndarray) -> list[dict]:
    inside=static_newton&(xyz[:,2]>=40-1e-9)&(xyz[:,2]<=200+1e-9)
    ag=np.flatnonzero(inside&(types==2));w=np.flatnonzero(inside&(types==1))
    q=gaussian_ag_source_power(xyz[ag])
    rows=[]
    for dx,dy in ((0.,0.),(0.,2.5),(2.5,0.),(2.5,2.5)):
        o=np.array([dx,dy,40.])
        _,va,_=_raw_cic_volume_power(xyz[ag]-o,q,(40,40,32),(5.,)*3,AG_ATOMIC_VOLUME_A3)
        _,vw,_=_raw_cic_volume_power(xyz[w]-o,np.zeros(len(w)),(40,40,32),(5.,)*3,W_ATOMIC_VOLUME_A3)
        total=va+vw;overflow=np.maximum(total-125.,0.);mask=overflow>0
        share=np.divide(va,total,out=np.zeros_like(va),where=total>0)
        category=np.where(share>=.8,0,np.where(share<=.2,1,2))
        names={0:"Ag_rich",1:"W_rich",2:"mixed_interface"}
        sums={}
        for c,nm in names.items():
            selected=mask&(category==c)
            sums[nm]={"cells":int(selected.sum()),"overflow_A3":float(overflow[selected].sum())}
        flat=np.argsort(overflow.ravel())[-25:][::-1]
        top=[]
        for cell in flat:
            if overflow.ravel()[cell]<=0: break
            ix,iy,iz=np.unravel_index(int(cell),overflow.shape)
            top.append({"cell_ijk":[int(ix),int(iy),int(iz)],
                "center_xyz_A":[float(dx+(ix+.5)*5),float(dy+(iy+.5)*5),float(40+(iz+.5)*5)],
                "raw_fill":float(total.ravel()[cell]/125),"overflow_A3":float(overflow.ravel()[cell]),
                "Ag_volume_fraction_of_material":float(share.ravel()[cell]),
                "region":names[int(category.ravel()[cell])]})
        rows.append({"offset_xy_A":[dx,dy],"raw_overflow_A3":float(overflow.sum()),"cells_overfilled":int(mask.sum()),
                     "overflow_by_local_phase":sums,"top_25_cells":top})
    return rows


def potential_energy_rates() -> dict:
    result={}
    for name,file in (("first_0p5ps","v9_cold_nvt_nve.log"),("followup_0p5ps","v9_cold_followup.log")):
        cols,rows=read_thermo(V9/file);idx={x:cols.index(x) for x in ("Time","PotEng","c_T_mobile","c_T_Ag","c_T_W")}
        a=np.asarray(rows);a=a[a[:,idx["Time"]]<=.50001]
        def fit(lo,hi):
            m=(a[:,idx["Time"]]>=lo)&(a[:,idx["Time"]]<=hi)
            rate=float(np.polyfit(a[m,idx["Time"]],a[m,idx["PotEng"]],1)[0])
            return {"interval_ps":[lo,hi],"potential_energy_rate_eV_ps":rate,
                    "temperature_start_end_K":[float(a[m,idx["c_T_mobile"]][0]),float(a[m,idx["c_T_mobile"]][-1])]}
        result[name]={"first_0p2ps":fit(0,.2),"middle_0p2ps":fit(.15,.35),"last_0p2ps":fit(.3,.5)}
    return result


def temperature_crosscheck(final: dict) -> dict:
    header,rows=read_thermo(CROSSCHECK_LOG)
    idx={name:header.index(name) for name in ("Step","Atoms","KinEng","c_T_mobile","c_T_W","c_T_Ag")}
    row=max(rows,key=lambda x:x[idx["Step"]])
    lmp={"step":int(row[idx["Step"]]),"atoms":int(row[idx["Atoms"]]),
        "kinetic_energy_eV":float(row[idx["KinEng"]]),
        "all_mobile_K":float(row[idx["c_T_mobile"]]),"W_mobile_K":float(row[idx["c_T_W"]]),
        "Ag_mobile_K":float(row[idx["c_T_Ag"]])}
    py=final["temperature"]["groups"]
    comp={"all_mobile_K":py["all_mobile"]["temperature_K_LAMMPS_compute_temp_convention_no_velocity_bias"],
        "W_mobile_K":py["W_mobile"]["temperature_K_LAMMPS_compute_temp_convention_no_velocity_bias"],
        "Ag_mobile_K":py["Ag_mobile"]["temperature_K_LAMMPS_compute_temp_convention_no_velocity_bias"]}
    delta={k:float(comp[k]-lmp[k]) for k in comp}
    return {"log":CROSSCHECK_LOG.name,"run_type":"LAMMPS run 0 only; no atom integration",
        "LAMMPS_values":lmp,"Python_recomputed_values":comp,"Python_minus_LAMMPS_K":delta,
        "agreement_within_0p01K":bool(max(abs(x) for x in delta.values())<.01)}


def main() -> None:
    data_ids,data_types,data_xyz=read_lammps_data_atoms(DATA)
    fixed_ids=data_ids[data_xyz[:,2]<20]
    pre,pre_ag,pre_w,pre_hist,xyz_pre,types_pre=state_report(PRE,"frozen V6/V7 restart before V9 cold continuation",fixed_ids)
    final,fin_ag,fin_w,fin_hist,xyz_final,types_final=state_report(FINAL,"V9 final NVE restart after 1 ps total NVT and 0.1 ps total NVE",fixed_ids)
    if not (np.array_equal(types_pre,types_final) and np.array_equal(data_ids,read_lammps_data_atoms(DATA)[0])):
        raise RuntimeError("pre/final element IDs differ")
    moved_counts={"Ag_L1_cell_count_change":int(np.abs(fin_ag-pre_ag).sum()),
                  "W_L1_cell_count_change":int(np.abs(fin_w-pre_w).sum()),
                  "Ag_equivalent_atoms_reassigned_between_10A_cells":float(.5*np.abs(fin_ag-pre_ag).sum()),
                  "W_equivalent_atoms_reassigned_between_10A_cells":float(.5*np.abs(fin_w-pre_w).sum()),
                  "Ag_fraction_of_counts_reassigned":float(.5*np.abs(fin_ag-pre_ag).sum()/fin_ag.sum()),
                  "W_fraction_of_counts_reassigned":float(.5*np.abs(fin_w-pre_w).sum()/fin_w.sum()),
                  "Ag_cells_changed":int(np.count_nonzero(fin_ag!=pre_ag)),
                  "W_cells_changed":int(np.count_nonzero(fin_w!=pre_w)),
                  "density_grid_10A_origin_A":[0,0,0]}
    with PAIR_CSV.open("w",newline="",encoding="utf-8") as f:
        wri=csv.writer(f,lineterminator="\n");wri.writerow(["r_center_A","W_W_directed_count_pre","Ag_W_directed_count_pre","Ag_Ag_directed_count_pre",
                                        "W_W_directed_count_final","Ag_W_directed_count_final","Ag_Ag_directed_count_final"])
        centers=np.arange(len(pre_hist[0]))*.025+.0125
        for row in zip(centers,*pre_hist,*fin_hist): wri.writerow(row)
    static_newton_ids=data_ids[(data_xyz[:,2]>=40.0-1e-9)&(data_xyz[:,2]<=300.0+1e-9)]
    static_newton=np.isin(data_ids,static_newton_ids)
    overfill=overfill_by_phase(xyz_final,types_final,static_newton)
    newton_current_z=xyz_final[static_newton,2]
    newton_z_counts={"static_Newton_atoms":int(static_newton.sum()),
        "currently_in_40_to_200A_inclusive":int(np.count_nonzero((newton_current_z>=40.-1e-9)&(newton_current_z<=200.+1e-9))),
        "current_below_40A":int(np.count_nonzero(newton_current_z<40.-1e-9)),
        "current_above_200A":int(np.count_nonzero(newton_current_z>200.+1e-9)),
        "current_exactly_at_200A_within_tolerance":int(np.count_nonzero(np.abs(newton_current_z-200.)<=1e-9)),
        "mapping_gate_note":"Mapping uses the frozen 40..200 A region; displaced Newton-group atoms outside it are excluded rather than clipped into edge cells."}
    with OVERFILL_CSV.open("w",newline="",encoding="utf-8") as f:
        wri=csv.writer(f,lineterminator="\n");wri.writerow(["offset_x_A","offset_y_A","cell_i","cell_j","cell_k","center_x_A","center_y_A","center_z_A","raw_fill","overflow_A3","Ag_fraction","region"])
        for ph in overfill:
            for cell in ph["top_25_cells"]:
                wri.writerow(ph["offset_xy_A"]+cell["cell_ijk"]+cell["center_xyz_A"]+[cell["raw_fill"],cell["overflow_A3"],cell["Ag_volume_fraction_of_material"],cell["region"]])
    result={
        "title":"Ag-W V10 local density, Voronoi, neighborhood, temperature and force diagnostics",
        "frozen_v9_commit":"579c855f25f8286a3cdb535ea7336be95962f994",
        "snapshots":{"pre_v9_cold":pre,"v9_final":final},
        "LAMMPS_temperature_crosscheck":temperature_crosscheck(final),
        "pre_to_final_10A_density_change":moved_counts,
        "v9_final_static_newton_z_location_counts":newton_z_counts,
        "overfilled_5A_cells_by_phase_final_state":overfill,
        "cold_relaxation_energy_rates":potential_energy_rates(),
        "voronoi_caveat":"Use only z=20..190 A interior cells for material-volume comparison. LAMMPS Voronoi cells partition space, including interstitial regions; no Voronoi cell volume is substituted for the fixed atomic-volume constants in the TTM mapper.",
        "decision":"Diagnostics only. No box, atoms, IDs or atomic volumes changed; no TTM step/source executed.",
        "pair_distance_histogram_csv":PAIR_CSV.name,
        "top_overfill_cells_csv":OVERFILL_CSV.name,
    }
    OUT.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"output":str(OUT),"density_change":moved_counts,
        "final_mobile_force_max_eV_A":final["forces"]["mobile_max_eV_A"],
        "final_overfill_phase_min":min(overfill,key=lambda x:x["raw_overflow_A3"])["offset_xy_A"],
        "pre_temperature":pre["temperature"]["groups"].get("all_mobile"),
        "final_temperature":final["temperature"]["groups"].get("all_mobile")},indent=2))


if __name__=="__main__": main()
