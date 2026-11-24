#!/usr/bin/env python3
"""Diagnostic-only transfer-distance sensitivity on the frozen V9 final state."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent
V9 = MODEL / "interface_reconstruction_v9_hybrid"
TTM = MODEL / "ttm"
sys.path.insert(0, str(TTM))
from agw_ttm_reference import (  # noqa: E402
    AG_ATOMIC_VOLUME_A3, W_ATOMIC_VOLUME_A3, _raw_cic_volume_power,
    gaussian_ag_source_power, joint_capacity_maxflow,
)
from diagnose_agw_ttm_mapping_v7 import read_dump, read_lammps_data_atoms  # noqa: E402

DATA = MODEL / "interface_reconstruction_v5/efs_static_v1_fixed_bottom/efs_minimized.data.gz"
RESTART = V9 / "v9_nve2_end.dump.gz"
OUT = HERE / "v10_transfer_distance_sensitivity.json"
BIN_EDGES = np.array([0., 10., 20., 30., 40., 60., 100., 150.])


def profile(power: np.ndarray, h: float, dx: float, dy: float) -> dict:
    nx, ny, nz = power.shape
    xx = (dx + (np.arange(nx) + .5) * h) % 200.
    yy = (dy + (np.arange(ny) + .5) * h) % 200.
    x, y = np.meshgrid(xx, yy, indexing="ij")
    x = (x - 100. + 100.) % 200. - 100.
    y = (y - 100. + 100.) % 200. - 100.
    r = np.broadcast_to(np.sqrt(x*x + y*y)[:, :, None], power.shape)
    z = 40. + (np.arange(nz) + .5) * h
    depth = 200. - z[None, None, :]
    total = float(power.sum())
    bins = np.digitize(r.ravel(), BIN_EDGES[1:-1])
    radial = np.bincount(bins, weights=power.ravel(), minlength=len(BIN_EDGES)-1)
    return {
        "power_eV_ps": total,
        "radial_bins_A": BIN_EDGES.tolist(),
        "radial_power_fraction": (radial / total).tolist() if total else None,
        "radial_RMS_A": float(np.sqrt(np.sum(power*r*r)/total)) if total else None,
        "mean_depth_A": float(np.sum(power*depth)/total) if total else None,
    }


def apply_transfer(raw_ag: np.ndarray, raw_w: np.ndarray, raw_source: np.ndarray,
                   flow: dict, h: float) -> tuple[np.ndarray, np.ndarray, dict]:
    ag, w, source = raw_ag.ravel().copy(), raw_w.ravel().copy(), raw_source.ravel().copy()
    nx, ny, nz = raw_ag.shape
    dist, moved = [], []
    for donor, targets in flow.get("transfer_fractions", {}).items():
        donor = int(donor)
        ix, iy, iz = donor // (ny*nz), (donor//nz) % ny, donor % nz
        for target, frac in targets:
            target = int(target)
            tx, ty, tz = target // (ny*nz), (target//nz) % ny, target % nz
            dx = min((tx-ix) % nx, (ix-tx) % nx) * h
            dy = min((ty-iy) % ny, (iy-ty) % ny) * h
            d = float(np.sqrt(dx*dx + dy*dy + ((tz-iz)*h)**2))
            moved.append((donor, target, float(frac)))
            dist.append((d, float((raw_ag.ravel()[donor]+raw_w.ravel()[donor])*frac)))
            va = raw_ag.ravel()[donor] * frac
            vw = raw_w.ravel()[donor] * frac
            vs = raw_source.ravel()[donor] * frac
            ag[donor] -= va; ag[target] += va
            w[donor] -= vw; w[target] += vw
            source[donor] -= vs; source[target] += vs
    td = np.asarray(dist, dtype=float)
    distance_info = {
        "actual_transfer_edges": len(moved),
        "actual_transfer_volume_A3": float(td[:,1].sum()) if len(td) else 0.,
        "actual_mean_distance_A_volume_weighted": float(np.average(td[:,0], weights=td[:,1])) if len(td) and td[:,1].sum() else None,
        "actual_max_distance_A": float(td[:,0].max()) if len(td) else None,
        "Ag_volume_conservation_error_A3": float(ag.sum()-raw_ag.sum()),
        "W_volume_conservation_error_A3": float(w.sum()-raw_w.sum()),
        "source_power_conservation_error_eV_ps": float(source.sum()-raw_source.sum()),
    }
    return source.reshape(raw_source.shape), ag.reshape(raw_ag.shape), w.reshape(raw_w.shape), distance_info


def diagnostic_partial_fractions(raw: np.ndarray, spacing: tuple[float,float,float], maxdist: float,
                                 capacity_resolution: float=1e-12) -> tuple[dict,int]:
    """Reconstruct the same direct donor/receiver graph to expose a partial flow.

    The production helper intentionally omits transfer fractions when its gate
    fails. This V10-only duplicate is diagnostic; it uses the same integer
    capacities, geometric edges and Dinic max-flow, and its flow value is
    checked against the production helper before use.
    """
    raw=np.asarray(raw,dtype=float); sx,sy,sz=spacing
    flat=raw.ravel(); cell_volume=float(np.prod(spacing))
    excess=np.maximum(flat-cell_volume,0.); free=np.maximum(cell_volume-flat,0.)
    threshold=max(cell_volume*1e-12,1e-14)
    donors=np.flatnonzero(excess>threshold); receivers=np.flatnonzero(free>threshold)
    demand=float(excess.sum())
    scale=min(1./capacity_resolution,1e9/max(demand,float(free.max(initial=0.)),1e-30))
    resolution=1./scale
    demand_units=np.floor(excess[donors]*scale).astype(np.int64)
    free_units=np.floor(free[receivers]*scale).astype(np.int64)
    use_d=demand_units>0; use_r=free_units>0
    donors,demand_units=donors[use_d],demand_units[use_d]
    receivers,free_units=receivers[use_r],free_units[use_r]
    nx,ny,nz=raw.shape; ncell=raw.size
    dnode=np.full(ncell,-1,dtype=np.int64); rnode=np.full(ncell,-1,dtype=np.int64)
    dnode[donors]=np.arange(1,len(donors)+1,dtype=np.int64)
    roffset=1+len(donors); rnode[receivers]=roffset+np.arange(len(receivers),dtype=np.int64)
    sink=roffset+len(receivers); demand_total=int(demand_units.sum(dtype=np.int64))
    rows=[]; cols=[]; caps=[]
    for donor,units in zip(donors,demand_units):
        rows.append(0);cols.append(int(dnode[donor]));caps.append(int(units))
    rx=min(nx//2,int(np.ceil(maxdist/sx))); ry=min(ny//2,int(np.ceil(maxdist/sy)))
    rz=min(nz-1,int(np.ceil(maxdist/sz)))
    for donor in donors:
        ix=int(donor//(ny*nz));iy=int((donor//nz)%ny);iz=int(donor%nz);seen=set()
        for dx in range(-rx,rx+1):
            tx=(ix+dx)%nx; ddx=min(abs(dx),nx-abs(dx))*sx
            for dy in range(-ry,ry+1):
                ty=(iy+dy)%ny;ddy=min(abs(dy),ny-abs(dy))*sy
                for dz in range(-rz,rz+1):
                    tz=iz+dz
                    if tz<0 or tz>=nz or (dx==0 and dy==0 and dz==0): continue
                    if np.sqrt(ddx*ddx+ddy*ddy+(dz*sz)**2)>maxdist+1e-12: continue
                    target=int(np.ravel_multi_index((tx,ty,tz),raw.shape));rn=int(rnode[target])
                    if rn<0 or target in seen: continue
                    seen.add(target);rows.append(int(dnode[donor]));cols.append(rn);caps.append(demand_total)
    for receiver,units in zip(receivers,free_units):
        rows.append(int(rnode[receiver]));cols.append(sink);caps.append(int(units))
    graph=coo_matrix((np.asarray(caps,dtype=np.int64),(rows,cols)),shape=(sink+1,sink+1)).tocsr()
    solved=maximum_flow(graph,0,sink,method="dinic")
    transfer={}
    for donor in donors:
        node=int(dnode[donor]);start,end=solved.flow.indptr[node],solved.flow.indptr[node+1]
        for dest,units in zip(solved.flow.indices[start:end],solved.flow.data[start:end]):
            if units<=0 or dest<roffset or dest>=sink: continue
            target=int(receivers[int(dest)-roffset])
            transfer.setdefault(int(donor),[]).append((target,float(units)*resolution/float(flat[donor])))
    return transfer,int(solved.flow_value)


def main() -> None:
    data_ids, data_types, data_xyz = read_lammps_data_atoms(DATA)
    ids, types, xyz = read_dump(RESTART)
    if not (np.array_equal(data_ids, ids) and np.array_equal(data_types, types)):
        raise RuntimeError("Frozen V9 state changed atom IDs/types")
    static_newton = (data_xyz[:,2] >= 40.-1e-9) & (data_xyz[:,2] <= 300.+1e-9)
    inside = static_newton & (xyz[:,2] >= 40.-1e-9) & (xyz[:,2] <= 200.+1e-9)
    ia = np.flatnonzero(inside & (types==2)); iw = np.flatnonzero(inside & (types==1))
    h, dx, dy = 10., 5., 5.
    origin = np.array([dx,dy,40.])
    q = gaussian_ag_source_power(xyz[ia])
    _, raw_ag, raw_source = _raw_cic_volume_power(xyz[ia]-origin,q,(20,20,16),(h,)*3,AG_ATOMIC_VOLUME_A3)
    _, raw_w, _ = _raw_cic_volume_power(xyz[iw]-origin,np.zeros(len(iw)),(20,20,16),(h,)*3,W_ATOMIC_VOLUME_A3)
    raw_total = raw_ag+raw_w
    raw_profile = profile(raw_source,h,dx,dy)
    rows=[]
    for maxdist in (20.,25.,30.,40.):
        flow=joint_capacity_maxflow(raw_total,(h,)*3,maxdist)
        partial,partial_units=diagnostic_partial_fractions(raw_total,(h,)*3,maxdist)
        if partial_units!=int(flow.get("flow_units",0)):
            raise RuntimeError(f"Partial diagnostic flow differs from production at {maxdist:g} A")
        flow_for_profile=dict(flow)
        flow_for_profile["transfer_fractions"]=(flow["transfer_fractions"] if flow["status"]=="PASS" else partial)
        mapped_source,mapped_ag,mapped_w,distinfo=apply_transfer(raw_ag,raw_w,raw_source,flow_for_profile,h)
        mapped_total=mapped_ag+mapped_w
        mapped_profile=profile(mapped_source,h,dx,dy)
        fill=mapped_total/1000.
        raw_fraction=np.divide(raw_ag,raw_total,out=np.zeros_like(raw_ag),where=raw_total>0)
        mapped_fraction=np.divide(mapped_ag,mapped_total,out=np.zeros_like(mapped_ag),where=mapped_total>0)
        occupied=raw_total>0
        fraction_delta=np.abs(mapped_fraction[occupied]-raw_fraction[occupied])
        radial_delta=[float(a-b) for a,b in zip(mapped_profile["radial_power_fraction"],raw_profile["radial_power_fraction"])]
        rows.append({
            "configured_max_transfer_distance_A":maxdist,
            "solver_status":flow["status"],
            "mathematical_assignment_complete":flow["status"]=="PASS",
            "overflow_demand_A3":flow["overflow_demand_A3"],
            "locally_unassignable_overflow_A3":flow["locally_unassignable_overflow_A3"],
            "global_free_capacity_A3":flow["global_free_capacity_A3"],
            "candidate_edges":flow["candidate_edges"],
            "partial_flow_diagnostic_used":flow["status"]!="PASS",
            "diagnostic_flow_units_match_production":True,
            "diagnostic_partial_assignment_only_not_a_valid_mapping":flow["status"]!="PASS",
            "actual_or_partial_transfer_diagnostic":distinfo,
            "max_fill_after_complete_or_partial_diagnostic_assignment":float(fill.max()),
            "composition_spatial_change":{
                "cells_with_nonzero_raw_volume":int(occupied.sum()),
                "mean_absolute_Ag_fraction_change":float(fraction_delta.mean()) if len(fraction_delta) else 0.,
                "p95_absolute_Ag_fraction_change":float(np.quantile(fraction_delta,.95)) if len(fraction_delta) else 0.,
                "max_absolute_Ag_fraction_change":float(fraction_delta.max()) if len(fraction_delta) else 0.,
            },
            "source_profile_before":raw_profile,
            "source_profile_after_complete_or_partial_diagnostic_assignment":mapped_profile,
            "source_profile_radial_RMS_change_A":float(mapped_profile["radial_RMS_A"]-raw_profile["radial_RMS_A"]),
            "source_profile_radial_fraction_change":radial_delta,
            "source_mean_depth_change_A":float(mapped_profile["mean_depth_A"]-raw_profile["mean_depth_A"]),
            "interpretation":"partial max-flow is exposed diagnostically when production fails; no incomplete result is a valid mapping, and even a complete assignment at this radius is not physically approved",
        })
        print(json.dumps({"max_distance_A":maxdist,"status":flow["status"],
                          "unassigned_A3":flow["locally_unassignable_overflow_A3"],
                          "max_fill_after_assignment":float(fill.max()),
                          "partial_transfer_edges":distinfo["actual_transfer_edges"]}),flush=True)
    result={
        "title":"V10 transfer-distance sensitivity; diagnostic only",
        "frozen_v9_commit":"579c855f25f8286a3cdb535ea7336be95962f994",
        "input":{"data_sha256":hashlib.sha256(DATA.read_bytes()).hexdigest(),
                 "restart_sha256":hashlib.sha256(RESTART.read_bytes()).hexdigest(),
                 "atoms_total":int(len(ids)),"mapped_Ag":int(len(ia)),"mapped_W":int(len(iw)),
                 "grid_A":h,"offset_xy_A":[dx,dy],"atomic_volumes_changed":False,
                 "coordinates_or_box_changed":False},
        "raw_source_profile":raw_profile,
        "sensitivity_results":rows,
        "decision":"No mapping pass is declared from transfer-range growth. Review actual distance, spatial redistribution and source-profile changes; no TTM/LAMMPS step was run.",
        "ttm_or_laser_step_executed":False,
    }
    OUT.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__": main()
