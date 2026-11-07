#!/usr/bin/env python3
"""Offline manufactured plane/sphere/rough interface temperature tests.
These are reconstruction sensitivity tests, not a validated VOF interface solver.
"""
from pathlib import Path
import csv, numpy as np
OUT=Path(__file__).resolve().parent
TAB=[]
with open(OUT/'inputs/Ag_evaporation_solver_table.csv',newline='') as f:
    for r in csv.DictReader(f):
        try: TAB.append((float(r['temperature_K']),float(r['p_sat_Pa']),float(r['alpha_ev']) if r['alpha_ev'] else 0.1428))
        except Exception: pass
Ttab=np.array([r[0] for r in TAB]); Ptab=np.array([r[1] for r in TAB]); Atab=np.array([r[2] for r in TAB])
Rspec=77.0798309247
N=36; h=1.; axis=(np.arange(N)+0.5); X,Y,Z=np.meshgrid(axis,axis,axis,indexing='ij')
points=np.stack((X,Y,Z),axis=-1); g_liq=-22.; g_gas=-8.

def geometry(name):
    if name=='plane':
        n0=np.array([0.3,-0.2,0.9327379]); n0=n0/np.linalg.norm(n0); c=np.array([18.,18.,18.])
        d=np.sum((points-c)*n0,axis=-1); n=np.broadcast_to(n0,(N,N,N,3)); ts=1400.+0.*X; surf=points-d[...,None]*n
        return d,n,ts,surf
    if name=='sphere':
        c=np.array([18.,18.,18.]); v=points-c; r=np.linalg.norm(v,axis=-1); n=v/np.maximum(r[...,None],1e-12); d=r-10.; ts=1400.+0.*X; surf=points-d[...,None]*n
        return d,n,ts,surf
    # Rough graph interface, liquid side below the surface. Signed distance is
    # the local normal projection of vertical offset; surface temperature has
    # a weak tangential modulation to expose fit contamination near roughness.
    x,y,z=points[...,0],points[...,1],points[...,2]
    zs=18.+1.4*np.sin(2*np.pi*x/17.)*np.cos(2*np.pi*y/19.)
    zx=1.4*(2*np.pi/17.)*np.cos(2*np.pi*x/17.)*np.cos(2*np.pi*y/19.)
    zy=-1.4*(2*np.pi/19.)*np.sin(2*np.pi*x/17.)*np.sin(2*np.pi*y/19.)
    q=np.sqrt(1+zx*zx+zy*zy); n=np.stack((-zx/q,-zy/q,1/q),axis=-1)
    d=(z-zs)/q
    ts=1400.+25.*np.sin(2*np.pi*x/23.)
    surf=points-d[...,None]*n
    return d,n,ts,surf

def temp_field(d,ts):
    return ts + np.where(d<0,g_liq*d,g_gas*d)
def flux(T):
    T=np.clip(T,Ttab[0],Ttab[-1]); p=np.interp(T,Ttab,Ptab); a=np.interp(T,Ttab,Atab)
    return a*p/np.sqrt(2*np.pi*Rspec*T)

def test(name):
    d,n,ts,surf=geometry(name); temp=temp_field(d,ts)
    grads=np.stack([ (np.roll(temp,-1,axis=ax)-np.roll(temp,1,axis=ax))/(2*h) for ax in range(3)],axis=-1)
    sel=(np.abs(d)<=0.55)&(np.abs(X-18)<14)&(np.abs(Y-18)<14)&(np.abs(Z-18)<14)
    out=[]; rad_area=[]
    # Deterministic local plane fit: selected liquid cells within normal distance
    # 3.2h and tangent distance 1.5h of each interface centroid.
    inds=np.argwhere(sel)
    for ii,jj,kk in inds:
        p=points[ii,jj,kk]; nn=n[ii,jj,kk]; ss=surf[ii,jj,kk]; exact=float(ts[ii,jj,kk])
        tc=float(temp[ii,jj,kk]); centered=tc+float(np.dot(grads[ii,jj,kk],ss-p))
        # small 9x9x9 search around point, matching a local one-sided normal fit
        i0=max(0,ii-4); i1=min(N,ii+5); j0=max(0,jj-4); j1=min(N,jj+5); k0=max(0,kk-4); k1=min(N,kk+5)
        cand=points[i0:i1,j0:j1,k0:k1].reshape(-1,3); dd=d[i0:i1,j0:j1,k0:k1].ravel(); tt=temp[i0:i1,j0:j1,k0:k1].ravel()
        vv=cand-ss; coord=vv@nn
        tangent=np.linalg.norm(vv-coord[:,None]*nn,axis=1)
        use=(coord<0)&(coord>=-3.2*h)&(tangent<=1.5*h)&(dd<0)
        if np.count_nonzero(use)<4: continue
        slope,intercept=np.polyfit(coord[use],tt[use],1)
        one=float(intercept)
        w=h*h*np.linalg.norm(nn)
        jcent=float(flux(centered)); jone=float(flux(one))
        out.append([name,ii,jj,kk,float(d[ii,jj,kk]),exact,centered,one,centered-exact,one-exact,int(np.count_nonzero(use)),w,jcent,jone])
    if not out: raise RuntimeError(f'no fitted interfaces: {name}')
    a=np.array(out,dtype=object); weights=np.asarray(a[:,11],float); ec=np.asarray(a[:,8],float); eo=np.asarray(a[:,9],float); jc=np.asarray(a[:,12],float); jo=np.asarray(a[:,13],float)
    total_area=float(weights.sum()); ratec=float(np.sum(weights*jc)); rateo=float(np.sum(weights*jo))
    summary=[name,len(out),np.average(np.abs(ec),weights=weights),np.max(np.abs(ec)),np.sqrt(np.average(ec**2,weights=weights)),np.average(np.abs(eo),weights=weights),np.max(np.abs(eo)),np.sqrt(np.average(eo**2,weights=weights)),total_area,ratec,rateo,100*(rateo/ratec-1.) if ratec else np.nan]
    return out,summary
allrows=[]; summaries=[]
for nm in ('plane','sphere','rough'):
    rr,ss=test(nm);allrows.extend(rr);summaries.append(ss)
fields=['geometry','i','j','k','signed_distance_cells','T_exact_surface_K','T_centered_K','T_liquid_sided_fit_K','centered_error_K','one_sided_error_K','fit_sample_count','interface_area_weight_m2','HK_flux_centered_kg_m2_s','HK_flux_one_sided_kg_m2_s']
with open(OUT/'p1p_manufactured_interface_temperature.csv','w',newline='') as f:
    w=csv.writer(f);w.writerow(fields);w.writerows(allrows)
with open(OUT/'p1p_manufactured_interface_summary.csv','w',newline='') as f:
    w=csv.writer(f);w.writerow(['geometry','interface_cells_fit','centered_MAE_K','centered_max_abs_K','centered_RMS_K','one_sided_MAE_K','one_sided_max_abs_K','one_sided_RMS_K','area_weight_sum_m2','total_HK_rate_centered_kg_s','total_HK_rate_one_sided_kg_s','one_sided_rate_change_pct']);w.writerows(summaries)
print('summary rows:')
for r in summaries: print(r)
