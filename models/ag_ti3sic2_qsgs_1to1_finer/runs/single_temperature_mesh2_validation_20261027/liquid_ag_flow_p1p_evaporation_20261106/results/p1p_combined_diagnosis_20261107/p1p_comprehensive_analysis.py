#!/usr/bin/env python3
"""Reproducible post-processing for archived Ag70 P1P checkpoint and histories."""
from pathlib import Path
import csv, struct, hashlib, math
import numpy as np

BASE=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
RUN=BASE/'results/p1p_real_smoke_4mw_3x_3ps_20261106'
REPLAY=Path('/tmp/p1p-combined/replay')
DX=3.01857484616085e-10
RHO_AG=10490.0
AG_ATOM=107.8682*1.66053906660e-27
FIELDS=['sample_vof','drag','ag_fraction','skeleton_fraction','liquid_fraction','temperature_K','enthalpy_density_J_m3','u_x_m_s','u_y_m_s','u_z_m_s','pressure_Pa','curvature_m_inv','surface_tension_N_m','sigma_force','marangoni_x','marangoni_y','marangoni_z','source_W_m3','drag_s_inv','material_code','pore_fraction']

def rows(path):
    with open(path,newline='') as f: return list(csv.DictReader(f))
def write_csv(path, header, data):
    with open(path,'w',newline='') as f:
        w=csv.writer(f); w.writerow(header); w.writerows(data)
def read_snap(path, expected_magic=b'P1HSNAP1'):
    b=Path(path).read_bytes(); fmt='@8s6I5d'; hsize=struct.calcsize(fmt)
    magic,ver,nx,ny,nz,nf,flags,t,dx,x0,y0,z0=struct.unpack(fmt,b[:hsize])
    if magic!=expected_magic: raise ValueError((magic,expected_magic))
    a=np.frombuffer(b,dtype=np.float32,offset=hsize).copy().reshape(nf,nx,ny,nz)
    return dict(magic=magic.decode(),version=ver,nx=nx,ny=ny,nz=nz,nfields=nf,flags=flags,time_ps=t,dx=dx,x0=x0,y0=y0,z0=z0),a

hist=rows(RUN/'p1h_history.csv')
flux=rows(RUN/'p1p_evaporation_flux.csv')
trans=rows(RUN/'p1j_transport_validation.csv')
arc=rows(RUN/'p1h_arc_pressure_history.csv')
# Preserve every history row in a transparent energy ledger.
energy=[]
for r in hist:
    E=float(r['enthalpy_J']); Ein=float(r['input_energy_J']); Es=float(r['bottom_sink_energy_J']); Er=float(r['energy_residual_J'])
    energy.append([r['time_ps'],r['source_power_W'],r['input_minus_bottom_mW'],E,Ein,Es,Ein-Es,E-float(hist[0]['enthalpy_J']),Er,
                   r['Tmax_K'],r['Tmax_skeleton_K'],r['liquid_Ag_volume_nm3'],r['pressure_residual'],r['thermal_residual'],r['dt_fs']])
write_csv(OUT/'P1P_ENERGY_LEDGER.csv',['time_ps','source_power_W','input_minus_bottom_mW','stored_enthalpy_J','cumulative_input_J','cumulative_bottom_sink_J','net_input_minus_sink_J','delta_enthalpy_J','energy_balance_residual_J','Tmax_K','Tmax_skeleton_K','liquid_Ag_volume_nm3','pressure_residual','thermal_residual','dt_fs'],energy)
# Mass ledger using the solver's reference density and logged evaporated mass.
init_vol=float(hist[0]['Ag_volume_nm3']); init_mass=init_vol*1e-27*RHO_AG
mass=[]
for r in hist:
    t=float(r['time_ps']); frow=min(flux,key=lambda q:abs(float(q['time_ps'])-t))
    xrow=min(trans,key=lambda q:abs(float(q['time_ps'])-t))
    V=float(r['Ag_volume_nm3']); M=V*1e-27*RHO_AG
    escaped=float(frow['escaped_mass_cumulative_kg'])
    residual=100*(M+escaped-init_mass)/init_mass
    mass.append([r['time_ps'],V,M,float(r['Ag_volume_drift_pct']),escaped,escaped/AG_ATOM,M+escaped,residual,
                 r['liquid_Ag_volume_nm3'],xrow['Ag_in_near_gas_cells'],xrow['partition_error_max']])
write_csv(OUT/'P1P_MASS_LEDGER.csv',['time_ps','Ag_volume_nm3','Ag_component_mass_kg','solver_logged_volume_drift_pct','evaporated_mass_kg','evaporated_Ag_atom_mass_equivalents','Ag_plus_evaporated_mass_kg','recomputed_mass_residual_pct','liquid_Ag_volume_nm3','Ag_in_near_gas_cells','component_partition_error_max'],mass)
# Save selected checkpoint fields in a compressed, portable archive.
sp=REPLAY/'p1h_marangoni_on_restart_loaded_t003.012805_ps.bin'
npz_path=OUT/'p1p_checkpoint_fields_t3p012805.npz'
if sp.exists():
    header,A=read_snap(sp)
    if header['nfields']!=len(FIELDS): raise ValueError('unexpected field count')
    np.savez_compressed(npz_path,**{name:A[i] for i,name in enumerate(FIELDS)},header=np.array([header['nx'],header['ny'],header['nz'],header['time_ps'],header['dx'],header['x0'],header['y0'],header['z0']]))
else:
    saved=np.load(npz_path); header={'nx':int(saved['header'][0]),'ny':int(saved['header'][1]),'nz':int(saved['header'][2]),'time_ps':float(saved['header'][3]),'dx':float(saved['header'][4]),'x0':float(saved['header'][5]),'y0':float(saved['header'][6]),'z0':float(saved['header'][7])}; A=np.stack([saved[name] for name in FIELDS])
# Point-in-time geometry and velocity/mass distribution.
nx,ny,nz=header['nx'],header['ny'],header['nz']; dx=header['dx']; dv=dx**3
x=(header['x0']+(np.arange(nx)+0.5)*dx)*1e9
y=(header['y0']+(np.arange(ny)+0.5)*dx)*1e9
z=(header['z0']+(np.arange(nz)+0.5)*dx)*1e9
X,Y=np.meshgrid(x,y,indexing='ij')
cx,cy=9.58832,9.65943951
L=nx*dx*1e9
dxw=(X-cx+L/2)%L-L/2; dyw=(Y-cy+L/2)%L-L/2
R=np.sqrt(dxw[:,:,None]**2+dyw[:,:,None]**2)
ag=A[2]; sc=A[3]; liq=A[4]; temp=A[5]; uxyz=A[7:10]; f=A[0]; mat=A[19].round().astype(np.uint8)
R3=np.broadcast_to(R,ag.shape)
ranges=[(0,3),(3,5),(5,8),(8,100)]
rad=[]
for lo,hi in ranges:
    m=(R3>=lo)&(R3<hi)
    rad.append([lo,hi,float(np.sum(ag[m])*dv*RHO_AG),float(np.sum(ag[m]*liq[m])*dv*RHO_AG),float(np.sum(ag[m]*sc[m])*dv*RHO_AG),float(np.sum(ag[m]*(mat[m]==2))*dv*RHO_AG),float(np.sum(ag[m]*(mat[m]==0))*dv*RHO_AG)])
write_csv(OUT/'P1P_CHECKPOINT_RADIAL_MASS.csv',['radius_min_nm','radius_max_nm','Ag_component_mass_kg','liquid_Ag_component_mass_kg','Ag_in_current_skeleton_marker_kg','Ag_in_original_Ti3SiC2_cells_kg','Ag_in_original_gas_cells_kg'],rad)
# Compare to the actual initial snapshot from the A/B/C common 300 K initialization, not a raw material-code count.
initial_path=OUT/'abc_short/A/p1h_marangoni_on_initial_t000.000000_ps.bin'
if initial_path.exists():
    ih,IA=read_snap(initial_path)
    ag_initial=IA[2]
    change=[]
    for lo,hi in ranges:
        m=(R3>=lo)&(R3<hi)
        m0=m[:,:,:ag_initial.shape[2]]
        mi=float(np.sum(ag_initial[m0])*dv*RHO_AG)
        mf=float(np.sum(ag[m])*dv*RHO_AG)
        ml=float(np.sum((ag*liq)[m])*dv*RHO_AG)
        change.append([lo,hi,mi,mf,mf-mi,100.*(mf-mi)/mi if mi else float('nan'),ml])
    write_csv(OUT/'P1P_CHECKPOINT_RADIAL_MASS_CHANGE.csv',['radius_min_nm','radius_max_nm','initial_Ag_mass_kg','checkpoint_Ag_mass_kg','delta_Ag_mass_kg','relative_change_pct','checkpoint_liquid_Ag_mass_kg'],change)
# Checkpoint summary and area-weighted velocity, center surface displacement from logged field.
liquid_weight=ag*liq
liquid_vol=float(np.sum(liquid_weight)*dv/1e-27)
liquid_speed=np.sqrt(np.sum(uxyz**2,axis=0))
wt=float(np.sum(liquid_weight))
mean_speed=float(np.sum(liquid_weight*liquid_speed)/wt) if wt>0 else float('nan')
radial_u=np.zeros_like(liquid_speed)
rr=np.broadcast_to(R,(nx,ny,nz))
radial_u=(uxyz[0]*dxw[:,:,None]+uxyz[1]*dyw[:,:,None])/(rr+1e-30)
mean_ur=float(np.sum(liquid_weight*radial_u)/wt) if wt>0 else float('nan')
# Composite surface height from f>0.5, using the exact initial snapshot and
# the solver's top-face convention (cell-center z + dx/2).
h=np.full((nx,ny),np.nan); h0map=np.full((nx,ny),np.nan)
for i in range(nx):
  for j in range(ny):
    ids=np.where(f[i,j,:]>0.5)[0]
    if len(ids): h[i,j]=z[ids[-1]]+0.5*dx*1e9
initial_path=OUT/'abc_short/A/p1h_marangoni_on_initial_t000.000000_ps.bin'
if initial_path.exists():
  ih,IA=read_snap(initial_path); fi=IA[0]
  for i in range(nx):
    for j in range(ny):
      ids=np.where(fi[i,j,:]>0.5)[0]
      if len(ids): h0map[i,j]=(ih['z0']+(ids[-1]+1)*ih['dx'])*1e9
finite=np.isfinite(h)&np.isfinite(h0map); center=(R[:,:,0]<3.); source=(R[:,:,0]<=5.)
dh=h-h0map; depression=-dh
mean_center_h=float(np.nanmean(dh[center&finite])) if np.any(center&finite) else float('nan')
mean_source_depression=float(np.nanmean(depression[source&finite])) if np.any(source&finite) else float('nan')
mean_outer_h=float(np.nanmean(dh[(R[:,:,0]>=5)&finite]))
write_csv(OUT/'p1p_checkpoint_surface_height.csv',['i','j','x_nm','y_nm','initial_height_nm','final_height_nm','delta_height_nm','depression_positive_down_nm'],[[i,j,x[i],y[j],h0map[i,j],h[i,j],dh[i,j],depression[i,j]] for i in range(nx) for j in range(ny) if finite[i,j]])
# Cluster counts using existing solver diagnostics, preserve all points.
peaks=rows(RUN/'p1h_peak_pool_dimensions.csv')
peak_rows=[[r['time_ps'],r['Tmax_skeleton_K'],r['components'],r['largest_cells'],r['largest_liquid_volume_nm3'],r.get('x_size_nm',''),r.get('y_size_nm',''),r.get('z_size_nm','')] for r in peaks]
write_csv(OUT/'P1P_POOL_CONNECTIVITY_HISTORY.csv',['time_ps','Tmax_skeleton_K','component_count_6_face_periodic_xy','largest_cells','largest_liquid_cluster_nm3','largest_x_size_nm','largest_y_size_nm','largest_z_size_nm'],peak_rows)
# Summary table for report.
last=hist[-1]; first=hist[0]; lastflux=flux[-1]
summary={
 'run_time_ps':float(last['time_ps']),'max_T_K':max(float(r['Tmax_K']) for r in hist),'max_skeleton_K':max(float(r['Tmax_skeleton_K']) for r in hist),
 'max_liquid_Ag_nm3':max(float(r['liquid_Ag_volume_nm3']) for r in hist),'liquid_Ag_at_checkpoint_nm3':liquid_vol,
 'max_ag_volume_drift_pct':max(abs(float(r['Ag_volume_drift_pct'])) for r in hist),'max_energy_residual_J':max(abs(float(r['energy_residual_J'])) for r in hist),
 'initial_enthalpy_J':float(first['enthalpy_J']),'final_enthalpy_J':float(last['enthalpy_J']),'input_J':float(last['input_energy_J']),'sink_J':float(last['bottom_sink_energy_J']),'final_residual_J':float(last['energy_residual_J']),
 'evap_kg':float(lastflux['escaped_mass_cumulative_kg']),'evap_atom_equivalents':float(lastflux['escaped_mass_cumulative_kg'])/AG_ATOM,
 'mean_liquid_speed_m_s':mean_speed,'mean_radial_liquid_velocity_m_s':mean_ur,'checkpoint_max_velocity_m_s':float(last['max_velocity_m_s']),
 'checkpoint_mass_kg':init_mass,'center_mean_surface_height_change_nm':mean_center_h,'source_radius_le_5nm_mean_depression_nm':mean_source_depression,'outer_mean_surface_height_change_nm':mean_outer_h,
 'snapshot_fields':len(FIELDS),'snapshot_dx_nm':dx*1e9,
}
with open(OUT/'P1P_CHECKPOINT_SUMMARY.csv','w',newline='') as f:
 w=csv.writer(f); w.writerow(['metric','value']); w.writerows(summary.items())
print('HEADER',header)
for k,v in summary.items(): print(k,v)
print('mass initial volume',init_vol,'checkpoint sum vol_nm3',float(np.sum(ag)*dv/1e-27),'liquid vol from field',liquid_vol)
