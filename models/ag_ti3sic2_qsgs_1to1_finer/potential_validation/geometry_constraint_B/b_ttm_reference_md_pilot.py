import json, math, os, sys
from pathlib import Path
import numpy as np
sys.path.extend(['/workspace/.venvs/lammps/lib/python3.12/site-packages','/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer'])
from mpi4py import MPI
from lammps import lammps
from ag_ttm_reference import AgParameters, AgTTMGrid, build_ag_grid_fields, lattice_temperature_by_cell, rescale_lattice_kinetic, KB_EV_K, M_VV2E_METAL
root=Path('/workspace/-/models/ag_ti3sic2_qsgs_1to1_finer'); bdir=root/'potential_validation/geometry_constraint_B'; candidate=root/'potential_validation/relaxed_backbone_candidate'
os.chdir(root)
report=json.loads((candidate/'geometry_report.json').read_text()); lx,ly,lz=report['material_dimensions_A']; nsteps=int(os.environ.get('AG_TTM_PILOT_STEPS', '100')); dt_ps=0.0005
if nsteps <= 0:
    raise ValueError('AG_TTM_PILOT_STEPS must be positive')
lmp=lammps(cmdargs=['-screen','none','-log','none','-sf','omp'],comm=MPI.COMM_WORLD)
commands=['units metal','atom_style atomic','boundary p p f','atom_modify sort 0 0',f'read_data {candidate}/ag_ti3sic2_qsgs.data.gz',f'read_dump {bdir}/ag_minimized_100.dump.gz 100 x y z box no replace yes',f'include {bdir}/potential.geometry_wca_fixed_skeleton.inc','neighbor 2.0 bin','neigh_modify delay 0 every 1 check yes','region r_bottom block INF INF INF INF 0.0 20.0 units box','group bottom region r_bottom','group Ag type 4','group skeleton type 1 2 3','group mobile_Ag subtract Ag bottom','fix freeze_bottom bottom setforce 0.0 0.0 0.0','fix freeze_skeleton skeleton setforce 0.0 0.0 0.0','mass 4 107.8682','velocity mobile_Ag create 300.0 190427 mom yes rot no dist gaussian','timestep 0.0005','fix integrate mobile_Ag nve','thermo 1','thermo_style custom step pe ke etotal','thermo_modify flush yes','run 0']
for s in commands: lmp.command(s)
x=lmp.numpy.extract_atom('x'); v=lmp.numpy.extract_atom('v'); typ=lmp.numpy.extract_atom('type'); nlocal=int(lmp.extract_setting('nlocal'))
idx=np.flatnonzero((typ[:nlocal]==4)&(x[:nlocal,2]>=40.0)&(x[:nlocal,2]<=lz)); xyz=np.array(x[:nlocal][idx],copy=True); masses=np.full(len(idx),107.8682)
spacing_nom=10.0; lengths=np.array([lx,ly,lz]); shape=np.ceil(lengths/spacing_nom).astype(int); spacing=lengths/shape
r2=(xyz[:,0]-lx/2)**2+(xyz[:,1]-ly/2)**2; q=(5.0/20.0)*np.exp(-r2/(30.0**2)-(lz-xyz[:,2])/20.0)
fields=build_ag_grid_fields(xyz,q,tuple(shape),tuple(spacing),4.09**3/4.0)
params=AgParameters.from_csv(root/'ttm_data/Ag_TTM_hybrid_figure_fit_300_20000K.csv',allow_low_temperature_extension=True)
grid=AgTTMGrid(params,fields.fraction_ag,tuple(spacing),fields.source_power_eV_ps,initial_temperature_K=300.0)
e0=float(grid.electron_energy_eV.sum()); et0=float(lmp.get_thermo('etotal')); source_cum=0.0; transfer_cum=0.0
print(f'init atoms={lmp.get_natoms()} newton_ag={len(idx)} grid={tuple(int(i) for i in shape)} source_power={fields.source_power_eV_ps.sum():.9f} eV/ps',flush=True)
for step in range(1,nsteps+1):
    lmp.command('run 1 pre no post no')
    et_md=float(lmp.get_thermo('etotal'))
    vnow=np.array(v[:nlocal][idx],copy=True)
    tl=lattice_temperature_by_cell(vnow,masses,fields.atom_map,smoothing_radius_cells=1)
    acc=grid.advance(dt_ps,tl,heat_on=True)
    vnew,applied=rescale_lattice_kinetic(vnow,masses,fields.atom_map,acc.lattice_energy_per_cell_eV)
    v[:nlocal][idx]=vnew
    source_cum+=acc.source_to_electron_eV
    step_transfer=float(applied.sum()); transfer_cum+=step_transfer
    tl_after=lattice_temperature_by_cell(vnew,masses,fields.atom_map,smoothing_radius_cells=1)
    active_cells=fields.fraction_ag.ravel()>0.0
    hot_cells=fields.source_power_eV_ps.ravel()>=0.01*fields.source_power_eV_ps.max()
    tl_after_flat=tl_after.ravel()
    hot_lattice_mean=float(np.average(tl_after_flat[hot_cells],weights=fields.source_power_eV_ps.ravel()[hot_cells]))
    # et_md already contains all prior-step lattice transfers; add only the
    # current step's velocity handoff when checking the post-TTM energy.
    combined=et_md+step_transfer+float(grid.electron_energy_eV.sum())
    resid=combined-(et0+e0)-source_cum
    if step%20==0 or step==1 or step==nsteps:
        vcur=np.array(v[:nlocal][idx],copy=True)
        ke=0.5*M_VV2E_METAL*np.sum(masses[:,None]*vcur**2)
        temp=2*ke/(3*len(idx)*KB_EV_K)
        active=fields.fraction_ag>0
        print(f'step={step} time_ps={step*dt_ps:.5f} T_newtonAg={temp:.3f}K Tl_source_weighted={hot_lattice_mean:.3f}K Tl_active_max={np.nanmax(tl_after_flat[active_cells]):.3f}K Te_minmax={grid.temperature_K[active].min():.3f},{grid.temperature_K.max():.3f}K source={source_cum:.6f}eV e2l={transfer_cum:.6f}eV substeps_last={acc.substeps} ledger_resid={resid:.6e}eV',flush=True)
lmp.command('run 0 pre no post no'); etfinal=float(lmp.get_thermo('etotal')); efinal=float(grid.electron_energy_eV.sum())
print(f'final_after_run0 cumulative_residual={(etfinal+efinal)-(et0+e0)-source_cum:.6e} eV; total_electron_to_lattice={transfer_cum:.6f} eV',flush=True)
lmp.close(); MPI.Finalize()
