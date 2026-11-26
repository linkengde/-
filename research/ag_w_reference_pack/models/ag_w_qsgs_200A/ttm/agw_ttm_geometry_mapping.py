"""Geometry-consistent Ag/W mapping from supplied phase-volume grids.

Input Ag/W phase fractions define material occupancy. CIC is used only for
atom kinetic-energy handoff and Ag source deposition.
"""
from __future__ import annotations
import numpy as np
from agw_ttm_reference import (
    AG_ATOMIC_VOLUME_A3, AgWMappingError, AtomGridMap,
    cloud_in_cell_atom_map, gaussian_ag_source_power,
)

def _phase_map(raw: AtomGridMap, fraction: np.ndarray, phase: str):
    ids=np.asarray(raw.cell_ids,dtype=np.int64)
    weights=np.asarray(raw.weights,dtype=float).copy()
    flat=np.asarray(fraction,dtype=float).ravel()
    if ids.shape!=weights.shape or ids.ndim!=2 or np.any(ids<0) or np.any(ids>=flat.size):
        raise AgWMappingError(f"{phase} CIC map is malformed")
    weights *= flat[ids] > 1.0e-12
    rows=weights.sum(axis=1)
    bad=np.flatnonzero(rows<=1.0e-14)
    if len(bad):
        raise AgWMappingError(f"{len(bad)} {phase} atoms have no CIC support in their own phase cells")
    weights /= rows[:,None]
    mapped=AtomGridMap(ids,weights,raw.grid_shape)
    support=np.bincount(ids.ravel(),weights=weights.ravel(),minlength=flat.size)
    active=flat>1.0e-12
    missing=active&(support<=1.0e-14)
    return mapped,{
        "atoms_without_phase_cell_support":int(len(bad)),
        "active_phase_cells":int(active.sum()),
        "active_phase_cells_without_atom_support":int(missing.sum()),
        "max_atom_row_sum_error":float(np.max(np.abs(weights.sum(axis=1)-1.0))),
        "minimum_active_cell_support_weight":float(support[active].min()) if np.any(active) else 0.0,
    }

def _map_power(mapping: AtomGridMap, per_atom_power: np.ndarray) -> np.ndarray:
    power=np.asarray(per_atom_power,dtype=float)
    if power.shape!=(mapping.cell_ids.shape[0],) or np.any(~np.isfinite(power)) or np.any(power<0):
        raise ValueError("per-atom power must be finite, non-negative, and match the map")
    return np.bincount(mapping.cell_ids.ravel(),
        weights=(power[:,None]*mapping.weights).ravel(),
        minlength=int(np.prod(mapping.grid_shape))).reshape(mapping.grid_shape)

def map_agw_newton_layer_geometry(
    xyz_absolute_A, types, fraction_ag, fraction_w, grid_shape, spacing_A,
    source_kwargs=None, newton_mask=None, grid_offset_xy_A=(0.0,0.0), box_xy_A=200.0,
):
    """Map atoms to supplied Ag/W geometry fractions and phase-aware CIC."""
    xyz=np.asarray(xyz_absolute_A,dtype=float); typ=np.asarray(types,dtype=int)
    shape=tuple(int(v) for v in grid_shape); spacing=tuple(float(v) for v in spacing_A)
    fa=np.asarray(fraction_ag,dtype=float); fw=np.asarray(fraction_w,dtype=float)
    if xyz.shape!=(len(typ),3) or np.any(~np.isfinite(xyz)) or np.any(~np.isin(typ,(1,2))):
        raise ValueError("Expected finite coordinates and type 1=W/type 2=Ag")
    if fa.shape!=shape or fw.shape!=shape or len(shape)!=3 or len(spacing)!=3:
        raise ValueError("Geometry fractions must match the requested 3D grid")
    if (np.any(~np.isfinite(fa)) or np.any(~np.isfinite(fw)) or np.any(fa<0) or np.any(fw<0)
        or np.any(fa>1) or np.any(fw>1) or np.any(fa+fw>1+1e-8)):
        raise AgWMappingError("Geometry fractions exceed finite-volume cell capacity")
    offset=np.asarray(grid_offset_xy_A,dtype=float)
    if offset.shape!=(2,) or np.any(~np.isfinite(offset)) or np.any(np.abs(offset)>1e-10):
        raise AgWMappingError("V11 geometry fractions require x/y grid origin (0,0)")
    candidate=np.ones(len(xyz),dtype=bool) if newton_mask is None else np.asarray(newton_mask,dtype=bool)
    if candidate.shape!=(len(xyz),): raise ValueError("newton_mask must match atom count")
    inside=(xyz[:,2]>=40-1e-9)&(xyz[:,2]<=200+1e-9)
    idx_ag=np.flatnonzero(candidate&inside&(typ==2)); idx_w=np.flatnonzero(candidate&inside&(typ==1))
    if not len(idx_ag) or not len(idx_w): raise AgWMappingError("Both phases must contribute Newton atoms")
    origin=np.array([0.0,0.0,40.0])
    raw_ag=cloud_in_cell_atom_map(xyz[idx_ag]-origin,shape,spacing)
    raw_w=cloud_in_cell_atom_map(xyz[idx_w]-origin,shape,spacing)
    map_ag,diag_ag=_phase_map(raw_ag,fa,"Ag")
    map_w,diag_w=_phase_map(raw_w,fw,"W")
    kwargs=dict(source_kwargs or {})
    kwargs.setdefault("center_xy_A",(box_xy_A/2,box_xy_A/2))
    kwargs.setdefault("periodic_box_xy_A",box_xy_A)
    q=gaussian_ag_source_power(xyz[idx_ag],**kwargs)
    atom_power=q*AG_ATOMIC_VOLUME_A3
    source=_map_power(map_ag,atom_power)
    expected=float(atom_power.sum()); error=float(source.sum()-expected)
    if abs(error)>max(1e-9,abs(expected)*1e-12):
        raise AgWMappingError("Phase-aware CIC did not conserve Ag source power")
    if np.any(source[fa<=1e-12]>1e-12):
        raise AgWMappingError("Ag source leaked into cells without Ag volume")
    cellvol=float(np.prod(spacing)); total=fa+fw
    nx,ny,nz=shape
    cx=(np.arange(nx)+.5)*spacing[0]; cy=(np.arange(ny)+.5)*spacing[1]
    cz=40+(np.arange(nz)+.5)*spacing[2]
    x,y,z=np.meshgrid(cx,cy,cz,indexing="ij")
    r2=(x-box_xy_A/2)**2+(y-box_xy_A/2)**2
    depth=np.maximum(0,200-z); ps=float(source.sum())
    diag={
      "mapping_method":"input phase-volume fractions + phase-restricted CIC",
      "atomic_volume_capacity_flow_used":False,
      "max_total_phase_fill":float(total.max()),
      "mean_void_fraction":float(np.mean(1-total)),
      "Ag_geometry_volume_A3":float(fa.sum()*cellvol),
      "W_geometry_volume_A3":float(fw.sum()*cellvol),
      "Ag_phase_cic":diag_ag,"W_phase_cic":diag_w,
      "Ag_source_power_eV_ps":ps,"expected_Ag_source_power_eV_ps":expected,
      "source_power_conservation_error_eV_ps":error,
      "source_radial_RMS_A":float(np.sqrt(np.sum(source*r2)/ps)) if ps else 0.0,
      "source_mean_depth_A":float(np.sum(source*depth)/ps) if ps else 0.0,
      "Ag_atom_map_max_row_sum_error":diag_ag["max_atom_row_sum_error"],
      "W_atom_map_max_row_sum_error":diag_w["max_atom_row_sum_error"],
    }
    return {
      "fraction_ag":fa.copy(),"fraction_w":fw.copy(),"source_power_ag_eV_ps":source,
      "atom_map_ag":map_ag,"atom_map_w":map_w,"atom_indices_ag":idx_ag,"atom_indices_w":idx_w,
      "candidate_newton_atoms":int(candidate.sum()),"grid_offset_xy_A":offset.tolist(),
      "box_xy_A":float(box_xy_A),
      "excluded_outside_ttm_region":{
        "below_40A":int(np.count_nonzero(candidate&(xyz[:,2]<40-1e-9))),
        "above_200A":int(np.count_nonzero(candidate&(xyz[:,2]>200+1e-9))),
        "below_40A_Ag":int(np.count_nonzero(candidate&(typ==2)&(xyz[:,2]<40-1e-9))),
        "below_40A_W":int(np.count_nonzero(candidate&(typ==1)&(xyz[:,2]<40-1e-9))),
        "above_200A_Ag":int(np.count_nonzero(candidate&(typ==2)&(xyz[:,2]>200+1e-9))),
        "above_200A_W":int(np.count_nonzero(candidate&(typ==1)&(xyz[:,2]>200+1e-9))),
      },
      "mapped_ag_volume_A3":float(fa.sum()*cellvol),
      "mapped_w_volume_A3":float(fw.sum()*cellvol),
      "max_total_fill_fraction":float(total.max()),
      "source_power_on_state_eV_ps":ps,"mapping_diagnostics":diag,
    }
