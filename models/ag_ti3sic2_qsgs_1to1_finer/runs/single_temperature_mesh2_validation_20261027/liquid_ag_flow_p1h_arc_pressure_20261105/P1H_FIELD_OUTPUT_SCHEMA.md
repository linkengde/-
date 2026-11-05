# P1H field snapshot schema

P1H writes the inherited thermal/VOF frame as `p1h_*.bin` and an additional
pressure frame as `p1h_pressure_*.bin`. Both are atomic little-endian binary
files. The common header is packed as `<8s6I5d` (72 bytes): magic, version,
`nx, ny, nz, nfields, case_flag`, then `time_ps, dx_m, x0_m, y0_m, z0_m`.
Float32 arrays follow in field-major C order `(field, x, y, z)`.

The main frame has magic `P1HSNAP1`, version 1, and 20 fields:

| Index | Name | Unit | Meaning |
| --- | --- | --- | --- |
| 0 | `vof` | 1 | Sample versus gas VOF fraction (`f=1` sample) |
| 1 | `level_set_m` | m | CLSVOF signed-distance field |
| 2 | `ag_marker` | 1 | Ag occupancy fraction |
| 3 | `skeleton_marker` | 1 | Ti3SiC2 occupancy fraction |
| 4 | `liquid_fraction` | 1 | Ag liquid fraction from the enthalpy model |
| 5 | `temperature_K` | K | Mixture temperature |
| 6 | `enthalpy_density_J_m3` | J/m^3 | Stored volumetric enthalpy |
| 7 | `velocity_x_m_s` | m/s | x velocity |
| 8 | `velocity_y_m_s` | m/s | y velocity |
| 9 | `velocity_z_m_s` | m/s | z velocity |
| 10 | `pressure_Pa` | Pa | Hydrodynamic pressure |
| 11 | `curvature_m-1` | 1/m | Interface curvature estimate |
| 12 | `surface_tension_candidate_N_m` | N/m | Current model's Ag surface-tension candidate at local temperature |
| 13 | `vof_weighted_surface_tension_N_m` | N/m | Ag/liquid-weighted surface-tension field used by the interface model |
| 14 | `marangoni_acceleration_candidate_x_m_s2` | m/s^2 | x component of projected Marangoni acceleration |
| 15 | `marangoni_acceleration_candidate_y_m_s2` | m/s^2 | y component of projected Marangoni acceleration |
| 16 | `marangoni_acceleration_candidate_z_m_s2` | m/s^2 | z component of projected Marangoni acceleration |
| 17 | `source_power_density_W_m3` | W/m^3 | Normalized volumetric heating source before the synchronized off-ramp scale |
| 18 | `skeleton_and_mushy_drag_s-1` | 1/s | Local Brinkman/mushy damping rate |
| 19 | `material_id` | category | 0 gas, 1 Ag, 2 Ti3SiC2 (stored as float32) |

The pressure sidecar has magic `P1HPRS01`, version 1, and four fields:

| Index | Name | Unit | Meaning |
| --- | --- | --- | --- |
| 0 | `arc_pressure_Pa` | Pa | Configured spatial pressure profile before surface masking |
| 1 | `arc_pressure_mask` | 1 | Ag surface coverage × liquid fraction × upper-facing normal weight |
| 2 | `interface_normal_pressure_Pa` | Pa | Pressure potential after the interface mask |
| 3 | `total_normal_stress_potential_Pa` | Pa | Effective capillary potential plus masked arc pressure |

The pressure history CSV records time, shared heat-source ramp scale, configured
pressure peak, maximum effective interface pressure, integrated z force in N,
enabled state, and profile name. Positive pressure on the tracked upper surface
should have a negative z force under the `f=1` sample / `f=0` gas convention.

The field frame is fixed-grid data sampled from the octree calculation at the
64^3 P1G/P1H reference lattice; it is a visualization/diagnostic handoff, not a
lossless representation of every adaptive cell in the solver.
Consequently, material volumes should be taken from the solver's conservative
time-series CSV, not recomputed by summing projected HDF5 fractions. Undefined
curvature samples use the finite Basilisk `nodata` sentinel (`1e30` in this
snapshot); consumers should mask non-interface cells before interpreting
curvature values.
