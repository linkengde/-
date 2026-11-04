# P1 final decision

## Outcome

**P1 numerical verification: pass for the tested dimensionless benchmarks.** The automatic-curvature path produced the correct hot-to-cold flow for negative dσ/dT, reversed it for positive dσ/dT, matched the planar Stokes solution to 0.7805% at H/Δx=64, conserved VOF mass to better than 2×10⁻⁹%, and showed timestep sensitivity below 0.001%. On a curved static interface, the automatic-curvature normal force and projected Marangoni tangential force converged to their analytic resultants with a 0.0137% net-force residual at R/Δx=64.

The controlled numerical benchmark is sufficiently separated from the P0 fine-grid parasitic-flow level to support continued method development. The exact-curvature/wavelet-AMR failure remains isolated and unchanged.

## Not validated by P1

- Direct reproduction of Basilisk's official dynamic `src/test/marangoni.c` example: that source was absent from the installed mirror, and this run instead used independent analytical Stokes and static force-closure benchmarks.
- A solved temperature field or full variable-σ transient drop migration.
- Dimensional Ag flow speed, Ag viscosity coupling, a Ti₃SiC₂-constrained channel, or an Ag70 melt pool. Atmosphere-specific literature slopes were included only as normalized sign-sensitivity cases.
- Ag surface-tension extrapolation above the cited measurement ranges, including 2277 K.
- Production heat/flow coupling, free-surface deformation, or permanent crater formation.

## Next allowed research step

Keep development on an isolated branch. First acquire and archive the official Basilisk Marangoni benchmark source/version, then reproduce its reference case with automatic curvature. Separately, obtain a justified Ag viscosity and thermal-gradient scale before converting the atmosphere-specific surface-tension slopes into dimensional velocities. Do not start the Ag70 melt-pool flow calculation from this P1 result alone.
