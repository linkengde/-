# P0 final decision

## A. Tested automatic-curvature VOF path

**Pass for the requested standalone static-drop verification scope.** The height-function curvature sphere and 72° sessile-drop tests passed at R/Δx = 16, 32 and 64 with pressure error below 0.51%, contact-angle error below 0.016°, and mass drift below 0.000281%. The R64 automatic-curvature sphere remained stable for the existing 91-step test to t=0.05. A matched R64 AMR short control also completed two pressure projections without warnings, with residual/tolerance ratios 0.2550 and 0.02694 and mass drift +1.21×10⁻⁶%. This does not establish that an integrated production heating/flow model or an Ag melt-pool calculation has been validated; neither was run here.

This is a dimensionless solver validation; it is not yet an Ag-specific or melt-pool flow validation.

## B. Exact-curvature auxiliary path

**Not passed at R64 on the wavelet-AMR topology.** A one-time wavelet mesh update to 681,472 cells causes failure even when the grid is then frozen. On the unmodified 732,264-cell geometrically refined mesh, the same exact-curvature case passes two projections. The same wavelet topology is stable with automatic curvature. This points to an interaction between prescribed curvature and the coarse/fine AMR force-pressure discretization, not a blanket VOF failure. Resetting the analytic Laplace pressure after wavelet adaptation leaves the first-projection failure unchanged, so simple transfer of the initial pressure jump is not the explanation. The local cell-level mismatch is not yet isolated, so no production or numerical correction was made.

The 512³ full-domain uniform R64 control was resource-limited before initialization at the explicit 30-GiB guard; it provides no pass/fail evidence. No cgroup OOM occurred.

## Gate outcome

P0 is **partially complete**: the automatic-curvature VOF path is validated for the tested static-drop scope, while the exact-curvature R64 auxiliary branch remains an open numerical issue. This task did not start P1, alter the production enthalpy solver, or run heating or melt-pool flow.

P1 should remain a separate controlled authorization step. If it is later approved, its baseline must use the tested automatic-curvature path and must not depend on the failing exact-curvature auxiliary branch. The limitation and the unresolved AMR diagnosis should accompany any P1 result.

The run-by-run evidence is in `R64_EXACT_CURVATURE_DIAGNOSIS.md`, `basilisk_results/r64_fixed_vs_amr.csv`, `basilisk_results/r64_pressure_projection.csv`, and `basilisk_results/r64_curvature_force_balance.csv`.
