# R/Δx=64 dynamic attempts

## Completed VOF-curvature tests

With the pinned Basilisk snapshot and four Open MPI ranks, both dynamic R/Δx=64 cases completed 91 steps to t=0.05 without pressure-convergence or VOF-CFL warnings:

- VOF-curvature sphere: pressure error +0.0309716%, `Umax=8.62060e-5`, mass drift +3.73149e-6%.
- 72° sessile drop: pressure error +0.0202520%, `Umax=2.82242e-4`, mass drift −1.73489e-6%, fitted angle 71.998499°.

The full per-case stdout/stderr logs are `sphere_vof_R64_dt0p5.*.log` and `sessile72_vof_R64_dt0p5.*.log` under `basilisk_cases/`.

A VOF-curvature sphere at R/Δx=32 was also run serially and with two MPI ranks. Pressure errors were +0.125732645% and +0.125734802%, respectively; maximum speeds differed by 0.8%. This supports using MPI for the finer-grid cases. The comparison logs are `mpi_R32_sphere_vof.*.log`.

## Analytical-curvature failure

The exact-curvature sphere remains unstable at R/Δx=64:

1. The original serial full run reached 732,264 leaf cells, then reported pressure non-convergence and VOF-CFL violations before any valid `RESULT`.
2. An MPI full run (four ranks, factor 0.5) reproduced the pressure failure (residual about 8.8×10¹²), VOF-CFL violation, and floating-point exception.
3. An MPI short run with time-step factor 0.25 still diverged.
4. Replacing the initial pressure step `p = 2` for `f > 0.5` with `p = 2f` did not fix it.
5. A diagnostic run with viscosity μ=10 also diverged; this was not treated as a material case.

These probes produce nonphysical velocities and pressure values. They are preserved for diagnosis only and must not be interpreted as validation results. Their logs are `mpi_sphere_exact_R64_*.*.log` under `basilisk_cases/diagnostics/`.

The VOF-curvature sphere and contact-angle cases succeed at this resolution, so the remaining issue is specific to the prescribed-curvature forcing / pressure-projection combination, not a general inability to run the R/Δx=64 grid. Overall P0 remains partial because the requested exact-curvature three-grid comparison is incomplete.

## Targeted capillary time-step and projection diagnostics

Follow-up diagnostic runs at capillary timestep factors 0.50, 0.35, and 0.20 all respected the Basilisk capillary bound
`dt_cap = sqrt(rho_m*dx^3/(pi*sigma)) = 0.00110193278037`
for `rho_m=sigma=1`, `dx/R=1/64`. Each first pressure projection converged, while each second projection reached 100 iterations and failed. The post-projection residual/tolerance ratios were 1.7736e10, 8.2990e9, and 3.3414e9; maximum speeds reached 1.45e23, 9.94e22, and 6.30e22. Before failure, `f` remained in [0,1], no nonfinite cells were counted, and mass drift remained below 5e-6 percent. The huge velocities then caused VOF CFL warnings and, in two probes, an MPI floating-point exception. Factor 0.50 returned process code zero despite an invalid state, so it is still classified as failed.

Code audit found one `iforce.h` and no duplicate `tension.h` force. The evidence points to a discrete force/pressure balance issue specific to imposing constant analytical curvature on the adaptive VOF interface. This is a high-confidence failure category, but the triggering discrete detail is not isolated enough to claim a fully proven root cause. The numerical data and raw logs are in `basilisk_results/vof_r64_exact_timestep_diagnostics.csv` and the `exact_R64_factor*.log` files. See `P0_FINAL_VALIDATION_REPORT.md` for the consolidated final decision.
