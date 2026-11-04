# P1b dynamic Marangoni droplet validation

## Decision

The pinned Basilisk CLSVOF reference run passes the project's preliminary dynamic-benchmark checks. The separate project automatic-curvature standard-VOF comparison does **not** reproduce the theoretical migration speed within tolerance and does not converge with the tested grid sequence. Therefore this validates the reference case and exposes a limitation in the current project VOF route; it does not clear that route for Ag melt-pool predictions.

## Pinned reference and run

The official `src/test/marangoni.c` was run from the user-provided Linux Basilisk snapshot `comphy-lab/basilisk-C v2026-07-20`, commit `9bb4f30aaf195feaa7a3a613468a1c888f5a4b9e`. The archive SHA256 is `4686294e8f649a424ca437a7fa94bfac9984c32bea5e9f33a27a5783cc692829`; the exact source and selected dependency hashes are recorded in `results/p1b_official/p1b_source_manifest.json`. This is a downstream-patched snapshot based on Basilisk/Darcs, not an untouched upstream release.

A harness-only stop event was appended so the recurring adaptation/logging events terminate after the source's `t=3*t0` measurement. The source body and physical settings were otherwise unchanged. Reproduction steps are in `results/p1b_official/BUILD_AND_RUN.md`; raw stdout, final velocity output, source, patch and the generated reference image are archived in `results/p1b_official/`.

The case is dimensionless and axisymmetric, with `R=1`, `rho=mu=1` in both phases, `Re=0.066`, source setting `Ca=0.66`, `Gamma_T=0.066`, `Gamma_0=0.1`, and theoretical `U_drop=-0.0088`. The source comments that `Ca=0.66` is ten times the value in Al Saud et al. (2018), chosen to accelerate computation, and that results are also good at `Ca=0.066`; the latter was not run here. The imposed linear surface-tension field is not a solved temperature field.

## Official CLSVOF results

| Cells per radius | Domain N | Measured `u_drop/U` | Theory `U_drop/U` | Relative speed error | Max logged mass drift |
|---:|---:|---:|---:|---:|---:|
| 8 | 128 | 0.128711 | 0.133333 | 3.467% | 0.05941% |
| 16 | 256 | 0.131629 | 0.133333 | 1.278% | 0.02714% |
| 32 | 512 | 0.132880 | 0.133333 | 0.340% | 0.01262% |

The speed error decreases at each refinement. At 32 cells per radius the error is below the provisional 5% goal; mass drift is below 0.1% at all three resolutions. This is a pass for the reference setup and these tested quantities, not a general validation of thermal Marangoni physics.

## Project automatic-curvature VOF comparison

A separate isolated adaptation of the same dimensionless case used `axi.h`, centered Navier–Stokes, standard `two-phase.h`/`tension.h`, automatic curvature, and a CSF representation of the spatially varying normal and tangential surface stresses. It did not use the reference CLSVOF plus `integral.h` stress discretization. No Ag geometry or Ag material parameters were used.

| Cells per radius | Domain N | Measured `u_drop/U` | Theory `U_drop/U` | Relative speed error | Max logged mass drift |
|---:|---:|---:|---:|---:|---:|
| 8 | 128 | 0.127917 | 0.133333 | 4.062% | 0.01918% |
| 16 | 256 | 0.120751 | 0.133333 | 9.437% | 0.00423% |
| 32 | 512 | 0.120003 | 0.133333 | 9.998% | 0.00423% |

Mass conservation is good, but the velocity error increases and then remains near 10%; the sequence does not meet the speed-error or grid-convergence checks. A preliminary tangential-only implementation produced the wrong migration direction because it omitted the variable normal `sigma*kappa*n` contribution. That incomplete run was stopped and is kept only in `results/p1b_autocurvature/diagnostics_tangential_only_incomplete/`; it is not included in the results table. The full-stress implementation still fails quantitative agreement, so the discrepancy must be resolved before using this project route for dynamic interface prediction.

The combined grid table and mass data are `p1b_velocity_convergence.csv` and `p1b_mass_conservation.csv`. Full logs and source for the second method are under `results/p1b_autocurvature/`.

## Scope

No thermal conduction equation, Ag-specific dimensional property, Ti₃SiC₂ skeleton, nanoscale pore, Ag70 geometry, or melt-pool flow was simulated. These results do not establish Ag velocities, melt-pool deformation, or crater formation. The existing P0/P1 controlled tests remain valid only for their previously reported configurations.
