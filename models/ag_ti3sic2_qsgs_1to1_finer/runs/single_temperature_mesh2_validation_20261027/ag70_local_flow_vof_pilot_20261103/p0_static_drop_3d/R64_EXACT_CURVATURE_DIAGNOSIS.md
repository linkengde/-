# R64 exact-curvature VOF diagnosis

## Finding

The R64 failure is specific to the combination of **prescribed constant curvature and the wavelet-coarsened AMR mesh**, rather than a general failure of the 3-D VOF solver. On the initial geometrically refined octree frozen at 732,264 cells, the exact-curvature case completed both pressure projections. Applying one wavelet adaptation produced a 681,472-cell coarse/fine topology; with that topology frozen before flow, the first pressure projection failed. Applying the same wavelet change after the first step instead allowed the first projection and failed the second. The automatic height-function curvature control passed both projections on the corresponding wavelet AMR path.

This isolates the failure category to exact-curvature force/pressure projection on the wavelet-AMR topology. Resetting the analytic Laplace pressure after the initialization-time AMR pass did **not** remove the failure, so simple transfer of the initial pressure jump is not the sole cause. The cell-local trigger is not yet proven: global f, pressure, curvature, and force-balance summaries change very little across the first adaptation, which points toward a localized coarse/fine prolongation or face-gradient balance issue. No unsupported numerical repair has been applied.

## Test setup and scope

All completed comparisons use the same pinned Basilisk mirror build, MPI-4, R=1, domain length 8R, ρ=μ=σ=1, `TOLERANCE=1e-4`, `CFL=0.5`, and the R/Δx=64 capillary step factor 0.5 (`dt_cap=0.00110193278037`, first step 0.000550966390183). Each projection test ends after two nominal steps to t=0.001 or stops when the pressure projection diverges. These are dimensionless solver tests, not Ag-material or melt-pool calculations.

The exact-curvature mode supplies φ=σκ=2 on mixed VOF cells through one `iforce.h`. No `tension.h` force is active at the same time. The automatic mode uses the Basilisk height-function curvature route. The pressure/force residual in the CSVs is an audit metric computed on interface-crossing faces:

`r = ∇p − φ∇f`, `r_rms = sqrt(sum(r² w)/sum(w))`,

with `w=|Δf|Δ²`. It is not the multigrid residual. The pressure projection residual and tolerance are recorded separately.

## Fixed mesh versus AMR

| Case | Mesh at projection | Result | Main evidence |
|---|---|---|---|
| D1 full-domain uniform R64 | 512³ = 134,217,728 cells requested | Resource-limited before solver initialization; no numerical result | Guard stopped at about 30 GiB; cgroup OOM counters remained zero |
| D2 geometric interface-refined octree, frozen | 732,264 cells, no wavelet coarsening | Both projections converged | 2 and 1 iterations; residual/tolerance 0.2524 and 0.02087; final Umax 9.59×10⁻⁵; mass drift +1.26×10⁻⁶% |
| D2 wavelet AMR once in initialization, then frozen | 732,264 → 681,472 cells | First projection diverged | 100 iterations; residual/tolerance 2.6702×10¹⁰; Umax 1.7822×10²³ |
| D2 same pre-adapted wavelet grid, reset initial Laplace pressure, then frozen | 681,472 cells | First projection diverged identically | Same 100 iterations, residual/tolerance 2.6702×10¹⁰, Umax 1.7822×10²³ |
| D2 wavelet AMR after first step, then frozen | 732,264 → 681,472 cells | First projection passed; second diverged | second residual/tolerance 1.7736×10¹⁰; Umax 1.4525×10²³ |
| D3 dynamic wavelet AMR, exact curvature | 732,264 → 681,472 → 678,112 cells | First projection passed; second diverged | Same second-projection blow-up as the one-update case |
| D3 dynamic wavelet AMR, automatic curvature | 732,264 → 681,472 → 678,112 cells | Both projections converged | 2 and 1 iterations; residual/tolerance 0.2550 and 0.02694; mass drift +1.21×10⁻⁶% |

The full uniform R64 run did not reach an `INIT` record. An initial per-process address-space-limited attempt ended with SIGSEGV; a second attempt was stopped by the explicit cgroup memory guard at a recorded peak of 32,243,068,928 bytes. The cgroup `oom` and `oom_kill` counters remained zero. Both are **resource-limited attempts**, not numerical failures or successes; see `basilisk_cases/diagnostics/r64_exact_uniform.status.json`.

The static wavelet-prefreeze test is the most discriminating result: ongoing adaptation is unnecessary for failure. The one-time wavelet topology alone is enough. Resetting the initial pressure on that topology gives the same first-projection divergence, so a stale transferred Laplace jump by itself is ruled out. The passing frozen geometric grid differs in topology and cell distribution, so its pass does not establish that all fixed meshes are stable.

## Pressure projection sequence

For the passing exact-curvature geometric frozen mesh, projection residual/tolerance falls from 0.25239 on step 0 to 0.02087 on step 1. The pressure field remains finite, f stays within [0,1], and mass drift remains under 1.3×10⁻⁶%.

On the wavelet-prefrozen exact mesh, the first projection reaches 100 iterations: residual before projection is 4.7663×10²⁸, residual afterward is 8.7961×10¹², tolerance is 329.42, and residual/tolerance is 2.6702×10¹⁰. The velocity reaches 1.7822×10²³ while f remains in [0,1] and no nonfinite cells have yet been counted. In the one-update and dynamic tests, the same pressure blow-up appears on step 1, with residual/tolerance 1.7736×10¹⁰. The resulting pressure magnitude is about 10²⁴ and velocity about 10²³; VOF CFL warnings follow. A zero process exit code in the factor-0.5 historical probe does not make this invalid field a passing run.

Pressure projection diagnostics are in `basilisk_results/r64_pressure_projection.csv`; raw logs are under `basilisk_cases/diagnostics/r64_*.stderr.log`. Historical factor-0.35 and factor-0.20 probes are included in the pressure table. They also fail on the second projection, consistent with the earlier conclusion that the capillary timestep bound is not the primary cause.

## Force and field-transfer observations

The wavelet step changes cell count from 732,264 to 681,472, but the logged global f bounds, number of mixed cells, curvature statistics, and maximum/minimum pressure remain nearly unchanged. For the exact case, the measured interface-weighted pressure-force residual after the first converged projection is 0.16048 RMS (relative 0.00225); after the second projection on the frozen geometric mesh it falls to 0.05397 (relative 0.000755). On the wavelet-prefrozen mesh the first projection instead creates a residual around 6.919×10²³ RMS (relative 9.683×10²¹). The automatic-curvature path on dynamic AMR remains finite, with relative residual RMS 0.01757 after step 0 and 0.01549 after step 1.

The key observation is that coarse global summaries do not reveal the instability before it occurs, and resetting the initial pressure does not cure it. That is consistent with, but does not prove, a localized coarse/fine mismatch between the prescribed surface-force potential, pressure prolongation, and the projection gradient. The next numerical debugging step should inspect face-local values on the coarsest/fine interface around the pressure-solver failure and compare their `f`, φ, p, and prolongation stencils against the height-function case.

The full force-balance series, including values before and after AMR events, is in `basilisk_results/r64_curvature_force_balance.csv`. The mesh snapshots include cell counts, pressure extrema/mean, f range, mixed-cell count and curvature mean/standard deviation/extrema. For compatibility with the first diagnostic build, the raw MESH log field named `kappa_rms` is `statsf(kappa).stddev` (standard deviation), not the root-mean-square curvature.

## Explicit decision boundary

- **Automatic-curvature VOF path:** passed the R64 standalone static-drop controls, including the existing 91-step `t=0.05` sphere and 72° sessile-drop cases, and the new same-topology two-projection AMR control. This supports the narrow P0 claim that the tested automatic-curvature VOF static-drop path is operational. It does not validate an integrated production heating/flow model or an Ag melt pool.
- **Exact-curvature auxiliary path:** has an independent instability on the wavelet-AMR topology. It remains unvalidated at R64 and should not be used for a result that depends on this test path until the coarse/fine force-pressure mismatch is fixed or its domain of validity is explicitly restricted.
- **Overall P0:** remains partial because the requested exact-curvature R64 diagnostic is unresolved and the full uniform-grid control could not be completed within the memory guard. This does not imply that every VOF path fails.
- **P1:** not run in this task. See `P0_FINAL_DECISION.md` for the final gate statement.

## Reproduction and artifacts

Build `basilisk_r64_diagnosis.c` with the pinned source and MPI setup described in `BASILISK_BUILD.md`, then run from this directory:

```sh
mpirun -np 4 /tmp/basilisk-r64-diag-mpi 9 sphere exact 0.5 0.001 1.0 freeze
mpirun -np 4 /tmp/basilisk-r64-diag-mpi 9 sphere exact 0.5 0.001 1.0 prefreeze
mpirun -np 4 /tmp/basilisk-r64-diag-mpi 9 sphere exact 0.5 0.001 1.0 prefreeze_resetp
mpirun -np 4 /tmp/basilisk-r64-diag-mpi 9 sphere exact 0.5 0.001 1.0 once
mpirun -np 4 /tmp/basilisk-r64-diag-mpi 9 sphere exact 0.5 0.001 1.0 dynamic
mpirun -np 4 /tmp/basilisk-r64-diag-mpi 9 sphere auto 0.5 0.001 1.0 dynamic
```

The uniform mode was attempted under an external memory monitor and is documented in `basilisk_cases/diagnostics/r64_exact_uniform.status.json`. Generate the three numerical CSVs from preserved stdout/stderr logs with:

```sh
python3 basilisk_results/build_r64_diagnosis_tables.py
```
