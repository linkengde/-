# P0 VOF static-drop validation report

## Decision

P0 is **partially validated, not fully passed**. The VOF-height-function sphere and the 72° sessile drop completed dynamic tests at all three resolutions, R/Δx = 16, 32, and 64. Their pressure error, contact-angle error, mass drift, and parasitic velocity improve consistently with refinement. The separate analytical-curvature sphere test passed at R/Δx = 16 and 32, but its pressure projection diverges at R/Δx = 64, including MPI and smaller-time-step attempts. Therefore the full requested matrix is not passed, and this report does not advance to P1.

No heating, production enthalpy-solver modification, Ag melt-pool prediction, or Marangoni calculation was performed.

## Solver and setup

The independent solver is Basilisk C/qcc on an adaptive 3-D octree grid, using its centered incompressible Navier–Stokes solver in the Stokes limit, VOF interface tracking, height-function curvature, and height-function contact-angle boundary. The analytical-curvature sphere supplies κ = 2/R to isolate the pressure/force balance; the VOF-curvature sphere calculates κ from the interface. All values are dimensionless verification data, not simulations with Ag properties:

- Liquid and gas density: 1; dynamic viscosity: 1; surface tension: 1; radius: 1.
- Gravity and imposed flow: zero; no evaporation or thermal gradient.
- Sessile-drop initial shape: analytical spherical cap with a 72° contact angle on z = Z0. Basilisk calls this plane `back`; the boundary is set on `h.t[back]` and `h.r[back]`.
- Domain length: 8R. R/Δx = 16, 32, and 64 correspond to levels 7, 8, and 9.
- Dynamic runs end at t = 0.05 with capillary time-step factor 0.5 and `TOLERANCE=1e-4`. The R/Δx=32 time-step repeats use factor 0.25 at the same end time.
- R/Δx=16 and 32 cases ran in serial. The successful R/Δx=64 dynamic cases ran on four Open MPI ranks. A same-case R/Δx=32 comparison showed the MPI and serial pressure error differed by 0.0000022 percentage points; maximum speed differed by about 0.8%.

The Basilisk source snapshot, compiler and MPI provenance are recorded in [vof_solver_version.json](vof_solver_version.json); build and reproduction instructions are in [BASILISK_BUILD.md](BASILISK_BUILD.md). The upstream Basilisk source endpoint returned HTTP 403 in this environment, so the run used the explicitly pinned 2019 third-party GitHub mirror snapshot. No official release tag or provenance signature was verified.

## Dynamic results

Pressure error is `100 × (measured Δp − 2σ/R)/(2σ/R)`. For these dimensionless properties, the reported capillary number is Ca = μUmax/σ = Umax. Mass drift is the final-minus-initial VOF volume relative to the initial VOF volume.

| Dynamic case | R/Δx | Pressure error | Maximum speed / Ca | Mass drift | Measured contact angle |
|---|---:|---:|---:|---:|---:|
| Exact-curvature sphere | 16 | −0.000200% | 2.44×10⁻⁶ | −0.000281% | — |
| VOF-curvature sphere | 16 | +0.502878% | 8.31×10⁻⁴ | −0.000280% | — |
| 72° sessile drop | 16 | +0.503560% | 1.09×10⁻³ | −0.000051% | 71.9840° |
| Exact-curvature sphere | 32 | −0.0000031% | 4.94×10⁻⁷ | +0.0000164% | — |
| VOF-curvature sphere | 32 | +0.125733% | 2.20×10⁻⁴ | +0.0000182% | — |
| 72° sessile drop | 32 | +0.108620% | 4.83×10⁻⁴ | +0.0000076% | 71.9949° |
| VOF-curvature sphere | 64 | +0.030972% | 8.62×10⁻⁵ | +0.0000037% | — |
| 72° sessile drop | 64 | +0.020252% | 2.82×10⁻⁴ | −0.0000017% | 71.9985° |

The VOF-curvature sphere and sessile-drop pressure error and maximum speed decrease at every completed grid refinement. The sessile contact-angle error decreases from −0.0160° to −0.00509° to −0.00150°. Every completed dynamic case has mass drift far below 0.1%.

The R/Δx=64 VOF-curvature sphere and sessile runs each completed 91 time steps. Their reported cell counts at the end state were 678,112 and 243,888, respectively. Neither run reported a pressure-convergence or VOF-CFL warning.

## Analytical-curvature R/Δx=64 failure

The analytical-curvature sphere did not produce a valid R/Δx=64 dynamic result. Its pressure projection fails to converge, velocity then becomes nonphysical, and the run reports VOF-CFL violations. The failure was reproduced in serial and MPI. It persisted in short MPI diagnostics with a smaller time-step factor (0.25), a smoothed initial pressure `p=2f`, and increased diagnostic viscosity μ=10. These altered runs are failure probes only; their invalid result lines are not included in the CSV summaries. Logs and details are preserved in `basilisk_cases/diagnostics/R64_DYNAMIC_ATTEMPTS.md`.

R/Δx=64 static geometry checks did complete. Before any hydrodynamic step, the sphere volume error was −0.0183% and mean curvature error +0.3870%; the 72° cap volume error was −0.0229%, curvature error +0.2825%, and fitted angle error −0.00115°. The initialized pressure jump equals 2 because it is prescribed as the initial condition; this is not counted as a dynamic validation result.

## Time-step sensitivity

At R/Δx=32, reducing the capillary time-step factor from 0.5 to 0.25 changed VOF-sphere pressure error from +0.12573% to +0.13130%, and sessile-drop pressure error from +0.10862% to +0.11303%. Contact-angle error changed from −0.005094° to −0.005218°. Maximum speed changed by about 11% for the sphere and 4% for the sessile drop, while mass drift stayed below 0.00002%. The pressure and angle conclusions are insensitive to this time-step change; the speed variation remains part of the numerical uncertainty.

## Threshold assessment

| Check | Result |
|---|---|
| Laplace pressure, VOF-computed curvature | Pass dynamically at R/Δx = 16, 32, and 64; errors decrease with refinement |
| Laplace pressure, prescribed analytical curvature | Pass dynamically at 16 and 32; R/Δx=64 projection diverges |
| 72° sessile contact angle | Pass dynamically at all three resolutions; finest error is −0.00150° |
| Mass conservation <0.1% drift | Pass for every completed dynamic case |
| Grid trend | Established for VOF-curvature sphere and sessile cases; analytical-curvature three-grid trend remains incomplete |
| Time-step sensitivity | Pressure and angle pass at R/Δx=32; speed sensitivity is reported above |
| Overall P0 | Partial; do not call the complete requested matrix passed |

## Files and reproduction

The Basilisk CSV summaries and run configuration are in `basilisk_results/`:

- `vof_static_drop_summary.csv`
- `vof_contact_angle_summary.csv`
- `vof_mass_conservation.csv`
- `vof_grid_convergence.csv`
- `vof_timestep_sensitivity.csv`
- `basilisk_run_configuration.json`

The primary matrix is reproduced with `python3 run_basilisk_matrix.py` after building qcc as described in `BASILISK_BUILD.md`. To reproduce the R/Δx=64 VOF-curvature sphere and sessile-drop runs, install MPI development headers, set `BASILISK_SRC` to the pinned Basilisk `src` directory, then run `P0_MPI_NP=4 python3 run_basilisk_r64_mpi.py`. This supplemental runner appends/upserts the two dynamic rows and stores their per-case logs. R/Δx=64 analytical-curvature failure probes are diagnostic only and are documented separately.
