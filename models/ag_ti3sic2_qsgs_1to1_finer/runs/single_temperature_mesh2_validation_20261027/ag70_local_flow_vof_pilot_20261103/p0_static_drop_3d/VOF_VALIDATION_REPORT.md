# P0 VOF static-drop validation report

## Decision

P0 is **partially validated, not fully passed**. Dynamic checks on the R/Δx=16 and 32 meshes meet the requested pressure, contact-angle, mass-conservation, and time-step criteria. R/Δx=64 passed static geometry and curvature checks, but its transient pressure/velocity projection did not finish within the available run limit. The t=0 R/Δx=64 measurements are not treated as dynamic pressure-balance or parasitic-current results. Do not advance to P1 on the claim that the full P0 matrix passed.

No heating, production enthalpy-solver modification, Ag melt-pool prediction, or P1 Marangoni calculation was performed.

## Solver and setup

The independent solver is Basilisk C/qcc on an adaptive three-dimensional octree grid, with the centered incompressible Navier–Stokes solver in the Stokes limit, VOF interface tracking, height-function curvature, and the height-function contact-angle boundary. The exact-curvature mode prescribes κ=2/R on a free sphere to isolate force/pressure balance. The automatic mode calculates curvature from the VOF interface. All cases are dimensionless verification cases, not physical Ag simulations:

- Liquid and gas density: 1; dynamic viscosity: 1; surface tension: 1; radius: 1.
- Gravity and external flow: zero; evaporation and thermal gradients: absent.
- Sessile drop: analytical spherical cap with 72° contact angle on z=Z0. In Basilisk's 3-D boundary names, this plane is `back`; the condition is applied to `h.t[back]` and `h.r[back]`.
- Domain length: 8R. R/Δx=16, 32, and 64 correspond to octree levels 7, 8, and 9.
- Dynamic cases: t=0.05 dimensionless time, capillary time-step factor 0.5. The R/Δx=32 time-step repeats use factor 0.25 at the same end time. `TOLERANCE=1e-4`.
- The initial pressure is set to the Laplace-compatible jump 2σ/R, and velocity is zero. Reported pressure errors are measured after the dynamic run, not inferred from the initial value.

The exact source snapshot, compiler, provenance limits, and build command are recorded in [vof_solver_version.json](vof_solver_version.json) and [BASILISK_BUILD.md](BASILISK_BUILD.md). The existing 2-D OpenFOAM pilot remains in its original files; it is not combined with these Basilisk results.

## Dynamic results

Pressure error is `100 × (measured Δp − 2σ/R)/(2σ/R)`. Parasitic velocity is the maximum velocity magnitude, and with the dimensionless properties above its reported capillary number is Ca=μUmax/σ. Initial volume error is against the analytical sphere or spherical-cap volume; mass drift is the final-minus-initial VOF volume relative to the initial VOF volume.

| Dynamic case | Pressure error | Maximum speed / Ca | Mass drift | Contact-angle error |
|---|---:|---:|---:|---:|
| Exact-curvature sphere, R/Δx=16 | −0.00020% | 2.44×10⁻⁶ | −0.000281% | — |
| VOF-curvature sphere, R/Δx=16 | +0.5029% | 8.31×10⁻⁴ | −0.000280% | — |
| 72° sessile drop, R/Δx=16 | +0.5036% | 1.09×10⁻³ | −0.000051% | −0.0160° |
| Exact-curvature sphere, R/Δx=32 | −0.0000031% | 4.94×10⁻⁷ | +0.0000164% | — |
| VOF-curvature sphere, R/Δx=32 | +0.1257% | 2.20×10⁻⁴ | +0.0000182% | — |
| 72° sessile drop, R/Δx=32 | +0.1086% | 4.83×10⁻⁴ | +0.0000076% | −0.0051° |

Across the dynamic R/Δx=16 and 32 cases, pressure error is below 1%, contact-angle error is below 2°, and mass drift is far below 0.1%. The VOF-curvature sphere and sessile-drop maximum speed both decrease from R/Δx=16 to 32. These two dynamic mesh levels show a consistent improvement, but do not establish convergence through R/Δx=64.

## Time-step sensitivity

At R/Δx=32, halving the capillary time-step factor from 0.5 to 0.25 changed the sphere pressure error from +0.12573% to +0.13130%, and the sessile-drop error from +0.10862% to +0.11303%. Contact-angle error changed from −0.005094° to −0.005218°. Maximum speed changed by about 11% for the sphere and 4% for the sessile drop; the absolute speeds stayed below 5.1×10⁻⁴. Both runs retained mass drift below 0.00002%. Pressure and contact-angle conclusions are insensitive to this step change; the small parasitic-speed change should remain visible in the results.

## R/Δx=64 static geometry results and dynamic limitation

At t=0, before any hydrodynamic time step, the VOF geometry and height-function curvature give:

| Static geometry check | Initial volume error | Curvature mean error | Fitted angle error |
|---|---:|---:|---:|
| Sphere | −0.0183% | +0.3870% | — |
| 72° spherical cap | −0.0229% | +0.2825% | −0.00115° |

These values verify the initialized geometry, static contact-angle boundary, and curvature field only. The pressure jump at t=0 is 2 because the analytical Laplace pressure was prescribed as the initial condition; it is not counted as a dynamic VOF pressure-balance result.

The R/Δx=64 sphere mesh contained 732,264 leaf cells. A dynamic exact-curvature run reached pressure/VOF warnings and terminated before a result. Follow-up short attempts with a Laplace-compatible initial pressure, looser `TOLERANCE=1e-4`, and even zero viscosity did not complete one step within 180 seconds. Logs and attempt details are preserved in `basilisk_cases/diagnostics/`. This leaves the fine-grid dynamic pressure error, parasitic velocity, and time-step sensitivity unverified. It is the reason the overall decision remains partial.

## Threshold assessment

| Check | Result |
|---|---|
| 3-D Laplace balance, prescribed curvature | Pass at dynamic R/Δx=16 and 32 |
| 3-D Laplace balance, VOF curvature | Pass at dynamic R/Δx=16 and 32; R/Δx=64 dynamic result unavailable |
| 72° static contact angle | Pass dynamically at R/Δx=16 and 32, and geometrically at R/Δx=64 |
| Mass conservation <0.1% drift | Pass for every completed dynamic case |
| Grid trend | Dynamic speed/error improve between 16 and 32; full 16/32/64 dynamic convergence not established |
| Time-step sensitivity | Pass for pressure/contact angle at R/Δx=32; velocity changes are reported above |
| Overall P0 | Partial; do not call the complete three-grid dynamic matrix passed |

## Files and reproduction

The machine-readable Basilisk results are in `basilisk_results/`:

- `vof_static_drop_summary.csv`
- `vof_contact_angle_summary.csv`
- `vof_mass_conservation.csv`
- `vof_grid_convergence.csv`
- `vof_timestep_sensitivity.csv`
- `basilisk_run_configuration.json`
- `compile.log`

Each completed case has matching stdout/stderr logs in `basilisk_cases/`. Run `python3 run_basilisk_matrix.py` after building the pinned qcc source as described in `BASILISK_BUILD.md`. The R/Δx=64 static rows use `run_mode=geometry_only`; the unavailable dynamic cases are not silently filled or extrapolated.
