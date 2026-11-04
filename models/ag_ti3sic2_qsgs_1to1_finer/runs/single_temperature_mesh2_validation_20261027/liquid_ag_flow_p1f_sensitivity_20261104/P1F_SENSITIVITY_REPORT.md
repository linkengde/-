# P1F numerical sensitivity study

## Decision

The prescribed P1F time-step and thermal-diffusivity sensitivity runs completed. The DT-2/DT-3 fitted center-of-mass speeds differ by 0.1943%, all runs conserve liquid volume within 0.1%, and the zero-temperature-gradient case has negligible net centroid translation. The four dimensionless thermal-diffusivity cases also complete and show less than 0.21% spread in fitted speed.

This is a pass of the requested sensitivity checks, not validation of physical Ag flow. A separate zero-gradient static-drop check still exhibits substantial local velocity and large raw pressure-Poisson residual values. The centroid criterion therefore passes while local static force balance remains unresolved. Do not use these runs to claim an Ag melt-pool velocity or proceed directly to a physical flow prediction.

## Scope and reproducibility

The runs use the dimensionless axisymmetric two-phase CLSVOF prototype from P1E, with an advected and diffused temperature tracer, temperature-dependent surface tension, and R16 wavelet AMR. The tests contain no Ag-specific material properties, Ti₃SiC₂ skeleton, latent heat, phase change, or experimental melt-pool geometry.

The fixed benchmark parameters are `R=1`, `rho1=rho2=1`, `mu1=mu2=1`, `Re=0.066`, `Ca=0.066`, `Gamma0=1`, and `GammaT=0.066`. Thus `Uref=0.066`, `t0=mu/GammaT`, and the analytical translation speed is `-0.1333333333 Uref`. The tracer starts from `T=gradT*x`, with `gradT=1` except for the zero-gradient control. Surface tension is `sigma=Gamma0+GammaT*T`.

The thermal Peclet number is `PeT=Uref*R/alphaT`. P1F-1 fixes `PeT=0.01` (`alphaT=6.6`) and changes only the Basilisk maximum timestep cap. The baseline cap, `0.00506658613`, is the largest logged P1E R16 timestep; the actual accepted step can be smaller. The three caps are 1×, 1/2×, and 1/4× this baseline. P1F-2 fixes the cap at 1× and changes only `PeT`.

For each run, the analyzer fits `xcm` against `t/t0` over the same `0.30–0.50 t0` interval. Since `Uref*t0=R=1`, the fitted slope is `U/Uref`. This avoids comparing interval velocities sampled at different iteration counts. Reproduce with:

```bash
export BASILISK_SRC=/path/to/pinned/basilisk/src
bash run_p1f_validation.sh
```

The script builds all three drivers, writes to a new timestamped `results/reproduction_*` directory, refuses to overwrite a nonempty output directory, runs P1F-1, gates P1F-2 on the P1F-1 decision, and includes the supplemental static-pressure initial-condition checks. The pinned Basilisk source revision used for the recorded runs is `9bb4f30aaf195feaa7a3a613468a1c888f5a4b9e`.

## P1F-1: timestep sensitivity

| Case | DT cap | Fitted speed (`U/Uref`) | Error vs analytical speed | Max volume drift | Fit R² |
|---|---:|---:|---:|---:|---:|
| DT-1 | 0.00506658613 | -0.130587553 | 2.0593% | 0.002117% | 0.999999989 |
| DT-2 | 0.002533293065 | -0.130826502 | 1.8801% | 0.000739% | 0.999999990 |
| DT-3 | 0.0012666465325 | -0.131081241 | 1.6891% | 0.000187% | 0.999999778 |
| Zero gradient, DT-3 | 0.0012666465325 | +0.000232410 | — | 0.004754% | 0.0651 |

DT-2 and DT-3 differ by 0.194337% in fitted speed, below the requested 1% limit. All listed cases have exit code zero and maximum volume drift below 0.1%. The zero-gradient fitted centroid speed is `0.1743%` of the magnitude of the analytical translation speed. For an explicit numerical interpretation of “no significant directional migration,” this report uses a 1%-of-reference-speed limit (`0.0013333 Uref`); this operational threshold is an analysis choice, not an independently supplied physical tolerance.

The analytical-speed error decreases from DT-1 to DT-3, but this is time-step sensitivity at a fixed grid, not a spatial grid-convergence result.

### Static-drop caveat

Although the zero-gradient centroid barely translates, it is not locally quiescent. Its maximum local speed is `0.06401` (all times) and `0.02243` after `0.30 t0`; the maximum raw `mgp.resa` in the fit window is `39.41`. Initializing the inside pressure to either sign of the nominal Laplace jump did not remove the local currents:

| Initial inside pressure | Fitted centroid speed (`U/Uref`) | Max local speed | Max local speed after `0.30 t0` | Max volume drift | Max raw `mgp.resa` in fit window |
|---|---:|---:|---:|---:|---:|
| 0 | +0.000232410 | 0.06401 | 0.02243 | 0.004754% | 39.41 |
| `+2 Gamma0/R` | +0.000258137 | 0.06699 | 0.06699 | 0.015729% | 40.33 |
| `-2 Gamma0/R` | +0.000551372 | 0.07502 | 0.02230 | 0.004971% | 39.99 |

Neither pressure sign produced a quiet static solution. The raw `mgp.resa` values are recorded as solver diagnostics; they are not normalized residuals and are not interpreted here as a dimensionless error. This local force-balance issue needs a targeted numerical investigation before treating a computed local velocity field as physical.

## P1F-2: thermal-diffusivity sensitivity

| `PeT` | Dimensionless `alphaT` | Fitted speed (`U/Uref`) | Speed change vs `PeT=0.01` | Interface temperature RMS departure | Tangential interface gradient RMS | Max volume drift |
|---:|---:|---:|---:|---:|---:|---:|
| 0.005 | 13.2 | -0.130665739 | -0.05987% | 0.00067197 | 0.817393 | 0.002200% |
| 0.010 | 6.6 | -0.130587553 | 0.00000% | 0.00095121 | 0.817391 | 0.002117% |
| 0.020 | 3.3 | -0.130536847 | +0.03883% | 0.00228297 | 0.817387 | 0.002107% |
| 0.050 | 1.32 | -0.130399300 | +0.14416% | 0.00652604 | 0.817553 | 0.002152% |

All four cases completed with exit code zero and mass/volume drift below 0.1%. The fitted-speed range from `PeT=0.005` to `0.05` is approximately 0.204% of the reference speed. In contrast, interface-temperature RMS departure increases by about 9.7× across that range. The benchmark speed is therefore weakly sensitive to this parameter over the tested dimensionless interval, while the interface temperature field changes appreciably. These `alphaT` values are dimensionless benchmark settings, not estimates of silver thermal diffusivity.

## Files and final status

The C drivers and shared implementation are in `results/p1f_*.[ch]`; `run_p1f_validation.sh` reproduces the runs. The raw per-step histories, stdout/stderr logs, and exit codes are retained alongside `p1f_timestep_summary.csv`, `p1f_diffusivity_summary.csv`, and `p1f_static_balance_summary.csv`.

**P1F sensitivity criteria: passed. Static local force-balance quality: unresolved.** The next numerical task should diagnose and reduce zero-gradient static-drop local currents and inspect the pressure-projection balance before introducing calibrated Ag properties. No physical melt-pool flow is predicted by this report.

## Separate P1E long-run record

P1F leaves the earlier P1E long-run artifacts untouched. The P1E report describes an earlier R16 run as ending near `0.478 t0` with exit code 139, while the archived `thermal_R16_Pe001.exit_code.txt` currently contains `124` and its stdout reaches approximately `0.502 t0`. This is an archive/provenance discrepancy; this P1F work does not overwrite either record or claim to diagnose the earlier termination. Reconcile that discrepancy as a separate task before relying on the long-run failure classification.
