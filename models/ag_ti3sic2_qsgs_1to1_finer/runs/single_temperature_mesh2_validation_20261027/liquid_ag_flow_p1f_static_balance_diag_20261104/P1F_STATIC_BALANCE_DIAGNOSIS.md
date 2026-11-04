# P1F static-drop force-balance diagnosis

## Decision

Repeated dynamic AMR is the main condition associated with the zero-gradient local-velocity bursts in this benchmark. R16 fixed and initially-adapted/frozen meshes remain quiet in the burst window, while dynamic AMR produces large liquid velocities. The same contrast appears at R32. The exact cause inside mesh transfer and interface-force reconstruction is not isolated yet, so this is a diagnosis of the failing path, not a production solver repair.

For the short R16 benchmark, adapting once at startup and then freezing the mesh is a viable controlled workaround: it retains the dynamic Marangoni migration speed and mass conservation through `0.5 t0`. It has only been tested over a centroid displacement of about one R16 cell and is not approved for longer interface motion or a real Ag melt pool.

## Method

The diagnostic reuses the P1F dimensionless axisymmetric CLSVOF and `integral.h` surface-tension implementation. It sets a constant temperature and constant surface tension for static tests, starts at zero velocity, and measures liquid-only maximum and volume-weighted RMS velocity, liquid centroid, pressure jump, meridional curvature, volume drift, cell count, and both Basilisk pressure projections. The liquid is the `1-f` phase for the signed-distance initialization used here.

The pressure solver reports an absolute residual. For comparison against Basilisk's `project()` stopping tolerance, the report uses `mgp.resa * dt^2 / TOLERANCE`; for the half-step projection it uses `mgpf.resa * (dt/2)^2 / TOLERANCE`. A value below one means the configured tolerance was reached. The historical raw `mgp.resa` spikes alone were not evidence of solver divergence.

Runs use R16 or R32, identical `0.20 t0` windows unless stated otherwise, and `DT=0.25` of the P1E maximum timestep cap. R16 fixed uniform, dynamic AMR, and frozen-after-first-AMR runs use the same physical timestep and initialization. The R32 frozen run ends at `0.10 t0`, spanning the R32 dynamic-AMR spike near `0.073 t0`.

## Static-drop results

| Case | Liquid peak speed | Time | Speed in the matched burst window | Other evidence |
|---|---:|---:|---:|---|
| R16 fixed uniform, zero gradient | 0.00242 | 0.00133 `t0` | 0.000194 at 0.173–0.176 `t0` | Volume drift `1.91e-5%`; mean pressure jump 2.00316; mean meridional curvature 0.99981 |
| R16 AMR once, then frozen, zero gradient | 0.00245 | 0.00133 `t0` | 0.000204 at 0.173–0.176 `t0` | Volume drift `1.91e-5%`; mean pressure jump 2.00317; mean curvature 0.99956 |
| R16 dynamic AMR, zero gradient | 0.06274 | 0.17401 `t0` | 0.06274 | Pressure jump ranges 1.9169–2.0862 and mean curvature error reaches 2.06% in that window |
| R16 dynamic AMR, half timestep | 0.03652 | 0.16376 `t0` | 0.000687 at 0.173–0.176 `t0` | Halving `DT` moves and reduces the burst, but does not remove bursts over the run |
| R32 dynamic AMR, zero gradient | 0.45335 | 0.07264 `t0` | 0.45335 at 0.065–0.080 `t0` | Pressure jump ranges 1.7138–2.2249; both pressure projections reach about 0.997 of their configured residual tolerance |
| R32 AMR once, then frozen, zero gradient | 0.00157 | startup | 0.000106 at 0.065–0.080 `t0` | Volume drift `6.50e-5%`; mean pressure jump 2.00078; mean curvature error below 0.04% |

The R16 dynamic-AMR peak is about 320 times the fixed-grid maximum in the same time window. At R32, the dynamic-AMR burst exceeds the frozen-grid value by more than three orders of magnitude in the matched window. The R16 dynamic run's maximum pressure residual was still below its configured tolerance; the R32 run approached that tolerance. Pressure initialization with either sign of the Laplace jump did not remove the R16 dynamic-AMR bursts.

Reducing the pressure tolerance tenfold did not resolve the issue: the R16 dynamic-AMR peak remained 0.07735 and the R32 peak remained 0.39691, with the burst shifting in time. Refining on the interface and temperature without including velocity was worse, reaching 0.17545 on R16. These controls point to repeated adaptive-grid/interface-field updates as the principal issue, rather than pressure-jump sign or a simple pressure tolerance setting. They do not identify whether velocity prolongation, level-set/VOF synchronization, or the discrete capillary force transfer is the exact trigger.

## Frozen-mesh dynamic control

With R16 AMR applied once at startup and then frozen, the temperature-gradient Marangoni run completes to `0.5 t0` with exit code zero. A least-squares fit of `xcm(t)` over `0.30–0.50 t0` gives:

- fitted speed: `-0.130558 Uref`;
- analytical speed: `-0.133333 Uref`;
- relative error: `2.0813%`;
- fit `R2`: `0.999999881`;
- maximum liquid-volume drift: `0.0001908%`.

This shows that the P1F dynamic test can meet its speed and conservation checks without repeated remeshing. During the interval, the droplet centroid moves approximately one R16 cell, so the result validates only this short controlled translation. Longer movement could leave the initially-refined band and requires a new test.

## P1E record discrepancy

P1F has not modified the earlier P1E artifacts. The P1E report says the earlier R16 run ended near `0.478 t0` with exit code 139. The archived `thermal_R16_Pe001.exit_code.txt` says 124, and its stdout reaches about `0.502 t0`. These records do not identify the same termination consistently. The raw P1E artifacts are retained; reconcile their provenance separately before classifying that run.

## Next numerical step

Keep repeated dynamic AMR disabled for this benchmark until a transfer repair passes the static drop. Diagnose the remesh transition by logging pre/post-remesh velocity, face flux divergence, `f`/`d` consistency, curvature, surface-force impulse, and the two pressure projections. Test a divergence-preserving velocity transfer and synchronized VOF/level-set reconstruction as separate changes. After a repair, repeat the R16 and R32 static controls and the frozen-mesh dynamic control before using physical Ag properties.

All values here are dimensionless benchmark results; none are Ag property or melt-pool predictions.
