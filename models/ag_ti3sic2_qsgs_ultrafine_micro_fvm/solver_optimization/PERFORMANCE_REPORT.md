# PERFORMANCE_REPORT

## Measured platform and method

Cloud allocation: 4 CPU cores (`cpu.max=400000 100000`), 32 GiB memory cgroup limit, and about 51 GiB free workspace disk during this stage. Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8, and VTK 9.3.1 were available. Numba was initially absent, then installed temporarily under `/tmp` as version 0.67.0; PyAMG 5.3.0 was also installed temporarily for the pilot. `nvidia-smi` was unavailable and no GPU was used. The code path uses SciPy sparse linear interfaces and FFT plus the optional Numba kernels; no commercial or external CFD solver was used.

The solver remains the existing single-temperature, implicit enthalpy FVM; the optimized CPU path uses Numba fused material property, face-conductance, and seven-point stencil kernels. Time integration and source, geometry, candidate material data, 1 K phase regularization, nonlinear tolerance, and latent heat were not changed. The EW-safe strategy changed only the inexact PCG solve tolerance per Newton iteration. FFT remains the default preconditioner.

## Results

| Case | Status | Wall time | CPU time | Peak RSS | Iterations / retries | Energy residual |
|---|---|---:|---:|---:|---|---:|
|50³ pure Ag Numba EW-safe, 4 threads, dt-factor 0.25|complete|8.50 s|17.21 s|283.1 MiB|120 steps / 0 retries|1.47×10⁻⁹ normalized|
|50³ QGS reference, fixed, dt-factor 0.25|complete|reference total wall not reliable; solver step sum 62.91 s|not retained in old summary|204.5 MiB|458 Newton / 16,291 PCG / 0 retries|9.51×10⁻¹⁰ normalized|
|50³ QGS Numba fixed, 1 thread|complete|54.42 s|54.15 s|288.1 MiB|458 / 16,293 / 0|9.51×10⁻¹⁰|
|50³ QGS Numba fixed, 2 threads|complete|48.80 s|82.14 s|286.5 MiB|458 / 16,293 / 0|9.51×10⁻¹⁰|
|50³ QGS Numba fixed, 4 threads|complete|42.38 s|114.95 s|288.3 MiB|458 / 16,293 / 0|9.51×10⁻¹⁰|
|50³ QGS Numba EW-safe, 4 threads|complete|25.89 s|61.27 s|287.8 MiB|510 / 6,217 / 0|9.98×10⁻¹⁰|
|100³ pure Ag Numba EW-safe, 4 threads|complete|69.55 s|125.17 s|501.1 MiB|63 steps / 255 Newton / 1,949 PCG / 3 retries|5.10×10⁻¹¹ normalized|
|100³ QGS Numba EW, 4 threads|complete|508.33 s|1,033.92 s|508.3 MiB|1,666 / 12,282 / 58|3.61×10⁻⁹|
|100³ QGS Numba EW-safe, 4 threads|complete|553.43 s|1,157.06 s|512.3 MiB|1,256 / 12,320 / 24|1.43×10⁻⁹|
|150³ QGS Numba EW-safe, 4 threads|16-step limit|172.06 s|310.00 s|1,093.6 MiB|105 / 498 / 5|2.24×10⁻¹⁰|

At 50³, the accepted-step solver time improved 3.94× under Numba + EW-safe relative to the original fixed-tolerance implementation. This satisfies the 3× solver-time target. The 100³ full numerical cycle completed in 9.22 minutes, meeting the 15-minute wall target under the same seed/source setup. There is no completed whole-cycle baseline for the old 100³ code, so a 100³ acceleration factor is not claimed. The old run stopped after more than 60 minutes without final output.

Pure-Ag controls completed in 8.50 s (50³) and 69.55 s (100³) using the same candidate heat pulse, 1 ns heating and 2 μs cooling. The 100³ pure-Ag run had no retry and a normalized energy residual of 5.10×10⁻¹¹. This provides a material/geometry control, not an Ag/Ti₃SiC₂ result.

The 100³ EW-safe trajectory traded about 9% more wall time than EW for fewer retries and a lower normalized energy residual. The local Ag-weighted hot liquid fractions differed by about 0.56%; exact time/space convergence is still required before physical interpretation.

The 150³ short stage reached the onset of melting, but its step size contracted to 2.34 ps and it used five retries in only 16 accepted steps. It did not complete the 1 ns heating period or cooling. No full-cycle runtime is extrapolated from these first steps, and the full 150³ cycle was not launched.

## Conclusion

- 50³ QGS short cycle: **performance GO**, numerical regression agreement demonstrated for the tested schedule.
- 100³ QGS short cycle: **conditional performance GO**; completed within 15 minutes, but temporal/spatial convergence and source/material validity remain open.
- 150³ QGS full cycle: **NO-GO this stage**; only the 16-step limit probe completed. Continue with improved hierarchy reuse/restart I/O profiling and a measured melt-band schedule before budgeting a full run.

Every heat input in these benchmarks is a numerical volumetric pulse, not a physical arc boundary.
