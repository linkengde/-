# R/Δx=64 dynamic run attempts

R/Δx=64 static geometry and curvature completed successfully, but transient pressure/velocity projection did not.

1. The first full exact-curvature sphere case (`sphere_exact_R64_dt0p5.stderr.log`) reached 732,264 leaf cells, then reported failure to converge pressure, a VOF CFL violation, and exited without a `RESULT` record.
2. A retry with the analytical Laplace-pressure initial field and `TOLERANCE=1e-6` did not reach a result in about 2 minutes and was stopped for diagnosis. Its captured log is `debug_R64_pressure_init.err`.
3. With `TOLERANCE=1e-4`, the same R/Δx=64 short test with μ=1 exceeded a 180 s timeout without completing a step (`debug_R64_tol1e4.err`).
4. With `TOLERANCE=1e-4` and μ=0 to skip the implicit viscous solve, the same short test also exceeded 180 s (`debug_R64_mu0.err`).

The remaining logs preserve initialization counts and termination status. These are resource/runtime failures or solver non-convergence, not successful dynamic tests. No R/Δx=64 dynamic pressure, velocity, contact-angle, or time-step conclusion is claimed.
