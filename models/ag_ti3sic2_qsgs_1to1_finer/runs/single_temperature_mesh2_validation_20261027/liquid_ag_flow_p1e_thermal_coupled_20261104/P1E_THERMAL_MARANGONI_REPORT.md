# P1E: evolving-temperature Marangoni prototype

## Decision

This run directory is on branch `codex/ag70-p1e-thermal-marangoni`, based on the P1D validation commit `743ea84e31c41af3eb1975bb7cd0471e68845dc3`.

The CLSVOF + `integral.h` route now advances a dimensionless temperature tracer by advection and diffusion, and refreshes surface tension from that evolving field. It passes the controlled half-`t0` checks below at R8 and R16. This supports continuing with this solver route for further numerical development. It does not validate a real Ag/Ti₃SiC₂ melt pool, and it is not a production heat-transfer model.

## Model and scope

- Basilisk `two-phase-clsvof.h` and `integral.h`; pinned source release `v2026-07-20`, commit `9bb4f30aaf195feaa7a3a613468a1c888f5a4b9e`.
- Axisymmetric equal-density/equal-viscosity test fluids, with the same controlled Marangoni translation benchmark as P1D: $Re=Ca=0.066$, $U_{ref}=0.066$, $t_0=1/0.066$, and reference migration speed $-0.133333U_{ref}$.
- Dimensionless temperature starts as $T=x$, is included as an advected tracer, and diffuses implicitly with equal diffusivity in both fluids. The numerical thermal Peclet number is $Pe_T=0.01$, so $\alpha_T=U_{ref}R/Pe_T=6.6$. Axial ends use fixed linear-gradient temperatures; radial boundaries use zero normal gradient.
- The interface coefficient is σ=Γ₀+ΓₜT with Γ₀=1 and Γₜ=0.066. It is refreshed after each tracer advection/diffusion step and passed to `integral.h`.
- The temperature scalar is appended by copying Basilisk's statically initialized tracer list; the initial in-place append caused an invalid realloc and was replaced with `list_concat`. The corrected source builds cleanly and the logged checks finish without solver warnings.
- The runs stop at $0.5t_0$. No Ag-specific thermal properties, Ti₃SiC₂ skeleton, experimental dimensions, or Ag70 geometry are used. No production solver or heating case was changed.

## Results

| Case | $R/\Delta x$ | Mesh | Speed / $U_{ref}$ | Reference error | Maximum logged mass drift | Interface temperature RMS departure |
|---|---:|---|---:|---:|---:|---:|
| Positive gradient | 8 | uniform | –0.127309 | 4.518% | 0.003817% | 0.0015874 |
| Positive gradient | 8 | AMR | –0.127201 | 4.600% | 0.003794% | 0.0015906 |
| Positive gradient | 16 | AMR | –0.130556 | 2.083% | 0.002374% | 0.0015281 |
| Zero gradient | 8 | AMR | −8.67×10⁻⁶ | not applicable | 0.002006% | 0 |

The R8 fixed-grid and AMR speeds differ by about 0.086%. Refining from R8 to R16 lowers the reference speed error from about 4.6% to 2.08%; this is a promising two-level trend, not a full grid-convergence result. The maximum logged mass drift is 0.00382% across the gradient cases and 0.00201% in the zero-gradient control, all below the 0.1% limit (exact values are in `results/p1e_thermal_coupled_summary.csv`). In the zero-gradient control, the centroid drift is negligible; the peak dimensionless velocity $2.71\times10^{-4}$ is residual capillary/interfacial numerical flow. Under the positive gradient, peak speed is $2.15\times10^{-2}$, roughly 80 times larger.

The temperature field remains close to the initial linear profile along the interface (RMS deviation about $1.5\times10^{-3}$). This is a numerical consistency check for this chosen diffusivity, not a material calibration.

## Earlier incomplete run

An earlier R16 run targeted $3t_0$ without an explicit short horizon. It ended before its target, at about $0.478t_0$, with exit code 139 and no final report. The cause was not established. Its raw log is retained. The controlled R16 $0.5t_0$ run completed with exit code 0, so the earlier result is not treated as a numerical pass or as evidence of physical instability.

## Limits and next step

This is a two-resolution, short-horizon numerical demonstration. It does not establish full grid convergence, time-step convergence, physical silver diffusivity, temperature-dependent density/viscosity, phase change, solid skeleton constraints, or free-surface deformation in a composite. Do not interpret its velocity or temperature as an Ag70 prediction.

Next, retain the evolving-temperature tracer but add a time-step sensitivity check and a temperature-diffusion sensitivity sweep; then add calibrated Ag properties and only later test a reconstructed local melt-pool geometry. Real Ag70 flow and crater predictions remain out of scope until interfacial properties and the approximately 1.55 nm throat continuum-limit issue are resolved.

## Reproduction

See `run_p1e_validation.sh`. The detailed logs and CSV summary are in `results/`.
