# NONLINEAR_CONVERGENCE

## Policies

The original reference solve uses Newton with a line search and matrix-free PCG preconditioned by an FFT approximation to a constant-coefficient Laplacian. Its default relative PCG and Newton tolerances remain `1e-8`, with 400 PCG and 48 Newton iterations available per accepted trial. No tolerance was loosened for these comparisons.

- `fixed`: reference PCG tolerance on every Newton update.
- `ew`: inexact Newton forcing up to 0.5; retained as a faster but more permissive comparison.
- `ew_safe`: first solve uses 0.05, then a residual-ratio forcing term capped at 0.1 and floored at 0.005 (or the original configured tolerance). Nonlinear acceptance remains `1e-8`.
- `amg`: optional smoothed-aggregation PyAMG preconditioner pilot. The production benchmark remains FFT by default.

When a trial fails, its reason and attempted/reduced dt are flushed to `retry_log.jsonl`; the same accepted state remains recoverable. The driver halves dt and retries. It does not clip temperature or liquid fraction to force acceptance.

## Measured comparisons

| Run | Accepted steps | Newton updates | PCG iterations | Rejected trials | Wall | Normalized energy residual |
|---|---:|---:|---:|---:|---:|---:|
| 50³ reference, fixed | 120 | 458 | 16,291 | 0 | step sum 62.91 s* | 9.51×10⁻¹⁰ |
| 50³ Numba, fixed, 4T | 120 | 458 | 16,293 | 0 | 42.38 s whole run | 9.51×10⁻¹⁰ |
| 50³ Numba, EW-safe, 4T | 120 | 510 | 6,217 | 0 | 25.89 s whole run | 9.98×10⁻¹⁰ |
| 50³ Numba, EW, 4T | 121 | 586 | 6,130 | 1 | 26.18 s whole run | 5.51×10⁻⁹ |
| 100³ Numba, EW, 4T | 128 | 1,666 | 12,282 | 58 | 508.33 s | 3.61×10⁻⁹ |
| 100³ Numba, EW-safe, 4T | 88 | 1,256 | 12,320 | 24 | 553.43 s | 1.43×10⁻⁹ |
| 150³ Numba, EW-safe, 4T | 16 | 105 | 498 | 5 | 172.06 s | 2.24×10⁻¹⁰ |

\*The reference 50³ run's old driver summary did not preserve a reliable whole-process wall value; the per-step sum was persisted and is used for this comparison. The Numba whole-run times include JIT warm-up and checkpoint/log I/O. The like-for-like solver-step sum is 62.91 s vs 15.98 s (3.94×).

The 100³ `ew_safe` cycle reduced Newton updates by about 24.6% and rejected trials from 58 to 24 relative to `ew`, while its PCG count was similar. Retry reasons were 16 Newton-iteration-limit failures and 8 line-search failures for `ew_safe`; the original `ew` run had 57 line-search failures and 1 Newton-iteration-limit failure. EW-safe ran about 9% slower wall-to-wall and changed hot-endpoint Ag-weighted liquid fraction from 3.31147×10⁻⁴ to 3.33011×10⁻⁴ (about 0.56% relative). This is a stability/energy-accounting tradeoff, not a universal speed win. Both runs resolved local Ag liquid fraction 1.0 and returned the system to zero liquid after cooling.

## AMG pilot

PyAMG 5.3.0 was installed under `/tmp` for a single 50³ preheated-state numerical-pulse step. Both cases were warmed up before timing. FFT used 3 Newton updates and 32 PCG iterations, taking 0.096 s. AMG used 2 Newton updates and 2 PCG iterations, but its hierarchy setup took 0.420 s and the step took 0.519 s. Maximum temperature differed by 4.62×10⁻¹¹ K and bulk liquid fraction by 4.48×10⁻¹² relative; the AMG normalized energy residual was 2.66×10⁻¹².

This confirms that AMG can reduce iterations and preserve this test result, but hierarchy assembly/setup dominated the one-step wall time. No full-cycle, 100³, or 150³ AMG result exists. It is not selected for production yet; reuse of a hierarchy and a high-resolution setup-cost study would be needed.

## Approaches not implemented in this stage

An enthalpy-primary-variable or active-set phase solver was not added. For this model, that requires a monotone enthalpy-to-temperature inversion in every residual evaluation plus a Jacobian treatment for temperature-dependent conductivity and mixed-cell latent heat. A quick variable substitution would not be a verified active-set method. The existing 1 K mushy-zone constitutive curve remains the common reference; changing its width would confound solver accuracy and performance.

The current dt controller is convergence-driven: it halves dt after a failed Newton/line-search trial and grows the next accepted dt by 1.5 up to the phase nominal cap. It does not yet estimate local time-discretization error or directly limit liquid-fraction change. Energy residual is measured at every accepted step and in the run summary, but is not an acceptance controller. Consequently the 100³ under-15-minute result is a performance result, not a time-converged solution claim.
