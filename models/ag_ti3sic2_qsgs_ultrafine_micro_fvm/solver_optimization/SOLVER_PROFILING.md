# SOLVER_PROFILING

## What was measured

Profiling ran on 4 allocated CPU cores, a 32 GiB memory limit, and no GPU. The production source case is the existing seed 1105 QGS geometry and localized numerical volumetric pulse. The 50³ cycle used `dt-factor=0.25`; 100³ and the 150³ short stage used `dt-factor=1`. The 100³ reference wall-time run had a 15-minute cap. All successful steps flush a JSONL record; checkpoints are atomically replaced after each accepted step.

The profiler separates property evaluation, face conductance, face flux/residual work, matrix-vector products, forward/inverse FFT calls, PCG, Newton line search, and per-Newton wall time. These timings are nested: do not add them together. Process CPU time, wall time, and peak resident memory are logged separately.

## 100³ complete EW-safe profile

The Numba / `ew_safe` cycle completed 88 steps in 553.43 s wall, 1157.06 CPU-s, and 512.3 MiB peak RSS. Average process CPU use was 2.09 cores. The accepted-step solver-time sum was 307.01 s; the remaining wall includes JIT startup, QGS remapping, per-step logging, compressed checkpoint writes, and other driver overhead.

| Profile component | Accumulated time (s) | Interpretation |
|---|---:|---|
| PCG calls | 185.25 | Dominant measured solve component |
| Newton line searches | 96.03 | Includes residual evaluation; strongly affected by the narrow 1 K mushy band |
| FFT forward transforms | 33.94 | Part of the FFT preconditioner and included in PCG wall |
| FFT inverse transforms | 36.75 | Part of the FFT preconditioner and included in PCG wall |
| Matrix-vector calls | 40.63 | Repeated stencil applications, overlapping with PCG |
| Face conductance | 25.00 | Conductivity/interface face coefficients |
| Material property updates | 20.01 | Fused Numba material fields |
| Face flux work | 16.00 | Residual/energy flux evaluations |
| Preconditioner spectrum setup | 2.77 | FFT spectrum assembly |

These inclusive timings overlap. Line search and PCG account for the main algorithmic cost. The successful run also recorded 24 rejected trials, mostly Newton-limit or line-search failures in the melt band. `retry_log.jsonl` retains each timestep, trial dt, reduced dt, and failure reason.

## Resolution and thread observations

- 50³ QGS solver-step time: reference 62.91 s; Numba + `ew_safe` 15.98 s, a 3.94× reduction on the same 120-step cycle and `dt-factor=0.25`.
- On the fixed-tolerance Numba 50³ cycle, 1/2/4 threads took 44.39/39.02/32.69 s of summed step-solver time. Four threads were fastest, but speedup was 1.36× over one thread, not linear. Whole-process CPU totals were 54.15/82.14/114.95 CPU-s and process wall 54.42/48.80/42.38 s; this suggests useful parallelism with overhead and imperfect scaling.
- 100³ `ew_safe` reached the end of the 1 ns heat and 2 μs cooling windows in 553 s. Heating took about 389 s, cooling about 164 s. Cooling advances much faster because most cells are near 300 K and the source is off.
- 150³ took 172 s wall for 16 accepted steps, reaching 0.715 ns and approaching melt onset. Solver step time summed to 54.63 s; checkpoint compression and flushes contribute significant additional wall time. Peak RSS was 1.09 GiB.

The actual 150³ step sequence became nonlinear near melt onset: its latest dt decreased to 2.34 ps after five retries. The test ended at the planned short-step limit and does not support a reliable full-cycle runtime extrapolation.
