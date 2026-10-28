# 2× short-run profiling

The cProfile run continued the supplied 0.35 ps checkpoint to 0.36 ps. It took 290.12 s under profiler instrumentation (baseline unprofiled wall time: 291.06 s), with peak RSS 784,272 KiB. The main cumulative costs were:

| Function | Calls | Cumulative time |
| --- | ---: | ---: |
| QSGSEnthalpyGrid.run | 1 | 286.64 s |
| advance | 990 | 284.27 s |
| diffusion_power_density_W_m3 | 990 | 128.59 s |
| state_from_enthalpy_density | 992 | 100.04 s |
| conductivity_W_m_K | 991 | 36.35 s |
| Ag conductivity_W_m_K | 991 | 33.96 s |
| numpy.isclose | 1000 | 16.98 s |
| numpy.roll | 7924 | 16.84 s |
| numpy.interp | 992 | 15.94 s |

The isclose scan was removed from the common Case A path by precomputing eligible indices. The same interval was rerun from the same checkpoint; all state fields and the energy ledger matched the unmodified baseline bit-for-bit. The observed wall time decreased from 287.850 s to 271.030 s, while peak RSS changed from 781,708 to 781,580 KiB.

Raw cProfile output: baseline_0p01ps.cprof; human-readable call table: baseline_0p01ps.top_functions.txt. Measurements and both baseline/profiled result summaries are in the adjacent JSON files and run directories.
