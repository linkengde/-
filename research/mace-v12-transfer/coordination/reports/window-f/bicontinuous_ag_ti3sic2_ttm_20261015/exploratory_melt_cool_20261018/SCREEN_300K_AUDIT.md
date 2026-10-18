# 300 K baseline screen: stopped after 0.2 ps

The test was run on `post_relax.data`, an independent copy of the user-supplied 39,486-atom model. The original model and potential files remain unchanged. Calibration-v0 mapping was preserved. The run used NVT at 300 K, 1 fs timestep, and was interrupted after the 200 fs frame; no heat pulse or TTM was applied.

## Observed force spike

| time | temperature from LAMMPS | Fmax (eV/Å) | mean force | atoms F>5 | atoms F>10 |
|---:|---:|---:|---:|---:|---:|
| 0.000 ps | 300.0 K | 2.450 | 0.051 | 0 | 0 |
| 0.100 ps | 199.0 K | 66.491 | 0.755 | 5 | 4 |
| 0.200 ps | 263.0 K | 67.010 | 0.861 | 20 | 5 |

At 0.1 ps, the maximum force was 66.49 eV/Å on C ID 19219; at 0.2 ps the LAMMPS log reports Fmax=67.01 eV/Å. The 0.1 ps neighborhood contains Ti ID 25627 at 1.809 Å, Ti ID 25628 at 1.914 Å, and Ag ID 38116 at 1.981 Å from that carbon. The same C ID had already reached 50.83 eV/Å near minimization iteration 475. Its source-region label is a geometric `interior` candidate but it is also an initial Ag-interface candidate (nearest Ag 2.851 Å); it is not evidence of a bulk-only force-field test. This recurring local response is the immediate stop condition.

## Skeleton and contact summary

- step 0: TSC graph at 2.70 Å has 1 component(s); Ag graph at 3.05 Å has 21 component(s), largest fraction 0.99708.
- step 100: TSC graph at 2.70 Å has 3 component(s); Ag graph at 3.05 Å has 25 component(s), largest fraction 0.99692.
- step 200: TSC graph at 2.70 Å has 2 component(s); Ag graph at 3.05 Å has 26 component(s), largest fraction 0.99643.

The test was intentionally stopped before equilibration because the repeated localized force spikes are inconsistent with proceeding to a driven melt-pool calculation under this calibration-v0 mapping. The partial trajectory frames (0, 100, 200 fs) are provided for OVITO inspection; they are not a completed 300 K stability test. No heating, cooling, TTM, DFT, or model/potential modification occurred.
