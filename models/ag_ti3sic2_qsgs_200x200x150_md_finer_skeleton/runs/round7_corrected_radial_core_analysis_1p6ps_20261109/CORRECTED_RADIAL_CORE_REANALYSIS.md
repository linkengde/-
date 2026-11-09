# Corrected Round7 spatial-mask and melt-candidate counts

The earlier 1.600 ps screen's source mask calculated `dy` but omitted it from the predicate. As a result, its 6,448-atom subset was an x–z band spanning y. This separate analysis reuses the same five saved trajectory snapshots and adds the missing periodic y contribution while preserving the prior dimensionless exponent cutoff of 2:

```text
core = initial_newton_layer & Ag
       & ((dx/30 Å)^2 + (dy/30 Å)^2 + depth/20 Å <= 2)
```

That mask follows the Round7 TTM spatial source form `exp(-((dx²+dy²)/R² + depth/λ))`, with R=30 Å and λ=20 Å. Its cutoff corresponds to source weight ≥ exp(-2) of the central surface value. It applies to initial positions; all three indicators use the same selected Ag IDs across the five snapshots.

## Exact thresholds and joint expression

- Q6: five-frame arithmetic mean of Ag-only Steinhardt q6, cutoff radius 4.0 Å, `< 0.4081069282833649`.
- PTM: five-frame fraction of non-FCC classifications `>= 0.8` (at least 4/5 snapshots); FCC is the sole rejected PTM class.
- Displacement: single-origin drift-corrected squared net displacement over 0.60–1.60 ps `>= 0.5171679303519052 Å²`.
- Joint expression: `core & (q6_mean < 0.4081069282833649) & (nonfcc_fraction >= 0.8) & (msd_i >= 0.5171679303519052)`.

Thresholds are unchanged from the prior screen and were inherited from Round5 references; they were not recalibrated on Round7. The displacement algorithm and temporal sampling are also unchanged.

## Counts from the corrected subset

| Screen | Passing atoms |
|---|---:|
| Corrected spatial subset | 1,787 |
| Q6 | 902 |
| PTM non-FCC fraction | 1,773 |
| Displacement | 1,748 |
| Q6 ∩ PTM | 895 |
| Q6 ∩ displacement | 881 |
| PTM ∩ displacement | 1,734 |
| All three | 874 |
| Fail at least one | 913 |

Failure patterns among the 913 that fail: Q6 only 860; displacement only 21; Q6 and displacement 18; Q6 and PTM 7; PTM only 7; all other combinations 0. The stored joint flags equal the recomputed AND for all 1,787 rows.

## Why PTM is nearly always non-FCC

The PTM breakdown for these same 1,787 Ag atoms is:

| Frame (ps) | FCC | HCP | BCC | Other |
|---:|---:|---:|---:|---:|
| 0.60 | 58 | 0 | 0 | 1,729 |
| 0.85 | 15 | 0 | 0 | 1,772 |
| 1.10 | 3 | 0 | 0 | 1,784 |
| 1.35 | 0 | 0 | 0 | 1,787 |
| 1.60 | 0 | 0 | 0 | 1,787 |

Thus the non-FCC persistence pass is almost entirely PTM **Other**, not conversion to HCP or BCC. The script does not save RMSD, so the data cannot distinguish thermal disorder from insufficient Ag neighbors or a local environment rejected by PTM's default normalized RMSD cutoff (0.1). At 1.60 ps, the Ag-only neighbor count within 4.0 Å has mean 9.61 among the 874 joint candidates; 669/874 have fewer than 12 neighbors, 408/874 fewer than 10, and 157/874 fewer than 8. The noncandidate subset is even more under-coordinated (mean 5.97). These counts show that coordination is a material confounder to check, not proof that PTM Other candidates are solid or liquid. Full per-frame class and coordination diagnostics are in `ptm_coordination_diagnostic.json` and can be regenerated with `diagnose_ptm_coordination.py`.

## Endpoint clusters

At 1.60 ps, the 874 candidate Ag atoms form 41 connected components using ≤3.5 Å Ag–Ag edges, periodic x/y, and open z. There is no minimum-size filter. The largest component contains 693 atoms (79.3% of candidates), with periodicity-aware dimensions 77.262 × 72.421 × 26.633 Å and convex-hull envelope 52,326.7 Å³. These are geometry statistics of the selected candidate atoms, not a validated liquid volume. The same potential and threshold limitations described in the neighboring Round7 audit still apply.

## Controlled comparison to the legacy spatial mask

Both columns use the same 0.60–1.60 ps frames, Q6/PTM/displacement thresholds, and 3.5 Å clustering. Only the spatial prefilter changes; the corrected mask adds the missing periodic y term.

| Statistic | Legacy x–z mask | Corrected radial mask |
|---|---:|---:|
| Selected Newton-layer Ag | 6,448 | 1,787 |
| Q6 passing | 3,078 | 902 |
| PTM non-FCC fraction passing | 6,418 | 1,773 |
| Displacement passing | 6,304 | 1,748 |
| Q6 ∩ PTM | 3,067 | 895 |
| Q6 ∩ displacement | 3,000 | 881 |
| PTM ∩ displacement | 6,276 | 1,734 |
| Joint candidates | 2,989 | 874 |
| Candidate components | 148 | 41 |
| Largest component | 1,683 | 693 |
| Largest-component dimensions (Å) | 84.32 × 146.96 × 40.21 | 77.26 × 72.42 × 26.63 |
| Largest-component convex-hull envelope (Å³) | 201,755.8 | 52,326.7 |

The largest component's y span falls from 146.96 Å to 72.42 Å under the radial mask, removing the spurious box-scale extension associated with omitting y. The corrected set is still a temporary candidate classification, not a validated melt pool.
