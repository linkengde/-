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

## Endpoint clusters

At 1.60 ps, the 874 candidate Ag atoms form 41 connected components using ≤3.5 Å Ag–Ag edges, periodic x/y, and open z. There is no minimum-size filter. The largest component contains 693 atoms (79.3% of candidates), with periodicity-aware dimensions 77.262 × 72.421 × 26.633 Å and convex-hull envelope 52,326.7 Å³. These are geometry statistics of the selected candidate atoms, not a validated liquid volume. The same potential and threshold limitations described in the neighboring Round7 audit still apply.
