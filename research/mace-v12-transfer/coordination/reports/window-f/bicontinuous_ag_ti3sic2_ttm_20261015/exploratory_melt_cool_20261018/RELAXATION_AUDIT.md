# 39,486-atom exploratory short-relaxation audit

## Identity and numerical checks

- Atom count: 39,486; IDs and types preserved in final data: **True**. PBC: `f p p`; cell: 110.000000 × 98.085134 × 100.920749 Å.
- Finite force vectors: **True**. Initial Fmax/mean: 3.8003/1.2371 eV/Å. Final Fmax/P95/P99/mean: 2.4501/0.1619/0.3836/0.0514 eV/Å.
- Maximum minimum-image displacement: 5.5234 Å (atom 39188); P95/P99: 2.7670/3.7391 Å; atoms displaced >1 Å: 8140.

## Connectivity at archived geometry cutoffs

- Ag, cutoff 3.05 Å: before 1 component(s), largest 1.000000; after 21 component(s), largest 0.997083.
- TSC, cutoff 3.30 Å: before 1 component(s), largest 1.000000; after 1 component(s), largest 1.000000.

## Pair-distance minima

| pair | before (Å) | after (Å) | count after <2.0 Å |
|---|---:|---:|---:|
| Ti-C | 2.0916 | 1.7292 | 1668 |
| Ti-Si | 2.6840 | 2.3215 | 0 |
| Ag-Ti | 2.5008 | 2.0271 | 0 |
| Ag-Si | 2.5016 | 2.1657 | 0 |
| Ag-C | 2.5091 | 1.9060 | 324 |
| Ag-Ag | 2.8834 | 2.2610 | 0 |

## Decision

The run reached its 500-iteration cap; it did not meet the requested 0.5 eV/Å Fmax tolerance. The minimization trace includes a large transient force spike. Review `MINIMIZATION_STEP_SUMMARY.csv` and `HIGH_FORCE_POST_RELAX.csv`; this result alone does not authorize heating. The next decision is based on whether the final high-force atoms and short contacts are isolated surface/interface rearrangements and whether both networks remain connected.
