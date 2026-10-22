# Cross-facet synthesis of frozen TSCBOP force diagnostics

Date: 2026-10-22
Scope: offline aggregation of already completed DFT labels and LAMMPS `run 0` forces for two prismatic Ti₃SiC₂ slab families. No DFT, LAMMPS, relaxation, fitting, or data-role changes were performed.

## Finding

The unchanged TSCBOP shows a large absolute force mismatch on both sampled prismatic facets. The full-slab vector RMSE is 1.9152 eV/Å on the 96-atom `(10-10)` slab and 1.9452 eV/Å on the 120-atom `(11-20)` slab. Ti and C dominate those errors; Si is much lower in both. A paired terminal-Ti displacement is tracked more closely in *differential force response* (vector RMSE 0.01348 and 0.01117 eV/Å, respectively) than in absolute force, but that measures only one prescribed displacement mode per slab.

This consistency across two orientations supports a **same-source, cross-facet diagnostic**: the large absolute Ti/C force mismatch is not confined to one sampled prismatic orientation. It does not establish a universal surface error or independent-parent generalization. Both slabs descend from COD 9009647 and share its crystallographic source; neither slab is an independent source-family blind test.

## Provenance and locked model

Both comparisons use the same original TSCBOP parameter file, SHA-256 `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`, with LAMMPS `pair_style tersoff`, element mapping `Ti Si C`, `units metal`, and `run 0`. The source CIF is `Ti3SiC2_COD_9009647.cif`, SHA-256 `6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`.

| Facet family | Atoms / composition | Reference input SHA-256 | Response input SHA-256 | DFT k-grid | Geometry status |
|---|---|---|---|---|---|
| `(10-10)`, 8 layers | 96; Ti₄₈Si₁₆C₃₂ | `60096a122270e57ac2fb1cbbc2c0d1b087eab12127bcee087f30d244dc8900e1` | `710b0abdd20f0d59103ebfe061ad65c834bae7a4b2f5a300453ba1527e63e030` | MP 12×2×1 | Frozen, unrelaxed slab; paired terminal-Ti response |
| `(11-20)`, 10 layers | 120; Ti₆₀Si₂₀C₄₀ | `fd5ff5c1cae1490bcbde9ffd08f329c7ec8cac687555f29568d725f4230cc5c3` | `6f91573bcb4a0feb7e9f05329f6eca422d3e7aac723763e961b52c7b83ddf9b6` | MP 7×2×1 | Frozen, unrelaxed slab; inversion-paired terminal-Ti response |

Both used PBE/PAW, PW 500 eV, Fermi–Dirac 0.10 eV, fixed coordinates and cell, and the same reported SCF thresholds/software. The different meshes were selected to maintain reciprocal-space sampling for different in-plane cells; this is not proof of cross-facet k-point convergence. The cell sizes, layer counts, terminations, and atom counts differ, so their force metrics are descriptive comparisons rather than a controlled one-variable experiment.

## Force comparison

`CROSS_FACET_FORCE_SUMMARY.csv` contains the copied per-facet, per-element metrics and the paired-response metrics. Absolute errors are `F_TSCBOP − F_DFT`; response error is `(F_TSCBOP,RESP−F_TSCBOP,REF) − (F_DFT,RESP−F_DFT,REF)`. Vector RMSE is the square root of the mean per-atom squared vector error. The two geometries in each facet are a paired perturbation family, not independent structures.

| Facet | Geometry | All-atom vector RMSE (eV/Å) | Ti (eV/Å) | Si (eV/Å) | C (eV/Å) |
|---|---|---:|---:|---:|---:|
| `(10-10)` | REF | 1.9152 | 2.1284 | 0.0801 | 2.0507 |
| `(10-10)` | Ti response | 1.9130 | 2.1276 | 0.0805 | 2.0459 |
| `(11-20)` | REF | 1.9452 | 2.1551 | 0.2165 | 2.0884 |
| `(11-20)` | Ti response | 1.9448 | 2.1552 | 0.2166 | 2.0871 |

The differential-response vector RMSE is 0.01348 eV/Å across 96 atoms and 0.04101 eV/Å within 3.3 Å of the displaced Ti atoms on `(10-10)`. On `(11-20)`, it is 0.01117 eV/Å across 120 atoms and 0.02678 eV/Å in the corresponding 3.3 Å neighborhood. Those values indicate that the tested finite response is much closer than the absolute force field on these two specific geometries. They do not validate other response modes, relaxed surfaces, finite-temperature configurations, or elastic response.

The local `(11-20)` environment audit also found high errors on near-surface Ti and subsurface C without an obvious anomalous short contact or simple low-coordination explanation. This is consistent with a broader same-geometry-family mismatch, but it does not identify a particular potential parameter as the cause.

## Decision and next step

- `TSCBOP_ABSOLUTE_FORCE_ACCURACY_ON_SAMPLED_PRISMATIC_SLABS=POOR`.
- `SINGLE_TERMINAL_TI_DIFFERENTIAL_RESPONSE_ON_TWO_SAME_SOURCE_FACETS=REASONABLY_CLOSE_WITHIN_TESTED_MODES_ONLY`.
- `INDEPENDENT_CRYSTALLOGRAPHIC_PARENT_VALIDATION=UNAVAILABLE`.
- No parameter fitting, candidate potential, production-readiness, or TRAIN/DEV/TEST reassignment is supported by this synthesis.

The useful next evidence is a genuinely new, fully traceable Ti₃SiC₂ parent or surface dataset, with its source and data-role policy frozen before any force scoring. A third facet derived from the same COD parent may broaden geometry coverage but cannot resolve the independent-source gap. Until a new source is available, preserve these two facet pairs as diagnostic evidence and do not call them an independent blind validation.

## Source reports

- `(10-10)` comparison summary SHA-256: `97d72f5f4d0d1e3eba9ccbc4df5ccfd37f06dbd27724ae6f3c4bdfb8efb04e08` (`analysis/DEV_TSCBOP_RUN0_PAIR_20261021/SUMMARY.json` in `tsc_surface_dev_prismatic_20261020`).
- `(11-20)` comparison summary SHA-256: `9355edb4d63d5eaf0048b010316685acf35623139ad87cf916e8b02a85b88aa8` (`analysis/TSCBOP_RUN0_PAIR_20261022/SUMMARY.json`).
- Source audit: `tsc_dev_parent_source_followup_20261021/INDEPENDENT_DEV_SOURCE_DECISION.md`.
