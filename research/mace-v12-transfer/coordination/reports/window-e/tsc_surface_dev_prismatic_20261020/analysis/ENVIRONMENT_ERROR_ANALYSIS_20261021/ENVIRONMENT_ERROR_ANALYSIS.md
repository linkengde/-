# TSCBOP force residual by local environment: prismatic slab pair

Date: 2026-10-21
Scope: post-process the already archived 96-atom fixed-geometry DFT/TSCBOP `run 0` pair only. No new scientific calculation was run, no potential or structure was changed, and data roles remain unchanged.

## Input integrity and method

The parent candidate manifest and the earlier pair-analysis manifest both pass their SHA-256 checks. The REF input hash is `60096a122270e57ac2fb1cbbc2c0d1b087eab12127bcee087f30d244dc8900e1`; its force-comparison CSV hash is `c0476f953a877dfa430fac7d7987f751b25693338dbbd85b41c8c08ff79438e6`. The displaced-point force CSV hash is `5cc7c95ea38e12a360dcbeb492aec1bed1c7a7b539ebb2e2f36a14a5241c6eab`. Source IDs, element symbols, and each CSV coordinate were rechecked against the corresponding REF or displaced input before grouping.

For each atom, the residual is `F_TSCBOP - F_DFT`; vector RMSE is `sqrt(mean_i(|residual_i|^2))`. Depth is the minimum Cartesian z distance to either outermost atomic plane of the REF slab. Descriptive depth bins are `<2.5 Å`, `2.5–5 Å`, and `≥5 Å`. First-shell coordination counts use explicit geometric cutoffs Ti–C `<2.50 Å` and Ti–Si `<3.10 Å`; these counts are descriptive and are not an independently fitted coordination classifier. All atoms remain under PBC with the archived 18 Å vacuum direction.

## Absolute force-error pattern

| Group | N | Vector RMSE (eV/Å) | P95 error norm (eV/Å) | Maximum (eV/Å) |
|---|---:|---:|---:|---:|
| All | 96 | 1.9152 | 2.7762 | 2.8630 |
| Ti | 48 | 2.1284 | 2.8629 | 2.8630 |
| C | 32 | 2.0507 | 2.2589 | 2.2589 |
| Si | 16 | 0.0801 | 0.1964 | 0.1964 |
| Outer, `<2.5 Å` from either surface | 24 | 1.5544 | 2.7131 | 2.7132 |
| Subsurface, `2.5–5 Å` | 24 | 2.0494 | 2.8630 | 2.8630 |
| Inner, `≥5 Å` | 48 | 2.0070 | 2.7762 | 2.7762 |

Ti and C dominate the absolute residual. The error is not confined to the outermost planes: the inner half of this finite slab still has a 2.007 eV/Å vector RMSE. Because the slab is only eight layers and shares the COD 9009647 parent lineage, this is evidence about this geometry, not a clean bulk-versus-surface attribution.

The geometric coordination grouping shows elevated Ti residuals for the most common Ti–C count group (`n=2`, 28 Ti; 2.7395 eV/Å), while the `n=4` group (12 Ti) is lower (0.2570 eV/Å). C atoms with four Ti neighbors (28 atoms) have 2.1758 eV/Å RMSE; the two observed C coordination groups are not balanced independent samples. These correlations are descriptive only: species, depth, local bonding, and slab periodicity covary in this one structure, so they do not isolate a unique Tersoff parameter or prove a causal coordination failure.

## Differential displacement response

For the existing symmetric terminal-Ti ±0.02 Å pair, the differential error is `(F_TSCBOP,RESP−F_TSCBOP,REF) − (F_DFT,RESP−F_DFT,REF)`. Its vector RMSE is 0.01348 eV/Å over 96 atoms and 0.04101 eV/Å for the nine atoms within 3.3 Å of either displaced Ti. By reference-depth bins, response-error vector RMSE is 0.02474 eV/Å at the outer layer, 0.01054 eV/Å in the subsurface bin, and 0.00140 eV/Å in the inner bin.

Thus this one local force-response mode is reproduced much more closely than the absolute forces. It does not cancel or repair the large absolute Ti/C residual and does not establish response accuracy for other displacement modes, terminations, or parent families.

## Decision

`TSCBOP_ABSOLUTE_FORCE_ACCURACY_ON_THIS_SLAB=POOR` and `TSCBOP_SINGLE_Ti_RESPONSE_MODE=REASONABLY_CLOSE_WITHIN_THIS_GEOMETRY_ONLY`. The result supports continued use only as a diagnostic comparison, not a quantitative accuracy or production claim. It is not an independent crystallographic-source DEV/TEST result, and no TRAIN/DEV/TEST role has been assigned. No fitting is initiated by this analysis.
