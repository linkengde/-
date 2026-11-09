# Round7 1.600 ps corrected radial source-core reanalysis

This is a separate post-processing reanalysis of the existing Round7 TTM-MD trajectory. It uses the same five snapshots (0.60, 0.85, 1.10, 1.35, 1.60 ps), the same Q6/PTM/displacement thresholds, and the same 3.5 Å cluster rule as the prior screen. No MD steps were rerun and the earlier results remain preserved.

The prior spatial mask omitted its y term. This reanalysis includes the source's actual periodic radial coordinates and uses the same dimensionless exponent cutoff, `((dx/30 Å)^2 + (dy/30 Å)^2 + depth/20 Å) ≤ 2`. The Round7 TTM source weight is proportional to `exp(-((dx²+dy²)/30² + depth/20))`, so this selects the initial Newton-layer Ag atoms whose source weight is at least `exp(-2)` of its center/surface value. The region is fixed from initial coordinates. This is a transparent correction of the old mask, not a new thermal simulation or an independently calibrated melting threshold.

## Reanalysis results

- Corrected source-core subset: 1,787 initial Newton-layer Ag atoms.
- Joint Q6 + PTM + displacement candidates: 874.
- Components at the 1.60 ps endpoint: 41; largest component: 693 atoms (79.3% of candidates).
- Largest component periodicity-aware dimensions: 77.26 × 72.42 × 26.63 Å; convex-hull envelope: 52,326.7 Å³.
- The candidate flag was recomputed from the three exact conditions for every CSV row; zero mismatches.

The per-criterion counts and fail-pattern breakdown are in `CORRECTED_RADIAL_CORE_REANALYSIS.md`. This remains exploratory: the thresholds come from Round5, the Round7 potential is not independently validated for liquid Ag/interface behavior, and five snapshots are coarse temporal sampling. Call these joint-screen candidates, not confirmed liquid atoms.

## Files

- `round7_1p6ps_joint_candidate_3d.png` — corrected-mask visualization.
- `round7_1p6ps_joint_candidates.dump` — all joint candidates at 1.60 ps.
- `round7_1p6ps_largest_candidate_cluster.dump` — largest component at 1.60 ps.
- `round7_1p6ps_joint_screen.json` — machine-readable settings and cluster summary.
- `round7_1p6ps_source_core_ag.csv` — per-atom values for the 1,787-atom corrected subset.
- `reanalyze_corrected_radial_core.py` — analysis script with the radial y term restored.
