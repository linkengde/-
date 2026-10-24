# Liquid Ag density anchor and geometry copies — 2026-10-24

Abdullaev et al. (2020) provide a directly retrieved, CC BY 3.0 primary paper on pure-silver density. Its liquid fit is valid from the 1234.93 K melting point through 1631 K:

```text
rho_m(T) = 9327.1 - 0.9956 * (T - 1234.93) kg m^-3
```

The fit gives 9.312096308 g/cm³ at 1250 K and 8.963636308 g/cm³ at 1600 K. At 1250 K, Table 1 reports density errors of 0.25% at the melting point and 0.26% at 1300 K; the paper does not state a coverage/confidence convention for those percentages. They are neighboring tabulated error estimates, not a 95% interval or a separately measured 1250 K point. At 1600 K the table reports 0.34%, with the same unspecified coverage convention.

The folder `candidates/` contains 12 separate geometry-only copies of the original 1250 K and 1600 K liquid snapshots, at density ratios 0.97, 1.00, and 1.03 to the paper's temperature-specific fitted density, with two seed lineages each. Their cells and positions were scaled isotropically while preserving fractional coordinates. Original candidates remain unchanged. The ±3% ratios are a sensitivity grid, not source uncertainty.

The 1050 K snapshots are intentionally omitted from this rescaled set: 1050 K is below the melting point and outside the paper's stated liquid-fit range. They remain only in the earlier geometry-preparation package with their model-derived density annotation.

All 12 copies pass atom-count, finite-coordinate, periodic minimum-distance, target-density, fractional-coordinate-preservation, unique-file-hash, and no-DFT-label checks in `rescale_candidates.py` and `rescaled_candidate_manifest.json`. They remain outside datasets and train/dev/test splits. No DFT was run or authorized by this geometry-only step; independent solid-Ag source/role review, Ti₃SiC₂ source/role policy, frozen DFT method/energy target, and V20 owner completion remain separate gates.

## Reproduction

From the repository root:

```bash
/workspace/.venvs/mace/bin/python research/mace-v12-transfer/coordination/reports/window-c/abdullaev2020_liquid_ag_density_rescaled_20261024/rescale_candidates.py
```

The source PDF is not redistributed here. Retrieved IOP PDF SHA256: `cd37e2b91be0369ebe1b9d72ffb7377a28dfc74171ed0dd3cd22892f95e23c3a`.
