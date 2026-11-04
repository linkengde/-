# P1b/P2 final decision

## Outcome

**Not ready for physical Ag70 melt-pool flow prediction.** The official dynamic benchmark is reproduced, but the project's automatic-curvature VOF comparison fails quantitative velocity convergence. The Ag viscosity record is now uncertainty-audited, yet its measurements stop at 1673 K and wall/interface momentum transfer remains unknown.

| Gate | Evidence | Status |
|---|---|---|
| P0 automatic-curvature static VOF tests | Prior controlled curvature, contact-angle and mass results | Passed for those tests only |
| P1 planar Marangoni/static force checks | Prior controlled dimensionless tests | Passed for those tests only |
| P1b official CLSVOF dynamic droplet | 0.340% fine-grid speed error; max logged mass drift 0.0594% | Passed for reference setup |
| P1b project auto-curvature VOF dynamic comparison | 9.44–10.00% errors on finer grids; no convergence | Not passed |
| Liquid Ag viscosity | Nine NIST points, 1273–1673 K, 95% expanded uncertainty; all publisher-fit residuals within U95 | Audited candidate in measured range only |
| Ag density | Abdullaev experimental correlation, 1234.93–1631 K | Candidate within range; reported error convention remains unspecified |
| Ag surface tension | Two atmosphere-specific Morohoshi fits with distinct ranges | Candidates only in their separate ranges |
| Ag/Ti₃SiC₂ wetting | 72° static datum under cited specimen/atmosphere/temperature | Condition-specific static value only |
| Ag/Ti₃SiC₂ friction/slip | No calibrated dynamic wall data | Unresolved |
| Nanopore continuum validity | ~1.55 nm bottleneck, ~0.317 nm voxels, no scale separation or wall law | Quantitative pore-scale NS not justified |

The next technical gate is to resolve the automatic-curvature dynamic stress discretization and repeat P1b. Separately, obtain or derive a defensible Ag/Ti₃SiC₂ momentum-transfer law. No Ag70 geometry, thermal field, or real melt-pool flow was used in this phase.
