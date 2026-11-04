# P1b/P2 final decision

## Outcome

**Not ready for a physical Ag70 melt-pool flow calculation.** The literature audit improves the property basis, but the official dynamic Marangoni reproduction remains unrun and important material/interface gates remain open.

## Decision table

| Gate | Current evidence | Status |
|---|---|---|
| P1 controlled Marangoni implementation | Planar Stokes sign, speed, timestep and mass tests; static curved-force closure | Passed for those dimensionless cases only |
| P1b official dynamic droplet migration | User confirms the official page exists and identifies its axisymmetric CLSVOF/theory/grid-convergence content; source retrieval from this runner returned HTTP 403 | Not run; no decision |
| Liquid Ag viscosity | Nine primary-paper points and fit, 1273–1673 K; user reports the NIST JSON’s nine values match and expanded-uncertainty fields are present | Candidate inside stated range; numeric NIST uncertainty fields still need local extraction; not validated above it |
| Liquid Ag density | Abdullaev experimental correlation, 1234.93–1631 K | Candidate within range; reported error convention is unspecified |
| Ag surface tension | Two separate Morohoshi atmosphere-specific fits: 1219–1734 K and 1350–1521 K | Candidate only within each separate range; do not combine or extrapolate |
| Ag/Ti₃SiC₂ static wetting | 72° under the cited polished-substrate/vacuum/970 °C condition | Condition-specific static datum only |
| Ag/Ti₃SiC₂ friction and slip | No calibrated data | Unresolved |
| Nanopore continuum validity | ~1.55 nm geometric z bottleneck, ~0.317 nm voxel; no scale separation or wall law | Quantitative pore-scale NS not justified |

## Permitted next work

1. Review and publish the saved restricted-egress draft if appropriate; it adds only `basilisk.fr` and `trc.nist.gov` to the existing custom allowlist. The draft is not active until published.
2. Retrieve and hash the official Basilisk source and revision, then finish P1b with the official method and a carefully matched project-method comparison.
3. Retrieve the NIST ThermoML JSON, preserve its SHA256, and extract the expanded-uncertainty values with their units and coverage convention. The nine-value equality has been externally checked in the current chat as reported by the user, but the raw JSON and uncertainty values are not yet in the repository.
4. Seek condition-matched Ag/Ti₃SiC₂ momentum-transfer data or develop a separately validated interface-calibration study. Do not treat the static contact angle as a friction/slip measurement.
5. Revisit the continuum model only with a documented scale-bridging basis and a justified driving condition. If pore-scale Navier–Stokes remains unsupported, develop a validated coarse-grained closure before flow predictions.

No production solver or original Ag70 geometry was changed. No melt pool was regenerated, no real Ag70 flow was run, and no 2000–2277 K Ag transport parameter was represented as experimentally validated.
