# TSC69 force outlier and data-role review

## Scope

Read-only review of the converged fixed-geometry result for `TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_REF`. No geometry, DFT label, or role assignment was changed, and no new calculation was run.

## Maximum-force environment

The largest force is on source atom ID 63 (Ti): `|F| = 2.625585 eV/Å`, with vector `[0.987880, 1.729197, 1.711043] eV/Å`. It lies 1.582 Å inward from the nearest outermost slab plane along the audited surface normal. Its inversion-related Ti atom (source ID 3) has `|F| = 2.625482 eV/Å` and the opposite vector to within `1.1×10⁻⁴ eV/Å`. This paired response is consistent with the audited symmetric Ti/Ti slab and does not look like an isolated atom-mapping or serialization error.

Nearest neighbors of Ti ID 63 are C ID 69 and C ID 60 at 2.088203 Å, Si IDs 65 and 56 at 2.681042 Å, and Ti ID 49 at 2.969608 Å. These distances are not anomalously shorter than the input geometry's audited global minimum of 2.088203 Å (source IDs 52 and 67). The slab was generated as an as-cleaved, unrelaxed neutral-atom surface; the high forces are therefore consistent with an unrelaxed surface response. This check cannot establish that the termination is mechanically stable or that the coordinates are near a relaxed geometry.

## Data role decision

- Current role remains `UNASSIGNED_SAME_SOURCE_FACET_DIAGNOSTIC`.
- The candidate manifest identifies source family `surface_pyramidal_10m11_cod9009647_symm_ti_69atom`, derived from COD 9009647.
- It is a new, previously unscored facet geometry but shares its crystallographic source root with other data. It is **not** an independent-source DEV or TEST family.
- This result has not been assigned to TRAIN, DEV, or TEST. Do not use it for fitting or claim independent validation based on this point.

## Conclusion

`FORCE_OUTLIER_GEOMETRY_QA=CONSISTENT_WITH_SYMMETRIC_UNRELAXED_SURFACE; NO_SHORT_CONTACT_OR_MAPPING_ANOMALY_FOUND`.

`DATA_ROLE=UNASSIGNED`; no role change is justified by this single fixed-geometry point. The 2.626 eV/Å force is retained as a diagnostic label and is not treated as a convergence failure, since the SCF itself converged.

Evidence: `results/TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_REF/RESULT.json`, `labelled.extxyz`, the geometry candidate manifest, and the corresponding SHA-256 manifests.
