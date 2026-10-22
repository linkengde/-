# Ti₃SiC₂ independent parent source scan

Date: 2026-10-22
Scope: read-only source search to determine whether a genuinely new, traceable Ti₃SiC₂ crystallographic parent is available for independent validation. No structure was generated or changed; no data roles were changed; no DFT, MD, relaxation, or fitting was run.

## Repository evidence

Scanned all CIF paths and their contents in these relevant remote branch snapshots:

| Branch | Commit | Ti₃SiC₂ CIF records found |
|---|---|---|
| `cloud-c-three-paper-crystallography-interface-audit-20261011` | `e999abccf69f717f29d255377b006d4d4a258266` | Three records listed below |
| `cloud-c-v19-minimal-dev-sources-20261011` | `b55c576dd83a248e69caee8c6fbb44a4d7cccb59` | COD 9009647 only |
| `cloud-c-tisc-dft-extension-plan-20261014` | `0a8585e03522b88075e5f62770a8b5409031bf69` | COD 9009647 only |
| `cloud-f-tscbop-refit` | `a3061d93d7d0bc9fc2d871fdb4ffad77f93cf889` | COD 9009647 only |
| `cloud-e-b1-dft` | `88c5407ad717af0a1c410fb505c405a23e559246` | COD 9009647 only |

The three source records present in the crystallography audit branch are:

| Record | SHA-256 | Source meaning |
|---|---|---|
| `periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif` | `6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3` | COD 9009647, citing Kisi et al. (1998) |
| `coordination/reports/window-c/v19_three_paper_crystallography_interface_audit_20261011/Ti3SiC2_jeitschko_nowotny_1967_reconstructed.cif` | `23ee5c4715f30b31dc41843fec88d6888a169b98a4d627628694876ca920f79a` | Reconstructed from Jeitschko–Nowotny (1967) |
| `coordination/reports/window-c/v19_three_paper_crystallography_interface_audit_20261011/Ti3SiC2_rawn_2000_reconstructed.cif` | `39db016144bf735f12dd20646060e049598fc32f7eaffa5cd9ad423d3ae02ce0` | Reconstructed from Rawn et al. (2000) |

The C crystallography audit (`cloud-c-three-paper-crystallography-interface-audit-20261011`, commit above) reports that these records share the Ti₃SiC₂ `P6₃/mmc` prototype, Wyckoff sites, composition and layer sequence. Rawn is nearly coincident with COD; the Jeitschko reconstruction differs more in its C internal coordinate but remains in the same crystallographic parent family. These reconstructions have no DFT energy/force labels. They therefore do not close the strict independent-parent blind-validation gap.

The GitHub branch scan covered complete CIF listings in the five named snapshots. Repeated occurrences of the COD file across branches are the same bytes and SHA; they are not independent sources. The current E branch contains no additional Ti₃SiC₂ CIF beyond that COD file.

## External source retrieval attempt

Read-only requests were attempted against the Crystallography Open Database search endpoint and the Materials Project, OQMD, AFLOW and NOMAD OPTIMADE/API endpoints for reduced composition Ti₃SiC₂. All requests failed before an HTTP response with `curl: (56) CONNECT tunnel failed, response 403`. No external record was downloaded or used. This is an access restriction of the current environment, not evidence that these databases lack a record.

## Decision

`NEW_INDEPENDENT_CRYSTALLOGRAPHIC_PARENT_FOUND=NO`; `STRICT_INDEPENDENT_DEV_READY=NO`; `SAME_SOURCE_CROSS_FACET_DIAGNOSTIC=AVAILABLE`; `DATA_ROLES_CHANGED=NO`.

Do not launch additional DFT on the COD, Jeitschko, or Rawn reconstructions solely to label them independent. The completed `(10-10)` and `(11-20)` slabs remain useful same-source, cross-facet diagnostics. To proceed with strict independent-parent validation, the next required input is a retrievable, provenance-complete structure source not derived from those three records, or an explicitly revised validation criterion that permits a geometry-family holdout while retaining its same-source limitation.

No candidate was promoted to TRAIN, DEV, or TEST. No TEST structure or label was accessed during this scan.
