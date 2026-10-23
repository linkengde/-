# Ti3SiC2 surface-family source and freeze audit

Audit date: 2026-10-23 (UTC)
Working branch: `codex/tsc-surface-dev-prismatic-tangential-20261022`
Audited baseline: `aab06ebe78d5f1abf57e4ff2ec24321fbf251419`

## Decision

`INDEPENDENT_CRYSTALLOGRAPHIC_DEV_PARENT = NOT_AVAILABLE`
`INDEPENDENT_CRYSTALLOGRAPHIC_TEST_PARENT = NOT_AVAILABLE`
`EXISTING_GEOMETRY_FAMILY_CANDIDATES = ALREADY_SCORED`
`NEW_UNSCORED_SAME_SOURCE_FACET_CANDIDATE = AVAILABLE`
`TRAIN_DEV_TEST_ROLES_CHANGED = NO`
`NEW_DFT_OR_MODEL_RUN = NO`

The repository contains one Ti3SiC2 crystallographic CIF: COD 9009647, SHA-256 `6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`. The prior source audit also reviewed the reconstructed Jeitschko–Nowotny (1967) and Rawn (2000) records. They share the same P6_3/mmc prototype, Wyckoff sites, composition, and layer sequence; Rawn is nearly coincident with COD. Those records have no DFT force labels and do not provide an independent crystallographic parent family for a blind test. The Zhang et al. (2026) interface source supplies scalar interface-spacing/adhesion values, not a complete atomic structure with per-atom DFT forces.

The source audit's earlier COD and Materials Project requests were blocked by the environment's 403 tunnel response. A read-only network probe on 2026-10-23 again returned `CONNECT tunnel failed, response 403`; no external file was downloaded. Consequently, a strict source-independent DEV and TEST split cannot be frozen from currently verified sources.

## Existing geometry families and leakage status

| Family | Existing fixed-geometry labels | Source root | Prior use | Eligible as blind DEV/TEST now? |
|---|---:|---|---|---|
| Basal Ti-terminated `(0001)` | Reference + paired normal response | COD 9009647 | C-reviewed TRAIN force-label candidates | No; retain with its TRAIN family |
| Basal Si-terminated `(0001)` | Reference + paired normal response | COD 9009647 | C-reviewed TRAIN force-label candidates | No; retain with its TRAIN family |
| Prismatic `(10-10)` | 96-atom reference + Ti normal response; both SCF-converged | COD 9009647 | Previously compared to TSCBOP | No; already scored, same source root |
| Prismatic `(11-20)` | 120-atom reference + Ti normal response + Ti tangential response; all SCF-converged | COD 9009647 | Previously compared to TSCBOP | No; already scored, same source root |

The prismatic slabs are separate geometric families from the basal slabs and from each other. They can support a **facet-held-out diagnostic** if that weaker claim is explicitly adopted, but they are neither an independent crystallographic-source test nor blind now because their results have already been examined. Their source/geometry hashes and prior status are indexed in `SURFACE_FAMILY_SOURCE_REGISTER.csv`.

After the source scan, a separate geometry-only candidate was generated from the same verified COD CIF: a symmetric Ti-terminated pyramidal `(10-11)` slab with 69 atoms plus one inversion-paired 0.02 A terminal-Ti response. Its exact inputs and geometry checks are in [the candidate package](../tsc_surface_geometry_dev_candidate_10m11_20261023/GEOMETRY_AUDIT.md). The pair is not yet scored and has no DFT labels. It is suitable only as a **same-source facet-held-out DEV candidate**, subject to later method/resource review; it cannot satisfy strict source-independent DEV or TEST.

## Freeze outcome and next valid action

No candidate is promoted to DEV or TEST in this audit. Existing roles remain unchanged, and no TEST coordinates or labels were accessed.

For strict independent-parent validation, the minimum missing item is a traceable Ti3SiC2 parent structure from an independently accessible source, with source metadata and license/provenance, followed by two disjoint surface-family descendants: one DEV family and one evaluator-controlled TEST family. The evaluator must retain the TEST structure identity and labels outside E's workspace.

If the project accepts a geometry-family holdout instead of an independent crystallographic-source claim, the new `(10-11)` pair is available as a DEV candidate. Its reference and deterministic response descendants must stay together as one family. A second, different facet must be generated and scored only by the independent evaluator for TEST. Both results must be described as same-source facet generalization, not source-independent validation. The existing `(10-10)` and `(11-20)` families are ineligible for a new blind split because they have already been scored.

The source-audit package itself only organizes existing evidence. The separate `(10-11)` package contains geometry only. No data role was changed, no TEST identity or label was accessed, and no DFT/MD/TTM or fitting was run; no potential was modified. A fresh process check found no GPAW/MPI job in this Codex workspace; ASE and GPAW are not installed in its base environment. ASE and Pymatgen were installed only in `/tmp/e-surface-ase-venv` for deterministic geometry construction, not for DFT.

## Evidence checked

- `tsc_dev_source_audit_20261020/DEV_SOURCE_AUDIT.md` and `.json`; their local SHA manifest verified.
- `tsc_dev_parent_source_followup_20261021/INDEPENDENT_DEV_SOURCE_DECISION.md` and `SOURCE_CANDIDATE_REGISTER.csv`; both SHA entries verified.
- Basal label qualification `tsc_surface_train_label_qualification_20261020/SURFACE_LABEL_QA.json`.
- Surface input/result audits for prismatic `(10-10)`, prismatic `(11-20)`, and its tangential response.
- Current repository file scan for Ti3SiC2 CIFs, returning only the COD 9009647 CIF.
- Read-only external-source probe; blocked by the network tunnel (HTTP status unavailable, curl exit due to 403).

The exact inspected file digests are recorded in `SURFACE_FAMILY_SOURCE_REGISTER.csv`. No source files or prior reports were modified.
