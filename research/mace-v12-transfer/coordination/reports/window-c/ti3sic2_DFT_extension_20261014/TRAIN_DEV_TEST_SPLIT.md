# Stage C Ti₃SiC₂ split proposal

## Decision

This is a proposed parent-family split, not a frozen role manifest. E6 reports `role_freeze.json` as empty. No existing label role is changed by this package.

| Set | Planned content | Parent-family rule | Current state |
|---|---|---|---|
| TRAIN | A volume/shear (4), B Si ±0.02 Å (2), C clean Ti/Si surfaces and local response (4) | All ten geometries descend from the same public 48-atom Ti₃SiC₂ baseline proposal (`157c1d…c3e4`). Keep the entire family together. | `TRAIN_CANDIDATE`, pending a separate role freeze. Existing numerical-control labels remain unchanged; if these new labels are admitted, the whole ancestry must be reviewed as one group. |
| DEV | Independent parent reference + deterministic response (2) | Must have a source root and geometry lineage disjoint from the baseline-derived TRAIN family and every previously scored family. | No qualifying source/coordinates identified in the reviewed E6 bundle; both rows are source-pending. |
| TEST | Independent parent reference + deterministic response (2) | Must be disjoint from TRAIN and DEV. An independent evaluator must hold the source geometry and labels; C must not inspect them before the evaluator releases a final score. | No qualifying source identified. No test structure or label was opened or enumerated. |

The 10 generated geometries in `inputs/` have valid atom counts, composition, finite cell/coordinates, TTT PBC and minimum periodic distance ≥1.5 Å. That establishes geometry preflight only. It does not establish a DFT label, role eligibility, relaxed surface stability or a production fit dataset.

## Lineage details

- The four A structures are byte-for-byte copies of the registered volume/shear proposal files, all sharing the 48-atom baseline ancestry.
- The two B structures displace only `local_new_id=7` (Si, `spacegroup_kinds=2`) by ±0.02 Å in z at fixed cell. They share the same baseline family.
- The four C structures are symmetric Ti- and Si-terminated basal slabs cleaved from that same baseline and its c-periodic images. Each reference/response pair stays grouped; the top-layer response is one atom moved outward by +0.02 Å.
- The 4 D rows in `DFT_STRUCTURE_LIST.csv` are planning slots only. Their parent structures are not fabricated from TRAIN geometries. Source IDs, cell, atom count, geometry SHA and k mesh remain blank/TBD until an independently sourced family is reviewed.

## Role safeguards

All roles here are proposals. The existing five numerical-control labels and all E-window archives keep their current state. No labels are copied into a training set. Do not split k-point variants, perturbations, slab images or other descendants from one parent across TRAIN/DEV/TEST. A separate manifest freeze must resolve whether the baseline-derived family can enter TRAIN without changing the existing records' roles. If not, the ten new labels remain on hold or the full family is handled under a separately approved role migration.

