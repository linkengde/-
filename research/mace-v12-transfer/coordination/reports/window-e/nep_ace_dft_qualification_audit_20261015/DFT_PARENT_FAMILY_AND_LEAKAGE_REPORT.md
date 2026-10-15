# DFT parent-family and split leakage audit

## Deduplication and family limits

The source scan contains 554 frames in 247 files, but only 120 unique ordered geometries at 1e-6 Å coordinate tolerance. This is a geometry count, not an independent-parent count. Repeated files across V9–V18 often represent the same geometry copied into successive training releases. Seven k-mesh entries for the 48-atom Ti₃SiC₂ baseline collapse to one geometry. The Ag FCC baseline/volume/shear/displacement cases and the Ti/C displacement cases are related controls derived from a small number of parent crystals. The three 77-atom gaps are a single parent family and have no DFT labels.

The CSV exposes `parent_structure_id` and `parent_lineage_evidence`. Where explicit parent metadata was absent, the ID is explicitly an inferred coarse family or unresolved marker. It must not be treated as a verified parent UUID. The 135-atom full-interface result and Ag-only/TSC-only subsets have explicit common-parent mapping in the E handoff; these are correlated components of one interface calculation.

## Recorded training/validation/test history

Across archived split files, 60 exact geometries have a V9–V18 split role recorded. The role-history occurrences sum to 258 TRAIN assignments, 18 VALID assignments and 20 TEST assignments across successive versions; these are repeated version assignments, not counts of independent structures. V9 records 18 train/2 test and no validation. V11–V16 record 23–35 train, 2 validation and 3 test each. V17 records 43 train/3 validation and no test; V18 records 30 train/3 validation and no test. The same three exact development geometries occur in both V17 and V18. The V11–V16 test files also reuse a small set of exact holdout geometries.

The split files do not provide a robust, independently held-out parent-family test for the requested four-element solid/liquid/interface deployment. Some source records lack parent IDs, so exact-frame non-overlap cannot rule out related slab, registry, or acquisition-family leakage. Training role history from a prior MACE release is provenance only; it does not grant permission to move a frozen test family into training.

## Required split policy

Keep current TRAIN/DEV/TEST roles frozen. For a future NEP/ACE dataset, assign whole parent structures and whole physical scenarios before conversion:

- TRAIN: bulk controls, selected displacement families, one clean-surface family, selected interface registry families, and decorrelated Ag-liquid snapshots.
- DEV: distinct surface termination/registry and distinct liquid-contact parent families, never adjacent frames from a TRAIN trajectory.
- TEST: a newly generated independent interface registry/termination and independent liquid-solid contact or coexistence parent family, sealed from training and hyperparameter decisions.

Every derived configuration must inherit the parent-family ID and role. Different k meshes or smearing settings for one geometry, snapshots from one short trajectory, and gap points from one scan are grouped and kept in one split.

`INDEPENDENT_TEST_READY = FAIL`. The V17/V18 development set is reused, V17/V18 have no test set, and current parent-family lineage is incomplete.
