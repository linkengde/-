# V17 dataset integration recipe for A review

## Fastest archive-supported rank restoration

The 29 directly supported training rows have exact rank3. One exact matched relabel can provide a 30-frame rank4 pool. All13 eligible V16 frame choices are 0–8 and10–13; frame9 does not lift rank. Full identities, formulas and declared settings are in the rank report. Choose based on scientific coverage and cost, not rank alone; no evidence here selects a best force-learning candidate. Recovering an original valid archive is preferable to a new calculation if available. No jobs are authorized or launched by this recipe.

## Alternatives

A: core29 plus one successful exact-geometry relabel ->30 rows/rank4/one job/zero partial rows. B: core29 plus one provisional partial ->30 rows/rank4/zero jobs/one unresolved row, or all inherited ->43 rows/rank4/14 partial. C: relabel all14 partial exact geometries ->43 rows/rank4/14 jobs; relabeling only the minimum subset has the same bookkeeping outcome as A. Intermediate added rows must each justify coverage or provenance benefit. All counts assume replacement of historical labels for selected geometries, never double counting old and new labels.

## Builder acceptance gates

1. Freeze A-approved training membership and hashes. Keep the29 directly supported row identities and exact copied REF fields; map any replacement to its original frame and new archive without pretending it recovers original-run settings.
2. Require each new archive input geometry/order/species/cell/PBC identity, converged SCF, finite labels, original log and settings evidence, energy-field identity and exact force copy. Preserve historical absent IDs; do not invent them.
3. Use a declared native extrapolated target; retain free energy separately. Common target does not imply forces are exact derivatives of extrapolated energy at finite smearing. Prior paired-response evidence motivates preserving both and reviewing derivative consistency; it does not prove archives are wrong. Unary bulk Gamma adequacy must be reviewed before adopting interface numerical settings. Distinguish XC agreement from kpoint/smearing/version agreement.
4. Recompute exact composition rank, row counts, all source hashes, replacement/no-duplicate ledger and role isolation. Keep development and withheld labels outside training and model selection; compare permitted inputs under frozen algorithms before integration. Do not call shared-parent tests morphology-independent.
5. Publish explicit provenance tiers and uncertainty, then require A review before production build/training. Algebraic rank and provenance certainty are acceptance checks, not model force gates or permission for MD/TTM.

## Scope

Only previously published training audit metadata was reused; no development/withheld candidates or outputs, scores, model inference, DFT, dataset or role edits. Source artifacts were SHA256 verified. Reproduce with analyze.py. Options JSON/CSV contains exact counts and prerequisites.
