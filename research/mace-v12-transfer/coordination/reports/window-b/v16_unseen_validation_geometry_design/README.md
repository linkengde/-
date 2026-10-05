# V16 unseen-validation geometry design

Authoritative proposal bundle: [`final_candidates/screening_report.md`](final_candidates/screening_report.md).

## Proposed geometries

| Interface | Ag registry translation (Å) | Marked Ag–X distance (Å) | Exact / near duplicates | Candidate contact screen |
|---|---:|---:|---:|---|
| AgC | (+0.522, −0.671, 0) | 1.989 | 0 / 0 | pass |
| AgSi | (+0.633, +0.403, 0) | 3.005 | 0 / 0 | pass |
| AgTi | (+0.346, −0.776, 0) | 3.101 | 0 / 0 | pass |

These are geometry-only proposals with preserved IDs/order/cell/PBC. Each translates the Ag sublattice rigidly and changes the local environment of framework atoms identified in the historical force-localization report. The framework itself remains fixed, so the candidates do not test framework motion or independent morphology.

Dedupe covered 59 geometry/source inventory entries, 81 hashed source files and all six V16 acquisition inputs. It compared same-composition structures across existing train/validation/test frames, reserved V13–V15 input geometries, planned V16 inputs and earlier B geometry proposals. The nearest existing geometries fail the prior near-duplicate rule because their Ag registry RMS exceeds 0.15 Å; full distances and IDs are in the candidate manifest.

Contact checks are species-specific and empirical: each element-pair minimum is compared with the smallest same-composition distance among existing split geometries and official input geometries, with a 0.05 Å allowance. The 1.741 Å C–Ti framework contact remains intact and is not treated as an all-element cutoff. Per-pair minima, empirical floor values and source paths are recorded in `final_candidates/candidate_manifest.json`.

The geometries inherit small-cluster motifs and remain unassigned proposals. A must freeze and assign their role before any V16 training; no DFT label was read or generated for them. The six planned acquisitions and all scored V15 probes remain distinct from these proposal geometries.

`final_candidates/SHA256SUMS.txt` checks the bundle. The root `SHA256SUMS.txt` covers this README, generator, final bundle and retained attempt diagnostics.
