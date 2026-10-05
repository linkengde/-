# Superseded design-run artifacts

These generated attempts are retained for audit. Do not use their candidate geometries as current proposals. The authoritative candidates and screening are in `../final_candidates/`.

- `candidate_set_01`: extxyz coordinate round-trip tolerance was too strict; corrected to the writer's 8-decimal coordinate precision.
- `candidate_set_02`: report-generation field mismatch; retained the failure record.
- `candidate_set_03`: contact-audit field type mismatch; retained the failure record.
- `candidate_set_04`: first successful geometry pass, but its empirical contact baseline included unapproved proposal geometries; superseded by the baseline limited to existing splits and official input geometries.
- `candidate_set_05`: same geometry selection as final, before near-duplicate comparison was extended across differing vacuum-cell sizes.

None of these attempts accessed DFT outputs or performed calculations. `final_candidates/` is the only candidate set proposed for A review.
