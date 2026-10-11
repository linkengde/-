# AG64 step8000 role decision

**Record ID:** `AG64_STEP8000_ROLE_DECISION`

**Decision:** **B — HOLD**

**Scope:** Offline review only. The XYZ, QA output, thresholds, and source data were not changed. No DFT or training was started.

The package identities pass: the 64-atom Ag geometry, cell/PBC, potential SHA, candidate file SHA, QA SHA, and manifest member hashes agree with the supplied handoff commit `38ff5afb5c3cbe14ca0889212910f0ca6558e991`. `STRUCTURE_PASS=F` is the deliberate manual-C gate state, not a failed physical diagnostic. The two final windows (steps 7000 and 8000) are marked valid by the package's automatic QA, and their reported metrics satisfy the thresholds embedded in `QA_report.json`.

The thresholds are now captured verbatim in `AG64_STEP8000_QA_FREEZE_RECORD.json`; they are an audit record, not a retroactive preregistration. The earlier C release configuration explicitly says `thresholds_frozen_by_C=false`, and the D package does not bind its QA result to a frozen-threshold record or the executable QA/generator code. The QA report also omits the executed sampling-configuration SHA. The candidate-producing script and its commit cannot be identified from the D handoff commit: that commit adds only `MANIFEST.json`, `QA_report.json`, `SHA256SUMS.txt`, `candidate_structure.xyz`, and `thermal.csv`.

The manifest declares checkpoint manifest SHA `44c5b966f52cbface56ed7a1fc135adc421288e91e76263c5ed53382e0230fa8`, but the checkpoint bytes are absent from this handoff, so that declared checkpoint hash could not be independently rehashed here. The QA report's input checkpoint at step 4000 has state/array hashes, but no bound sampling-config SHA or independently verifiable checkpoint manifest in the handoff.

These are provenance and gate-freeze gaps, not evidence that the final structure failed the listed screening metrics. They are sufficient to keep the formal data release on HOLD. The report describes a 1600 K hot-stage snapshot; even after the blockers are resolved, it must not be represented as an equilibrated 1250 K liquid.

## Required to reconsider

1. Supply the exact candidate-generation/QA script path and immutable commit (or script SHA), plus the exact sampling configuration SHA used for this trajectory.
2. Bind the thresholds to the versioned pre-run configuration or provide contemporaneous evidence that the listed values were fixed before these windows were scored. A retrospective record alone cannot establish this.
3. Provide the referenced checkpoint archive or an independently retrievable manifest and members whose SHA-256 values can be rechecked. Clarify whether the step-4000 QA input checkpoint or the step-8000 released checkpoint is the intended source identity.
4. Re-run only the offline provenance and threshold-binding audit after those records arrive; do not regenerate or edit the frozen XYZ.

If the hold is cleared, this source family may be considered for a **fixed-geometry hot-stage DEV single point** only, after a separate DFT release. Keep the family out of TRAIN and TEST; a single hot-stage point does not establish equilibrium liquid behavior or validate melting/solidification physics.
