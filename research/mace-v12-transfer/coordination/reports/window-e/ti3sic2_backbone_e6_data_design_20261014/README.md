# Stage E6 — minimum dataset and constrained-fit design

This directory records a design decision based on the E2–E5 archived audits. It contains no fitted parameters, changed potential, newly generated structure, or new scientific calculation. Current TRAIN/DEV/TEST role freeze remains unchanged and empty.

Read in this order:

1. `E6_DATA_DESIGN.md` — evidence and gaps;
2. `PARAMETER_SELECTION.md` — candidate 2/3/4 scalar tied blocks and measured local conditioning;
3. `MINIMUM_TRAINING_SET.md` — family-level train/validation split proposal;
4. `DFT_REQUIREMENTS.md` — reuse versus new DFT labels;
5. `FIT_READINESS_DECISION.md` — release gate.

`fit_candidate_conditions.csv` contains recomputed, column-normalized force-Jacobian submatrix diagnostics using the E5 per-component matrix. These are current-geometry local diagnostics only. The file does not predict conditioning after proposed new data are added.
