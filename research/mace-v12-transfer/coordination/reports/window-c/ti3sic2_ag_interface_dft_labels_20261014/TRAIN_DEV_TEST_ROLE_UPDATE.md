# TRAIN / DEV / TEST role update — proposal only

This file records the intended use of the Stage C-D1 candidates; it does not change any existing record's role.

## Group A — Ti₃SiC₂ backbone

The six new volume/shear/Si-displacement geometries descend from `TSC_KISI_COD9009647_BASELINE_DERIVED_FAMILY`. Keep this entire family together. Proposed role after labels and an explicit role freeze: `TRAIN_CANDIDATE_FOR_TSCBOP_CONSTRAINED_EVALUATION`.

The five existing baseline/Ti4f/C4f DFT labels retain `NUMERICAL_CONTROL/DIAGNOSTIC` status; this task neither copies them into a fit set nor relabels them. The earlier E5 interface-derived 114-atom TSC slab remains force-diagnostic only. These records have been used in prior diagnostics and cannot serve as independent DEV or TEST.

## Group B — Ag/Ti₃SiC₂ interface

The nine gaps descend from one 217-atom Si-terminated interface parent family. Proposed role after DFT and separate approval: `TRAIN_CANDIDATE_FOR_FUTURE_UNIFIED_AG_TI_SI_C_MODEL`. These data cannot constrain the three-element Ti–Si–C TSCBOP parameterization. Keep the full gap scan in one family and do not mix it into the Group A TSCBOP fit.

## DEV and TEST

- `DEV`: still `SOURCE_PENDING`; E7 did not identify a genuinely independent Ti₃SiC₂ source family, and no reproducible independent Ag/Ti₃SiC₂ coordinate family with labels was established here.
- `TEST`: remains evaluator-controlled. No sealed test geometry or label was opened, copied, or enumerated.

Do not split perturbation signs, interface gaps, k-point variants, or descendants of one parent between TRAIN and DEV/TEST. This candidate set is not a validated fitting dataset: no labels exist and no independent DEV/TEST is populated.
