# Stage E7 data integration and pre-fit audit

This package imports and audits Stage C's ten Ti₃SiC₂ **geometry-only** candidates, reviews E5/E6 fitting evidence, and records TRAIN/DEV/TEST readiness by source family. It does not contain newly completed Stage C DFT labels. No fitting, DFT, MD, TTM, or potential modification was performed.

Start with `DFT_DATASET_READY_REPORT.md`. Machine-readable records are `DFT_LABEL_INVENTORY.csv`, `DATASET_ROLE_LEDGER.csv`, `STRUCTURE_AUDIT.csv`, and `PARAMETER_ENVIRONMENT_ACTIVATION.csv`. `audit_dataset.py` and `audit_parameter_environment.py` regenerate the geometry audit tables; `audit_existing_labels.py` verifies the five previously archived B bulk energy/force records, all using the existing ASE 3.29.0 environment. E5 diagnostics are included as unmodified CSV reference copies under `e5_reference/` with upstream SHA records in `PROVENANCE.json`.

The outcome is `START_TRIAL_FIT=NO`: all ten Stage C files lack DFT labels; existing numerical-control roles are unchanged; independent DEV and evaluator-held TEST parent families are absent.
