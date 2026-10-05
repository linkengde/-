# Entry draft readiness

Completed static Python compilation (without executing evaluation or importing Torch/MACE), Bash syntax validation, and default safe preflight. Dataset is the A-integrated 29/2/3 screening dataset; actual composition rank4; all finite labels, pinned hashes, historical valid/test byte identity and frozen INPUT isolation pass.

Actual training, inference, blind archive runtime checks, epoch/checkpoint inspection and evaluation were not run. No holdout calculations/summary/label files or A live checkpoints were opened. Readiness means statically reviewed drafts with safe input checks, not demonstrated runtime correctness or model reliability. A must review and run the runtime checks after DFT ends.

Training retains historical V14 model/defaults and seed45/80 epochs/energy100/force1000/CPU/batch1/lr0.0001/weight_decay5e-7/E0s estimated. Archives must pass before training; trainer uses only official train/valid/test, not fresh registry probes. Frozen historical test remains regression evidence and is not for checkpoint selection. Architecture/default compatibility must be reviewed against A's actual runtime.

Evaluation requires a selected model hash and A-reviewed 80-epoch completion/selection evidence before reading blind labels. It reports per-structure and role/interface aggregate energy/vector/separation metrics. Provisional gates are 10 meV/atom, 0.05 eV/Å and 0.10 eV/Å; status is FAIL on a failed probe or UNDETERMINED if narrow gates are all met. No automatic PASS, broad transfer claim or MD/TTM authorization.

Limitations carried forward: 14 partially documented and 12 explicitly declared inherited train rows lack complete original-run revalidation. Three residual archives were verified in dataset construction. Physical-method consistency remains UNKNOWN. New registry probes retain parent small-cluster motif correlations. The draft deliberately does not infer missing reference settings or substitute free energy.

See README for exact intended commands, selection schema and reviewer checklist. Verify hashes from repository root:

```bash
sha256sum -c research/mace-v12-transfer/coordination/reports/window-b/v15_training_evaluation_entry_draft/SHA256SUMS.txt
```
