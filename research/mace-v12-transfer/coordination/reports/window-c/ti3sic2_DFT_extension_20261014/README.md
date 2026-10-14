# Ti₃SiC₂ supplemental DFT preparation — Stage C

This folder contains a 14-point proposed data plan, ten source-backed candidate geometries, their geometry validation and SHA-256 provenance. The four independent DEV/TEST geometry slots are intentionally uncreated because a distinct, auditable source family has not been identified.

The structure builder uses ASE 3.29.0 and never creates a calculator. No DFT, potential fitting, MD or TTM was run. Read `FINAL_DECISION.md` before using any input. `START_DFT=NO`; none of the structures or data roles is execution-authorized or frozen for training.

Requested deliverables are `DFT_STRUCTURE_LIST.csv`, `TRAIN_DEV_TEST_SPLIT.md`, `DFT_INPUT_AUDIT.md`, `RESOURCE_ESTIMATE.md`, `FINAL_DECISION.md`, and `SHA256SUMS.txt`. Extra provenance and builder files make the source and geometry checks reproducible.

