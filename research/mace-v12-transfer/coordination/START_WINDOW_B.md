# Window B current assignment

V14 training/evaluation is complete and screened SCREEN_FAIL. The v14 per-atom force localization report is complete and published. Window A is running three residual-targeted Ag-C/Ag-Si/Ag-Ti DFT labels for a possible v15 cycle.

## Run now: audit fresh v15 blind-holdout design

This is a read-only geometry audit. Do not start DFT, run MACE inference/training, edit any dataset or holdout input, or run MD/TTM.

First sync the repository and confirm the owner:

    cd /workspace/-
    git fetch origin main
    git merge --ff-only origin/main
    python3 research/mace-v12-transfer/coordination/sync_tasks.py status --remote

Then inspect:
- coordination/reports/window-b/v14_force_localization.md and .json;
- mace_periodic_v14_interface_energy/results/v14_validation_assessment.json;
- pbe_interface_v15_targeted_acquisition/input_manifest.json and its input geometries;
- v14 training/validation/test data and all existing v13/v14 holdout/acquisition geometries.

Publish coordination/reports/window-b/v15_holdout_design_audit.md. For Ag-C, Ag-Si and Ag-Ti separately, document the source geometry and marked pair, candidate-to-training/holdout structural overlap evidence, contact distances and nearest-neighbor environment, and whether the v15 acquisition structures cover the localized residuals. Recommend one fresh blind geometry per interface, selected from a distinct available motif/registry if one exists. If the repository has no defensible independent candidate, report the missing source data instead of fabricating coordinates. Proposed candidates must stay out of v15 training and validation until their blind score is recorded.

Before analysis, publish a B progress event:

    python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_holdout_design_audit --state running --iteration 0 --note "Auditing v15 candidate geometry independence and fresh blind holdout coverage."

After writing the report, publish it and mark the task complete:

    python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_holdout_design_audit --state completed --iteration 1 --note "Published read-only v15 blind-holdout geometry design audit."
    python3 research/mace-v12-transfer/coordination/sync_tasks.py publish window-b

If fetch/fast-forward fails, the B owner differs from c5035b48-f43b-4da0-b8e4-e2862817f86a, or any overlap/geometry check is ambiguous, stop and report the exact issue. Do not change the owner or start a calculation.

## Task ownership

Window B remains owned by instance c5035b48-f43b-4da0-b8e4-e2862817f86a. B may publish only its task/progress and coordination/reports/window-b/. Window A owns the v15 DFT inputs/results, main work log and global SHA256 manifest. The three active A DFT candidates remain excluded from every v14 score and any future blind set.
