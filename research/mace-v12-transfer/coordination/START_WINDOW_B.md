# Window B current assignment

The four B DFT labels, v13 residual report, and the single v14 training/evaluation run are complete. V14 screened `SCREEN_FAIL`; do not retrain it or use it for production MD/TTM.

## Run now: localize v14 holdout force errors

Window A is running three residual-targeted Ag-C/Ag-Si/Ag-Ti DFT candidates for a possible v15 cycle. While those run, B should diagnose which atoms and neighbors dominate the v14 holdout force errors. This is post-hoc analysis only: do not retrain, edit the v14 dataset/holdouts, or start DFT.

Send this command to Window B:

```bash
cd /workspace/-
git fetch origin main
git merge --ff-only origin/main
bash research/mace-v12-transfer/coordination/reports/window-b/run_v14_force_localization.sh
```

The script verifies the B owner and the exact v14 checkpoint/reference hashes, analyzes the three scored blind holdouts, and publishes a per-atom JSON/Markdown report under `coordination/reports/window-b/`.

## After the report

Wait for Window A to publish and verify the targeted DFT labels and issue a new explicit assignment before B starts another model cycle. Any v15 training must use newly reserved blind geometries. Do not run production MD or TTM before the validation gates pass.

## Task ownership

Window B remains owned by instance `c5035b48-f43b-4da0-b8e4-e2862817f86a`. If the claim fails or the instance ID differs, stop and report it; do not change the owner. B may publish only its task/progress, `coordination/reports/window-b/`, its assigned acquisition results, and the assigned v14 model-cycle artifacts. A owns the main work log and global SHA256 manifest.
