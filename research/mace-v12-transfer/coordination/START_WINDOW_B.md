# Window B current assignment

Window B's Ag-C species-aware geometry review is complete. A is running the four-rank Ag-Ti residual PW-PBE label. B has a separate, read-only audit task that can proceed without using DFT/MACE resources or touching A's work.

## Assignment: audit v15 dataset integrity and label provenance

Goal: give Window A an independent, evidence-based checklist for building v15 after the current DFT labels finish. This assignment is an audit only. Do not edit scripts, datasets, labels, manifests, or holdout inputs.

First synchronize and verify this is the registered B machine:

```bash
set -e
cd /workspace/-
git fetch origin main
git merge --ff-only origin/main
python3 research/mace-v12-transfer/coordination/sync_tasks.py identity
python3 research/mace-v12-transfer/coordination/sync_tasks.py status --remote
```

The B identity must be `c5035b48-f43b-4da0-b8e4-e2862817f86a`. Confirm `v15_dataset_integrity_audit` is assigned in `window-b.json`. If the identity differs, the branch cannot fast-forward, or the assignment is absent, stop and report the issue.

Publish the start event:

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_dataset_integrity_audit --state running --iteration 0 --note "Auditing existing v15 data provenance, split leakage, and builder acceptance checks; read-only, no calculations."
```

Review existing repository files only, including:

- v12, v13, and v14 dataset builders, manifests, splits, and training-frame audits;
- the v15 residual acquisition input manifest and the already published Ag-C and Ag-Si results;
- the three frozen v15 registry holdout inputs and their manifest; their labels are still pending until A publishes them;
- `V15_FORCE_SEPARATION_PLAN.md` and the V14 validation assessment.

Record evidence for:

1. Each existing training/validation/test frame's label types, source method, configuration tag, and file provenance; distinguish force-only LCAO labels from PW-PBE energy-and-force labels.
2. Exact and near-duplicate risks across train/valid/test and reserved holdouts, comparing atom IDs, elements, cell/PBC, and same-ID coordinates. Give the threshold used and identify any split that must not count as independent evidence.
3. The four legacy force-only LCAO training frames: their identities, what provenance is present, what is missing, and a conservative v15 inclusion rule. Do not infer undocumented GPAW/basis settings.
4. Expected v15 frame and energy/force label counts if verified residual labels enter training, while all three v15 registry holdouts remain blind and excluded from training, validation, checkpoint selection, and hyperparameter selection.
5. The elemental-composition matrix rank and checks a v15 builder should enforce, plus hash and byte-preservation checks for frozen splits.
6. A concise pass/fail/unknown recommendation for each audit item, with paths and hashes as evidence.

Do not inspect or depend on A's `/tmp` checkpoint or live process. Treat Ag-Ti residual and all three new registry labels as pending until their verified compact archives appear on `main`. Do not start DFT, run MACE inference/training, run MD/TTM, or perform any additional calculation. Do not edit data/build scripts. Publish only:

```text
research/mace-v12-transfer/coordination/reports/window-b/v15_dataset_integrity_audit.md
research/mace-v12-transfer/coordination/reports/window-b/v15_dataset_integrity_audit.json
research/mace-v12-transfer/coordination/reports/window-b/v15_dataset_integrity_audit_SHA256SUMS.txt
```

Include the analyzed `main` commit and input SHA256 values. Before publishing, fetch `main` again; if relevant inputs changed, refresh the audit against the latest committed files. Then finish:

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_dataset_integrity_audit --state completed --iteration 1 --note "Published read-only v15 dataset-integrity and label-provenance audit with evidence hashes."
python3 research/mace-v12-transfer/coordination/sync_tasks.py publish window-b
```

Do not claim that v15 is ready to train. Window A owns the eventual dataset build, training, and scoring.
