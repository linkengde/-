# Window B: inherited-method provenance trace

Task: `v15_inherited_method_provenance_trace`.
Registered owner: `c5035b48-f43b-4da0-b8e4-e2862817f86a`.

## Start

```bash
set -e
cd /workspace/-
git fetch origin main
git merge --ff-only origin/main
python3 research/mace-v12-transfer/coordination/sync_tasks.py claim window-b
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_inherited_method_provenance_trace --state running --iteration 0 --note "Tracing 15 unknown inherited training methods from historical evidence; no calculations."
```

If ownership or synchronization fails, stop and report the exact conflict.

## Work

Use the 15 `unknown_inherited_method` training records in
`reports/window-b/v15_builder_draft/generated_data/dataset_manifest.json`.
Read the prior dataset-integrity audit first; do not repeat its completed checks.

Trace each record through existing V9-V14 datasets/manifests, original archived
inputs/outputs, calculation logs and scripts, and relevant Git history. Search
committed repository evidence only. Do not inspect A's live files, checkpoints,
V15 blind outputs, labels or summaries. No new DFT or model execution.

For each record publish:
- config type, source split/index/hash and persistent atom IDs;
- traced original calculation paths/commit IDs and evidence hashes;
- method/basis or PW cutoff, XC, k-points, PBC/cell, convergence and software
  version only where evidenced; mark missing fields unknown;
- geometry and energy/force identity to original output, with comparison
  tolerances stated; do not identify provenance merely by similar geometry,
  a filename, a config tag or a global manifest;
- separate evidence-backed recovery, partial recovery and unresolved status;
- potential reference-method conflicts and missing artifacts needed to resolve
  them. Where evidence is absent, list exact existing frame geometries for
  possible matched PW-PBE relabeling; do not start relabeling.

Produce a compact per-frame JSON/CSV ledger, readable report and SHA256 list
under `coordination/reports/window-b/v15_inherited_method_provenance_trace/`.
Do not duplicate models or large historical datasets. Do not modify the builder
draft, production data/scripts, model settings, A task, main worklog or global
manifest. No DFT/MACE inference/training/MD/TTM. Unknowns must remain unknown;
this task need not recover all 15 to finish honestly. A decides integration.

## Finish

Publish report artifacts, then update completion through the coordination script:

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_inherited_method_provenance_trace --state completed --iteration 1 --note "Published per-frame inherited-method trace and unresolved-evidence remediation ledger; no calculations."
python3 research/mace-v12-transfer/coordination/sync_tasks.py publish window-b
```
