# Window B: V15 training and evaluation entry drafts

Job `v15_training_evaluation_entry_draft`; owner remains
`c5035b48-f43b-4da0-b8e4-e2862817f86a`.

```bash
set -e
cd /workspace/-
git fetch origin main
git merge --ff-only origin/main
python3 research/mace-v12-transfer/coordination/sync_tasks.py claim window-b
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_training_evaluation_entry_draft --state running --iteration 0 --note "Preparing isolated V15 training/evaluation entry scripts; no training, inference or blind-label access."
```

Stop on synchronization/owner conflict. A owns active DFT, production integration,
actual training and evaluation. Write only under
`coordination/reports/window-b/v15_training_evaluation_entry_draft/`.

Read A's integrated `mace_periodic_v15_interface_energy/data/dataset_manifest.json`,
README and builder, plus the historical V14 training/evaluation scripts.
Do not repeat completed provenance audits or alter official data/settings.

Deliver an isolated `run_v15_training.sh` draft and an evaluation Python draft
with explicit root/output arguments independent of script location. Training must
retain the V14 foundation model, architecture/defaults, seed45, max80 epochs,
energy100/force1000, CPU and batch1. Do not tune based on blind outcomes. Preserve
nonempty-output and existing-model refusal; do not use blind labels for checkpoint
selection. Add preconditions for manifest counts/rank/hashes and all three blind
archives being complete/verified before A starts actual training. Preserve pipefail
and real exit status; require A to inspect epoch/checkpoint completion evidence.

Evaluation draft should use historical evaluation logic with V15 model paths,
reporting per-structure and aggregate energy/force-vector errors and marked-pair
separation errors against provisional gates (10 meV/atom, 0.05 eV/A vector RMSE,
0.10 eV/A separation). Distinguish historical regression scores from fresh local
registry probes and keep provenance/motif-independence limitations explicit.
Read blind labels only when A explicitly executes evaluation after model/checkpoint
selection; this B task must not execute such reads. Missing model/reference fields,
nonfinite outputs, wrong IDs/geometry, incomplete archives must fail clearly.
No automatic model PASS from training completion; distinguish FAIL and UNDETERMINED.

Provide README with intended commands and prerequisites, code-review checklist,
SHA256 list and readiness report. Static compilation/shell syntax and input checks
may be performed if useful; no MACE inference/training, DFT, MD or TTM. Do not read
V15 holdout calculation files/summaries/labels or A live checkpoints while preparing.
Do not modify production scripts, A task, main worklog or global manifest.

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_training_evaluation_entry_draft --state completed --iteration 1 --note "Published isolated V15 training/evaluation entry drafts and readiness notes; actual execution remains with A."
python3 research/mace-v12-transfer/coordination/sync_tasks.py publish window-b
```

## Continue automatically

After this task, follow `B_CONTINUOUS_WORK.md` and the ordered B task queue.
Do not stop merely to request the next per-task command from the user.
