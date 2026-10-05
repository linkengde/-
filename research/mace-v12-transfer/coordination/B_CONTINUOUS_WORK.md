# Standing A/B collaboration instruction

User authorizes continuing the project and B assistance. A maintains B tasks in
main; user does not need to issue a command for every task. B's registered owner
and separate environment must remain unchanged.

## B loop while the session remains active

1. Finish current assigned work, verify and publish compact deliverables and status.
2. Preserve local edits; synchronize main only from a clean tracked state with
   `git fetch origin main` and `git merge --ff-only origin/main`. On a conflict,
   report it; no reset, force push or owner change.
3. Read `tasks/window-b.json`, this document and current task instructions. Select
   the first not-completed assigned job in `queued_work_order`. Start it through
   `sync_tasks.py progress`; no additional user confirmation is needed.
4. Completed jobs are never rerun simply because START_WINDOW_B still describes
   them. Publish each stage before moving on. If no assigned work is ready, publish
   waiting status and check main periodically while the session is active; use
   short waits and keep progress messages. Do not invent DFT/training assignments.
5. Pause on calculation failure, verification failure, owner conflict or unusable
   environment with exact recovery instructions. Also respect true dependencies;
   do not bypass A review or use blind labels while waiting.

A reviews B outputs, integrates production files, and replenishes this queue.
B cannot start A jobs or alter production settings. This is an agent workflow,
not a background daemon: Git commits cannot awaken a session that has ended.

## Assigned queue

### v15_training_evaluation_entry_draft

Use the detailed assignment in START_WINDOW_B.md. Finish and publish it first.
No actual training/inference or V15 blind-output reads.

### v15_acquisition_batch_priority

Write only `reports/window-b/v15_acquisition_batch_priority/`. Use existing
wide-span design and candidate manifest, historical V14 localization, and
provenance reports; do not repeat candidate generation or read V15 blind labels.
Recommend a small initial future training-acquisition batch from the 14 passing
proposals: normally two per interface, one wider distance endpoint and one
registry proposal where available. Review all 11 candidate near-pairs to avoid
redundancy. If an interface has no distinct approved registry, say so rather than
promote a flagged contact. Provide exact candidate IDs/hashes, novelty by axis,
species minimum distances, parent lineage, priority and backup choice. Screening
PASS is geometry triage, not physical validation. Include input hashes, JSON and
report. Do not calculate, change candidates or create new blind tests. A decides
which receive DFT after its current screening.

### v15_reference_convention_review

Write only `reports/window-b/v15_reference_convention_review/`. Reuse existing
historical provenance ledger; investigate reference-energy conventions in committed
training frames and archived acquisition scripts/logs, never V15 blind results.
Tabulate what REF_energy corresponds to (native energy/free_energy where evidenced),
smearing/k-point declarations, and which energy convention current residual
acquisition scripts actually request. Distinguish original evidence from declarations.
Quantify available paired native energy/free_energy offsets and explain whether
consistent bookkeeping is established or still unknown. Do not infer that offsets
alone explain force errors. Produce recommendations for future matched relabeling
and preservation of convention, with explicit unresolved fields. No dataset edits,
new DFT, model inference/training or MD/TTM. Publish compact report/JSON/hashes.

## Stage progress and publication

Replace JOB below with the actual assigned name. Begin with state running,
iteration 0. Completion requires actual finished artifacts and hashes:

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job JOB --state completed --iteration 1 --note "Completed assigned report/entry stage; syncing queue for next work."
python3 research/mace-v12-transfer/coordination/sync_tasks.py publish window-b
```

Progress `iteration` for these report jobs is a stage counter, not SCF iteration.
Do not print credentials or upload large checkpoints. Existing B allowed paths
remain in force. Keep main worklog/global manifest and A task untouched.
