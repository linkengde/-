# Window B: provenance-aware builder revision

Job `v15_provenance_aware_builder_revision`. Owner remains
`c5035b48-f43b-4da0-b8e4-e2862817f86a`.

```bash
set -e
cd /workspace/-
git fetch origin main
git merge --ff-only origin/main
python3 research/mace-v12-transfer/coordination/sync_tasks.py claim window-b
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_provenance_aware_builder_revision --state running --iteration 0 --note "Updating isolated provenance ledger and exercising explicit unresolved-row exclusion; no calculations."
```

Stop on owner mismatch or synchronization conflict. Read your completed builder
and inherited-method trace, plus A review in WORK_LOG_2026-10.md.

## Scope

Create a new copy of the builder under
`coordination/reports/window-b/v15_provenance_aware_builder_revision/`.
Do not change the original published draft or production files. Permit output
only in child directories of this new report directory, with nonempty-directory
refusal. Keep all existing archive, hash, finite-label, geometry, and blind-input
isolation checks. Never inspect fresh V15 holdout labels, summaries or outputs.

1. Reconcile the provenance ledger with your trace. Preserve raw `source_method`
and `dft_method` separately, recorded declarations versus original-run evidence,
and partial/unresolved status. Do not invent settings. Expose available k-point,
cutoff, smearing and energy-convention declarations with precise evidence paths.
Maintain original serialized REF labels; never substitute free_energy.
2. Add an explicit exclusion option for the unresolved
`AgTi_2p3871A_periodic_PW_PBE_force_only` row. Match it by config plus source
frame/hash, not config text alone; require exactly one match. Do not automatically
exclude the other 14 partially documented frames.
3. Exercise two isolated scenarios: baseline 30/2/3 and exclusion 29/2/3.
Check each composition rank 4/4, finite labels, retained frame geometry/label
identity, historical valid/test byte identity, and fresh blind INPUT isolation.
Record the exclusion ledger and exact commands. No training or model inference.
4. Report evidence-tier counts and composition rank for each scenario. Explain
that excluding the one unresolved row reduces a provenance risk but does not
fully verify the other inherited rows. Recommend which scenario is suitable for
an explicitly limited numerical screening run; A decides production integration.
5. Provide script, README, two manifests, compact comparison report and hashes.
Avoid republishing duplicate old outputs unless needed for the scenario evidence.

No DFT, MACE, MD or TTM. No production dataset/builders, A task, main log, global
manifest or training settings edits. Do not use blind label results for choices.
If a check fails, retain diagnostics and report; never bypass it.

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_provenance_aware_builder_revision --state completed --iteration 1 --note "Published provenance-aware isolated builder and checked 30/29-frame scenarios; production decision remains with A."
python3 research/mace-v12-transfer/coordination/sync_tasks.py publish window-b
```
