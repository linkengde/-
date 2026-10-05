# Window B: V15 builder draft and provenance audit

Job: `v15_builder_draft_and_provenance`. Owner must remain
`c5035b48-f43b-4da0-b8e4-e2862817f86a`.
A continues V15 blind DFT; B prepares isolated dataset engineering in parallel.

## Start

```bash
set -e
cd /workspace/-
git fetch origin main
git merge --ff-only origin/main
python3 research/mace-v12-transfer/coordination/sync_tasks.py claim window-b
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_builder_draft_and_provenance --state running --iteration 0 --note "Preparing isolated V15 builder and provenance audit; no DFT or training."
```

Stop if identity/ownership conflicts or fast-forward fails. Read the existing
V15 dataset-integrity audit, V15_FORCE_SEPARATION_PLAN.md, V14 builder and dataset
manifest before implementation.

## Deliverables and boundaries

Work only under:
`research/mace-v12-transfer/coordination/reports/window-b/v15_builder_draft/`.
Write a reusable builder draft with explicit repository-root and output-directory
arguments; resolve paths independently of where the draft lives. Refuse an
existing nonempty output directory. Never write production data or overwrite
historical files.

1. Start from V14 training data. Explicitly exclude the four undocumented
force-only LCAO frames identified in your audit. Record config types and reasons;
do not fabricate missing energies or per-frame provenance.
2. Add the three completed AgC/AgSi/AgTi residual-shell V15 training acquisitions
from `pbe_interface_v15_targeted_acquisition/` only after verifying archive PASS,
convergence, source hashes, input/output geometry identity, and finite labels.
Do not add the 33 unlabeled wide-span proposals.
3. Preserve historical valid/test bytes and document the historical AgTi
training/test geometry overlap even when PBC differs. Report independent-geometry
status separately from complete-input identity. Do not silently resplit history.
4. Read frozen V15 registry holdout INPUTS only to check geometry exclusion.
Never read holdout outputs, energies, forces, summaries, or active A files.
Check candidate train/valid/test overlap by IDs/elements/cell/coordinates, with
PBC differences separately reported. Fail if new V15 blind geometry enters
train/valid; record inherited history overlap without calling it independent.
5. Produce dataset manifest with frame/label counts, source and output SHA256,
energy composition matrix rank (required 4/4), excluded-frame ledger and
per-frame provenance status. Unknown inherited provenance must remain explicitly
unknown; explain whether it blocks physical-method consistency claims.
6. Exercise the draft on committed training sources, writing generated data only
inside this report directory. Provide the exact command and a verification report.
Keep expected training size (31 - 4 + 3 = 30) as a diagnostic; investigate any
mismatch rather than bypassing validation. Document actual rank/counts/results.
7. Include source script, README, report, manifest, and SHA256 list covering
all deliverables. No new dependency reinstall if existing ASE/numpy suffice.

Do not change production builders, task A, frozen inputs, training settings,
models, main worklog or global manifest. Do not run DFT, MACE inference/training,
MD or TTM. A reviews and integrates the draft; this does not authorize training
or establish V15 validation PASS.

## Publish

Commit report artifacts and then publish completion through the existing script.
If checks fail, preserve diagnostics and mark failed with the exact blocker.

```bash
python3 research/mace-v12-transfer/coordination/sync_tasks.py progress window-b --job v15_builder_draft_and_provenance --state completed --iteration 1 --note "Published isolated V15 builder draft, generated-data checks and provenance evidence; A integration pending."
python3 research/mace-v12-transfer/coordination/sync_tasks.py publish window-b
```
