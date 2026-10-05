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

## next_acquisition_execution_pack

Prepare only under `reports/window-b/next_acquisition_execution_pack/`. Use the
existing six selected candidates in batch_recommendation.json, preserving their
exact bytes/hashes, IDs, cell/PBC and label absence. Rescreen these INPUT geometries
against current committed train/valid/test, existing acquisition INPUTS, and frozen
blind INPUTS only; record exact/near overlap and source inventory hashes. Do not
read V15 blind labels, summaries, calculations or live A files. Do not regenerate
candidates or change their geometry.

Produce an input manifest and review-only DFT runner/archive drafts adapted from
the existing targeted-acquisition workflow, with explicit root/output arguments
and nonempty/active-job guards. Keep PW-PBE500/Gamma/Fermi0.1 settings for these
same-motif cells, original energy convention and four MPI ranks. Record both
native energy and free energy when supported, without changing the chosen REF
convention. Include geometry/input/output hash, finite-label and convergence
checks. Do not import or execute GPAW during pack preparation.

Propose balanced future ownership: B three wider-distance candidates (AgC, AgSi,
AgTi); A three registry candidates. This proposal does NOT claim jobs or change
owner cards. Planned labels must be distinct from all existing IDs. Mark every
planned calculation `blocked_pending_A_review_and_V15_screen`, with launch disabled
until A publishes an explicit authorization/owner assignment. Do not create a
self-authorization flag or executable shortcut that starts six jobs immediately.
Actual production job integration/claims stay with A.

Provide README, six input files, source manifest, duplicate/screening report,
runner/archive drafts, readiness checklist and SHA256 list. No calculation,
training/inference, MD/TTM, production edits or old artifact overwrite. Preserve
all failure diagnostics. Publish this stage and continue the standing queue;
if no assigned jobs remain, report waiting rather than repeating completed work.

## v15_metric_implementation_review

While A trains V15, independently review the integrated evaluate_v15.py and
preflight.py under mace_periodic_v15_interface_energy. Write only
reports/window-b/v15_metric_implementation_review/. Do not edit production files.
Check energy units/atom normalization, atom-weighted vector versus component
RMSE, maximum-error definitions, separating-force sign and marked ID mapping,
coordinate/PBC conventions, historical-versus-fresh aggregation, gate handling,
selection-record checks and fail-fast behavior. Exercise pure mathematical
examples with synthetic arrays if needed; no model inference, reference label
reads, training or GPAW. Report definite defects with reproduction and proposed
patch separately from limitations. Never read V15 blind outputs, summaries or
models. Publish report/compact checks/hashes and completion, then sync queue.

## v15_force_error_localization — post-evaluation authorization

A has finished 80 epochs, frozen model selection and scored all V15 probes once.
This specific job now permits reading the published V15 selected model, selection
record, completed evaluation and its already-scored reference labels, and CPU
inference. Earlier preparation-only restrictions continue for other tasks.
Write only reports/window-b/v15_force_error_localization/. Verify model/reference
hashes. Reuse your historical localization method for the three fresh V15 probes:
per-atom vector and Cartesian errors, species breakdown, contact-shell versus
outside-shell RMSE with stated geometric cutoff, largest-error IDs/neighbors, and
marked-pair projected contributions. Confirm the aggregate metrics reproduce A's
published values within justified numerical tolerance. Compare V14/V15 patterns
without treating different probes as a matched direct model comparison.

Report AgSi's high vector RMSE but small marked-pair separation error distinctly.
Suggest which existing wider-distance/registry candidates address observed local
coverage gaps and which do not. V15 probe labels are no longer blind after scoring:
if reused for future training, their future role must be acquisition/history, not
new independent validation. Do not generate/tune a new model, inspect other unseen
labels, start DFT/MD/TTM, modify production files, or claim validation PASS.
Publish script, compact CSV/JSON/report/hash and completion, then sync queue.

## v16_distance_dft_acquisition — actual DFT authorized

A reviewed V15 FAIL/localization and assigned three complementary wider-distance
training labels in pbe_interface_v16_parallel_acquisition/input_manifest.json.
Run these only on the registered B instance. Verify no B DFT/training is active,
reuse existing GPAW MPI/ScaLAPACK environment and four actual cores. Execute:

```bash
bash research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v16_parallel_acquisition/run_assigned_labels.sh
```

Set GPAW_PYTHON if the known venv path differs; verify software/feature support
before launch, never silently switch calculator settings. The runner claims B,
publishes per-label progress and verified compact archives and skips completed
labels. Existing unarchived runs must be inspected, never duplicated. Keep large
state.gpw local. Stop on SCF/verification/ownership/environment failure with
diagnostics. No A registry labels, dataset edits, training, MD/TTM. Complete job
state only after all three assigned archives pass and publish; then sync queue.
These are future training acquisitions, not new blind tests.

## v16_unseen_validation_geometry_design

After assigned acquisitions finish (or alongside running DFT only if it does not
compete for MPI/CPU resources), prepare label-free proposals under
reports/window-b/v16_unseen_validation_geometry_design/. Need a withheld validation
set distinct from the six training acquisitions and all scored V15 probes. Read
committed geometries and historical localization only. Never use unseen DFT labels.

Propose one candidate per AgC/AgSi/AgTi with a new registry/environment combination;
prioritize transverse Ag/Si response and framework-neighbor coverage based on
observed errors. Screen species-aware contacts and periodic geometry; preserve
original IDs/order/cell/PBC where applicable, no energy/force fields. Compare exact
and near geometries against all current train/valid/test and planned acquisition
INPUTS. Record construction axes, parent lineage, species minima, ID displacements
and every source hash. A freezes accepted candidates before V16 training; DFT
label calculation requires separate explicit owner assignment. Do not claim motif
independence if inherited from existing small-cluster parents. Explain remaining
morphology/thermal coverage limits. No candidate receives a training role silently.
If no suitable distinct inputs can be produced, publish the geometric blocker and
required new parent source, not a false independence claim. No new DFT/model
training/inference/MD/TTM under this design job. Publish report/manifest/hashes.
