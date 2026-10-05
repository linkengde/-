# V15 training/evaluation entry drafts — A review required

B prepared these isolated scripts; B did not train, infer, run DFT/MD/TTM, read fresh blind label files or inspect A live checkpoints. All output must be a new child of this report directory. A may review/adapt these drafts for production integration; do not silently repoint their safeguards. Explicit repository/output arguments make invocation independent of working directory and script location.

## Safe input check (executed by B)

```bash
/workspace/.venvs/gpaw-mpi/bin/python -B /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_training_evaluation_entry_draft/preflight.py --repo-root /workspace/-
```

Default mode never opens holdout calculation files. It verifies pinned A dataset/manifest/builder/foundation bytes, 29/2/3 counts, finite REF labels, actual composition rank4, historical valid/test byte identity, and frozen INPUT hashes/label absence/geometry isolation. Changed pins fail; A must review changed artifacts rather than regenerate hashes to bypass failure. `input_pins.json` includes the V14 foundation file SHA256; no foundation override or architecture tuning is exposed. V14 defaults, seed45, CPU, batch1, max80, lr0.0001, decay5e-7, energy100/force1000 and estimated E0s are retained.

## Actual training — A only, after review

Wait for all DFT to end and all three blind archives to be complete and PASS. Runtime preflight checks no active MPI launcher, archive source hashes, PASS checks, convergence and declared archive member hashes. It reads archive metadata and computes file hashes; it does not parse blind extxyz labels or use them for model selection. B did not run this mode. MPI process absence alone is insufficient: A must also confirm its queue has ended.

```bash
bash /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_training_evaluation_entry_draft/run_v15_training.sh \
 --repo-root /workspace/- \
 --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_training_evaluation_entry_draft/a_training_01
```

Runtime defaults: `/workspace/.venvs/mace-v12/bin/python` and `mace_run_train`; A can set `MACE_PYTHON`/`MACE_RUN_TRAIN` to its reviewed runtime. The launcher preserves pipefail and real trainer exit status, refuses nonempty output/existing checkpoints, and never retries automatically. Failure outputs remain for diagnosis. Exit0 plus a checkpoint is not training completion evidence or model PASS. A must inspect the full epoch log for the intended 80 epochs (zero-based logs may end at epoch79), checkpoint selection and actual effective settings, including version/default architecture. Frozen historical test is not a checkpoint-selection set; selection uses training/validation only. Blind data are never passed to the trainer.

## Freeze selection before actual evaluation — A only

A writes a reviewed selection JSON outside the evaluation output, **before any blind evaluation**, with these fields (replace placeholders with actual evidence; this is a schema example, not fabricated completion):

```json
{
 "reviewed_by": "window-a",
 "model_path": "/absolute/path/to/selected.model",
 "model_sha256": "ACTUAL_SHA256",
 "training_exit_code": 0,
 "completed_epochs": 80,
 "selected_epoch": 71,
 "selection_reason": "ACTUAL_TRAIN_VALID_SELECTION_REASON",
 "selection_uses_only_training_validation": true,
 "blind_labels_used_for_selection": false,
 "training_stdout": {"path": "/absolute/path/to/train.stdout", "sha256": "ACTUAL_SHA256"},
 "epoch_completion_evidence": {"path": "/absolute/path/to/epoch_log", "sha256": "ACTUAL_SHA256"}
}
```

The example epoch71 is illustrative, not a proposed V15 checkpoint. A's reviewed completion declaration is required; the script does not guess completion from a filename. The selection record and evidence hashes bind a reviewed frozen selection but do not independently prove A's declarations. Do not select or retune using blind outcomes.

```bash
/workspace/.venvs/mace-v12/bin/python -B /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_training_evaluation_entry_draft/evaluate_v15.py \
 --repo-root /workspace/- \
 --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_training_evaluation_entry_draft/a_evaluation_01 \
 --model /absolute/path/to/selected.model \
 --selection-record /absolute/path/to/a_reviewed_selection.json
```

Only this explicit evaluation invocation parses blind reference labels, after checking selection/evidence/model hashes. It revalidates all archives and exact ordered species/IDs/marked-pair/cell/PBC/coordinates. Missing fields/files, incomplete archives, nonfinite labels/predictions and mismatched geometry fail. Reference keys are explicit; free_energy is never substituted. CPU float64 follows V14 evaluator logic. Force-vector RMSE is sqrt(mean over atoms of squared 3-vector error), not component RMSE. Separating force uses `(F_X-F_Ag) dot unit(r_X-r_Ag)` with the historical direct-coordinate cluster convention. Historical regression and fresh local registry scores are reported separately, with per-structure rows, aggregate energy MAE, atom-weighted vector RMSE/max, and marked-pair errors.

Provisional fresh-probe gates: 10 meV/atom, 0.05 eV/Å force-vector RMSE and 0.10 eV/Å separation error. Any failed interface gives FAIL; all narrow gates met gives UNDETERMINED, never automatic model PASS. Historical correlation and inherited-method uncertainty remain explicit. No broad morphology/thermal transfer claim or MD/TTM authorization follows. `evaluation.json`, `per_structure.csv` and hashes are saved only to the chosen isolated output directory.

## Code-review checklist

- [x] Separate root/output arguments; new report children only, nonempty refusal.
- [x] Agreed V14 foundation/default architecture and fixed numerical parameters.
- [x] Dataset/foundation/builder/input hash pins; actual counts/rank/finite labels.
- [x] Historical split bytes; fresh INPUT leakage guard; blind data absent from trainer arguments.
- [x] Runtime archive completeness/PASS/hash gates; active MPI refusal.
- [x] pipefail; no automatic retrain; A epoch/checkpoint review required.
- [x] Frozen selection/evidence hashes checked before blind reference parsing.
- [x] Runtime reference/prediction finiteness and geometry/ID checks; explicit label fields.
- [x] Per-structure/aggregate vector and separation metrics; FAIL versus UNDETERMINED.
- [ ] A reviews its actual package versions/defaults and these scripts before running.
- [ ] A confirms finished DFT queue/all blind archives and reviews 80-epoch/checkpoint evidence.
- [ ] Actual training/evaluation runtime validation remains unexecuted by B.
