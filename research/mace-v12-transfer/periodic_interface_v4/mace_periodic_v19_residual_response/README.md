# V19 residual-response candidate

The A-approved dataset contains 36 training frames and three development frames, with no test split and no sealed labels opened. Its source-family limitation and method caveats are recorded in `data/dataset_manifest.json` and `entry_review_approval.json`.

`preflight_v19.py` passed before the fixed recipe in `run_v19_training.sh` was run. The CPU fit completed epochs 0–79. Epoch-79 checkpoint coverage and identity with the exported model were verified by `verify_and_evaluate_v19.py`; train/dev-only error analysis is reproducible with `diagnose_v19_generalization.py`.

The candidate **fails** the current force-vector and separation-force development gates. It is not qualified for production or impact simulation. Training fitted train36 substantially better than dev3, but the dev structures differ in chemistry, periodic boundary conditions, and k-point settings from the training data. The three development frames also share one COD 9009647 source family. These results diagnose a generalization gap; they do not establish a single cause.

The archived run and evaluation are in `runs/v19_train36_dev3_20261023_A/`. The summary, per-epoch metrics, diagnostic, logs, final epoch-79 checkpoint/model hashes, and next data steps are in `../../coordination/reports/window-a/v19_train36_dev3_epoch79_evaluation_20261023/README.md`.

Do not open test/sealed data. Resolve energy-key lineage and numerical-method/domain mismatch, build a source-diverse dataset with force coverage for Ti/C backbone sites, obtain independent audit, then rerun the frozen preflight and recipe against a versioned approved split.
