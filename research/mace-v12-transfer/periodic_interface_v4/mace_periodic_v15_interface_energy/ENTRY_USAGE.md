# Integrated training/evaluation entries

A reviewed and adapted B entry drafts to use new child run directories under
this V15 folder. Input pins and numerical settings are unchanged; changes and
source hashes are in entry_integration.json. Training is not started while
four-rank DFT is active. Use run_v15_training.sh with --repo-root and
--output-dir (e.g. this folder/training_01) after all blind archives pass.
Run preflight.py without --check-blind-archives for label-free readiness checks.
Before evaluate_v15.py, freeze the selected model and training-completion evidence
in a reviewed selection record; see the original B entry README for its schema.
All narrow gates passing yields UNDETERMINED, not physical model PASS.
