# V19 residual-response entry

Training36 candidate is verified; no independent development labels have been integrated. No training has started.

Run build_training_candidate.py --check to verify hashes and exact label preservation. preflight_v19.py requires official train36/dev3/test0 data, reviewed independent lineage, source hashes, geometry isolation, B audit and explicit A approval. It currently exits BLOCKED as expected. run_v19_training.sh retains V18 seed45/medium/batch1/80epochs/lr1e-4/wd5e-7/energy100/forces1000/final79 policy and records successful launcher completion explicitly; no test_file is used. Syntax and missing-data gate tested; training command has not executed. Final model must still undergo all80-epoch and selected79 tensor-identity verification before freezing/evaluation.
