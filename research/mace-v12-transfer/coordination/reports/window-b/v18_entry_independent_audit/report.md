# V18 independent entry audit

Train30/dev3/test0, exact rank4, candidate/dev byte identity, all pinned source hashes, finite training labels and eight intended signed diagnostics pass. Train/dev geometry isolation passes at published exact/near thresholds. Frame0 source/order/cell preserved; original1x4x2 versus newGamma kpoint-convergence caveat remains. No withheld structure or label was opened.

Default recipe matches seed45/medium/lr1e-4/wd5e-7/batch1/80epochs/energy100/force1000, four threads, no test_file. Runner requires final epoch79 export; finalizer checks all80 epochs and completion markers, freezes selection before development scoring; evaluator exports per-atom predictions and fixed gates. Syntax/AST checks passed without training/inference.

Limitations/recommendations:
- Default threads4, allowed runtime override1–4: A should record actual threads to preserve controlled comparison.
- Finalizer is safe in intended set-e runner flow; called standalone it sets training_exit_code0 from log evidence without independent launcher receipt.
- Evaluator verifies current dataset against manifest but does not bind selection dataset_manifest/train/dev hashes back to current inputs; standalone selection/data drift guard could be strengthened.
- Pair projection uses unwrapped Cartesian displacement, safe for current marked pairs but not a generic MIC implementation.

No production edits or model calculations. This entry check does not prove numerical/physical convergence, force accuracy or independent generalization.
