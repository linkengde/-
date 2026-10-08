# V18 clean30 controlled comparison

A reviewed B's V17 review and clean30 builder at main16bbe99. This cycle accepts the minimum fixed-recipe comparison:21 directly supported V16 frames,8 verified V17 signed diagnostic labels and1 matched frame0 relabel. The30-frame training proposal is copied byte-for-byte. The same3 V17 development frames are copied byte-for-byte; no test split is supplied. Composition rank4 is required.

Training keeps medium foundation, seed45, four CPU threads, batch1,80 epochs, lr1e-4, wd5e-7, energy100/force1000 and final epoch79 selection. All checkpoints are retained and the runner refuses a different exported epoch. Native extrapolated energy is used. Frame0 was relabelled at Gamma versus its historical1x4x2 declaration; k-point convergence is unresolved. Removing14 partial rows removes coverage too, so this comparison cannot prove those labels caused the error. Neither the six withheld labels nor thermal/pressure applicability are assessed.

Accepted B recommendations: fixed-recipe membership comparison, unchanged gates, explicit AgTi separation/max regression tracking, per-atom diagnostics in the next evaluation. Deferred: optimizer changes and new augmentation directions until the comparison and vector diagnosis are reviewed.

From the repository root:
```bash
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/build_v18_dataset.py --repo-root "$PWD"
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/build_v18_dataset.py --repo-root "$PWD" --check
bash research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/run_v18_training.sh --repo-root "$PWD" --output-dir research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01
```
The runner automatically freezes epoch79 and scores the three development frames on successful completion; it saves per-atom vectors for B's diagnosis. V18 remains a provisional numerical screen.

The successful finalizer also publishes only the selected model/epoch79 checkpoint, complete logs and development/per-atom evidence to main using ordinary Git and owner checks. Unrelated local changes, artifact checks or publication conflicts stop publication with results preserved. It does not launch another version or awaken an inactive B session.
