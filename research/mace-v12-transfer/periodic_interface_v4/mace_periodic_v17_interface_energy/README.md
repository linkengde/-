# V17 provisional data-only screening

This cycle preserves the V16 CPU training recipe and adds the eight archive-verified signed V17 diagnostic labels. It carries forward all 35 V16 training frames, including rows whose original calculation evidence is incomplete. The dataset and model are therefore marked **provisional screening**; a passing development score cannot establish a fully consistent PW-PBE training target or independent transferability.

`build_v17_dataset.py` writes 43 training frames and three V17 development-validation frames. It does not create a test split and does not open any V17 withheld-test output, label, summary, or score. The three development structures are reused for checkpoint validation/reporting and are not independent final tests. B has reconciled its aggregate `archive_manifest.json`; the builder requires it complete and also verifies all eight per-label archives.

The run uses the V16 optimization settings: MACE-MP-0b3-medium foundation, seed 45, CPU, batch size 1, 80 epochs, learning rate `1e-4`, weight decay `5e-7`, energy weight 100, and force weight 1000. This five-core instance uses four CPU training threads for throughput; `MACE_NUM_THREADS` can set 1–4. Thread count can alter floating-point reduction order slightly, so it is recorded with the run. There is no `--test_file` argument. The runner keeps epoch checkpoints and verifies that the exported model is the fixed final epoch 79; only then does A freeze a selection record and report development metrics. The six V17 withheld structures remain sealed.

From the repository root:

```bash
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/build_v17_dataset.py --repo-root "$PWD"
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/preflight_v17.py --repo-root "$PWD"
bash research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/run_v17_training.sh --repo-root "$PWD" --output-dir research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_03
```

A froze the final epoch79 model and completion evidence in `training_provisional_03/selection_record.json`; the completed run was scored with:

```bash
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/evaluate_v17_dev.py --repo-root "$PWD" --output-dir research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/evaluation_dev_01 --model research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_03/models/MACE_periodic_v17_provisional_screening.model --selection-record research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_03/selection_record.json
```

Run03 completed and is recorded under `training_provisional_03/`; development-only results are under `evaluation_dev_01/`. The current result is FAIL: energy passes all three structures, force-vector RMSE fails all three, and separation error passes only AgC. This is a numerical screening step only; it does not authorize MD, TTM, or physical claims. Any failed gate requires a targeted data/source or optimization diagnosis; it is not remedied by lowering the gate.
