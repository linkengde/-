# V17 provisional data-only screening

This cycle preserves the V16 CPU training recipe and adds the eight archive-verified signed V17 diagnostic labels. It carries forward all 35 V16 training frames, including rows whose original calculation evidence is incomplete. The dataset and model are therefore marked **provisional screening**; a passing development score cannot establish a fully consistent PW-PBE training target or independent transferability.

`build_v17_dataset.py` writes 43 training frames and three V17 development-validation frames. It does not create a test split and does not open any V17 withheld-test output, label, summary, or score. The three development structures are reused for checkpoint validation/reporting and are not independent final tests. The B-owned `partial_archive_manifest.json` currently lags one verified AgTi per-label archive; the builder verifies each of the eight per-label archives directly and records this aggregate-manifest discrepancy.

The run uses the V16 settings: MACE-MP-0b3-medium foundation, seed 45, CPU, batch size 1, 80 epochs, learning rate `1e-4`, weight decay `5e-7`, energy weight 100, and force weight 1000. There is no `--test_file` argument. The fixed final epoch is reported on the three development structures; the six V17 withheld structures remain sealed until A freezes a checkpoint and selection record.

From the repository root:

```bash
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/build_v17_dataset.py --repo-root "$PWD"
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/preflight_v17.py --repo-root "$PWD"
bash research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/run_v17_training.sh --repo-root "$PWD" --output-dir research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_01
```

After inspecting the full 80-epoch log and final checkpoint, A records the selected model hash and completion evidence in `selection_record.json`, then runs:

```bash
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/evaluate_v17_dev.py --repo-root "$PWD" --output-dir research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/evaluation_dev_01 --model research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_01/models/MACE_periodic_v17_provisional_screening.model --selection-record research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/training_provisional_01/selection_record.json
```

This is a numerical screening step only. It does not authorize MD, TTM, or physical claims. Any failed gate requires a targeted data/source or optimization diagnosis; it is not remedied by lowering the gate.
