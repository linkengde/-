# Reproduce V17 paired response analysis

This package reviews four matched positive/negative V17 diagnostic acquisition pairs with the frozen V16 model. It reads only the eight diagnostic archives and V14 training parent frames 28, 29, and 30. It does not access V17 withheld-test outputs or run DFT, training, MD, or TTM.

The recorded inference environment is Python 3.12.14, ASE 3.29.0, PyTorch 2.5.1 CPU, and mace-torch 0.3.16. Run from the repository root:

```bash
mkdir -p /tmp/window-b-mpl-cache
MPLCONFIGDIR=/tmp/window-b-mpl-cache /workspace/.venvs/mace-v14-assist/bin/python \
  research/mace-v12-transfer/coordination/reports/window-b/v17_paired_response_analysis/analysis.py \
  --repo-root /workspace/- \
  --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v17_paired_response_analysis/analysis_reproduction_01 \
  --model /workspace/-/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy/training_01/checkpoints/MACE_periodic_v16_interface_energy_run-45.model \
  --selection-record /workspace/-/research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy/selected_model_record.json
```

The script refuses a non-empty output directory. It rechecks the V16 model and selection record, V14 dataset hash, acquisition manifests, archive verification files, atom identities, structures, declared perturbations, and finite DFT/MACE values. It writes `analysis.json` plus five compact CSVs. The report documents how to interpret finite secants, native/free energies, force projections, the source-parent geometry-hash discrepancy, and limits on validation claims.
