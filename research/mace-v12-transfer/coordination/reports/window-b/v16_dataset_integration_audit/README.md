# V16 dataset integration audit

This is a read-only audit of the integrated V16 dataset, six training-source
archives, and the V16 training/preflight entry points. It reads only the three
frozen V16 holdout **input geometries**; it never opens holdout calculation
outputs or reads their reference energy/force values. It does not run DFT,
training, model inference, MD, or TTM.

Reproduce the structural and provenance checks with:

```bash
/workspace/.venvs/gpaw-mpi/bin/python \
  research/mace-v12-transfer/coordination/reports/window-b/v16_dataset_integration_audit/audit_v16_dataset.py \
  --repo-root /workspace/- \
  --output-json research/mace-v12-transfer/coordination/reports/window-b/v16_dataset_integration_audit/audit.json
```

See `report.md` for findings, `audit.json` for the machine-readable checks, and
`SHA256SUMS.txt` for artifact hashes.
