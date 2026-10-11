# D-only launch command (hash-pinned, not run by C)

D must first fetch the C integration package and verify `SHA256SUMS.txt`, then perform a live node/potential/source/output-root preflight. This command is only for the single task covered by `SIGNED_SEGMENT_EXECUTION_RELEASE.json`; it must not be edited to widen the scope.

Replace only the two source-file path placeholders with D's locally fresh-fetched and SHA-verified copies. Keep all other arguments and the output root exactly as signed.

```bash
python research/mace-v12-transfer/coordination/reports/window-c/d_ag64_archive_integration_20261011/production_segment_adapter.py \
  --release-config research/mace-v12-transfer/coordination/reports/window-c/d_ag64_archive_integration_20261011/SIGNED_SEGMENT_EXECUTION_RELEASE.json \
  --release-config-sha256 026be260ca65ab55b60c96a407c6181649e5f245bcf2d138993e246d2aa56d9b \
  --adapter-sha256 e855ad26cf215281e04b612163cb9db3c075ca8a6366686cafcec5becde264df \
  --backend-sha256 8fd001a9b1706fc8832e90a77c22783a22a77eab06f50a66bd1d1eeb97c00346 \
  --owner D \
  --task-id C-D-AGL-1250K-rho097-seed01 \
  --start-step 1000 \
  --stop-step 2000 \
  --source-archive <D_FRESH_VERIFIED_STEP1000_ARCHIVE> \
  --source-manifest <D_FRESH_VERIFIED_STEP1000_ARCHIVE_MANIFEST> \
  --potential /workspace/WMH2006_Ag-local-validation/Ag.eam.alloy \
  --output-root /workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot/segment-step1000-2000 \
  --execute
```

The adapter itself rejects an incorrect owner/task/boundary, release or code SHA mismatch, changed D branch head, incorrect source or potential SHA, occupied output root, and storage failure. No other task or step range is released.
