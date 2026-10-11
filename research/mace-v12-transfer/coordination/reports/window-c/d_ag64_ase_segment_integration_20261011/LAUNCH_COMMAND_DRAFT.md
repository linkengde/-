# Blocked startup command draft

This is a syntax draft for the single D task only. It has not been run or sent to D. The current lock is unsigned and has `remote_archive_verified=false`; the adapter must exit HOLD. Do not change either flag locally.

```bash
cd /workspace/-/research/mace-v12-transfer/coordination/reports/window-c/d_ag64_ase_segment_integration_20261011
python ase_segment_adapter.py \
  --release-lock C_SIGNED_RELEASE_LOCK.json \
  --release-lock-sha256 <C_PUBLISHED_LOCK_SHA256> \
  --execution-code-sha256 e0bba3702c54bbdcb39f5f7e5641e1e1969a2ea5af8faceee8a0d50316f440dd \
  --owner D \
  --task-id C-D-AGL-1250K-rho097-seed01 \
  --start-step 1000 \
  --stop-at-step 2000 \
  --source-archive <D_FRESH_FETCHED_STEP1000_TAR_GZ> \
  --source-archive-manifest <D_FRESH_FETCHED_STEP1000_ARCHIVE_MANIFEST_JSON> \
  --source-archive-sha256 b62a818e58a2cf423e4cb1b9ed0468ae73ad578a75c4b2cb9fb3f28baf75e7b1 \
  --potential /workspace/WMH2006_Ag-local-validation/Ag.eam.alloy \
  --potential-sha256 958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447 \
  --output-root /workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot/segment-step1000-2000 \
  --execute
```

The placeholders for D's fresh-fetched archive paths and the future C-published lock SHA are unresolved. The source archive must be freshly fetched and verified in D's own session first. The current CLI also deliberately lacks a verified GitHub push/fresh-fetch backend and returns `VERIFIED_GITHUB_ARCHIVE_ORCHESTRATION_NOT_IMPLEMENTED` even if a signed-shaped lock is supplied. Therefore this command is not presently executable, and the code SHA above must not be interpreted as release approval.

Before any future launch, C must add and test the remote archive backend, obtain D's actual current archive paths/branch head, verify D's live resources and potential hash, update the code and tests, republish a new code SHA and lock SHA, and separately sign only this task and exact 1000-step bound.
