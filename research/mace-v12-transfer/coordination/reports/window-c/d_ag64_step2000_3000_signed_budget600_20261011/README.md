# Signed manual segment release

Task C-D-AGL-1250K-rho097-seed01: start2000, end3000, maximum1000 added steps. Walltime600s per invocation. No automatic continuation; further execution requires manual review and a new release. C did not start MD.

CONFIG_SHA256: `cc330c682e0052d2442dee3d8a6945abedaa0e2673b5af860be058c2ba8ea229`

Source D branch head pinned to bc3ff0631db3b7d78dcc17556befa98c2f454480, checked before signing. Source state and archive pins, fixed potential SHA, coordinates, momenta, RNG and original scientific schedule unchanged. Still1600K hot stage; not final1250K sampling.

Budget-only new executor version preserves the earlier integration review and archive test evidence. Signed adapter/core validation and negative scope/parameter tests pass. Exact required runtime Python3.12.14, ASE3.29.0, NumPy1.26.4 is enforced before MD; D must perform fresh resource, potential and remote-channel preflight. Runtime digest records the required reference versions, not a new live D resource measurement.

600 seconds may stop before3000. Save the last fully completed step, archive and verify; do not restart from2000 after a partial segment. The MD-loop timeout does not bound slow individual force calls or remote backup duration. Scientific outputs remain separate and immutable; original archives untouched. Other tasks, DFT and training remain HOLD.

From this package directory using the matching existing Python environment:

```bash
python continuation_adapter.py --release-config SIGNED_CONTINUATION_CONFIG_STEP2000.json --release-config-sha256 cc330c682e0052d2442dee3d8a6945abedaa0e2673b5af860be058c2ba8ea229 --potential /workspace/WMH2006_Ag-local-validation/Ag.eam.alloy --output-root /workspace/sampler-work/cloud-d/wmh2006-ag64-sixpilot/segment-step2000-3000-budget600-v1 --preflight
```

After preflight PASS, D may replace `--preflight` with `--execute` for this segment only. This command was not executed by C.
