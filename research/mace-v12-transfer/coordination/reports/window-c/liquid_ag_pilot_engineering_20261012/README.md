# D rho100 engineering preflight

This package is bound to C commit `287b60af424f7bcc3a7b3598e524f5b543884ec9` and its `LIQUID_AG_PILOT_LOCK.json` SHA-256 `e1adaf12d542184c709d2eb2a5363cbfe9578c168ebc5e7ad260ee39a4faf09b`. The lock says `execution_authorized=false`; the package preserves that HOLD.

The engineering command verifies the C lock, Ag EAM identity, ASE 3.29.0 parser and units, periodic image offsets at cutoff greater than half the box, C's exact PCG64 packing reference, FixCom degrees of freedom, temperature logging, checkpoint SHA and complete three-stream RNG restoration, interrupted boundary saves, and storage/step/wall guards. It also tests the 3599/3600/over-3600-second gates, a completed integration step that crosses the wall limit, and exception-exit checkpoint saving between periodic boundaries. Its six integration steps use a temporary fixture under `/tmp`; they are not the pilot trajectory and are not copied into this release.

Run from the activated D environment:

```sh
source /workspace/onboarding/runtime.sh
cd /workspace/-/research/mace-v12-transfer/coordination/reports/window-c/liquid_ag_pilot_engineering_20261012
/workspace/.venvs/gpaw-mpi/bin/python -m unittest -v test_pilot_engine.py
/workspace/.venvs/gpaw-mpi/bin/python run_pilot.py \
  --mode engineering-preflight \
  --lock /tmp/LIQUID_AG_PILOT_LOCK.json \
  --lock-sha256 e1adaf12d542184c709d2eb2a5363cbfe9578c168ebc5e7ad260ee39a4faf09b \
  --potential-path /workspace/WMH2006_Ag-local-validation/Ag.eam.alloy
```

The formal runner requires a separate C-published authorization file and its externally supplied SHA-256. That authorization must bind the same task, lock, runner SHA, environment SHA, 1000-step ceiling, and 3600-second ceiling. The engineering command does not generate that file or set its authorization field. The current lock is rejected before formal output paths or coordinates are created.

The setfl metadata inspector reports `C_TAIL_REVIEW_REQUIRED` because the parser itself does not approve tail policy. The frozen C lock separately approves only the exact source-hash-bound one-grid-bin near-zero tail; this package checks that the parser's measured tail is within one bin and leaves the policy decision attributed to C.

No Ag potential, paper, formal structure, candidate data, or MD result is included in this package.
