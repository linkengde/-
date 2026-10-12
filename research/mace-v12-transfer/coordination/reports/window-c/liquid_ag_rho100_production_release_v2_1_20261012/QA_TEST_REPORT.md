# D production QA v2.1 integration report

Captured: 2026-10-12 UTC  
D owner branch: `cloud-d-b1-dft`  
QA lock: `848f8db72d83cea91dcf4cd42ee9fbbf49a80e7b5e3759cb4f535be8cc5873d2`

## Automated tests

`test_production_runner_v2_1.py`: **12 passed, 0 failed** under Python 3.12.14 and ASE 3.29.0. The tests cover per-frame historical Q6 comparison, report-only local Q6, the strict global-Q6 threshold and 80% window fraction, OVITO PTM FCC/HCP and periodicity, finite-cell RDF normalization and moving-block uncertainty, complete synthetic liquid assessment, insufficient/unresolved data HOLD behavior, 1000-step and 3600-second boundaries, automatic phase transitions and phase/global caps, archive-failure continuation blocking, step-1000 positions/momenta/RNG restoration, restored-boundary state/RNG continuity, and isolated local-Git archive recovery.

The source history contains 11 frames. The archived legacy implementation and the new production implementation agree per frame (maximum absolute difference recorded in `Q6_REGRESSION_REPORT.json`); the legacy-global mean is `0.08871920638262233` and the local-Q6 diagnostic mean is `0.3516369060394911`. The v2.1 gate uses global Q6 `< 0.20` in each final window. Local Q6 has no threshold or trend gate.

## Real GitHub archive integration

The production archive backend pushed a fabricated bundle to `cloud-d-b1-dft`, fetched the pushed commit in a new shallow checkout, restored all four members, compared each SHA and validated the synthetic checkpoint payload. All five archive flags passed. The bundle contains no real checkpoint, trajectory, potential, or simulation output. Test commit: `c4f954e0b00bce9a75e7cc120fe03625d88a2ba1`; source branch head: `57495a6b0a50dbd8dd23cdb02122f7e2eeb57ca9`.

## Live preflight and execution state

The raw `D_LIVE_PREFLIGHT.json` reports `PREPARED_HOLD`. The original step-1000 checkpoint SHA, positions/momenta and all three PCG64 streams restore; Ag EAM size/SHA and ASE alloy parser pass; QA v2.1 historical regression passes; environment, resource, output-path and real GitHub branch-head checks pass. The observed cgroup quota is 4 CPUs, memory limit 32 GiB, `/workspace` available space is about 28.3 GB, and free inode count exceeds two million.

The packaged config remains `execution_authorized=false`; no formal MD, DFT, AIMD, LAMMPS, TTM or MACE work was run. This report requests C's independent review and a separate release bound to the final D branch head and all listed SHAs.
