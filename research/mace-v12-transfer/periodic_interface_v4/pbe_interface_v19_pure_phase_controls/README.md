# Registered pure-phase numerical controls

A integrated the B draft, corrected production repository resolution and explicit unassigned-owner rejection; all 20 static tests passed before enabling. Only Ag_baseline_k2x2x2_v19 and Ag_baseline_k4x4x4_v19 are registered to new B. All remaining proposed records stay disabled and unassigned. Input bytes and method are unchanged. These are numerical controls, not independent validation or accepted training labels.

Run sequentially from repository root:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/run_queue.py Ag_baseline_k2x2x2_v19
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/run_queue.py Ag_baseline_k4x4x4_v19
```

Run the second only after the first verifies and publishes successfully. Existing outputs block repeated launches. Each queue publishes status, runs four MPI ranks with one thread each and automatically verifies compact archives before publishing; state.gpw stays local. Read mesh_plan.json for spacing and limitations. Equilibrium forces do not certify displaced/strained forces. Native extrapolated and free energies remain separate. Do not rerun disabled-fixture tests on this enabled production manifest; their preserved validation describes pre-enablement tests, not physical DFT readiness. No Ag DFT has been started by A.

## k6 follow-up registration

k2/k4 completed but energy budget failed. B is now assigned `Ag_baseline_k6x6x6_v19`; run production queue with this exact label, verify/publish then compare k4/k6. k8 remains disabled. This adds no new structure or target convention.

## Displaced-force mesh control

Baseline k4/k6 locally passes. B now runs Ag_displacement_p01_k4x4x4_v19 then Ag_displacement_p01_k6x6x6_v19 on identical +0.01Å x input. Verify and publish each, compare forces and within-mesh displacement response. Not stability or independent model validation. All other entries stay as previously registered/disabled.
