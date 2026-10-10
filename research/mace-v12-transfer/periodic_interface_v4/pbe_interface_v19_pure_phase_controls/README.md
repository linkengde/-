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

## Ti3SiC2 baseline registration

A registered Ti3SiC2_baseline_k4x4x2_v19 and Ti3SiC2_baseline_k6x6x2_v19, sequentially after its active slab dipole-off job verifies/publishes. Fixed48-atom CIF bulk, all other method settings unchanged. Compare native/free energies and vector differences; baseline symmetry does not validate displaced forces.

## Ag symmetric displacement/volume controls

B registered Ag_displacement_m01_k6x6x6_v19, Ag_volume_m02_k6x6x6_v19, Ag_volume_p02_k6x6x6_v19 in that order. Existing positive/baseline reused. Volume k6 remains provisional: no broad convergence, optimized lattice or thermal stability claim.

## Ag volume numerical follow-up

B registered Ag_volume_m02_k8x8x8_v19 then Ag_volume_p02_k8x8x8_v19 for matched k6/k8 comparison. k8 baseline remains disabled; no claimed k8 energy-volume curve.

## Ag matched k8 volume curve

B now registered Ag_baseline_k8x8x8_v19 to complete same-mesh baseline with completed +/-2% k8 volumes. Check baseline k6/k8 then assemble baseline-subtracted k8 curve; not an optimized lattice or thermal-stability claim.

## Ag shear next task

B registered Ag_shear_m005_k6x6x6_v19 then Ag_shear_p005_k6x6x6_v19, only after assigned k8 baseline and reports finish. Shear k6 convergence remains provisional. No stress/elastic or thermal validation is claimed.

## Ti3SiC2 k8 follow-up

A registered Ti3SiC2_baseline_k8x8x4_v19 after k4/k6 vector difference0.014596eV/Å exceeded0.01. Compare to completed k6, no old repeats.

## Ag shear mesh follow-up

B registered Ag_shear_m005_k8x8x8_v19 and Ag_shear_p005_k8x8x8_v19. Compare unchanged inputs against completed k6 and relative energies against existing k8 baseline. No elastic/thermal validation claim.


A completed the same-geometry 6×6×2 / 8×8×4 comparison: energy budgets pass, while force-vector RMS (0.01359 eV/Å) exceeds the 0.01 eV/Å budget. A has registered 8×8×2 to isolate the in-plane mesh increment against the 6×6×2 label; do not interpret the earlier mixed-mesh comparison as a direction-specific failure.

## Current Ti3SiC2 in-plane mesh isolation

A completed the verified 6×6×2 / 8×8×4 comparison. Energy budgets pass; force-vector RMS is 0.01359 eV/Å and fails the 0.01 eV/Å budget. A registered the same baseline at 8×8×2 to isolate the in-plane increment at fixed kz=2. The nine interface gap-scan records remain disabled pending slab vacuum/dipole controls and reference-target review.

The 8×8×4 local `state.gpw` restart file (7,361,378,540 bytes) was removed only after both compact archives passed `verify_result.py` and no GPAW process was active. Its path and SHA256 are recorded in the A mesh-comparison report. Compact archives remain on main; 8×8×2 is ready to launch.

## Ti3SiC2 reciprocal-z isolation control

B is assigned one 6×6×4 calculation on the same 48-atom baseline input used by A’s 8×8×2 job. Compared with the archived 6×6×2 baseline, this fixes the in-plane mesh and raises only kz; it complements A’s fixed-kz in-plane test. Compare native and free energies separately and all-atom force-vector RMS/max, with the unchanged gates of 2 meV/atom and 0.01 eV/Å. This is a static numerical control, not validation or thermal stability. Run only after verifying the B registration and no active B DFT job.

## Ti₃SiC₂ backbone force-response batch

B registered four source-audited, 48-atom fixed-cell z-displacements against the completed `Ti3SiC2_baseline_k10x10x4_v19` method. The canonical method hash is `702df71e586a61cd5785afc62653d7135a707032c91b72842b035b62198fd709`; PBE, 500 eV PW, 0.10 eV Fermi–Dirac smearing, 10×10×4 k-points, PBC TTT, GPAW 26.7.0, ASE 3.29.0, and gpaw-data 1.2.1 are unchanged. Native extrapolated energy and free energy remain separate. These are candidate backbone response data, not independent validation or training labels.

Inputs are copied byte-for-byte from the checked proposal report and owned by `window-b` / `5f182750-e2d3-428b-bddc-5ab37bf09e5e`. Run one label at a time in this order, verifying and publishing each archive before the next:

- `Ti3SiC2_Ti4f_id0003_z_m020A_k10x10x4_v19` — `inputs/Ti3SiC2_Ti4f_id0003_z_m020A_PROPOSAL.extxyz` — SHA-256 `5529a040f0f35957f826b78754910970741bff741474e321dc1646a5c0c12bbf`
- `Ti3SiC2_Ti4f_id0003_z_p020A_k10x10x4_v19` — `inputs/Ti3SiC2_Ti4f_id0003_z_p020A_PROPOSAL.extxyz` — SHA-256 `5829ebcbb9ae331f6cf1dd7dc5abca0787c07667a55778e3c1eea241938b227b`
- `Ti3SiC2_C4f_id0009_z_m020A_k10x10x4_v19` — `inputs/Ti3SiC2_C4f_id0009_z_m020A_PROPOSAL.extxyz` — SHA-256 `19f1a4011620f5c7f245d0a8c893964aa8ebd1eb647f8a9b6f4547a75cac586f`
- `Ti3SiC2_C4f_id0009_z_p020A_k10x10x4_v19` — `inputs/Ti3SiC2_C4f_id0009_z_p020A_PROPOSAL.extxyz` — SHA-256 `6308961f33ad5fae3f58930570a47044a236e062de243821c20577c78955ff72`

Stop at the first failed or unconverged run and preserve its checkpoint. Do not run interface DFT; source verification and A’s slab method review remain pending.
