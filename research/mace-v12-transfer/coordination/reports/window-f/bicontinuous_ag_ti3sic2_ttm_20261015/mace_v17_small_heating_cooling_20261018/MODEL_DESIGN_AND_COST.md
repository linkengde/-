# MACE V17 representative Ag–Ti₃SiC₂ model

## Structure and boundaries

This is an independent, layered interface slab used as a CPU-sized screening model. It is not copied from a historical melt-pool trajectory and it is not a geometric reproduction of the 39,486-atom bicontinuous architecture.

- Ti₃SiC₂: 1,200 atoms from the ideal COD 9009647 crystal, 10 × 10 basal repetitions and one conventional c-period, C-terminated (0001) surface.
- Ag: 726 atoms from ideal FCC Ag, (111) orientation, ABC stacking and six planes.
- Total: 1,926 atoms; composition is 66.6749 wt% Ag and 33.3251 wt% Ti₃SiC₂ (40:60 by atom count is not the mass ratio).
- Cell: `a=(30.575,0,0) Å`, `b=(-15.2875,26.4787,0) Å`, `c=(0,0,48.8579) Å`; XY periodic, Z nonperiodic, 18 Å top vacuum.
- A 2.70 Å nominal gap was set from the top C termination to the bottom Ag plane. Geometry-audit minimum cross-phase distance was 2.7048 Å.
- Ag's matched in-plane lattice parameter is 3.8907% compressive relative to ideal Ag(111). The match is a deliberate finite-cell coincidence, not a strain-free interface.
- The lowest Ti plane (100 atoms) is fixed as a mechanical support. All other Ti₃SiC₂ and Ag atoms can move. No large fixed buffer was used.
- Initial geometric connectivity at 3.25 Å for Ti₃SiC₂ and 3.10 Å for Ag was one component in each phase.

## Potential and resource settings

Only the frozen MACE V17 screening weight was used: `MACE_periodic_v17_provisional_screening.model`, SHA256 `559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3`. Runs used ASE with MACE 0.3.16, Torch 2.5.1+cpu, CPU float64 and four Torch threads. No EAM, TSCBOP or Morse term was added. The original 39,486-atom data and the V17 weight were not modified.

The environment has 4 CPU-equivalents, 32 GiB RAM, no CUDA device, and no installed OVITO application/module. `/workspace/.local/bin/lmp` is present but is an older build without ML-IAP. A separate ML-IAP/Python LAMMPS was built and smoke-tested; its CPU MACE bridge lacks the ghost-feature exchange required by MACE, while the Kokkos MACE path exits with an invalid-pointer abort. Those LAMMPS attempts produced no accepted force results. See `LAMMPS_MACE_ACCELERATION_AUDIT.md`. Extended XYZ trajectories are written for direct import into OVITO on a machine where OVITO is installed.

The existing 39,486-atom frozen V17 run-0 used about 29 GiB peak RSS in this same 32 GiB container. That static inference left too little headroom for a safe full-system trajectory with velocity, neighbor and integration buffers. Full-size MD was not attempted; this representative 1,926-atom slab is a layered-interface model, not the original bicontinuous architecture.

Measured resource costs:

| Operation | Work | Wall time | CPU/RAM |
|---|---:|---:|---:|
| Initial static force evaluation | 1,926 atoms, float64 | 16.1 s | peak RSS 17.41 GiB |
| FIRE screening relaxation | 60 steps | 935 s | peak RSS 17.64 GiB |
| 300 K early screen | 20 steps = 0.02 ps | about 1,602 s | about 5,843 aggregate CPU s; observed RSS about 17 GiB |
| Local heating, 1 fs | 90 steps = 0.09 ps | about 1,710 s | about 6,194 aggregate CPU s; observed RSS at least 16.2 GiB |
| Continued local heating, 0.2 fs | 100 steps = 0.02 ps | 1,789 s | 6,491 aggregate CPU s; peak RSS 17.65 GiB |
| Heat-off cooling, 0.2 fs | 100 steps = 0.02 ps | 1,799 s | 6,494 aggregate CPU s; peak RSS 17.67 GiB |

The three thermal segments together used about 19,179 aggregate CPU seconds (5.33 CPU-hours) over 5,298 wall seconds (88.3 minutes). These are measurements from this run, not a calibrated benchmark. The full 100-step 300 K baseline was stopped after 20 steps and is not a completed stability validation. Heating at 1 fs was stopped at 0.09 ps after a short-contact/high-force event, then continued from the saved state in a separate 0.2 fs branch; the short branch did not establish melting. A 0.02 ps cooling segment was completed after turning off the Ag heat source, but it is too short to establish a full cooldown.

The original heat-source window contained 92 Ag IDs. The 0.09 ps restart reselected 90 atoms from the moving circular window; the offline order/mobility analysis and final OVITO hot-region temperature use the fixed original 92 IDs. The cooling run's raw diagnostic window was also spatially reselected, but no Ag heat bath was applied during cooling; final metadata recomputes the fixed-ID temperature.

The target model composition is about 59.95 wt% Ag. The representative slab is 66.6749 wt% Ag and is layered rather than bicontinuous. Its trajectory supports interface-level exploration only; composition, morphology, and boundary conditions differ from the original model.
