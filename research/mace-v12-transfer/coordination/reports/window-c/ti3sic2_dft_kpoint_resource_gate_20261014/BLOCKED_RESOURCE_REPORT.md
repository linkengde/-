# Blocked Resource Report — Stage C-D k-point pilots

```text
RESOURCE_GATE_PASS=NO
KPOINT_CONVERGENCE_PASS=NO
DFT_LABEL_GENERATION_READY=NO
SCF_STARTED=0
```

The two pilot SCFs were not started. The live 32-GiB cgroup had approximately 5.62 GiB free, its recorded peak was only 73,728 B below the cap, and `memory.max` had been reached 15 times. A prior project resource assessment had already marked this allocation unsafe for GPAW; there is no successful memory benchmark for the 217-atom interface. Starting only the smaller bulk point would not satisfy the requested two-system convergence gate and would consume resources despite the interface blocker, so neither SCF was launched.

No convergence conclusion can be drawn from these unrun comparisons. In particular, this report does not claim that either proposed k-point pair passes or fails numerically.

## Deferred pilot definition

| Structure | Meshes | Required comparison after a future resource release |
|---|---|---|
| Existing `TSC_VOL_M02` Ti₃SiC₂ bulk structure | 2×2×2 and 4×4×4 | Native energy/atom, force-consistent free energy/atom, maximum atomic force norm, full stress tensor |
| Existing `AGTSC_SI_GAP_2p4A` interface structure | 1×1×1 and 2×2×1 | Same metrics; stress remains diagnostic because of slab vacuum |

Only the mesh may differ within each pair. Keep the original fixed geometries, species, atom order, cell, PBC, cutoff, occupations, functional, PAW setups, mixer, and convergence criteria. Do not launch the full 15-point set based on these pilots alone.

## Minimum condition to reopen

Run the pilot only in an allocation with a fresh, task-specific resource preflight and a demonstrated safe peak-memory reserve for both structures, including GPAW scratch/checkpoint and recovery storage. Keep the four-rank MPI configuration if that remains the approved input. Do not lower precision or change the meshes to fit the current container.
