# Memory Requirement Estimate — Stage C-D

## Known measurements and limits

| Evidence | Measurement | Interpretation |
|---|---:|---|
| Current cloud cgroup limit | 34,359,738,368 B (32 GiB) | Hard memory ceiling in the audited container |
| Current use at preflight | 28,327,940,096 B | About 6,031,798,272 B (5.62 GiB) headroom at that instant |
| Historical cgroup peak | 34,359,664,640 B | Only 73,728 B below the limit |
| `memory.events` | `max=15`, `oom=0`, `oom_kill=0`, `oom_group_kill=0` | Allocation pressure has reached the cgroup limit; no recorded OOM kill in this snapshot |
| Previous 394-atom interface attempt | Exceeded a 32-GiB environment during initialization | Confirms that 32 GiB was insufficient for that larger model, not the exact requirement for the 217-atom model |
| 217-atom interface peak | No successful matched benchmark | Unknown; cannot be inferred safely by scaling atom count |
| 48-atom Ti₃SiC₂ peak | No peak-memory measurement bound to the proposed pilot meshes | Runtime analogues exist, but they do not prove safe memory headroom |

## Assessment

The exact peak RAM requirement for either requested pilot is **unknown**. In particular, no successful measurement exists for the 217-atom interface at 500 eV with the proposed 1×1×1 and 2×2×1 meshes. The current 32-GiB cgroup is not an acceptable execution target: only about 5.62 GiB remained at the snapshot, the cgroup had previously reached its limit 15 times, and its historical peak nearly touched the hard ceiling.

Do not convert the 394-atom failure into an exact 217-atom RAM estimate, and do not infer that the 217-atom system will fit merely because it is smaller. GPAW plane-wave memory depends on cell volume, cutoff, k-point distribution, parallel layout, wavefunctions, eigensolver workspace, and checkpoint/scratch behavior.

## External resource planning

Before release, the external operator must provide a fresh task-specific preflight with the actual memory limit, current memory use, peak-memory monitoring, MPI layout, disk/scratch budget, and persistent checkpoint/restore space. The pilot must stop safely if memory reserve approaches the allocation limit; it must not reduce precision or change the approved method to fit.

No guaranteed minimum RAM is claimed here. A 64-GiB-class host may be used as a **provisional planning target** for arranging an initial measured pilot, but this is not a demonstrated sufficiency threshold or execution release. The 217-atom run must be monitored and its peak recorded. If the available host cannot provide an adequate reserve for both the calculation and recoverable outputs, keep `WAITING_FOR_COMPUTE_RESOURCE=YES` and run no SCF.

## Disk and time caveats

The prior audit observed about 22.07 GB free on `/workspace` overlay. This is not persistent storage proof and does not establish enough capacity for all GPAW scratch/checkpoint files plus archive recovery copies. The analogous 48-atom baseline ran for 3.26 h at 10×10×4 on four MPI ranks, and four related response points took 2.83–4.20 h; these are not wall-time predictions for the proposed meshes. No reliable wall-time or disk estimate exists for the 217-atom interface pilot without a measured run.
