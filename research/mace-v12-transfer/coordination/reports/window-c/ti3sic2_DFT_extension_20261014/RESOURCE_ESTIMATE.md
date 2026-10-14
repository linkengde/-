# Stage C resource estimate and preflight

## Current environment snapshot

Read-only preflight on 2026-10-14:

| Resource | Observation | Assessment for these DFT points |
|---|---:|---|
| CPU quota | 4 CPU-equivalent (`cpu.max=400000/100000`) | Matches historical four-rank runs, but does not imply physical-core exclusivity. |
| Memory limit | 34,359,738,368 B (32 GiB) | Hard cgroup cap. |
| Current cgroup use | 28,351,229,952 B; ~6.01 GB under the cap | Insufficient evidence of a safe DFT reserve in this live environment. |
| Cgroup high-water mark | 34,359,664,640 B, only 73,728 B below cap; `oom=0`, `oom_kill=0`, `max=15` | The environment has previously approached its limit; do not launch GPAW here. Much of current charge is file/shmem cache, but that does not certify headroom for GPAW's peak allocations. |
| `/workspace` disk available | 22,074,118,144 B | Enough for the small candidate structure files; not certified for large GPAW checkpoints plus concurrent scratch and recovery copies. |
| Free inodes | 2,018,535 | Sufficient for this small package. |
| Software | Python 3.12.14, GPAW 26.7.0, ASE 3.29.0, NumPy 2.5.3, SciPy 1.18.1 | Version preflight passes; it is not a DFT resource release. |

## Wall-time range (planning only)

The verified 48-atom Ti₃SiC₂ baseline at 10×10×4 took 3.26 h on four MPI ranks; the four Ti/C fixed-cell responses took 2.83–4.20 h. These are the closest measured analogues, not guarantees for new configurations.

| Batch | Points | Planning range, sequential on 4 ranks | Basis |
|---|---:|---:|---|
| A volume/shear + B Si response | 6 | ~18–30 wall-hours | 3–5 h/point extrapolated from the archived 48-atom cases. |
| C 100-atom clean slabs | 4 | ~32–96 wall-hours | ~8–24 h/point, an intentionally broad unbenchmarked estimate for 100 atoms, a 55.247 Å cell, and 10×10×1. Must be replaced by an approved measured pilot. |
| D independent DEV/TEST | 4 | Not estimable until parent cells and meshes exist | Atom count, cell and k mesh are unavailable. |
| Current 10 source-backed points | 10 | ~50–126 wall-hours | Sum of the planning ranges above, before retries. |
| Full 14-point plan | 14 | At least the above plus an unestimated D batch | Do not treat as a fixed budget or guarantee. |

Ranges assume serial point scheduling on the same four-rank allocation. Wall-clock time for distributed execution depends on actual independent nodes and queue limits. Slab memory has no matched archived benchmark, so time estimates do not establish memory feasibility.

## Resource gate

`RESOURCE_AVAILABLE=NO` for starting DFT in this environment. Although disk and software preflight pass for the planning files, the live memory headroom is ~6.01 GB and the recorded peak nearly equals the 32 GiB cap. There is no current per-task checkpoint/scratch budget or isolated DFT allocation attached to this request. Do not start even the 48-atom calculations here; do not lower cutoff, mesh, smearing or convergence criteria to fit. Before a future release, capture a fresh host/cgroup snapshot, reserve adequate memory and scratch for the exact point, confirm persistent result storage and checkpoint policy, and run only the separately approved pilot.

