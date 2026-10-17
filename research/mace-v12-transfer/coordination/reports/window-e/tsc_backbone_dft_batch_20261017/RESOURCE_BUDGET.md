# Resource and execution-budget review (offline)

## Current Codex Cloud allocation

Read-only checks on 2026-10-17 found a cgroup CPU quota of 4 cores (`cpu.max=400000 100000`), cgroup memory limit 34,359,738,368 bytes (32 GiB), and `/workspace` free space 16,659,927,040 bytes. The workspace had about 2.02 million free inodes. Host-reported available memory was about 30.36 GB at inspection. These are task-overlay resources, not persistent-volume guarantees.

The existing GPAW environment is Python 3.12.14, GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1; `gpaw info` reports MPI enabled, ScaLAPACK enabled, GPU disabled. The setup directory is `/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups`. A four-rank import test confirmed every rank saw GPAW 26.7.0, ASE 3.29.0, the same hostname and setup path. No DFT calculator was initialized.

PAW hashes were re-read from the active setup directory and match the archived batch values: Ti `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f`, Si `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0`, C `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983`.

## Estimated single-point resources

Each prepared slab has 100 atoms, a cell volume of 1,724.323 Å³, a 6×6×1 mesh and 500 eV cutoff. The audited 48-atom bulk run at 6×6×2 had volume 570.710 Å³, 384 valence electrons, 235 bands, 72 k-points, 4 MPI ranks, and GPAW log storage of 1,923,871,680 bytes for wave functions plus 597,628,224 bytes for projectors. Its converged 41-step result took 4,194.7 seconds.

For a rough basis-memory scale only, plane-wave count scales approximately with cell volume at fixed cutoff; occupied bands scale with valence electrons; and the k-point count halves. This gives a wave-function storage estimate of approximately `1.924 GB × (1724.323/570.710) × (816/384) × (36/72) = 6.18 GB`. GPAW's projector array scales approximately with atom/projector count, bands and k-points, giving `0.598 GB × (100/48) × (816/384) × (36/72) = 1.32 GB`. These are not peak-RSS estimates; they omit FFT, density, eigensolver, mixer, Python and allocator overhead. No 100-atom GPAW initialization or peak-RSS measurement was performed. **Memory fit within 32 GiB is plausible but not proven; C should treat this as conditional, not PASS.**

Using the same coarse work-size scaling, a 41-iteration surface run could approach 3.2× the 48-atom reference runtime (roughly 3.7 hours). SCF iteration counts and slab convergence are unknown. A 4-hour per-point cap would leave little margin and may time out. No platform-level job wall-clock allowance was found, so the exact wall-time limit remains unverified.

## Storage and interruption policy proposal

The four input extxyz files together are small (under 100 kB). Labelled outputs are expected to be small, while a GPAW checkpoint/scratch file is unmeasured. Reserve at least 10 GiB for one active task's temporary checkpoint and scratch; 16.66 GB is currently free, but only about 6.66 GB would remain as headroom. Run one task at a time. If a task completes and labels pass validation, do not retain a large recovery checkpoint. If interrupted, retain a temporary checkpoint only for that task and continue from it only if the reviewed runner/checkpoint method supports a validated restart. The checkpoint mechanism has not been tested for these 100-atom slabs. Do not clean historical files to make room.

Persistent storage is not verified beyond this task workspace. To preserve results, publish the completed report/labels to the E branch after C-approved execution; do not describe the task overlay itself as durable storage.

## Gate

`RESOURCE_REVIEW=CONDITIONAL_HOLD`: the geometry/method estimate is available, but there is no measured slab peak memory, exact checkpoint size, platform wall-time evidence, or validated interruption/restart procedure. C must explicitly accept the estimates and set a per-point time/space policy before execution.
