# Stage C-D7 — M77 2.40 Å DFT resource pilot preflight

## Outcome

```text
DFT_STARTED=NO
SCF_INITIALIZATION=BLOCKED_BEFORE_CALCULATOR_INITIALIZATION
KPOINT_DESCRIPTOR_CHECK=PASS
OUTPUT_SCHEMA_CHECK=PASS_STATIC_ONLY
RESOURCE_GATE=FAIL
DFT_EXECUTION_READY=NO
```

The requested reference input is the D6 77-atom, 2.40 Å fixed-geometry structure. This review did not run an SCF or create DFT outputs. The current memory gate fails, so the review stopped before GPAW calculator/SCF initialization and did not lower the cutoff, alter the k mesh, or change the structure.

## Bound input identity

| Item | Bound value |
|---|---|
| Task | `C-D6-GAP-2.40A` |
| Structure | `Ag_Ti3SiC2_M77_gap_2p40A.extxyz` |
| Structure SHA-256 | `1de937b99b2d1f5b4be3b9b45e4d13e4e634bd99531110a080abf7ce970319dd` |
| Config | `M77_SHARED_STRAIN_GAP_2.40A.json` |
| Config SHA-256 | `5ab320dec595b45a345d28e24f15ec227346f64b28e11e24eb39124ba914c629` |
| D6 runner SHA-256 | `542d1ddfd8af8bb6f324d23815e2fe6274c4e7fc0162ac4d9fc734ba380028d9` |
| Geometry | 77 atoms; Ag₂₈Ti₂₁Si₁₄C₁₄; fixed cell and coordinates; PBC TTT |
| Frozen candidate method | PBE/PAW, PW 500 eV, Fermi–Dirac 0.10 eV, 4×4×1, SCF tolerances 1e-5, Pulay beta 0.05 / nmaxold 8 / weight 100; maxiter 160; 4 MPI ranks / 1 OMP thread |
| Authorization | `execution_authorized=false`, `launch_enabled=false` |

No geometry or input file was modified. The signed execution gate remains absent.

## Current environment snapshot

| Check | Observed | Assessment |
|---|---:|---|
| CPU quota | cgroup `400000/100000` = 4 CPU equivalents; 5 CPUs visible | Passes the configured four-rank count as a quota, not proof of dedicated physical cores |
| cgroup memory limit | 34,359,738,368 B (32 GiB) | Hard ceiling |
| cgroup memory current | 28,333,223,936 B | Only 6,026,514,432 B (~5.62 GiB) headroom at preflight |
| cgroup memory peak | 34,359,664,640 B | 73,728 B below the limit historically |
| cgroup memory events | `max=15`, `oom=0`, `oom_kill=0` | Repeated allocation-limit events; no OOM kill in this snapshot |
| Host-visible MemAvailable | 16,984,948,736 B | Does not supersede the tighter cgroup limit |
| `/workspace` free space | 22,067,056,640 B | Temporary overlay; persistent scratch/restore is not established |
| Free inodes | 2,017,864 | Adequate for this package; does not clear DFT scratch |

The historical resource report also recorded a 32-GiB environment with about 5.62 GiB available, a near-limit peak, and 15 `memory.max` events; the 217-atom peak remains unmeasured. There is no measured M77 GPAW peak-memory value. With only about 5.62 GiB currently free under the cgroup cap, and prior memory-limit events, this host does not have defensible headroom for an SCF resource pilot. This is an evidence-based hold, not a claim that every 32-GiB host is insufficient.

## Software, PAW, and k-point checks

Read-only environment inspection reports Python 3.12.14, ASE 3.29.0, GPAW 26.7.0, `gpaw-data` 1.2.1, Open MPI 5.0.7. `gpaw info` reports MPI enabled, ScaLAPACK enabled, LibXC 5.2.3, and PAW data installed. The four setup files exist and match the D6 expected hashes: Ag `bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97`; Ti `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f`; Si `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0`; C `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983`.

The GPAW k-point descriptor accepted the frozen integer mesh `[4,4,1]` and materialized 16 full-zone points with zero extra offset. This verifies mesh construction only. It does not prove a GPAW SCF can initialize under the current memory limit, that four-rank k-point distribution succeeds, or that 4×4×1 is numerically converged. No MPI DFT job or SCF was started.

## SCF initialization and output completeness

SCF/calculator initialization was **not attempted** because the resource gate failed first. The D6 runner's static output-path audit confirms it is written to request native energy, force-consistent free energy, all atomic forces, full stress, SCF convergence/iteration count, elapsed time, and a result extxyz plus summary JSON. The stress is explicitly diagnostic for this vacuum-containing slab. This is a source/schema check only: there is no actual result archive or result-integrity test from a DFT run.

## Release decision

Do not start the 2.40 Å SCF pilot in this current container. Reconsider only after a fresh host-specific preflight demonstrates adequate reserved memory and persistent scratch/checkpoint/restore capacity, then separately authorize the single fixed-geometry point. Keep D6 input hashes unchanged. Do not reduce precision, cutoff, or the configured k mesh to fit this host.
