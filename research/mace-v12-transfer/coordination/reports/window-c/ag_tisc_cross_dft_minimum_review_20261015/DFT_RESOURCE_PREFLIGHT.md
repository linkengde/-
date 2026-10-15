# Current-resource and QE preflight — 2026-10-15

## Actual snapshot

Machine-readable values and timestamp are in `AUDIT_EVIDENCE.json`. The snapshot was collected during a read-only audit, not a DFT run.

| Item | Current observation | Decision |
|---|---|---|
| CPU | cgroup 400000/100000 = 4 CPU equivalents | Suitable count for a proposed 4-rank pilot; no claim of exclusive physical cores |
| Memory limit | 34,359,738,368 B = 32 GiB | Hardware ceiling |
| Memory used/headroom | 388,345,856 B / 33,971,392,512 B (about 31.64 GiB headroom) | Current headroom has improved; do not reuse yesterday's 5.62 GiB headroom as today's reading |
| Cgroup peak/events | 388,608,000 B peak; max=0, oom=0, oom_kill=0 | This snapshot has no recorded limit event; not a measured QE workload peak |
| Work-disk free | 22,066,921,472 B, overlay filesystem | Capacity/persistence for QE scratch plus archive/restore copies is not established |
| MPI | Open MPI 5.0.7, `/usr/bin/mpiexec` present | Launcher present; no QE MPI launch/compatibility run performed |
| QE `pw.x` | Not found in PATH or checked standard locations; no installed `quantum-espresso` package found | **BLOCKED** |
| QE UPF files/version/hash | Not supplied or verified | **BLOCKED** |
| GPAW | Existing 26.7.0 environment, ASE 3.29.0 and verified archived PAW files | Useful for historical verification; not an automatic substitute for requested QE |
| Exact new input / parent manifest | E handoff not present | **BLOCKED** |

No software was installed and no calculator, SCF, MD, fitting or force evaluation was launched. Current memory availability alone does not prove that a 135-atom or reduced QE model will fit.

## Cost evidence and planning budget

The six verified historical GPAW points took 0.707–1.834 hours each on four ranks, totaling 7.400 hours (29.601 core-hours). They contain 48–51 atoms in large cells. These measurements are **not** a QE or 135-atom wall-time benchmark; neither atom-count scaling nor code-name matching can produce a reliable absolute forecast.

For preliminary scheduling only, a future reviewed pilot could be given a **4-hour wall-time budget** on four ranks, with an earlier internal safe-stop margin and monitored memory. This is a proposed budget cap, not a prediction, scientific parameter change, or execution authorization. At that cap, the conditional 4–7 new-input first batch reserves 16–28 wall-hours (64–112 core-hours); E's reported ten-point batch would reserve 40 wall-hours (160 core-hours), excluding additional convergence controls. A point may time out without a usable label. Replace these caps with measured per-system timing and SCF iteration rates after the first valid pilot; do not promise completion within them.

No absolute peak-RAM, QE wavefunction/checkpoint size, or total scratch estimate is defensible without the actual UPF files, cutoffs, k-point layout, number of bands and geometry. The 77-atom/217-atom proxy in prior reports cannot certify a 135-atom QE memory requirement.

## Restart/storage conditions for a later release

Pin QE version/build, MPI layout, pseudopotential hashes, complete input and structure hashes, rank count and thread count. Preserve the full consistent QE `<prefix>.save` and any necessary wavefunction files, not only the final XML or a partial directory. Record atom IDs, cell/PBC, charge/spin, global calculation state, elapsed time and native units. A same-input restart must be demonstrated on non-production engineering data or an explicitly authorized pilot.

Use QE's supported safe wall-time stop (`max_seconds`/graceful exit mechanism for the verified version), leave time for checkpoint and archive completion, and validate the saved state before restart. An external kill alone is not a restart guarantee. Checkpoint/storage failure must stop subsequent calculations. Archive failed/unconverged diagnostics explicitly as incomplete, with no valid energy/force label claim.

Before execution, budget scratch + checkpoint + local archive + fresh remote restoration simultaneously. Ordinary Git's conservative 50 MiB single-file limit still applies. If a necessary QE file exceeds it, use a separately approved and tested durable storage channel; do not force-push large binaries or delete the only valid copy. None of those capacities or recovery checks is established by this preflight.
