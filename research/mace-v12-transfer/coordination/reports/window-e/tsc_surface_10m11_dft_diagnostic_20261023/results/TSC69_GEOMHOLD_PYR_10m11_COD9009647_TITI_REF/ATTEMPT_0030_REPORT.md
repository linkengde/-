# 69-atom Ti3SiC2 (10-11) fixed-geometry SCF status

Date: 2026-10-23 UTC
Structure: `TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_REF`
Role: unassigned same-source facet diagnostic; not an independent DEV/TEST.

## Result

`SCF_CONVERGED = NO`; no valid energy, force, or stress label was produced.
There is no `RESULT.json` or labelled structure. Values in the SCF table below
are intermediate iterations only and must not be used as DFT labels. No
geometry, cell, PBC, k-point, cutoff, smearing, functional, PAW, eigensolver,
mixer, or convergence threshold was changed.

The first invocation reached iteration 20, then resumed from that checkpoint
with the corrected GPAW 26.7.0 restart API and completed ten additional SCF
iterations. The runner's total iteration ceiling remains 160. The resumed
attempt was safely stopped after the iteration-30 checkpoint was saved and
independently reloaded; no further SCF was run.

| Resumed iteration | Global iteration | Energy (eV; intermediate) | `eigst` log column | `dens` log column |
|---:|---:|---:|---:|---:|
| 1 | 21 | -532.983030 | -2.92 | — |
| 2 | 22 | -541.784012 | -1.53 | -1.90 |
| 3 | 23 | -532.573731 | -2.09 | -1.38 |
| 4 | 24 | -532.362572 | -3.39 | -2.00 |
| 5 | 25 | -532.325580 | -2.97 | -2.10 |
| 6 | 26 | -532.286754 | -3.41 | -2.12 |
| 7 | 27 | -532.268314 | -3.36 | -2.18 |
| 8 | 28 | -532.262143 | -3.39 | -2.24 |
| 9 | 29 | -532.285058 | -3.41 | -2.33 |
| 10 | 30 | -532.330724 | -3.13 | -2.35 |

The `eigst` and `dens` columns are GPAW's printed log10 residuals. They
oscillate and remain well above the locked `1e-5` convergence criteria; the
energy also changes over the resumed sequence. No convergence marker appears.
Near-stationary energy values do not establish convergence.

## Current recovery point

- Checkpoint: `checkpoint_iter_0030.gpw`
- Size: 4,583,533,958 bytes
- SHA-256: `387a825db5e277fd03f1e6946b25125627ff5106626741f03ef08fc1e9a5fabe`
- Checkpoint input SHA-256: `4b84cfc8f76827c7f81d9401f60e0c161355d61adc57cfe9ba867a78c780f473`
- Geometry SHA-256: `fd3351b1ba798ea9e0567013740fba8651d801fa8b62199b89da2c4d5ef37669`
- Config SHA-256: `15ddc8f596281298656a5a21df67b2915f0aebd5ba3bb18316fd2935ee8729c2`
- Four-rank read-only reload: PASS; 69 atoms, stored total maxiter 160,
  remaining iteration cap 130; SCF was not started by the reload test.
- Disk free after checkpoint commit: 11,995,127,808 bytes.
- Workspace free space after stopping: 11,995,115,520 bytes.
- Cgroup `memory.peak` read after stopping: 34,359,717,888 bytes against a
  34,359,738,368-byte limit. This is a cgroup-wide high-water mark, so it
  cannot be attributed to this run; no live GPAW process remained at the
  measurement. The current calculation process had been observed around
  17 GB between SCF steps. Do not treat the cgroup peak as a per-process RSS
  measurement.

The checkpoint is retained locally for a possible same-parameter continuation.
It is 4.58 GB and is not stored in the Git commit; remote files contain the
checkpoint identity and checksum only. The checkpoint can be lost if this
temporary workspace expires, so the current remote commit is not a binary
checkpoint backup.

## Next gate

Do not start the paired Ti displacement response while the reference remains
unconverged. Before any continuation, preserve the same 160-step ceiling and
all locked settings, verify this checkpoint SHA and the latest pointer, and
ensure adequate disk headroom for atomic replacement. If a subsequent
unchanged-parameter continuation still oscillates without useful residual
decay, stop this facet point and report it unresolved instead of repeatedly
spending full walltime. No DFT role is assigned by this diagnostic.
