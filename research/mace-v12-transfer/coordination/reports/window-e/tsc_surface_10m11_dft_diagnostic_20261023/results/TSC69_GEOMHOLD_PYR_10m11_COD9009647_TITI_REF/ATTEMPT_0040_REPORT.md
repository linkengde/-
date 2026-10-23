# Ti3SiC2 (10-11) reference SCF continuation through iteration 40

Date: 2026-10-23 UTC
Task: `TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_REF`

## Outcome

`SCF_CONVERGED = NO`. The run has produced no valid DFT label: no final energy,
force, or stress is available. The only energy values below are intermediate
SCF iterations. The response geometry remains blocked until this reference
converges. No geometry or scientific input parameter was changed.

This continuation began from the verified iteration-30 checkpoint and stopped
after the iteration-40 checkpoint was committed and read-only reloaded. The
locked maximum remains 160 total SCF iterations; 120 iterations remain in the
checkpoint state. The local SCF loop restarted its printed iteration counter
at 1, while the global count is 31–40.

| Resumed iteration | Global iteration | Energy (eV; intermediate) | `eigst` log10 residual | `dens` log10 residual | Energy criterion |
|---:|---:|---:|---:|---:|---|
| 1 | 31 | -532.331988 | -3.82 | — | not marked |
| 2 | 32 | -532.409492 | -2.85 | -2.38 | not marked |
| 3 | 33 | -532.327565 | -3.44 | -2.23 | not marked |
| 4 | 34 | -532.257696 | -4.70 | -2.35 | not marked |
| 5 | 35 | -532.254742 | -4.74 | -2.65 | not marked |
| 6 | 36 | -532.258965 | -4.48 | -2.68 | passed this iteration (`c`) |
| 7 | 37 | -532.256271 | -4.35 | -2.72 | passed this iteration (`c`) |
| 8 | 38 | -532.263706 | -4.30 | -2.77 | not marked |
| 9 | 39 | -532.260919 | -4.31 | -2.82 | not marked |
| 10 | 40 | -532.262537 | -4.74 | -2.96 | passed this iteration (`c`) |

The energy criterion passed on a few individual iterations, but the eigenstate
and density criteria did not pass, so the full SCF did not converge. Across
iterations 30–40 the density residual improved gradually, while the eigenstate
residual oscillated. This is insufficient to treat the result as a label.

## Verified recovery point

- File: `checkpoint_iter_0040.gpw`
- Size: 4,583,533,956 bytes
- SHA-256: `12a40c58d5f4ba08909911ad4540c0ab11440b12d97ca3f69b399e650b88b788`
- Input SHA-256: `4b84cfc8f76827c7f81d9401f60e0c161355d61adc57cfe9ba867a78c780f473`
- Geometry SHA-256: `fd3351b1ba798ea9e0567013740fba8651d801fa8b62199b89da2c4d5ef37669`
- Config SHA-256: `15ddc8f596281298656a5a21df67b2915f0aebd5ba3bb18316fd2935ee8729c2`
- Four-rank read-only reload: PASS; 69 atoms, stored total maxiter 160,
  remaining cap 120; no SCF was started by the probe.
- Workspace free space after stop: 11,994,951,680 bytes.

The cgroup high-water mark was 34,359,717,888 bytes against a
34,359,738,368-byte limit; `memory.events` showed `max=185` and `oom_kill=0`.
This is a cgroup-wide high-water mark rather than a per-process RSS figure.
The active GPAW ranks had been observed around 17 GB between iterations. The
near-limit high-water mark is a resource warning, so no further continuation
is started in this report.

The 4.58 GB checkpoint remains local and reloadable, but is not included in
the remote Git commit. The remote report and manifest record its identity only.
No TRAIN/DEV/TEST role is assigned; this COD 9009647 facet is not an
independent-source blind validation family.
