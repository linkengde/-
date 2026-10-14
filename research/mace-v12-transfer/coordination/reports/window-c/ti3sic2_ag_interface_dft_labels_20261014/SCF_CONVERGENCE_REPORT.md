# SCF convergence report

No SCF calculation was started. There are therefore no iteration counts, convergence flags, total energies, forces, or stress tensors for these 15 structures.

| Group | Planned points | SCF started | Converged | State |
|---|---:|---:|---:|---|
| A — Ti₃SiC₂ bulk | 6 | 0 | 0 | `HOLD_RESOURCE_GATE` |
| B — Ag/Ti₃SiC₂ interface | 9 | 0 | 0 | `HOLD_RESOURCE_GATE;KPOINT_CONVERGENCE_PENDING` |

The full per-task statuses are listed in `STRUCTURE_RESULTS.csv`; `ENERGY_FORCE_STRESS_SUMMARY.csv` deliberately leaves label columns blank. A future run must preserve separate native total energy and force-consistent free energy, complete atomic forces, stress tensor, SCF convergence evidence, iteration count, logs, and hashes. An unconverged run must remain failed/incomplete and must not enter a fit dataset.
