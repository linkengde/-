# REGISTRY_PHASE1_RECOVERY_REPORT

Status: **PHASE1_RECOVERY_COMPLETE**. R1's isolated diagnostic copy converged in 105 SCF steps; the original R1 screening result remains unchanged and failed at its original 160-step limit.

Source batch: E branch commit `ab364698dca879284be61c66f0340b02a6aaccd7`.

## Scope and integrity

- R2/R4 were recovered from their existing `gpaw.log` files only. Neither was recalculated.
- R1's diagnostic uses a separate copy under `R1_diagnostic_case/`; the original R1–R4 `input.py` and `structure.xyz` remain unchanged. Their four input and structure SHA-256 pairs were rechecked against the recorded source hashes and all matched.
- The diagnostic retains PBE/PAW, 500 eV, 6×6×1, 0.10 eV Fermi–Dirac smearing, fixed cell and coordinates, and four MPI ranks / one OpenMP thread. Only `maxiter=200` and `Mixer beta=0.10` differ from the original R1 calculation.
- No 8×8 calculation, geometry optimization, registry edit, or data-role change was made.

## R2/R4 log recovery

| Registry | SCF | Steps | Total energy (eV) | Force-consistent free energy (eV) | Forces | Stress |
|---|---:|---:|---:|---:|---|---|
| R2 | Converged | 103 | -142.621455 | -142.792091 | 20 atoms recovered in original atom order | `STRESS_UNAVAILABLE` |
| R4 | Converged | 124 | -120.079642 | -120.215165 | 20 atoms recovered in original atom order | `STRESS_UNAVAILABLE` |

The GPAW logs contain the energy summary and full force table. Stress had been evaluated in memory before the result-reader exception, but no stress tensor was persisted in either a result file or log; no value is reconstructed or substituted. The original readout exception was an attempted `ASECalculator.scf` attribute access after calculation completion.

Complete recovered values and force vectors: [`R2_R4_LOG_RECOVERY.json`](R2_R4_LOG_RECOVERY.json).

## R1 diagnostic copy

Test: `maxiter=200`, `Mixer(beta=0.10, nmaxold=8, weight=100)` on the SHA-verified R1 copy. The test is isolated and does not replace the original R1 result.

**Outcome: CONVERGED in 105/200 steps.** The three GPAW convergence markers passed. At the final step, energy was `-120.155503 eV`, the eigensolver residual was `log10(error)=-6.86`, and the density residual was `log10(error)=-5.08`. The final three-cycle energy span was `4×10⁻⁶ eV`. The SCF history records large transient energies early in the run (minimum `-708032.776136 eV`); these are retained as observed log values, not treated as physical results. The run settled to the converged energy above.

The diagnostic result includes free energy `-120.294555 eV`, forces for all 20 atoms, stress, and `geometry_modified=false`: [`R1_diagnostic_case/run_maxiter200_mixer_beta010/result.json`](R1_diagnostic_case/run_maxiter200_mixer_beta010/result.json). Raw log: [`R1_diagnostic_case/run_maxiter200_mixer_beta010/gpaw.log`](R1_diagnostic_case/run_maxiter200_mixer_beta010/gpaw.log). The full energy/residual history is [`R1_DIAGNOSTIC_SCF_HISTORY.csv`](R1_DIAGNOSTIC_SCF_HISTORY.csv); machine-readable convergence state is [`R1_DIAGNOSTIC_CONVERGENCE_STATUS.json`](R1_DIAGNOSTIC_CONVERGENCE_STATUS.json).

Original R1/R3 SCF analysis: [`R1_R3_SCF_SUMMARY.json`](R1_R3_SCF_SUMMARY.json) and [`R1_R3_SCF_HISTORY.csv`](R1_R3_SCF_HISTORY.csv).

## Decision

Phase 1 is complete as a diagnostic/recovery report. R2/R4 energies and forces are recovered from existing logs; their stress is `STRESS_UNAVAILABLE`. The independent R1 copy converged with beta 0.10, but it does not replace the original screening result. No 8×8 validation was run or authorized by this report.
