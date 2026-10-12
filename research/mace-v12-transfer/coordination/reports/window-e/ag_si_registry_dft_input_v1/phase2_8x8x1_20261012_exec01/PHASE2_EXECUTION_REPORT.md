# Ag–Si Phase2 execution report

- Authorization SHA-256: `dd9eb403db1116e9f4524822741c00617f4239ab5caa1f9492844be836edf77f`
- Scope: R1 followed by R2 only; fixed-geometry PBE/PAW single points at 8×8×1; MPI 4; OMP 1.
- Role: `REGISTRY_SCREENING_ONLY`.
- Preflight bindings passed for authorization, parameter lock, runner, collector, and each input and structure SHA.
- Continue-on-error policy was applied. No retries or additional calculations were started.

| Registry | Status | Energy (eV) | Free energy (eV) | Forces | Stress | SCF steps |
|---|---|---:|---:|---|---|---:|
| R1 | SUCCESS | -188.92149225542852 | -189.16802546281826 | AVAILABLE | AVAILABLE | 152 |
| R2 | FAILED_SCF | unavailable | unavailable | NOT_AVAILABLE | NOT_AVAILABLE | 160 |

R2 reached the configured 160-step SCF limit and raised `KohnShamConvergenceError`. Its full GPAW log and traceback are preserved. R2 was not retried.

Detailed energies, force array, stress tensor, provenance SHAs, and collection statuses are in `REGISTRY_PHASE2_RESULTS.json`. Per-case manifests and SHA lists are under `R1/` and `R2/`; `EXECUTION_MANIFEST.json` records the authorized order and settings.
