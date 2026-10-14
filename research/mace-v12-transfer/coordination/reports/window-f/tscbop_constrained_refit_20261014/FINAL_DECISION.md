# Stage F decision: constrained TSCBOP pilot

## Decision

**HOLD before parameter sensitivity and fitting.** The available archive confirms a roughly 1.99 eV/Å TSCBOP force mismatch on a fixed 48-atom Ti₃SiC₂ geometry, but the current task environment cannot freshly reproduce it because ASE and LAMMPS are absent. The independent Stage E2 root-cause report named in the task was not found in the E commit tree. E's force-environment report and A's checksummed `check.py` provide useful prior implementation evidence, but the raw 114-atom surface labels are missing from this checkout, the audited DFT records contain no stress, and the authorized set has no independent validation geometry. These gaps block a defensible sensitivity matrix and constrained pilot fit.

The check of 48-atom records found no concrete atom ordering, element mapping, geometry, or force-unit defect in the archived scorer logic. That finding is limited to the stored 48-atom outputs and code assertions; it is not a replacement for the missing E2 report or a fresh LAMMPS reproduction. The original parameter file remains untouched. No candidate parameters were produced.

## Required status flags

| Flag | Status | Basis |
|---|---|---|
| `DATA_AUDIT_PASS` | **FAIL** | 48-atom archive integrity passes, but fitting inventory is incomplete: raw surface labels are absent, stress/strain labels are absent, and independent geometry is unavailable. |
| `BASELINE_REPRODUCED` | **PASS** | Checksummed archived A run reports 1.996540 and 1.987192 eV/Å on the same geometry at k4/k6; a fresh run in this container was not tested. |
| `PARAMETER_SENSITIVITY_COMPLETE` | **NOT_TESTED** | ASE/LAMMPS absent; no finite-difference parameter runs. |
| `PILOT_FIT_COMPLETED` | **NOT_TESTED** | Fit gate held; no candidate file or optimization. |
| `INDEPENDENT_VALIDATION_PASS` | **NOT_TESTED** | No independent validation family was available or evaluated. |
| `PRODUCTION_READY` | **FAIL** | The archived baseline has large Ti/C force errors; independent validation and elastic/strain response evidence are missing. Production use is not supported. |

## Actions and boundaries

The audit only read the specified E reports, the original TSCBOP file, the authorized 48-atom records and their manifests, and prior A diagnostic reports. It did not open sealed test data, change the source potential, modify any other window archive, install software, launch DFT, or run MD/TTM/training. No Ag EAM or Ag–Ti/Si/C cross term was changed.

The four existing volume/shear proposals and recovery of the raw 114-atom archive are recommendations for a later authorized stage. Until missing audit evidence and a usable static evaluator are supplied, Stage F remains stopped before fitting.
