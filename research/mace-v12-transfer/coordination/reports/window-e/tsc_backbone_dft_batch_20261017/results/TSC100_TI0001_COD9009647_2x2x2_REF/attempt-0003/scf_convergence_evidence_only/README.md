# Converged SCF evidence only — incomplete result

GPAW completed and converged the fixed-geometry SCF from the preserved step-20 checkpoint. The log records 20 additional iterations, final extrapolated energy -798.275237 eV, and `Converged in 20 steps`.

The worker then raised an assertion while checking the shared log for the convergence marker. It exited before calculating or saving forces and stress. This directory therefore contains convergence evidence, not a complete DFT label. Do not import it as force/stress training data or assign a data role.

The original checkpoint remains unchanged (SHA-256 in `SCF_EVIDENCE.json`). No further calculation was started after the original cumulative time budget.
