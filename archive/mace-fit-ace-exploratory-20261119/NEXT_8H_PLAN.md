# Next compute window: approximately 7–8 hours

The uploaded archive is the frozen baseline. The following work continues from it in the current cloud environment. It does not modify `main` or relabel data before the calculations converge.

| Order | Task | Budget | Completion gate |
|---|---|---:|---|
| 1 | Finish the V17 `Ti_gap_2p3` fixed-geometry GPAW PBE/PW 500 eV calculation at 1×4×2 k-points, Fermi width 0.1 eV. Current job is active. | Up to 2.5 h | SCF meets the archived energy, density, and eigenstate thresholds; otherwise preserve the partial log and do not use its label. |
| 2 | If task 1 converges, run the same structure at a denser 1×6×3 mesh and compare total energy and per-atom forces. | Up to 2.5 h | Establish mesh sensitivity; use the denser label only if the calculation converges. Keep both labels for audit. |
| 3 | Refit ACE on a copy of the candidate set with the verified label, comparing full force weight against the existing 0.05 sensitivity scale. | 0.5 h | Report train/development RMSE, each development frame, and Si-site proxy without changing the frozen baseline. |
| 4 | Refit MTP10 with the same candidate and compare force weights 0.05, 0.10, and 0.20. | 0.75 h | Preserve all variants and report whether the force-component gate is met. |
| 5 | Run source-group holdout validation for the best ACE and MTP candidates; keep all frames from a source family together. | 0.75 h | Report every fold and aggregate error; do not describe cross-validation as an independent blind test. |
| 6 | Check model loading on small representative LAMMPS cells, build a target-environment label shortlist from the force audit, hash the outputs, and commit the follow-up to the same archive branch. | 0.5 h | Compatibility checks and checksums pass; no million-atom integration unless smaller-cell validation supports it. |

The time is a compute budget, not a reason to weaken convergence or manufacture labels. If the 1×4×2 DFT calculation does not converge, skip the denser-grid and refit steps that depend on it, and use the remaining time for grouped validation and a documented shortlist of additional labels. The current data still lack an independent blind test set and registered Si hotspot-center IDs.
