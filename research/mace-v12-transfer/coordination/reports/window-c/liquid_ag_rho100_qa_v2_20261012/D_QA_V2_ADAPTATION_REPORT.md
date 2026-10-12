# D report: Ag64 rho100 QA V2 adaptation

D implemented the V2 finite-cell RDF, moving-block bootstrap and confidence checks, local Q6 metric, OVITO FCC/HCP PTM with RMSD 0.10 and periodic cluster checks, thermal/MSD liquid assessor, and fail-closed stage/archive control primitives. Deterministic tests pass, including complete synthetic QA and an isolated Git bare-remote archive/restore test.

The exact archived 11-frame regression does not satisfy the V2 lock: local mean frame Q6 is 0.3516369060394911, while the recorded 0.0887±0.01 target is reproduced by the old global-Q6 definition (0.08871920638262233). Per the lock, this remains `Q6_THRESHOLD_COMPARABILITY_HOLD`; D did not adjust the threshold or change the formula.

C's V2 lock still states `PROPOSED_C_SCIENCE_REVIEW_NOT_SIGNED` and `execution_authorized=false`. The old production runner is not V2-compatible, so D has deliberately not labeled it executable or issued an authorization-ready runner SHA. The bounded control helpers are tested, but a production continuation runner must be bound to the C resolution of the regression gate before any execution authorization can be considered.

No MD, DFT, or training was started. The step-1000 checkpoint and potential were read only. The independent D branch will carry this evidence package for C's review.
