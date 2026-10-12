# Ag64 rho100 Q6 adjudication (science lock v2.1.0)

**Decision:** the two reported values are different observables, not a failed regression. The archived coherent/global Steinhardt Q6 is `0.08871920638262233`; the V2 mean-local Q6 is `0.3516369060394911`. The former averages each atom's `q_6m` over atoms before taking the norm; the latter takes each atom's norm first and then averages. Equality is neither mathematically expected nor required.

The legacy global-Q6 implementation regression remains active on the exact archived 11-frame trajectory (`80961895ae1b552aec908c35511de421fab74d26026063cd9d0c3693f93dacd5`): expected mean `0.08871920638262233`, tolerance `0.01`. The legacy `global_Q6 < 0.20` condition remains a necessary joint liquid screen, requiring at least 80% of frames to pass in each of the final two consecutive non-overlapping 2 ps windows. It is not sufficient by itself to establish liquid behavior.

V2 local Q6 must be computed and reported independently. Its 11-frame reference mean is `0.3516369060394911` for implementation diagnostics only. It has no pass/fail threshold and must not be compared to `0.20` or required to equal global Q6. Do not fit a local-Q6 cutoff from these 11 frames.

A liquid decision is conjunctive: global Q6 screen + PTM FCC/HCP crystalline fraction and cluster limits + finite-cell RDF/coordination point and uncertainty checks + positive-diffusion MSD test + two-window thermal stationarity. Missing/unresolved/underpowered metrics return PENDING/HOLD. Local Q6 is report-only pending separate calibration. Existing non-Q6 V2 definitions are unchanged.

## Execution status

This science ruling does **not** authorize MD. The D commit `9f7b4fe5e740acca9f5764374f8c51c9d7bceefd` contains QA V2 analysis primitives and synthetic tests, but its own QA report says the copied production runner is not V2-compatible and is excluded. D must bind this v2.1 lock into the actual bounded automatic runner, add the distinct-Q6 regression and conjunctive-gate tests, exercise the actual runner's segment/archive path against real GitHub with a fresh fetch/restore, and submit the final code/config/environment/checkpoint/potential identities. C will perform one final review; execution remains HOLD until C signs a separate bounded authorization.
