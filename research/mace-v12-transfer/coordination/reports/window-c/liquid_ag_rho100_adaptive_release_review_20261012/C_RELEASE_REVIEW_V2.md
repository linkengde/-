# Ag64 rho100 adaptive continuation: C release review

## Current decision

`RELEASE_STATUS = BLOCKED_PENDING_D_PACKAGE_ARCHIVE_AND_REVIEW`

No continuation authorization is issued and no MD was started. The C strategy baseline is commit `53cc7bb8f33cf94568212b3b36a04c7f0a4a5709`. The D branch currently visible at GitHub is `cloud-d-b1-dft` HEAD `87ac1acf05a938f77180a0a3acfd4e76cf11c466`. Its tree contains the earlier rho100 pilot and step-3000 archive, but does not contain `followup_runner.py`, `liquid_criteria.py`, the claimed code-hash manifest, or the new synthetic integration/fault-test report. Therefore the supplied hashes cannot be checked against source bytes and the runner cannot be audited. This is a publication/evidence blocker, not a scientific failure.

The old step-1000 authorization covers only step 0–1000. It does not authorize this continuation. Once D publishes the exact package and all checks pass, C can issue one bounded authorization for the whole adaptive plan; internal 1000-step checkpoint/archive boundaries need not prompt the user, but each must pass remote restore/hash validation before the next chunk begins.

## Science QA V2 decision: RDF and sustained order

The prior strategy's phrases “broad first peak” and “no persistent sharp fcc/hcp multi-shell peaks” are qualitative, with no fixed estimator or uncertainty rule. A 64-atom cubic cell has `L = 10.7174553153 Å`; only `r < L/2 = 5.3587277 Å` has complete spherical shells under minimum-image PBC. RDF beyond that range cannot be used. A visual peak-width or higher-shell sharpness threshold is not calibrated against matched solid/liquid 64-atom references, so it must not decide PASS.

Version 2 retains reproducible first-shell diagnostics and locked location/coordination bounds, but removes peak width and higher-shell sharpness as binary gates. It instead requires both consecutive 2 ps windows to pass the explicit Q6/PTM order limits and requires positive diffusive MSD. Peak width and higher-shell features are still saved and reported. If the first peak/minimum is statistically unresolved, or the trajectory cannot support the required block uncertainty estimate, the result is HOLD, never default PASS.

`LIQUID_QA_SCIENCE_LOCK_V2.json` specifies the full definitions: finite-size-corrected ordered-pair RDF normalization, 0.05 Å histogram bins, 0.10 Å Gaussian smoothing for extremum detection only, peak/valley ranges, CN integration, moving-block bootstrap, the local q6 formula, PTM FCC/HCP classifier requirements, cluster connectivity, and unresolved-data behavior. This is a proposed science spec, not yet signed. D's `liquid_criteria.py` must be updated to match it and pass synthetic analytic/edge-case tests; then C must inspect the actual source and environment pins.

The old mean-Q6 threshold of 0.20 is only comparable if D reproduces the prior archived 11-frame diagnostic mean `0.0887 ± 0.01` with the specified 3.8 Å first-shell cutoff and exact q6 definition. Failure of that regression test holds the order gate for review; the threshold must not be tuned to force a pass.

## Required D package before review can continue

D must push an immutable package on its own branch containing the exact runner and QA source bytes, environment/config lock, SHA256SUMS, package README, and reports for mock-only end-to-end tests. Tests must show: valid automatic chunk progression; every boundary is ≤1000 steps; stage transitions occur only after the signed gate; no step beyond the global cap; injected numerical/time/interruption faults stop safely; archive upload and fresh-fetch restore preserve positions, velocities, RNG, step, scientific parameters, and task identity; and any remote archive/restore failure prevents the next chunk from starting. No real Ag MD belongs in this preauthorization package.

D must also include a fresh-remote-fetch receipt and the exact remote branch head. C will verify each file SHA from fetched Git objects, inspect the code paths that invoke the real ASE integrator and stop gate, and check the target step-1000 checkpoint (`c820437a…f34046`), EAM (`958f4623…c58447`), prior runner/environment, config, and archive backend bindings. Any code/spec change requires a new code/config hash and a repeat of affected tests.

## Execution boundary after review

Even after successful review, the only permitted run is the bounded rho100 task in the user's instruction: continue from the archived step-1000 state without reinitializing positions, velocities, or RNG; retain fixed EAM/cell/density/timestep/friction; stop at the locked global cap; checkpoint, archive, fresh-fetch, and verify every ≤1000-step segment; automatically progress only after a passing restore; stop on the 1250 K candidate-window terminal QA or any phase/global/resource/numerical/archive stop. No other trajectory, DFT, or MACE training is included.

`execution_authorized=false` until D's package is archived and C has independently completed the source, test, restore, and live-preflight review.
