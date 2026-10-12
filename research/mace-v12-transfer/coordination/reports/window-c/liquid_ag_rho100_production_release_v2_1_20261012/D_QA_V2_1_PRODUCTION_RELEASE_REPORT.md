# D Ag64 rho100 production runner release for C review

D integrated QA v2.1 into the actual bounded continuation runner and bound the package to C's exact adjudication commit `457a341ad7712704a607f2c5038d1aa41de6bae9` and lock SHA `848f8db72d83cea91dcf4cd42ee9fbbf49a80e7b5e3759cb4f535be8cc5873d2`. The existing source step-1000 checkpoint, WMH2006 EAM potential, frozen parameters, and RNG are unchanged.

The production gate now applies historical regression and strict global Q6 `< 0.20` to coherent/global Q6 only. At least 80% of frames in each final consecutive nonoverlapping 2 ps window must satisfy it. Mean local Q6 is logged as a report-only diagnostic and has no threshold or trend gate. The runner retains OVITO PTM FCC/HCP classification at RMSD 0.10, periodicity and crystalline cluster checks, QA v2 finite-cell RDF/coordination and moving-block confidence requirements, MSD and thermal/energy stationarity checks. Insufficient evidence remains PENDING; missing, nonfinite, unresolved or unavailable required evidence stops in HOLD.

Stage control is bound to the existing plan: step 1000 start; 1600 K heating; exact 10,000-step, 5 ps linear cooling updates; 1250 K equilibration capped at 10 ps; 5 ps candidate window; global step cap 60000; 1000-step and 3600-second per-chunk caps. At each boundary the runner writes checkpoint/log state, pushes an immutable archive, fresh-fetches and compares hashes, restores positions/momenta/RNG/phase state from the fetched checkpoint, and permits continuation only after every archive/restore flag passes. No initialization, reseeding, or potential change is performed.

## Review evidence

- QA science lock SHA: `848f8db72d83cea91dcf4cd42ee9fbbf49a80e7b5e3759cb4f535be8cc5873d2`
- Production config SHA: `9258af2d01ef64025b3674ab184da1f1aca216dc9787784a599b424ba244e663`
- Environment SHA: `e71e04b8edc989001a702ffb308b1b897d6b521d0ebfdcb41b43592ae2bff1fb`
- Historical Q6 regression: PASS, 11 frames, local statistic report-only
- Synthetic engineering suite: 12 passed, 0 failed
- Real GitHub synthetic archive upload/fresh-fetch/restore/SHA: PASS, commit `c4f954e0b00bce9a75e7cc120fe03625d88a2ba1`
- Live read-only preflight: `PREPARED_HOLD`; resource, environment, EAM and source checkpoint checks pass
- Formal MD steps executed in this release task: 0

The config is deliberately `execution_authorized=false`. This package is submitted to C for review. D will not start the production continuation until C publishes a separate authorization blob bound to the final D branch head and all required lock, runner, criteria, archive, test, environment, checkpoint, potential, and step-limit hashes.
