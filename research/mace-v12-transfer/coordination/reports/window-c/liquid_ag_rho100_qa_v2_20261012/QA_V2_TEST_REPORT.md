# Ag64 rho100 QA V2 engineering test report

Date: 2026-10-12 UTC
Owner: D
Task family: `C-D-AGL-1250K-rho100-followup01`

## Result

The V2 analysis and fail-closed control components pass the synthetic engineering suite. The archived 11-frame Q6 regression does **not** pass the comparator required by C's V2 lock, so the actual scientific release status is **HOLD**. No threshold, formula, potential, checkpoint, structure, velocity, or RNG state was changed. No MD, DFT, or training was run.

## Archived-data regression

The read-only analysis used the archived `trajectory.traj` (11 frames; SHA-256 `80961895ae1b552aec908c35511de421fab74d26026063cd9d0c3693f93dacd5`) and the unchanged step-1000 checkpoint (SHA-256 `c820437a2512c4b3d7dc9e60bb60ac3ced6588aa1f2a629d58cfdcabb24f34046`). The C V2 local definition gives mean frame Q6 **0.3516369060394911**. The archived global Steinhardt Q6 implementation gives **0.08871920638262233**, reproducing the documented legacy value but not the V2 local value. The V2 requirement is `0.0887 ± 0.01`; its local-Q6 comparator therefore reports `Q6_THRESHOLD_COMPARABILITY_HOLD`. The old `0.20` threshold is not reused against the V2 local metric.

## Test coverage

`test_qa_v2.py` ran 8 tests successfully:

- Recomputed both Q6 definitions over the exact 11 archived frames and verified the V2 regression returns HOLD without tuning.
- Ran actual OVITO PTM 3.16.1 on synthetic periodic FCC and HCP crystals at RMSD cutoff 0.10; both labels and periodic cell handling matched. A nonperiodic input returned a periodicity failure.
- Checked the locked `rmax=5.35 Å` against `L/2=5.3587276576 Å` and the finite-cell RDF normalization, first peak/valley, coordination number, 2,000 moving-block replicates, and confidence-bound path on synthetic counts.
- Ran the complete thermal/RDF/Q6/PTM/MSD assessor on deterministic synthetic frames; the synthetic comparator path passes, while a failed historical comparator returns HOLD.
- Verified incomplete sampling and unresolved RDF return pending/HOLD rather than PASS; verified missing OVITO/PTM returns `PTM_UNAVAILABLE`.
- Exercised 1000-step chunk boundaries and the 3599 s, 3600 s, and over-limit start-step checks. An already-finished segment cannot start another step.
- Tested archive commit, push, fresh clone/fetch, file restore, checkpoint payload comparison, SHA verification, and potential-file rejection against an isolated local bare Git remote. An archive restore failure prevents phase advancement.

The test emitted only upstream ASE/NumPy deprecation warnings while reading the archived trajectory; all assertions passed. The local bare-remote test is isolated synthetic infrastructure and does not substitute for the real D-branch GitHub publication/fresh-fetch check recorded by the package commit.

## Release gate and remaining blocker

C's file `LIQUID_QA_SCIENCE_LOCK_V2.json` remains the original unsigned proposal (`execution_authorized=false`). The V2 lock explicitly requires the historical Q6 regression to match and says the threshold is not comparable otherwise. D has not revised that lock or invented an alternative threshold. Consequently no execution runner or authorization-ready configuration is released, and automatic liquid continuation cannot safely transition into a signed run. The legacy runner remains outside this package because it consumes global-Q6/CNA and the old qualitative RDF policy; presenting it as a V2 runner would be misleading.

A C decision is required on the exact Q6 regression/threshold interpretation. Once C publishes a versioned resolution, D can bind the final V2 assessor and execution adapter to that exact lock and rerun the release-bound checks. Until then, `QA_V2_CONFIG.json` deliberately keeps `execution_authorized=false`.
