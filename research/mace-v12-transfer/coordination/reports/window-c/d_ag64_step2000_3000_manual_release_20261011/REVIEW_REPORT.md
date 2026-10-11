# Ag64 manual segment review

REVIEW_STATUS=APPROVED

This approves the scientific reasonableness of a manual review segment only. No signed configuration is issued; no execution_authorized=true configuration is published; no MD was started.

Task: C-D-AGL-1250K-rho097-seed01. Source step2000 → manual inspection at step3000, maximum1000 additional steps, 0.5ps.

## Source checkpoint

D source commit bc3ff0631db3b7d78dcc17556befa98c2f454480 was independently confirmed as its remote branch head. Source archive SHA c3a88a759e0bd258e74632aba0a68e4f1a765fcbe8324034d5bf127b5aa69965; checkpoint state SHA5d03864f07a98f177c3c0cbdcb850158ed89a64fa2a6ab7f879e8c7ab003be2a; arrays SHA0f7830ebeba5f9a83da77b15e153546b24fb3efaf0299de767d9b4fbe4639f78. Archive identity and actual ASE read-only restoration passed. Original coordinates, momenta, cell, PBC and PCG64 state remain the continuation source.

## Parameters

64 Ag atoms; rho9.03273341876g/cm³; cube10.826824658074Å; seed125009701; time step0.5fs; Langevin friction0.01fs⁻¹. Original schedule remains1600K hot stage until step20000, then cooling toward1250K. Steps2000–3000 remain1600K; this is not fixed1250K sampling. Original potential SHA958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447. No reinitialization is justified. A3600second walltime budget is a resource limit, not a scientific parameter change.

## Manual inspection point

Step3000 is a reasonable short diagnostic boundary for temperature/energy trend, RDF, Q6, coordination, minimum distance and MSD. It is not guaranteed to establish equilibrium or prove liquidity. Stop there; any further segment requires separate manual review.

## Risks and limits

Only1ps is currently completed; this next segment adds0.5ps. Temperature autocorrelation and64atom finite-size fluctuations limit conclusions. Short-time MSD is not a reliable diffusion coefficient; RDF/Q6 alone cannot prove melting. A wallclock interruption must preserve the last fully completed step and RNG before any restart. Archive failure must preserve the unique local output. Fresh D resource/runtime/remote-channel checks remain necessary before actual execution.

Supplementary non-scientific tests passed: source restore, synthetic interrupted/continuous state and RNG equality, timeout checkpoint, step3001 rejection, real backend upload/fresh-fetch/restore. These do not execute the real Ag potential or certify real MD thermalization. Code included here is a reviewed candidate only and still requires signed gates; it does not constitute a new manual-mode execution authorization.

All historical archives are untouched. Other five tasks, DFT and training remain HOLD.
