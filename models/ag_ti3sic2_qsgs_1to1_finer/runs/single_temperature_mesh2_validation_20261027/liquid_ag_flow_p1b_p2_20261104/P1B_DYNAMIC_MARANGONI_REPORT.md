# P1b dynamic Marangoni benchmark — status report

## Decision

**Not run; no P1b pass/fail decision is made.** The available P1 results establish controlled planar Stokes flow and static curved-interface force closure. They do not reproduce a translating droplet, solve a temperature field, or exercise a moving interface under a temperature-dependent surface tension.

## Baseline and what it establishes

Baseline commit: `66961178b383fa7c507bea6318f8672e1cb4c08e` (`codex/ag70-p1b-p2`). The inherited P1 report records an automatic-curvature planar Stokes benchmark with velocity error 0.7805% at H/Δx=64, reversal when the imposed sign of dσ/dT changes, maximum mass drift 1.53×10⁻⁹%, and a separate static-circle force-closure residual of 0.0137%. That implementation uses a flat interface, an imposed dimensionless tangential stress, and a fixed interface. It is not the Basilisk official dynamic droplet-migration case.

The inherited P1 configuration records that the local Basilisk source mirror did not contain `src/test/marangoni.c`; it also lacks a Darcs revision identifier. The exact executable and selected headers are hashed in the inherited `results/p1_solver_configuration.json`, but those hashes cannot identify the official dynamic benchmark source or its complete dependency revision.

## Source retrieval and reproduction status

The official reference target is `https://basilisk.fr/src/test/marangoni.c`. The current managed environment uses restricted egress and does not currently allow `basilisk.fr`. A configuration draft adding exactly `basilisk.fr` and `trc.nist.gov` was saved for review; the environment tool returned `requires_publish=true`, so that draft has **not** changed the running environment. No source was downloaded and no result has been fabricated.

The P1b work therefore remains blocked at source acquisition. After the environment configuration is reviewed and published, the required sequence is:

1. Retrieve the complete official test source and every included local dependency from the same Basilisk revision.
2. Record the source URL, revision identifier, SHA256 hashes, compiler, and run command.
3. Run the official method and compare the droplet migration speed, interface trajectory, and mass conservation against its stated reference/theory at the published resolutions.
4. Implement the same physical case in the project’s automatic-curvature path only after confirming equivalent nondimensional parameters and boundary conditions; report method differences separately.
5. Test at no fewer than three resolutions and preserve raw logs and CSV data.

The official test uses a different interface-tracking and surface-force formulation from the project’s P1 planar implementation. Results from those methods must not be presented as an apples-to-apples code comparison unless the geometry, stress law, nondimensional groups, and boundary conditions are matched.

## Scope limits

No Ag70 geometry, Ti₃SiC₂ skeleton, dimensional Ag flow, temperature solver, or melt-pool calculation was used. P1b remains a controlled numerical verification task only. The known P0 exact-curvature/wavelet-AMR defect remains unchanged.
