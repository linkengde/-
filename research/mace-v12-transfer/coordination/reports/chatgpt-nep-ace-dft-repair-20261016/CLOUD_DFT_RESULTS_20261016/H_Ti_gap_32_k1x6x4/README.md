# H 32-atom Γ-centered 1×6×4 single-point result

This directory contains the sole P1 task authorized by the cloud instruction. GPAW SCF converged in 49 steps; comparison against P0 fails all three predeclared k-point gates, so no follow-up calculation was started.

- Full result and 32×3 forces: `RESULT.json`
- Compact comparison record: `DFT_P1_RESULT.json`
- SCF log: `gpaw.log`
- Indexed coordinates and forces: `forces_indexed.csv`
- Method, thresholds, gate outcome: `SCF_AND_KPOINT_AUDIT.md`
- Environment and PAW evidence: `ENVIRONMENT_AND_PAW_SHA.json`, `preflight/`
- Checksums: `SHA256SUMS.txt`

Input SHA-256: `575686f5c9749c23c2258e7491633cd265367bad00432d0b9078feb89fe0f6b8`. P0 source result SHA-256: `0a6fe5669af67200beea1af442d078408e277c850b5e1464b8e9313b062999df`. Review commit: `3774b4c23abcb7351507349232c2960d04c0b44f`.
