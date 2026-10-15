# Existing-data reuse review — 2026-10-15

## Decision and evidence boundary

`REVIEW_STATUS=PARTIAL_BLOCKED_BY_UNAVAILABLE_HANDOFF`.

E's supplied handoff is local to another cloud task. In C, neither `AG_TISC_C_DFT_HANDOFF_20261015/` nor its tar.gz exists, and Git object `a0e94666eb6a7a355b9c1e5bcaba58c9868cd8c2` is absent. A fresh fetch of E's published branch resolves to `bdba4115569f74578e2fc2e2bee031cebf3f502d`; that tree does not contain the requested handoff or `DFT_REQUEST_TO_C.md`. This is not a checksum failure: C does not yet have the bytes to check.

E reports a 634,192-byte tar.gz with expected SHA-256 `4e35ce85e6a47295e07c6e4e688925bc343fee7cdd966d39e3619b9d00d7026f`, 208 checked handoff files, and 132 checked V19 files/16 reference cases. Those are sender-reported facts, not independently reproduced C checks. The 2,612,538-byte expanded size is also sender-reported. C needs the tar.gz attached here, or a reachable Git commit/path containing it. No upload of restricted third-party material is requested.

## Independently verified accessible historical data

C audited the six explicitly public residual-acquisition records in `research/mace-v12-transfer/exports/cross_potential_dft_reference_bundle_20261009/`. **Six records pass**, representing **three parent families**, not six independent validation families. The audit checks 38 bundle checksum entries, the original verifier's artifact hashes, input/result symbol order and atom IDs, positions, cell/PBC, all forces, energies, archived convergence evidence, and PAW fingerprints. It does not invoke a calculator. Full results and exact input/output hashes are in `AUDIT_EVIDENCE.json`; the reproducible read-only script is `audit_existing_data.py`.

| Label | Atoms | Native energy (eV) | Free energy (eV) | SCF iterations | Historical wall time (h) |
|---|---:|---:|---:|---:|---:|
| AgC_residual_mode_m_06_v19_acq | 48 | −255.335735269 | −256.492009845 | 58 | 0.791 |
| AgC_residual_mode_p_06_v19_acq | 48 | −254.228691059 | −255.383426980 | 58 | 1.733 |
| AgSi_residual_mode_m_03_v19_acq | 51 | −207.728837664 | −208.659987043 | 48 | 0.707 |
| AgSi_residual_mode_p_03_v19_acq | 51 | −207.701251175 | −208.632357852 | 47 | 1.560 |
| AgTi_residual_mode_m_03_v19_acq | 49 | −174.409515271 | −175.064442234 | 40 | 0.776 |
| AgTi_residual_mode_p_03_v19_acq | 49 | −174.204986873 | −174.857873270 | 40 | 1.834 |

All six use GPAW 26.7.0, PBE/PAW, 500 eV, integer 1×1×1 mesh, Fermi–Dirac 0.10 eV, energy/density/eigenstates tolerances 1e-5, Pulay beta 0.05/nmaxold 8/weight 100, and four MPI ranks. The logged Ag/Ti/Si/C PAW XML fingerprints match the currently archived setup XMLs; exact setup-file SHA-256 values are in the audit evidence. Full finite force arrays are present as `PW_PBE_forces`; native extrapolated energy and free energy are recorded separately. Stress labels are not established for these six records by this audit.

These are local residual-direction perturbations of small extracted geometries in large periodic boxes (28.14 Å for AgC, 27.86 Å for AgSi, 28 Å for AgTi), not complete periodic surface-termination models or uniform interface-distance scans. SCF convergence is verified, but numerical cutoff/k-point/box-size convergence for a new fitting accuracy target is not established by the six-record audit alone. They are reusable as archived force/response diagnostics without repeating their existing SCF. Inclusion in a new cross-potential objective remains conditional on scope, numerical error, and role review. Their original role is `training_acquisition_reused_development_family`; no role was changed. They cannot become independent TEST data.

These six may overlap E's reported six old points or V19's 16 cases; that identity cannot be asserted until the handoff is available. Do not add 6 + 16 as 22 distinct structures.

## 135-atom matched components

Five E analysis documents from the published branch pass their own SHA manifest. They cite a relaxed 135-atom geometry with SHA `eda2586c63c6c3c793ffce99417e7c04569941b0a2a916fda1b5402b7ba52d9f` and an initial geometry with SHA `1a3e341e1132bba9b349ab0fe6e3a3be11d1e1661671b96f6810a54f2e1f54ae`. E's new message additionally reports the full-interface, Ag-only, TSC-only matched labels and an original 2.0 Å full-system single point. Those raw files are not accessible to C here. Therefore **0 matched 135-atom triplets are independently verified in this review**, rather than zero triplets existing.

On receipt, compare full-system coordinates with the filtered phase coordinates by stable atom ID, without relaxation, reordering, recentering, changing cell/PBC, or changing strain. Verify software/pseudopotential identities, valence/spin/charge, cutoff, density cutoff where applicable, exact k mesh/offset, occupations, smearing, SCF criteria and actual convergence. Save native energies and units; do not combine GPAW and QE absolute energies.

Only a matched set supports `E_cross = E_full − E_Ag − E_TSC` and atom-mapped `F_cross,DFT = F_full,DFT − F_phase,DFT`. A relaxed isolated phase is not a matched reference for a frozen interface component. Translation-equivalent reference reuse requires a demonstrated numerical translation-error bound; different internal geometries cannot share a phase reference. No energy subtraction was performed in this review.

## Reuse counts

- Independently checked accessible labels: **6 energy/force records, 3 families**, reusable for their existing diagnostic scope.
- Independently checked 135-atom matched triplets: **0 verified here; source pending**.
- E's claimed handoff inventory and at-least-10-SCF plan: **pending receipt**, not rejected and not accepted without bytes.
- Confirmed avoidance: do not repeat the six existing converged diagnostic SCFs. Additional avoided full/phase calculations require exact structural/method matching; the count is not frozen.

The fixed classical components are not altered. E's verified analysis reports about 2 eV/Å bulk TSCBOP force mismatch. Matched DFT cross labels can isolate interfacial contributions, but fitting a cross term to compensate for that backbone error would not establish a reliable total potential. Preserve this limitation in subsequent validation.
