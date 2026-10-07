# V17 training pool provenance audit

## Result

21 of 35 V16 rows are directly archive verified; 14 remain partially supported. All eight V17 signed training diagnostics were already archive verified and retain native and free energies separately. Therefore the proposed 43-row pool has 29 directly supported rows and 14 partially supported rows. It cannot be certified as a completely original-run-proven homogeneous target.

## Evidence and method

Current V16 rows were compared exactly against their inherited source rows. Newly reconciled archives (rows 14–22) and residual/V16 archives were compared for ordered species, coordinates, cell/PBC, exact REF energy and force copies, converged summaries and original-log extrapolated energy (six-decimal log tolerance 5.1e-7 eV). Three parent archives and eight signed diagnostics reuse the completed hash-pinned parent audit. No development-validation or withheld output was read. Source declarations are retained separately from original evidence in ledger.json.

Directly supported interface labels use native extrapolated GPAW PBE energies, PW500/Gamma/0.1 eV smearing where recorded. Rows 0–13 retain their original declaration evidence and precise missing-source paths/hashes. Rows 0–9 have serialized native aliases, which do not independently prove the original getter; four elemental references have unproven energy convention. Older gap declarations use 1x4x2 rather than Gamma: common XC alone does not establish identical numerical settings. No evidence establishes that a free-energy label was substituted into the directly verified native target.

## Integration recommendation

Do not describe all 43 labels as fully certified. Preserve the 14 uncertain provenance flags and method differences; seek original logs/inputs/outputs or authorize matched relabeling before claiming complete consistency. Do not silently replace native energies with free energies, and keep both fields for future diagnostics. Full composition rank is 4; archive-supported-only rank is 3, so removing uncertain elemental anchors is not a neutral filter. A must review target convention and reference-anchor policy before construction. Provenance verification does not prove force accuracy, model quality, or thermodynamic transferability.

## Reproduction

Run inventory.py, then audit_pool.py with the existing GPAW Python environment. Both scripts operate only on the listed training evidence. SHA256SUMS.txt covers artifacts, not itself.
