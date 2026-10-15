# Minimum calculation matrix — provisional pending E handoff

## Not an execution release

`N_NEW_SCIENTIFIC_STRUCTURES=NOT_FROZEN`; `N_NEW_SCF=NOT_FROZEN`; `START_DFT=NO`.

E reports at least 10 new SCF points under specific assumptions in `C_MINIMUM_DFT_TASKS.md`. C does not yet have that file or the six source-audit files. Consequently C cannot verify what the ten points contain, their atom counts, whether component references are included, or their source-family assignments. The table below is C's conditional selection/order, not a reconstruction of the unavailable E request or a competing ten-point release.

## Task selection order

| Priority / ID group | Candidate calculation | Minimum purpose / count basis | Reuse and release condition |
|---|---|---|---|
| P0 / REUSE-OLD6 | Existing six public local perturbation labels | 6 records; 3 known parents; **0 new SCFs** | Retain existing roles, do not treat as normal-distance scans or independent TEST |
| P0 / REUSE-135 | Existing relaxed full interface + same-geometry Ag-only + TSC-only, plus initial 2.0 Å full point | Exact count/identity to be read from handoff; do not repeat valid results | All files, convergence, settings and matched component identity must pass |
| P1 / GAP-REF | One reviewed parent with compressed/reference/expanded normal separation | 3 full-interface geometries; if the reference is already valid, **2 new full geometries**; otherwise 3 | Rigidly translate complete Ag component only. Freeze actual distances from source; do not call a reference gap equilibrium without evidence |
| P1 / PHASE-REF | Exact same-cell projected phase-only references for GAP-REF | One identical TSC geometry can be reused across the gap scan; translated Ag references need individual results unless numerical translation invariance is validated | Do not subtract a different strain, relaxed phase, cell, cutoff, k mesh, or energy convention |
| P2 / REGISTRY | Second/third non-equivalent lateral registrations | Use E's proposed two registries × three gaps only after reading the actual definitions; remove overlap with P1 | Crystallographically specified translations; no random registry. New registry within the same parent is coverage, not independent validation |
| P2 / TERMINATION | Ti/Si/C-exposing surfaces or documented local contact environments | Add only the contact channels not genuinely sampled by accepted old labels/P1/P2 | Audit termination, stoichiometry, slab thickness, bottom face, vacuum, pair-distance distributions and Ag–Ti/Si/C coverage. A Si surface is not proof that all three pair channels are identifiable |
| P2 / REPULSION | Deterministic compressed geometry per missing contact channel | Reuse an accepted compressed point if it covers the needed range; otherwise freeze a new point | No arbitrary extreme overlap; geometric screen and force/numerical error budget before release |
| P3 / VALIDATION | Unfitted parents and registry/termination holdouts | At least one genuinely held-out parent before model claims; broader three-channel validation if required | No same-parent seed/rigid shift/relabeling as independent TEST. No sealed structure access or new role assignment here |

## Conditional first-batch count, not authorization

For one three-gap fixed-internal-geometry family, a strict coordinate-matched plan comprises three full-interface points, three Ag-only points at the actual coordinates, and one identical TSC-only point: **7 distinct SCF inputs**. If a valid reference-gap full/Ag/TSC triplet is already present, it removes three of these: **4 new SCF inputs (2 full interfaces + 2 Ag-only)**. Further reduction requires validated numerical translation invariance, not an assumption. If component internal coordinates also change, these counts must be recomputed.

This is only the smallest gap-response closure for one parent, **not** a complete three-pair Morse fitting dataset. It does not replace E's reported at-least-10-point scheme. The exact global minimum cannot be frozen before that scheme, existing-data overlap, terminal environments, and validation-source isolation are checked. No new coordinates or calculation tasks were created.

## Energy and force targets

Keep full-system native energy, free energy/smearing contribution, complete forces, and (where supported) stress. Derive interaction energies/forces only from validated matched components at the same geometry and method. Store full raw labels separately from derived cross residuals. Use consistent per-cell or per-interface-area normalization and specify how many interfaces exist; do not confuse an asymmetric vacuum slab with a two-interface periodic stack. Slab stress is vacuum-diluted unless a separate convention is defined.

Do not compare different atom counts, stoichiometries, strain partitions, slab thicknesses, or codes by directly subtracting their total energies. A 50–80-atom screening model would be a **new parent family**, not a strict energy-equivalent replacement for the 135-atom interface. Its terminal chemistry, Ag thickness and strain must be documented. No smaller model is selected or generated here.

## Numerical convergence checks before scientific labels

1. Read E's actual QE/GPAW inputs before pinning any new setting. Same PBE name does not make QE UPF and GPAW PAW absolute energy zeros or convergence tolerances equivalent.
2. Use an unchanged representative geometry. Compare the inherited baseline mesh against one denser in-plane mesh with the same offset convention. Test cutoff (and QE density cutoff) separately; do not silently change geometry or smearing.
3. Tighten SCF accuracy on that same input to establish energy/force noise. QE `conv_thr` in Ry is not interchangeable with GPAW's energy/density/eigenstate criteria of 1e-5.
4. Proposed review targets, **not frozen acceptance thresholds**: changes below 1 meV/atom and force-vector RMSE below 0.02 eV/Å with max force difference below 0.05 eV/Å. Interaction-energy error must also be small relative to the cross-potential signal and reported per interface area. Final tolerances must follow the fitting error budget and phase-subtraction uncertainty.
5. Distinguish new physical geometries from extra numerical-control SCFs on the same geometry. Neither the 7-input conditional example nor E's unverified ten-point count includes unspecified new convergence tests.

All old data roles remain unchanged. All new entries are proposals; no TRAIN/DEV/TEST labels are assigned by this matrix, and no fitting is started.
