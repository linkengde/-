# AG_TISC_C_DFT_HANDOFF_20261015

Purpose: transfer E's already completed, read-only Ag–Ti/Si/C cross-potential audit and its cited reusable DFT evidence to C. This package preserves the local source commit `a0e94666eb6a7a355b9c1e5bcaba58c9868cd8c2`; the source audit is not rewritten. No DFT, MD, or potential fitting was run while assembling this package.

Current decision from the audit: `CROSS_POTENTIAL_DATA_READY=NO`, `MORSE_FIT_FEASIBLE=NO`, `ADDITIONAL_DFT_REQUIRED=YES`. See `audit/DFT_REQUEST_TO_C.md` and the actionable execution queue in `C_MINIMUM_DFT_TASKS.md`.

Contents:
- `audit/`: exact six audited files, including original SHA256SUMS.
- `existing_dft/135atom_full_relaxation/`: converged full-system geometry/energy/free-energy/forces/stress and source logs.
- `existing_dft/135atom_components/`: converged Ag-only and Ti3SiC2-only results at that same final geometry, plus ID mapping, input subsets, settings, forces and logs.
- `existing_dft/135atom_original_2A_fixed_singlepoint/`: an already-converged full-system point on the original 2.000 A Ag–Si contact geometry. It has no matched phase-only references and is not presently a valid cross-interaction difference-of-differences point.
- `existing_dft/135atom_2x2x1_failed_attempt/`: failed attempt evidence only; not a DFT label.
- `existing_dft/v19_reference_bundle_20261009/`: the complete 132-file signed bundle: six local reference cases plus ten numerical controls; its embedded verifier previously passed.
- `existing_diagnostics/`: completed candidate02 hybrid-potential and Ag–C central-force lower-bound diagnostics cited by the audit.
- `SOURCE_INVENTORY.csv` and `SOURCE_VALIDATION.md`: source locations, byte counts, checksums, and validation results.
- `SHA256SUMS.txt`: checksums for every regular file in this handoff except the checksum file itself. Run `sha256sum -c SHA256SUMS.txt` from this directory.

Newly proposed registry/gap geometries are not present as frozen inputs. C must decide the model-form/acceptance gate, select and freeze exact structures and roles, then issue any DFT authorization. The supplied structures are existing references, not permission to launch calculations.

Large `.gpw` restart files (Ag component 305,563,401 bytes; TSC component partial checkpoint 711,450,360 bytes) were intentionally excluded to keep the handoff compact. Their source archive hashes were checked in place; final converged force/energy labels and full logs are included. Originals remain untouched in the repository.
