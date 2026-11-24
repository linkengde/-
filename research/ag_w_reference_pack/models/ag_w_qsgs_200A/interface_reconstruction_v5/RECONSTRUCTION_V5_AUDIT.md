# Ag-W QSGS Interface Reconstruction v5 audit

## Result

**Overall: FAIL.** No v5 structure passed the Dai EFS static gate, and no 300 K MD or TTM run was started. The TTM solver and its parameter files were not modified.

The best geometry-only candidate met the mass-ratio, estimated phase-density, and bicontinuity checks. It still had seven W-W pairs below the frozen distance screen. After a 10,000-iteration spring-free Dai EFS minimization with the intended fixed bottom layer, the structure developed substantially shorter contacts and atoms moved into the vacuum region. The minimizer stopped at its iteration cap.

## Fixed target and audit gates

- QSGS material: 200 x 200 x 200 A; box: 200 x 200 x 300 A.
- Boundary: `p p f`; W-BCC and Ag-FCC; target Ag:W mass ratio 1:1.
- Keep the original QSGS model and all v1-v4 records unchanged.
- Mass-ratio error <= 0.1%; each phase number-density error <= 1%.
- Largest connected component >= 98% in each phase; XY winding rank 2 and Z spanning in both phases.
- Frozen nearest-distance screens: Ag-W 2.54 A, W-W 2.60 A, Ag-Ag 2.48 A.
- No atoms outside z=0..200 A; spring-free Dai EFS static minimization must converge.
- For moving atoms: maximum force <= 0.1 eV/A and RMS force <= 0.02 eV/A.

The distance screens were carried forward from earlier local interface scans. They are used as fixed geometry checks in this audit, not described as universal experimental bounds.

## What was rebuilt

The original 80^3, 2.5 A QSGS morphology was regularized with a Gaussian filter and thresholded back to the original W voxel fraction of 0.3545. This reassigned 43,054 voxels (8.41%) while preserving the W voxel volume fraction. Separate W-BCC and Ag-FCC populations were regenerated. The new population contains 482,601 atoms: 178,555 W and 304,046 Ag, a +130 atom change from the source. Its Ag:W mass ratio is 0.999127, within the 0.1% limit; estimated phase number-density errors are -0.0757% for W and +0.0910% for Ag.

The signed-distance shift alone failed. A later temporary same-phase Morse plus cross-phase soft-core preconditioner produced the best candidate. It retained large connected components (W 99.327%, Ag 99.997%), XY winding rank 2, and Z spanning in both phases. The best candidate's minimum distances were W-W 2.5801 A, Ag-Ag 2.5281 A, and Ag-W 2.6134 A. Seven W-W pairs remained below the frozen 2.60 A screen, so this was not a geometry pass.

Temporary constraints were removed for the physical static check. The bottom z=0..20 A layer was fixed as intended for the later MD model. The spring-free EFS minimizer ran 10,000 iterations and stopped at the iteration cap:

| Static check | Result | Gate |
| --- | ---: | --- |
| W-W minimum distance | 2.4403 A | Fail: 10,302 pairs below 2.60 A |
| Ag-Ag minimum distance | 2.2945 A | Fail: 681 pairs below 2.48 A |
| Ag-W minimum distance | 2.4343 A | Fail: 797 pairs below 2.54 A |
| Atoms above z=200 A | 1,871 | Fail |
| Moving-atom maximum force | 0.3809 eV/A | Fail |
| Moving-atom RMS force | 0.00314 eV/A | Pass |
| Minimizer | 10,000-step cap | Fail: not converged |
| W/Ag mass ratio, density, connectivity | Within target | Pass |

The EFS-relaxed structure therefore fails both the strict static gate and the exploratory readiness gate. The severe short contacts, escaped surface atoms, and residual local force make a no-source 300 K pilot unsuitable. No laser heating, cooling, or TTM source remapping was performed.

## Trial ledger and files

- `precondition_trial_summary.json` records all geometry and LAMMPS precondition variants, including failed attempts.
- `exploratory_readiness.json` records the separate, less demanding check for a small 300 K no-source pilot and why it also failed.
- `reconstruct_interface_v5.py` regenerates the QSGS-based geometry candidates.
- `prepare_lattice_springs_v5.py` writes the temporary same-phase/cross-phase geometry precondition inputs. The temporary potential is never used for final static acceptance.
- `audit_static_v5.py` independently computes the frozen geometry, force, topology, density, atom-count, and minimizer gates.
- `trials/` contains the signed-normal geometry candidates and reports.
- `samephase_pair_precondition_screenbuffer_fast/` contains the best geometry-only candidate and its input, manifest, and LAMMPS log.
- `efs_static_v1_fixed_bottom/` contains the spring-free EFS input, final failed structure, dump, log, and `static_audit.json`.
- `efs_static_v1/` records the separate no-fixed-bottom run that lost one atom before output.
- `SHA256SUMS.txt` lists the key source, candidate, final failed output, and audit hashes.

The source `qsgs_morphology.npz`, original v1-v4 structures, and Ag-W TTM files remain unchanged. `static_audit.json` marks `md_start_gate: FAIL`; the failed EFS output is retained for inspection and is not a ready-to-run MD input.
