# V19 train36 energy and method lineage audit

## Result

This was a read-only review of the approved train36/dev3 records and their existing training manifests. It did not read test or sealed labels, alter dataset roles, launch DFT, or retrain a model.

Train36 has `REF_energy` for all 36 rows, but only 24 rows trace through current manifests to a verified archived DFT label. **Twelve inherited V16 rows remain unresolved at the original-source level:** the V16 manifest marks their method/energy lineage as inherited and not independently reverified, sets `independent_original_getter_verified=false`, and says original source settings are absent. Their serialized `source_method` string says Gamma PW-PBE 500 eV, but that declaration does not establish the actual run settings or whether `REF_energy` came from the intended GPAW energy getter.

| Train source group | Rows | Energy/method evidence |
|---|---:|---|
| Inherited V16 source rows | 12 | Original settings and energy getter not independently verified |
| V15 residual archive labels | 3 | Archived label evidence; source convention recorded as native extrapolated energy |
| V16 new acquisitions | 6 | Archive verified; native energy and free energy separately recorded |
| V17 targeted acquisitions | 8 | Archive verified; Gamma-point method metadata |
| V17 rank-recovery label | 1 | Archive verified; old source declared 1×4×2, replacement used Gamma; k-point convergence unresolved |
| V19 residual-mode acquisitions | 6 | Archive verified; native energy and free energy recorded separately |

The stored frame metadata reports 36/36 as fully periodic TTT and records Gamma-point PW-PBE 500 eV for all rows. The method strings split into 27 generic “Gamma single point” entries and 9 explicit “Gamma Fermi0.1” entries. The V16 audit still labels overall training method consistency **UNKNOWN**, because the inherited source details for 12 rows are unavailable. Train36 stores `REF_energy` in all rows, a `PW_PBE_energy_eV` alias in 13, and `free_energy_eV` in six; the original inherited rows do not retain their free-energy value in train.extxyz.

All dev3 rows instead use PBE/PW 500 eV, a 6×6×1 mesh, sigma 0.10 eV, xy dipole correction, and TTF slab PBC. Thus the training/development comparison changes k-point sampling, z periodicity, and slab geometry at once. The dev formula `C16Ag4Si8Ti24` also does not occur in train36. This prevents assigning the large energy/force gap to a single cause.

The training composition matrix has rank 4, but its singular values range from 152.759 to 0.397 (condition ratio ≈385). MACE estimated elemental offsets from train36, with corrections of −16.709 eV (C), +9.876 eV (Si), +6.428 eV (Ti), and +4.748 eV (Ag). The shared dev formula and near-constant −0.65 to −0.67 eV/atom energy errors make composition-baseline extrapolation plausible, but this is not proven as the sole cause.

## Next steps

1. Recover original GPAW source records for the 12 inherited frames. If source outputs cannot be recovered, preserve that uncertainty and decide which labels can safely remain in an energy-and-force training split; do not silently treat the method declaration as verified.
2. Keep dev3 out of fitting while repairing the candidate. Establish a consistent slab target method (including k-points, smearing, dipole, and vacuum) before comparing it with training labels.
3. Add source-diverse Ti₃SiC₂ backbone and Ag–Ti/Si/C interface labels with chemistry/composition closer to the dev structures. The existing pure Ti₃SiC₂ controls are COD 9009647-derived and retain control-only roles, so they cannot be used to claim source-independent validation against this dev3 family.
4. Find a second source-verified Ti₃SiC₂ parent. The local checkout contains two byte-different CIF files, but they both identify COD 9009647, the same article and crystallographic entry; after replacing the differing COD URL line, their full file bytes match. Unauthenticated Materials Project OPTIMADE and GitHub API requests both returned HTTP 403 through the network tunnel, so no qualified second parent was found in this pass. Do not synthesize an “independent” parent by perturbing this CIF.
5. Have B independently audit any new data entry. Then A can approve a versioned split, run the existing preflight, and only then repeat the fixed training recipe.

## Reproduction

Run from the repository root in the pinned MACE environment:

```bash
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/coordination/reports/window-a/v19_train36_energy_method_lineage_audit_20261023/audit_train36_energy_method_lineage.py
```

`audit.json` records the per-frame source category, energy keys, method metadata, PBC, input hashes, composition conditioning, and E0 offsets. The script does not read `test.extxyz` or any sealed label.
