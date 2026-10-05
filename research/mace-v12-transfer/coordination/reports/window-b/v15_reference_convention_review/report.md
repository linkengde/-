# Historical reference-energy convention review

B reviewed committed V11–V15 train/valid/test data, prior provenance ledger, the current residual acquisition runner, and three verified residual archives. V15 blind calculations/labels were not opened. No DFT, inference/training, MD/TTM or dataset changes occurred.

## Evidence result

**Bookkeeping is partially established, not uniformly verified.** Ten inherited V11 training rows preserve native energy/free_energy aliases. REF_energy equals native energy exactly for all ten. Fourteen inherited rows retain method declarations, but their original getter/convergence/setup evidence remains absent. Other explicit-method inherited rows were not newly certified. The three residual acquisitions have direct script/log/archive-to-REF linkage.

The residual runner explicitly calls `atoms.get_potential_energy()` without force_consistent, separately calls `atoms.get_forces()`, writes PW_PBE_energy_eV, and the reviewed builder maps that exact value to REF_energy. Archived energies match the log **Extrapolated** entry to its six-decimal precision, not **Free energy**. This archived numeric match supports the actual convention used; a method-name tag alone does not. No new calculation was run.

## Earliest inherited paired values

| Config | native energy − free_energy (eV/cell) | meV/atom | REF − native energy (eV) |
|---|---:|---:|---:|
| `Ti_gap_2p3` | 0.196650143 | 6.1453 | 0.0 |
| `Ti_gap_2p6` | 0.209564384 | 6.5489 | 0.0 |
| `Ti_gap_2p9` | 0.213123793 | 6.6601 | 0.0 |
| `Si_gap_2p6` | 0.120175979 | 6.0088 | 0.0 |
| `C_gap_2p6` | 0.241884875 | 12.0942 | 0.0 |
| `Si_gap_2p4` | 0.120360797 | 6.0180 | 0.0 |
| `C_gap_2p4` | 0.232513979 | 11.6257 | 0.0 |
| `Si_gap_2p8` | 0.120637511 | 6.0319 | 0.0 |
| `C_gap_2p8` | 0.244260604 | 12.2130 | 0.0 |
| `AgTi_2p55_periodic_PW_PBE_energy_force` | 0.635283401 | 12.9650 | 0.0 |

These aliases are duplicated values within a dataset, not independent original calculation outputs. V12–V15 copies are lineage evidence and must not be counted as repeated original runs. Serialization may drop native calculator results; absent pairs are unknown, not zero offset. The JSON retains every available version/split pair and all 29 current train rows.

## Verified residual values

| Label | REF eV/cell | log free energy eV/cell | approximate offset eV/cell | meV/atom |
|---|---:|---:|---:|---:|
| `AgC_residual_shell_v15_01` | -257.906079437 | -259.066578 | 1.160499 | 24.1771 |
| `AgSi_residual_shell_v15_01` | -208.578424129 | -209.495017 | 0.916593 | 17.9724 |
| `AgTi_residual_shell_v15_01` | -177.004343283 | -177.652938 | 0.648595 | 13.2366 |

Residual offsets use six-decimal log free energies, so their last decimals are approximate; do not treat rounded differences as full-precision calculator outputs. Logs and archive hashes are recorded in JSON.

## Smearing and sampling evidence

- Nine inherited gap frames declare PW-PBE 500 eV, k grid1x4x2 and smearing0.1 eV; smearing distribution is not recovered for those gap declarations.
- AgTi2.55 declares PW-PBE500, Gamma and Fermi0.1 eV. Four solid elemental references declare PW-PBE500 and Fermi0.1 eV but exact phase-specific k meshes are missing.
- Residual scripts and archived summaries record PW-PBE500, Gamma1x1x1, FermiDirac0.1 eV and GPAW version. Dataset/trace JSON separates raw source_method/dft_method declarations from original verification; absent settings remain absent.
- The excluded AgTi2.387 row has no direct original method or matched energy/force getter binding. This review does not reinstate it.

Different k meshes may be appropriate for different cells. They require numerical convergence review, not a blanket inconsistency verdict. Finite-smearing extrapolated energy and electronic free energy are distinct targets. Separately requested forces are not a proof that every stored force is the derivative of the chosen extrapolated target; that finite-smearing convention should be documented and checked for future references. Energy offsets alone do not establish a cause of vector force errors.

## Recommendations

1. Preserve existing REF labels and name their convention as native/extrapolated where evidenced. Never replace them silently with free_energy or apply an entropy offset inferred from another frame.
2. Future matched labels should record force_consistent argument, both full-precision native energy and free_energy, units, smearing distribution/width, XC/cutoff/k-point grid/offset, software/PAW versions, convergence criteria and input/output hashes. State which value feeds REF_energy and what thermodynamic convention forces represent.
3. Recover original historical logs/inputs first; otherwise A chooses exact-geometry matched relabeling with ordered IDs/cell/PBC and consistent target convention. Do not infer absent elemental k grids or use zero for missing offsets.
4. Confirm k-point/smearing convergence for actual cells; do not force every differently shaped periodic cell to Gamma merely to match a string. Keep method declarations separate from verified run evidence.
5. Keep future acquisition convention fixed for an A-reviewed numerical screening cycle; assess convention sensitivity with a separate documented design rather than tuning against blind results.

Input SHA256, per-row evidence and explicit unknown fields are in convention_review.json. SHA256SUMS.txt covers the compact deliverables and all accessed committed evidence.

## Additional inherited aliases first appearing in V12

Four explicitly declared historical training rows also preserve paired native aliases in V12. Their corresponding current V15 REF energies remain identical; original getter/convergence provenance is not newly certified. All 14 unique paired historical train records and the full version/split inventory are in JSON.

| Config | native energy − free_energy eV/cell | meV/atom | REF − native energy eV |
|---|---:|---:|---:|
| `AgSi_d2p2_periodic_PW_PBE_energy_force` | 0.966011144 | 18.9414 | 0.000000000 |
| `AgSi_d2p8_periodic_PW_PBE_energy_force` | 0.965766202 | 18.9366 | 0.000000000 |
| `AgC_d2p0_periodic_PW_PBE_energy_force` | 1.138065990 | 23.7097 | 0.000000000 |
| `AgC_d2p4_periodic_PW_PBE_energy_force` | 1.143411008 | 23.8211 | 0.000000000 |
