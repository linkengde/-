# Ag–Ti₃SiC₂ cross-potential audit

Date: 2026-10-15
Branch: `cloud-e-b1-dft`
Scope: retain the Ti₃SiC₂ TSCBOP and Ag EAM unchanged; assess only Ag–Ti, Ag–Si and Ag–C cross interactions. This audit is read-only with respect to all source structures, DFT archives and potential files.

## Decision context

The approved model role is a solid Ti₃SiC₂ heat-conducting skeleton with Ag as the phase that melts and resolidifies. Ti₃SiC₂ melting, impact, crater and splash behavior are out of scope. The previously reported TSCBOP force mismatch is retained as a limitation of the backbone; it is not a target for the Ag cross terms to cancel.

The accessible checkout was searched for reports explicitly named Stage71, Stage72 or Stage73 and for matching commit subjects. None are present in this local branch/ref set. Network access to `origin` was unavailable during this task, so remote-only Stage71–73 material could not be ruled out. Equivalent available evidence was reviewed instead: E's 135-atom component-force attribution, A's candidate02 static diagnosis/decomposition, C's v19 reference export, and the existing cross-potential reports. Their hashes and verification outcomes are recorded in `DFT_INTERFACE_DATA_INVENTORY.md`.

## Existing calibration-v0 implementation

The reviewed include is `research/mace-v12-transfer/coordination/reports/window-e/low_mismatch_217_dft_lite_20261013/minimal_135atom_model/cross_potential_validation_run0_20261014/potential.interface_calibration_v0.inc`, SHA-256 `34deced5e08174ba369290b87f5c3bfd2b284b272991adf302b2b0d28c12815f`. The independently archived copy in `window-b/current_heating_potential_handoff_20261009/` has the same SHA. It declares `units metal`, `pair_style hybrid eam tersoff morse 8.0`, atom types 1=Ti, 2=Si, 3=C, 4=Ag, with no `pair_modify shift yes` (unshifted truncated Morse). The Morse convention is `D0` in eV, `alpha` in Å⁻¹, `r0` in Å, cutoff in Å.

| Pair | D0 (eV) | alpha (Å⁻¹) | r0 (Å) | Cutoff (Å) |
|---|---:|---:|---:|---:|
| Ag–Ti | 0.94814194 | 2.00000000 | 2.33708864 | 8.0 |
| Ag–Si | 0.93944948 | 2.00000000 | 2.40407837 | 8.0 |
| Ag–C | 1.35077230 | 2.00000000 | 2.15255924 | 8.0 |

Ag–Ag is routed to `Ag_u3.eam` (SHA-256 `cb3d605cd27ec6efb50c4dd7462c2ab3d967e63775dc4400de9134315a44d807`); Ti/Si/C is routed to `TSCBOP.tersoff` (SHA-256 `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`). The run0 manifest records `ZBL_active=false`. No source potential was edited.

The include comment traces the Morse starting values to a single-plane HCP surrogate fitted to reported Ag(111)/Ti₃SiC₂(0001) adhesion targets. This is not a pair-resolved DFT force dataset, a matched distance curve, or a validation of the three Morse parameter triplets. The `window-b` manifest identifies this as candidate02 `calibration-v0` and explicitly says Stage69 identity is unconfirmed. The local checkout does not establish that this file is the historical Stage67 production package; no claim of identity is made. Stage67 files and results were not touched.

## Existing cross-potential checks

The candidate02 static DFT comparison already scored 16 cases: six signed local Ag–Ti/Si/C perturbation points plus ten numerical controls of one frozen AgSi slab. It used LAMMPS 22Jul2025 update4, `metal` units and run0. All six signed points came from reused training-family references, not an untouched validation family. Existing total-force-vector RMSEs are 1.870–2.607 eV/Å, and the marked-pair force errors are 2.618–5.606 eV/Å. For the paired local energy differences, calibration-v0 gives 0.5093 eV (Ag–C), 0.3212 eV (Ag–Si) and 0.3431 eV (Ag–Ti), versus DFT force-consistent free-energy differences of 1.1086, 0.02763 and 0.20657 eV respectively. These are local signed-mode comparisons, not normal-separation curves.

The candidate02 decomposition report found a 0.906–0.909 eV/Å irreducible vector-force lower bound for two Ag–C diagnostic structures after projecting the DFT-minus-intraphase residual onto all available Ag–C radial directions inside the 8 Å cutoff. This permits arbitrary independent radial force amplitudes per pair, more freedom than a shared Morse law. Therefore the current central Ag–C cross-force basis cannot meet the cited 0.05 eV/Å vector gate on those structures while the intraphase terms and pair basis remain fixed. This is a scoped model-form limitation for those inputs, not proof that every cross-potential form or physical configuration must fail. Ag–Ti and Ag–Si geometric lower bounds near zero are only necessary-condition results, not validation.

The E 135-atom component audit gives one matched frozen-geometry effective cross-force comparison: the existing Morse vs DFT effective cross-force vector RMSE is 1.3513 eV/Å over all atoms and 1.5178 eV/Å over atoms within 8 Å. By receiving atom species, RMSE is 2.1496 eV/Å for Ag, 0.3316 for Ti, 2.7422 for Si and 0.0946 for C. These species totals cannot uniquely assign each residual to one Ag–X pair because each atom's force combines every nearby cross pair. They also represent one Si-terminated interface geometry only.

## Audit conclusion

Existing work establishes that calibration-v0 is not quantitatively validated and has large local force/energy-response discrepancies. The 135-atom DFT-A/DFT-B phase references provide a valuable single-geometry cross-force and interaction-energy anchor, but the combined evidence has no clean, untouched family of matched energy/force validation data and no sufficient distance/registry/termination coverage to identify three independent Morse triplets. The Ag–C lower bound further cautions that fitting the existing radial form may not satisfy the prior strict force gate. Do not run an optimizer, select parameters from these reused structures, or label a candidate production-ready.

## Integrity and activity

The 132-file C v19 DFT reference bundle passed its embedded SHA/file-set and SCF checks. E's component-attribution and run0 archives, plus A's candidate02 diagnosis/decomposition archives, passed their local `SHA256SUMS.txt` checks. No DFT, LAMMPS, MD, minimization, TTM or fitting calculation was run for this audit. No Stage67 or calibration-v0 source file was changed.
